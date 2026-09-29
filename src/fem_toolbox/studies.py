"""Versioned, local stationary parameter studies with verified result reuse."""

from __future__ import annotations

import itertools
import json
import math
import os
import re
import tomllib
from collections.abc import Mapping
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
from tempfile import NamedTemporaryFile, TemporaryDirectory
from typing import Any, NoReturn

from fem_toolbox import __version__
from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.io import export_result_bundle, read_gmsh
from fem_toolbox.model import Model
from fem_toolbox.provenance import PROVENANCE_SCHEMA_VERSION

_STATIONARY = {
    "heat",
    "elasticity",
    "electrostatics",
    "piezoelectric",
    "thermopiezoelectric",
}
# Changing orchestration/default interpretation requires a revision bump.
_STUDY_REVISION = 1


def _fail(path: str, reason: str) -> NoReturn:
    raise ModelSchemaError(path=path, reason=reason)


def _json(value: object) -> str:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ModelSchemaError(
            path="study", reason="expected finite JSON data"
        ) from error


def _hash(value: object) -> str:
    return sha256(_json(value).encode()).hexdigest()


def _table(value: Any, keys: set[str], path: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) != keys:
        _fail(path, f"expected exactly {sorted(keys)}")
    return dict(value)


def _inside(root: Path, name: Any) -> Path:
    if not isinstance(name, str) or not name or Path(name).is_absolute():
        _fail("path", "expected a nonempty relative path")
    path = (root / name).resolve()
    if not path.is_relative_to(root):
        _fail("path", "path escapes input directory")
    return path


def _toml(path: Path) -> dict[str, Any]:
    try:
        with path.open("rb") as stream:
            return tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise ModelSchemaError(path=str(path), reason=str(error)) from error


def _parent(document: dict[str, Any], target: str) -> tuple[Any, str | int]:
    parts = target.split(".")
    node: Any = document
    try:
        for part in parts[:-1]:
            if isinstance(node, list) and not re.fullmatch(r"0|[1-9][0-9]*", part):
                raise KeyError(part)
            node = node[int(part)] if isinstance(node, list) else node[part]
        key: str | int = int(parts[-1]) if isinstance(node, list) else parts[-1]
        if isinstance(key, int) and key < 0:
            raise KeyError(key)
        node[key]
    except (KeyError, IndexError, TypeError, ValueError) as error:
        raise ModelSchemaError(
            path=target, reason="unknown parameter target"
        ) from error
    return node, key


def _target_unit(document: dict[str, Any], target: str) -> str:
    parent, key = _parent(document, target)
    parts = target.split(".")
    if not isinstance(parent, dict):
        _fail(target, "target an entire scalar/vector/tensor property, not a component")
    unit: Any = None
    if parts == ["mesh", "path"]:
        return "1"
    if parts == ["mesh", "length_scale_to_m"]:
        return "1"
    if parts[0] == "solver" and key in {
        "method",
        "block_scaling",
        "electric_thermal_strain",
    }:
        return "1"
    if parts[0] == "materials" and len(parts) == 4:
        if not isinstance(parent.get("units", {}), dict):
            _fail(target, "material units must be a table")
        if key in {"conductivity", "young", "stiffness", "permittivity"}:
            unit = parent.get("units", {}).get(key, parent.get("unit"))
        elif key in {
            "compliance",
            "stress_coupling",
            "strain_coupling",
            "permittivity_at_stress",
        }:
            unit = parent.get("units", {}).get(key)
        elif key in {"poisson", "axes"}:
            unit = "1"
        elif key in {"alpha", "principal_alpha"}:
            unit = parent.get("alpha_unit")
        elif key == "reference_temperature":
            unit = parent.get("reference_temperature_unit")
    if parts[0] in {"boundary", "loads"}:
        if key in {"value", "coefficient"}:
            unit = parent.get("unit")
        elif key == "ambient_temperature":
            unit = parent.get("ambient_unit")
        elif key in {"potential", "charge"}:
            unit = parent.get(f"{key}_unit")
    if not isinstance(unit, str) or unit not in {
        "1",
        "K",
        "m",
        "V",
        "C",
        "Pa",
        "1/Pa",
        "W/(m K)",
        "F/m",
        "C/m^2",
        "C/N",
        "1/K",
        "W/m^2",
        "W/(m^2 K)",
        "N/m^2",
        "W/m^3",
        "N/m^3",
        "C/m^3",
    }:
        _fail(target, "unsupported parameter property or non-SI source unit")
    return str(unit)


def _value(value: Any, template: Any, target: str) -> None:
    if isinstance(template, bool):
        valid = isinstance(value, bool)
    elif isinstance(template, str):
        valid = isinstance(value, str) and bool(value)
    elif isinstance(template, list):
        valid = isinstance(value, list) and len(value) == len(template)
        if valid:
            for v, t in zip(value, template, strict=True):
                _value(v, t, target)
    else:
        valid = type(value) in (int, float) and math.isfinite(value)
    if not valid:
        _fail(target, "value must match the finite scalar/array type and shape")


class Study:
    """A detached schema-1 manifest; input files are reread on each resolution.

    ``manifest`` has study_schema_version=1, a relative TOML ``model`` path,
    and ``parameters`` (name, target, SI unit, finite nonempty values). Parameters
    form a Cartesian product in declaration order (last parameter varies fastest).
    Targets address existing whole properties, e.g. boundary.dirichlet.0.value.
    ``root`` is the manifest directory; the base model must stay inside it.
    Schema/IO errors raise ModelSchemaError. Physical case errors are recorded by
    run(). No expressions, transient models, remote execution or new meshing.
    """

    def __init__(self, manifest: Mapping[str, object], *, root: str | Path) -> None:
        self._manifest: dict[str, Any] = json.loads(_json(dict(manifest)))
        self._root = Path(root).resolve()
        _table(self._manifest, {"study_schema_version", "model", "parameters"}, "study")
        version = self._manifest["study_schema_version"]
        if type(version) is not int or version != 1:
            _fail("study_schema_version", "expected integer 1")

    @classmethod
    def from_file(cls, path: str | Path) -> Study:
        """Read a UTF-8 .toml manifest with paths relative to its directory."""
        source = Path(path).resolve()
        if source.suffix.lower() != ".toml":
            _fail("study.file", "expected .toml")
        return cls(_toml(source), root=source.parent)

    def resolve(self) -> list[dict[str, Any]]:
        """Return detached deterministic case records, without solving any model.

        Case IDs hash named parameter values/units/targets; input fingerprints
        additionally cover complete model data, mesh bytes, scale, explicit
        solver settings, schema/toolbox/dependency versions and study revision.
        A changed base input invalidates reuse while retaining the case ID.
        """
        model_path = _inside(self._root, self._manifest["model"])
        if model_path.suffix.lower() != ".toml":
            _fail("model", "expected .toml")
        document = _toml(model_path)
        if (
            not isinstance(document.get("physics"), str)
            or document["physics"] not in _STATIONARY
        ):
            _fail("physics", "only stationary Model physics are supported")
        parameters = self._manifest["parameters"]
        if not isinstance(parameters, list):
            _fail("parameters", "expected array of parameter tables")
        names: set[str] = set()
        targets: set[str] = set()
        for p in parameters:
            _table(p, {"name", "target", "unit", "values"}, "parameters")
            name, target = p["name"], p["target"]
            if not isinstance(name, str) or not re.fullmatch(
                r"[A-Za-z][A-Za-z0-9_]*", name
            ):
                _fail("parameters.name", "expected an ASCII identifier")
            if not isinstance(target, str) or name in names or target in targets:
                _fail("parameters", "duplicate name/target or invalid target")
            names.add(name)
            targets.add(target)
            if p["unit"] != _target_unit(document, target):
                _fail(target, "parameter unit must match the SI model property")
            values = p["values"]
            if not isinstance(values, list) or not values:
                _fail(target, "expected a nonempty values array")
            parent, key = _parent(document, target)
            for v in values:
                _value(v, parent[key], target)
            if len({_json(v) for v in values}) != len(values):
                _fail(target, "duplicate values")
        versions = {
            "toolbox": __version__,
            "numpy": version("numpy"),
            "scipy": version("scipy"),
            "meshio": version("meshio"),
            "provenance_schema": PROVENANCE_SCHEMA_VERSION,
            "study_revision": _STUDY_REVISION,
        }
        cases = []
        for combination in itertools.product(*(p["values"] for p in parameters)):
            resolved = json.loads(_json(document))
            selected = {}
            for p, value in zip(parameters, combination, strict=True):
                parent, key = _parent(resolved, p["target"])
                parent[key] = value
                selected[p["name"]] = {
                    "value": value,
                    "unit": p["unit"],
                    "target": p["target"],
                }
            mesh = resolved.get("mesh", {})
            if not isinstance(mesh, dict):
                _fail("mesh", "expected a mesh table")
            mesh_path = _inside(model_path.parent, mesh.get("path"))
            try:
                mesh_hash = sha256(mesh_path.read_bytes()).hexdigest()
            except OSError as error:
                raise ModelSchemaError(path="mesh.path", reason=str(error)) from error
            # Locations are references, not numerical inputs. Mesh content is hashed.
            canonical = json.loads(_json(resolved))
            canonical["mesh"].pop("path", None)
            fingerprints = {
                "model": _hash(canonical),
                "mesh": mesh_hash,
                "solver": _hash(
                    {
                        "physics": resolved["physics"],
                        "settings": resolved.get("solver", {}),
                        "versions": versions,
                    }
                ),
            }
            fingerprints["input"] = _hash(
                {"fingerprints": fingerprints, "versions": versions}
            )
            cases.append(
                {
                    "case_id": _hash(selected),
                    "parameters": selected,
                    "document": resolved,
                    "fingerprints": fingerprints,
                    "versions": versions,
                }
            )
        return cases

    def run(
        self, output: str | Path, *, max_cases: int | None = None
    ) -> dict[str, Any]:
        """Run/retry cases serially and return the persisted schema-1 catalog.

        ``max_cases`` limits newly attempted cases, not cached successes. None
        runs all cases. Case exceptions are persisted and subsequent cases run;
        KeyboardInterrupt/SystemExit persist ``interrupted`` and are re-raised.
        Catalog/locking errors propagate. A lock prevents concurrent writers;
        after a hard process kill remove .study.lock only after checking no
        writer remains. Partial attempt directories are retained, never reused.
        """
        if max_cases is not None and (type(max_cases) is not int or max_cases < 0):
            _fail("max_cases", "expected a nonnegative integer or None")
        cases = self.resolve()
        root = Path(output).resolve()
        root.mkdir(parents=True, exist_ok=True)
        lock = root / ".study.lock"
        try:
            stream = lock.open("x")
        except FileExistsError as error:
            raise FileExistsError(
                f"Study locked: {lock}; verify no writer before recovery"
            ) from error
        try:
            with stream:
                stream.write(str(os.getpid()))
            return self._run_locked(root, cases, max_cases)
        finally:
            lock.unlink()

    def _run_locked(
        self, root: Path, cases: list[dict[str, Any]], limit: int | None
    ) -> dict[str, Any]:
        path = root / "catalog.json"
        previous = (
            _load_catalog(path) if path.exists() else {"cases": [], "retired_cases": []}
        )
        old = {c["case_id"]: c for c in previous["cases"]}
        ids = {c["case_id"] for c in cases}
        catalog: dict[str, Any] = {
            "study_schema_version": 1,
            "cases": [],
            "retired_cases": previous["retired_cases"]
            + [c for c in previous["cases"] if c["case_id"] not in ids],
        }
        for case in cases:
            saved = old.get(case["case_id"])
            current = dict(
                case,
                status="pending",
                attempts=0,
                error=None,
                result=None,
                attempt_directory=None,
                history=[],
            )
            if saved is not None:
                current.update(
                    {
                        k: saved[k]
                        for k in (
                            "status",
                            "attempts",
                            "error",
                            "result",
                            "attempt_directory",
                            "history",
                        )
                    }
                )
                reason = None
                if saved["fingerprints"] != case["fingerprints"]:
                    reason = "input fingerprint changed"
                elif saved["status"] == "succeeded" and not _valid_result(root, saved):
                    reason = "missing, incomplete or corrupt result bundle"
                if reason:
                    current["history"] = saved["history"] + [
                        {
                            "status": "invalid",
                            "reason": reason,
                            "fingerprints": saved["fingerprints"],
                            "result": saved["result"],
                            "attempts": saved["attempts"],
                        }
                    ]
                    current.update(status="invalid", error=reason, result=None)
                elif saved["status"] == "running":
                    current.update(
                        status="interrupted",
                        error="previous writer stopped before completion",
                    )
            catalog["cases"].append(current)
        _write_catalog(path, catalog)
        attempted = 0
        for case in catalog["cases"]:
            if case["status"] == "succeeded":
                continue
            if limit is not None and attempted >= limit:
                break
            attempted += 1
            # Never overwrite an orphan from a terminated process.
            while True:
                case["attempts"] += 1
                relative = f"cases/{case['case_id']}/attempt-{case['attempts']:06d}"
                directory = _inside(root, relative)
                try:
                    directory.mkdir(parents=True, exist_ok=False)
                except FileExistsError:
                    continue
                break
            case.update(
                status="running", error=None, result=None, attempt_directory=relative
            )
            _write_catalog(path, catalog)
            try:
                document = case["document"]
                model_path = _inside(self._root, self._manifest["model"])
                mesh_path = _inside(model_path.parent, document["mesh"]["path"])
                raw = mesh_path.read_bytes()
                if sha256(raw).hexdigest() != case["fingerprints"]["mesh"]:
                    _fail(
                        "mesh.path",
                        "mesh changed after study resolution; resolve again",
                    )
                with TemporaryDirectory() as temporary:
                    snapshot = Path(temporary) / "mesh.msh"
                    snapshot.write_bytes(raw)
                    mesh = read_gmsh(
                        snapshot,
                        length_scale_to_m=document["mesh"]["length_scale_to_m"],
                    )
                model = Model(mesh=mesh, document=document)
                run = model.solve()
                fields, manifest = export_result_bundle(
                    directory / "result", model.mesh, run
                )
                case["result"] = {
                    "manifest": str(manifest.relative_to(root)),
                    "fields": str(fields.relative_to(root)),
                    "manifest_sha256": sha256(manifest.read_bytes()).hexdigest(),
                    "fields_sha256": sha256(fields.read_bytes()).hexdigest(),
                    "model_input_sha256": run.provenance.input_sha256,
                }
                case["status"] = "succeeded"
            except (KeyboardInterrupt, SystemExit) as error:
                case.update(
                    status="interrupted", error=type(error).__name__, result=None
                )
                _write_catalog(path, catalog)
                raise
            except Exception as error:
                case.update(
                    status="failed",
                    error=f"{type(error).__name__}: {error}",
                    result=None,
                )
            _write_catalog(path, catalog)
        return catalog


def _write_catalog(path: Path, catalog: dict[str, Any]) -> None:
    """Replace a catalog only after its complete new contents have been synced."""
    temporary: Path | None = None
    try:
        with NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, delete=False
        ) as stream:
            temporary = Path(stream.name)
            stream.write(_json(catalog) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _load_catalog(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        _table(data, {"study_schema_version", "cases", "retired_cases"}, "catalog")
        if (
            type(data["study_schema_version"]) is not int
            or data["study_schema_version"] != 1
        ):
            raise ValueError("unsupported catalog schema")
        if not isinstance(data["cases"], list) or not isinstance(
            data["retired_cases"], list
        ):
            raise ValueError("invalid case lists")
        ids = []
        for c in data["cases"]:
            if not isinstance(c, dict) or not re.fullmatch(
                "[0-9a-f]{64}", c["case_id"]
            ):
                raise ValueError("invalid case ID")
            ids.append(c["case_id"])
            if c["status"] not in {
                "pending",
                "running",
                "interrupted",
                "succeeded",
                "failed",
                "invalid",
            }:
                raise ValueError("invalid case status")
            if (
                type(c["attempts"]) is not int
                or c["attempts"] < 0
                or not isinstance(c["history"], list)
            ):
                raise ValueError("invalid attempt history")
            for key in ("fingerprints", "result", "error", "attempt_directory"):
                c[key]
        if len(set(ids)) != len(ids):
            raise ValueError("duplicate case IDs")
        return dict(data)
    except (ValueError, KeyError, TypeError) as error:
        raise ModelSchemaError(
            path="catalog", reason=f"corrupt catalog: {error}"
        ) from error


def _valid_result(root: Path, case: dict[str, Any]) -> bool:
    try:
        result = case["result"]
        manifest = _inside(root, result["manifest"])
        fields = _inside(root, result["fields"])
        if sha256(manifest.read_bytes()).hexdigest() != result["manifest_sha256"]:
            return False
        if sha256(fields.read_bytes()).hexdigest() != result["fields_sha256"]:
            return False
        data = json.loads(manifest.read_text(encoding="utf-8"))
        return bool(
            data["bundle_schema_version"] in (1, 2)
            and data["files"]["fields"] == fields.name
            and data["files"]["fields_sha256"] == result["fields_sha256"]
            and data["provenance"]["input_sha256"] == result["model_input_sha256"]
        )
    except (OSError, ValueError, KeyError, TypeError, ModelSchemaError):
        return False
