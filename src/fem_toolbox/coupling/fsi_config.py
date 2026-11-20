"""Strict fsi_schema_version=1 TOML adapter and run artifacts (#93).

File I/O stays in this adapter; the coupling solver itself reads no files.
The TOML references an existing FluidDomain JSON and an existing elasticity
Model TOML (Gmsh mesh), both relative to the configuration directory. Every
coupling value is mandatory; there are no silent defaults.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NoReturn

import numpy as np

from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.rans import RansResult
from fem_toolbox.cfd.sa_model import _sa_constants
from fem_toolbox.coupling.fsi import (
    FsiCouplingError,
    FsiIteration,
    StationaryFsiProblem,
    StationaryFsiResult,
    solve_stationary_fsi,
)
from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.model import Model

__all__ = [
    "FsiConfig",
    "load_fsi_config",
    "run_fsi_config",
    "write_fsi_artifacts",
]

_COUPLING_KEYS = {
    "wall": str,
    "length_scale_m": float,
    "velocity_scale_m_s": float,
    "external_pressure_Pa": float,
    "relaxation": float,
    "max_coupling_iterations": int,
    "divergence_window": int,
    "span_origin_m": float,
    "preserve_inlet_flow": bool,
}
_COMMON_OPTIONS = {
    "force_N_m3",
    "pressure_reference_Pa",
    "convection_scheme",
    "velocity_relaxation",
    "pressure_relaxation",
    "max_iterations",
    "max_simple_steps",
    "max_unknowns",
}
_OPTIONS = {
    "laminar": frozenset(_COMMON_OPTIONS),
    "sa-noft2": frozenset(_COMMON_OPTIONS | {"sa_relaxation"}),
}
HISTORY_COLUMNS = (
    ("iteration", "iteration"),
    ("displacement_residual_m", "displacement_residual_m"),
    ("relative_displacement_residual", "relative_displacement_residual"),
    ("load_change_N", "load_change_N"),
    ("fluid_force_y_N", None),
    ("external_pressure_force_y_N", None),
    ("structure_reaction_y_N", None),
    ("transfer_relative_error", "transfer_relative_error"),
    ("maximum_displacement_ratio", "maximum_displacement_ratio"),
    ("maximum_displacement_gradient", "maximum_displacement_gradient"),
    ("minimum_area_ratio", "minimum_area_ratio"),
    ("fluid_iterations", "fluid_iterations"),
)


def _fail(path: str, reason: str) -> NoReturn:
    raise ModelSchemaError(path=path, reason=reason)


def _table(
    value: object,
    path: str,
    required: set[str],
    optional: frozenset[str] | set[str] = frozenset(),
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        _fail(path, "expected table")
    keys = {str(key) for key in value}
    if not required <= keys or not keys <= required | set(optional):
        _fail(
            path,
            f"required keys {sorted(required)}; optional {sorted(optional)}; got {sorted(keys)}",
        )
    return value


def _relative(root: Path, name: object, path: str, suffix: str) -> Path:
    if not isinstance(name, str) or not name or Path(name).is_absolute():
        _fail(path, "expected a relative file path")
    target = (root / str(name)).resolve()
    if not target.is_relative_to(root):
        _fail(path, "path escapes configuration directory")
    if target.suffix.lower() != suffix:
        _fail(path, f"expected {suffix} file")
    if not target.is_file():
        _fail(path, f"missing file {name}")
    return target


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _typed(value: object, kind: type, path: str) -> Any:
    if kind is float:
        if isinstance(value, bool) or not isinstance(value, int | float):
            _fail(path, "expected a real number")
        if not np.isfinite(value):
            _fail(path, "expected a finite number")
        return float(value)
    if type(value) is not kind:
        _fail(path, f"expected {kind.__name__}")
    return value


@dataclass(frozen=True, slots=True)
class FsiConfig:
    """Validated problem plus the exact input files and their SHA-256."""

    problem: StationaryFsiProblem
    source: Path
    document: Mapping[str, Any]
    input_files: Mapping[str, str]


def load_fsi_config(path: str | Path) -> FsiConfig:
    """Read a strict fsi_schema_version=1 TOML and build the coupling problem."""
    source = Path(path).resolve()
    if source.suffix.lower() != ".toml":
        _fail("file", "expected .toml")
    try:
        with source.open("rb") as stream:
            document = tomllib.load(stream)
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise ModelSchemaError(path="file", reason=str(error)) from error
    root = source.parent
    _table(
        document, "document", {"fsi_schema_version", "fluid", "structure", "coupling"}
    )
    version = document["fsi_schema_version"]
    if type(version) is not int or version != 1:
        _fail("fsi_schema_version", "expected integer 1")
    fluid = _table(
        document["fluid"],
        "fluid",
        {"domain", "model"},
        {"initial_nu_tilde", "options"},
    )
    model_name = fluid["model"]
    if model_name not in _OPTIONS:
        _fail("fluid.model", f"one of {sorted(_OPTIONS)} required")
    files: dict[str, Path] = {}
    files["fluid.domain"] = _relative(root, fluid["domain"], "fluid.domain", ".json")
    domain = FluidDomain.from_file(files["fluid.domain"])
    options = dict(fluid.get("options", {}))
    if not isinstance(fluid.get("options", {}), Mapping):
        _fail("fluid.options", "expected table")
    unknown = sorted(set(options) - _OPTIONS[model_name])
    if unknown:
        _fail("fluid.options", f"unsupported for {model_name}: {unknown}")
    nu_tilde: Any = None
    if model_name == "sa-noft2":
        if "initial_nu_tilde" not in fluid:
            _fail("fluid.initial_nu_tilde", "explicit SA start required")
        raw = fluid["initial_nu_tilde"]
        if isinstance(raw, str):
            files["fluid.initial_nu_tilde"] = _relative(
                root, raw, "fluid.initial_nu_tilde", ".json"
            )
            loaded = json.loads(files["fluid.initial_nu_tilde"].read_text("utf-8"))
            _table(loaded, "fluid.initial_nu_tilde", {"nu_tilde_m2_s"})
            nu_tilde = np.asarray(loaded["nu_tilde_m2_s"], dtype=float)
            if nu_tilde.shape != (len(domain.mesh.quad4),):
                _fail("fluid.initial_nu_tilde", "expected one value per fluid cell")
        else:
            nu_tilde = _typed(raw, float, "fluid.initial_nu_tilde")
    elif "initial_nu_tilde" in fluid:
        _fail("fluid.initial_nu_tilde", "only valid for sa-noft2")
    structure = _table(document["structure"], "structure", {"model", "facet_set"})
    files["structure.model"] = _relative(
        root, structure["model"], "structure.model", ".toml"
    )
    model = Model.from_file(files["structure.model"])
    mesh_name = tomllib.loads(files["structure.model"].read_text("utf-8"))["mesh"][
        "path"
    ]
    files["structure.mesh"] = (files["structure.model"].parent / mesh_name).resolve()
    facet_set = _typed(structure["facet_set"], str, "structure.facet_set")
    coupling = _table(document["coupling"], "coupling", set(_COUPLING_KEYS))
    values = {
        key: _typed(coupling[key], kind, f"coupling.{key}")
        for key, kind in _COUPLING_KEYS.items()
    }
    problem = StationaryFsiProblem(
        fluid_domain=domain,
        wall_name=values.pop("wall"),
        structure=model,
        facet_set=facet_set,
        fluid_model=model_name,
        initial_nu_tilde_m2_s=nu_tilde,
        fluid_options=options,
        **values,
    )
    hashes = {"config": _sha256(source)}
    for name, file in files.items():
        hashes[name] = _sha256(file)
    return FsiConfig(problem, source, document, hashes)


def _versions() -> dict[str, str]:
    """Versions of the code actually imported, not of installed metadata."""
    import scipy  # type: ignore[import-untyped]

    from fem_toolbox import __version__

    return {
        "fem_toolbox": __version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "python": platform.python_version(),
    }


def _json_ready(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _json_ready(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [_json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None if np.isnan(value) else ("inf" if value > 0 else "-inf")
    return value


def _write_json(path: Path, document: Any) -> None:
    with path.open("x", encoding="utf-8") as stream:
        stream.write(
            json.dumps(_json_ready(document), sort_keys=True, indent=2, allow_nan=False)
            + "\n"
        )


def _history_rows(history: tuple[FsiIteration, ...]) -> list[list[Any]]:
    rows = []
    for row in history:
        values = []
        for name, attribute in HISTORY_COLUMNS:
            if attribute is not None:
                values.append(getattr(row, attribute))
            elif name == "fluid_force_y_N":
                values.append(float(row.fluid_force_N[1]))
            elif name == "external_pressure_force_y_N":
                values.append(float(row.external_pressure_force_N[1]))
            else:
                values.append(float(row.structure_reaction_N[1]))
        rows.append(values)
    return rows


def write_fsi_artifacts(
    output: str | Path,
    config: FsiConfig,
    outcome: StationaryFsiResult | FsiCouplingError,
) -> Path:
    """Write report.json, history.csv and, on success, field JSON files.

    The output directory must not exist. All quantities carry SI unit suffixes.
    A failed run still records kind, reason, history and the last iterate.
    """
    target = Path(output)
    if target.exists():
        raise FileExistsError(f"output already exists: {target}")
    problem = config.problem
    report: dict[str, Any] = {
        "fsi_schema_version": 1,
        "configuration": config.document,
        "input_files_sha256": config.input_files,
        "versions": _versions(),
    }
    if problem.fluid_model == "sa-noft2":
        report["sa_noft2_constants"] = _sa_constants()
    history = outcome.history
    if isinstance(outcome, FsiCouplingError):
        report["status"] = "failed"
        report["failure"] = {"kind": outcome.kind, "reason": str(outcome)}
        report["last_displacement_m"] = outcome.displacement_m
    else:
        result = outcome
        flow = (
            result.fluid.flow if isinstance(result.fluid, RansResult) else result.fluid
        )
        structure = result.structure
        final = history[-1]
        report.update(
            status="converged",
            iterations=result.iterations,
            input_sha256=result.input_sha256,
            settings=dict(result.settings),
            mesh={
                "reference_mesh_sha256": result.mesh_motion.reference_mesh_sha256,
                "mesh_sha256": result.mesh_motion.mesh_sha256,
                "domain_sha256": result.mesh_motion.domain_sha256,
                "interface_sha256": result.mesh_motion.interface_sha256,
                "minimum_area_ratio": result.mesh_motion.minimum_area_ratio,
                "maximum_relative_laplace_residual": result.mesh_motion.maximum_relative_laplace_residual,
            },
            coupling={
                "displacement_residual_m": final.displacement_residual_m,
                "load_change_N": final.load_change_N,
                "fluid_force_N": final.fluid_force_N,
                "external_pressure_force_N": final.external_pressure_force_N,
                "structure_reaction_N": final.structure_reaction_N,
                "transfer_force_error_N": final.transfer_force_error_N,
                "transfer_moment_error_N_m": final.transfer_moment_error_N_m,
                "transfer_relative_error": final.transfer_relative_error,
                "maximum_displacement_ratio": final.maximum_displacement_ratio,
                "maximum_displacement_gradient": final.maximum_displacement_gradient,
            },
            fluid={
                "iterations": result.fluid.iterations,
                "inlet_mass_flow_kg_s": flow.inlet_mass_flow_kg_s,
                "outlet_mass_flow_kg_s": flow.outlet_mass_flow_kg_s,
                "global_mass_balance_kg_s": flow.global_mass_balance_kg_s,
                "global_momentum_balance_N": flow.global_momentum_balance_N,
                "power_balance_W": flow.power_balance_W,
                "maximum_momentum_residual_N_m3": float(
                    np.max(np.abs(flow.momentum_residual_N_m3))
                ),
                "maximum_mass_residual_kg_m3_s": float(
                    np.max(np.abs(flow.mass_residual_kg_m3_s))
                ),
                "provenance": result.fluid.provenance.to_dict(),
            },
            structure={
                "relative_residual": structure.relative_residual,
                "force_balance_N": structure.force_balance,
                "moment_balance_N_m": structure.moment_balance,
            },
        )
        if isinstance(result.fluid, RansResult):
            report["fluid"]["maximum_y_plus"] = float(np.max(result.fluid.wall_y_plus))
            report["fluid"]["maximum_sa_residual_m2_s2"] = float(
                np.max(np.abs(result.fluid.sa_residual_m2_s2))
            )
    target.mkdir(parents=True, exist_ok=False)
    _write_json(target / "report.json", report)
    with (target / "history.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow([name for name, _ in HISTORY_COLUMNS])
        writer.writerows(_history_rows(history))
    if isinstance(outcome, StationaryFsiResult):
        result = outcome
        domain = result.mesh_motion.domain
        flow = (
            result.fluid.flow if isinstance(result.fluid, RansResult) else result.fluid
        )
        fluid_fields: dict[str, Any] = {
            "points_m": domain.mesh.points,
            "quad4": domain.mesh.quad4,
            "cell_centres_m": domain.mesh.cell_centres_m,
            "mesh_displacement_m": result.mesh_motion.displacement_m,
            "velocity_m_s": flow.velocity_m_s,
            "pressure_Pa": flow.pressure_Pa,
            "mass_flux_kg_s": flow.mass_flux_kg_s,
            "wall_traction_Pa": flow.wall_traction_Pa,
            "wall_boundary_indices": domain.wall_quality.boundary_indices,
        }
        if isinstance(result.fluid, RansResult):
            fluid_fields.update(
                nu_tilde_m2_s=result.fluid.nu_tilde_m2_s,
                eddy_viscosity_m2_s=result.fluid.eddy_viscosity_m2_s,
                wall_y_plus=result.fluid.wall_y_plus,
            )
        mesh = problem.structure.elasticity_inputs()["mesh"]
        structure_fields = {
            "reference_points_m": mesh.points,
            "tetra4": mesh.tetra4,
            "geometry_displacement_m": result.geometry_displacement_m,
            "structure_displacement_m": result.structure_displacement_m,
            "interface_nodal_forces_N": result.loads.nodal_forces_N,
            "node_ids": result.structure.node_ids,
            "reactions_N": result.structure.reactions,
            "stress_Pa": result.structure.stress,
        }
        _write_json(target / "fluid_fields.json", fluid_fields)
        _write_json(target / "structure_fields.json", structure_fields)
    return target


def run_fsi_config(path: str | Path, output: str | Path) -> StationaryFsiResult:
    """Load, refuse an existing output, solve and always write artifacts.

    A coupling stop writes its failure artifacts and is then re-raised.
    """
    if Path(output).exists():
        raise FileExistsError(f"output already exists: {output}")
    config = load_fsi_config(path)
    try:
        result = solve_stationary_fsi(config.problem)
    except FsiCouplingError as error:
        write_fsi_artifacts(output, config, error)
        raise
    write_fsi_artifacts(output, config, result)
    return result
