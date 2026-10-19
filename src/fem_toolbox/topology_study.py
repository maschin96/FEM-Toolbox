"""Declarative local topology studies with immutable-input and artifact hashes."""

from __future__ import annotations

import hashlib
import json
import tomllib
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from fem_toolbox.errors import ValidationError
from fem_toolbox.io import write_vtu
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.topology import TopologyProblem
from fem_toolbox.topology_optimization import (
    OptimizationOptions,
    OptimizationResult,
    optimize_topology,
)


def _json(value: Any) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode()


def _keys(value: Any, allowed: set[str], name: str) -> dict[str, Any]:
    if not isinstance(value, dict) or set(value) - allowed:
        raise ValidationError(f"unknown or invalid {name} keys", field="configuration")
    return value


def problem_document(problem: TopologyProblem) -> dict[str, Any]:
    """Detached effective SI inputs; no hidden Model parsing or host paths."""
    return {
        "points_m": problem.mesh.points.tolist(),
        "tetra4": problem.mesh.tetra4.tolist(),
        "triangle3": problem.mesh.triangle3.tolist(),
        "stiffness_Pa": [m.stiffness.tolist() for m in problem.materials],
        "cases": [
            {
                "name": c.name,
                "conditions": [
                    {
                        "dofs": np.asarray(b.dofs).tolist(),
                        "values_m": np.asarray(b.values).tolist(),
                    }
                    for b in c.conditions
                ],
                "body_force_N_m3": np.asarray(c.body_force).tolist(),
                "tractions_N_m2": {
                    str(k): np.asarray(v).tolist()
                    for k, v in (c.tractions or {}).items()
                },
                "point_loads_N": {
                    str(k): np.asarray(v).tolist()
                    for k, v in (c.point_loads or {}).items()
                },
            }
            for c in problem.cases
        ],
        "weights": problem.weights.tolist(),
        "filter_radius_m": problem.filter_radius,
        "penalty": problem.penalty,
        "minimum_stiffness": problem.minimum_stiffness,
        "fixed_solid": problem.solid.tolist(),
        "fixed_void": problem.void.tolist(),
        "random_seed": None,
    }


class TopologyStudy:
    """One explicit local optimization; output directories must be new."""

    def __init__(
        self,
        problem: TopologyProblem,
        volume_fraction: float,
        *,
        options: OptimizationOptions | None = None,
    ) -> None:
        problem.initial_density(volume_fraction)
        self.problem = problem
        self.volume_fraction = float(volume_fraction)
        self.options = options or OptimizationOptions()
        self.source_inputs: dict[str, str] = {}

    @classmethod
    def from_file(cls, path: str | Path) -> TopologyStudy:
        """Load schema-1 TOML referencing existing static Model TOML files."""
        path = Path(path).resolve()
        root = path.parent
        raw = _keys(
            tomllib.loads(path.read_text()),
            {"schema_version", "models", "weights", "design", "optimizer"},
            "study",
        )
        if raw.get("schema_version") != 1 or isinstance(
            raw.get("schema_version"), bool
        ):
            raise ValidationError(
                "topology study schema_version must be 1", field="schema"
            )
        names = raw.get("models")
        if (
            not isinstance(names, list)
            or not names
            or not all(isinstance(n, str) for n in names)
        ):
            raise ValidationError("nonempty Model path list required", field="models")
        models = []
        inputs = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()}
        for name in names:
            candidate = (root / name).resolve()
            if Path(name).is_absolute() or not candidate.is_relative_to(root):
                raise ValidationError(
                    "Model path escapes study directory", field="models"
                )
            model = Model.from_file(candidate)
            models.append(model)
            inputs[name] = hashlib.sha256(candidate.read_bytes()).hexdigest()
            mesh_path = (
                candidate.parent / model.to_document()["mesh"]["path"]
            ).resolve()
            if not mesh_path.is_relative_to(root):
                raise ValidationError(
                    "mesh path escapes study directory", field="models"
                )
            inputs[str(mesh_path.relative_to(root))] = hashlib.sha256(
                mesh_path.read_bytes()
            ).hexdigest()
        design = _keys(
            raw.get("design"),
            {
                "volume_fraction",
                "filter_radius",
                "penalty",
                "minimum_stiffness",
                "fixed_solid",
                "fixed_void",
                "solid_sets",
                "void_sets",
            },
            "design",
        ).copy()
        if "volume_fraction" not in design or "filter_radius" not in design:
            raise ValidationError(
                "volume_fraction and filter_radius required", field="design"
            )
        volume = design.pop("volume_fraction")
        for source, target in (
            ("solid_sets", "fixed_solid"),
            ("void_sets", "fixed_void"),
        ):
            selected = design.pop(source, [])
            if not isinstance(selected, list) or not all(
                isinstance(n, str) and n in models[0].mesh.volume_sets for n in selected
            ):
                raise ValidationError("unknown design volume set", field=source)
            explicit = design.get(target, [])
            if not isinstance(explicit, list):
                raise ValidationError("design indices must be a list", field=target)
            design[target] = explicit + [
                int(i) for name in selected for i in models[0].mesh.volume_sets[name]
            ]
        optimizer = _keys(
            raw.get("optimizer", {}),
            {"max_iterations", "max_evaluations", "tolerance"},
            "optimizer",
        )
        study = cls(
            TopologyProblem.from_models(models, weights=raw.get("weights"), **design),
            volume,
            options=OptimizationOptions(**optimizer),
        )
        study.source_inputs = inputs
        return study

    def run(self, output: str | Path) -> OptimizationResult:
        """Optimize and save schema-1 JSON/VTU plus SHA-256 inventory.

        Budget/solver status is preserved; no stopped state becomes successful.
        Invalid initial states raise before an output directory is created.
        """
        output = Path(output)
        if output.exists():
            raise FileExistsError(output)
        result = optimize_topology(
            self.problem, self.volume_fraction, options=self.options
        )
        output.mkdir(parents=True, exist_ok=False)
        document = {
            "schema": 1,
            "problem": problem_document(self.problem),
            "volume_limit": self.volume_fraction,
            "optimizer": asdict(self.options),
            "source_inputs": self.source_inputs,
        }
        content = _json(document)
        provenance = RunProvenance.from_input_bytes(
            content,
            solver_name="fem_toolbox.topology_study.TopologyStudy.run",
            solver_parameters={
                "optimizer": "scipy.optimize.SLSQP",
                "filter": "volume_weighted_linear",
                "schema": 1,
            },
        )
        (output / "input.json").write_bytes(content)
        state = result.final
        rho = state.physical_density
        volumes = self.problem.volumes
        record = {
            "schema": 1,
            "status": result.status,
            "message": result.message,
            "units": {
                "compliance": "J",
                "volume": "m^3",
                "density": "1",
                "stress": "Pa",
            },
            "density": state.density.tolist(),
            "physical_density": rho.tolist(),
            "compliance_J": state.objective,
            "initial_compliance_J": result.initial_objective,
            "volume_fraction": state.volume_fraction,
            "volume_m3": float(volumes @ rho),
            "volume_limit": result.volume_limit,
            "evaluations": result.evaluations,
            "history": list(result.history),
            "provenance": provenance.to_dict(),
            "intermediate_cell_volume_fraction": float(
                volumes @ ((rho > 0.1) & (rho < 0.9)) / volumes.sum()
            ),
            "gray_measure": float(volumes @ (4 * rho * (1 - rho)) / volumes.sum()),
            "states": [
                {
                    "name": case.name,
                    "compliance_J": float(np.sum(s.applied_loads * s.displacement)),
                    "strain_energy_J": s.energy,
                    "relative_residual": s.relative_residual,
                    "force_balance_N": s.force_balance.tolist(),
                    "moment_balance_N_m": s.moment_balance.tolist(),
                }
                for case, s in zip(self.problem.cases, state.states, strict=True)
            ],
            "experimental_validation": "not_performed",
        }
        (output / "result.json").write_bytes(_json(record))
        volume_mesh = Mesh(
            points=self.problem.mesh.points,
            tetra4=self.problem.mesh.tetra4,
            volume_sets=self.problem.mesh.volume_sets,
        )
        write_vtu(
            output / "design.vtu",
            volume_mesh,
            cell_data={
                "density": {"tetra4": state.density},
                "physical_density": {"tetra4": rho},
            },
        )
        receipts = {
            name: hashlib.sha256((output / name).read_bytes()).hexdigest()
            for name in ("input.json", "result.json", "design.vtu")
        }
        (output / "SHA256.json").write_bytes(_json({"schema": 1, "files": receipts}))
        return result


def read_topology_result(path: str | Path) -> dict[str, Any]:
    """Read a schema-1 result only after verifying all bound artifact bytes."""
    root = Path(path).resolve()
    receipts = json.loads((root / "SHA256.json").read_text())
    if receipts.get("schema") != 1 or set(receipts.get("files", {})) != {
        "input.json",
        "result.json",
        "design.vtu",
    }:
        raise ValidationError("invalid topology artifact inventory", field="schema")
    for name, digest in receipts["files"].items():
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != digest:
            raise ValidationError(f"export hash mismatch: {name}", field="artifact")
    result: dict[str, Any] = json.loads((root / "result.json").read_text())
    if result.get("schema") != 1:
        raise ValidationError("unsupported topology result schema", field="schema")
    expected = hashlib.sha256((root / "input.json").read_bytes()).hexdigest()
    if result.get("provenance", {}).get("input_sha256") != expected:
        raise ValidationError("topology provenance hash mismatch", field="provenance")
    return result
