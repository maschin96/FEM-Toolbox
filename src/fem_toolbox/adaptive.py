"""Controlled stationary Tet4 heat adaptation with inherited problem data."""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
from typing import Literal

import numpy as np

from fem_toolbox.assembly import indices
from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.mesh.refinement import (
    Refinement,
    refine_tet4,
    validate_refinement,
    validate_tet4,
)
from fem_toolbox.physics.heat import Convection, ScalarField
from fem_toolbox.post import HeatIndicator, heat_indicator, mark_bulk
from fem_toolbox.solvers.heat import HeatResult, solve_heat, thermal_dofs

Status = Literal["tolerance", "dof_budget", "step_budget", "stagnation", "error"]


@dataclass(frozen=True)
class AdaptiveOptions:
    """Indicator tolerance sqrt(W K), active DOFs and number of solved meshes.

    Stagnation means eta fails to fall by min_reduction over stagnation_window
    successive refinements. Tolerance is heuristic, never a certified error bound.
    """

    tolerance: float = 1e-3
    max_steps: int = 10
    max_dofs: int = 100000
    fraction: float = 0.5
    min_reduction: float = 0.01
    stagnation_window: int = 3

    def __post_init__(self) -> None:
        for name in ("max_steps", "max_dofs", "stagnation_window"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValidationError("expected positive integer", field=name)
        for name in ("tolerance", "fraction", "min_reduction"):
            value = getattr(self, name)
            if isinstance(value, bool) or not np.isfinite(value):
                raise ValidationError("expected finite real value", field=name)
        if (
            self.tolerance < 0
            or not 0 < self.fraction <= 1
            or not 0 <= self.min_reduction < 1
        ):
            raise ValidationError(
                "invalid tolerance, fraction or reduction", field="options"
            )


@dataclass(frozen=True)
class HeatProblem:
    """Stationary SI data. Boundary keys are INPUT triangle row indices.

    Dirichlet values (K) may be scalar or point callable. Unlisted facets are
    insulated. Cellwise material and all facet maps inherit exactly through
    refinement parents; callable fields are reevaluated at physical points.
    Point-only constraints are deliberately excluded from this workflow.
    """

    material: ThermalMaterial | Sequence[ThermalMaterial]
    dirichlet: Mapping[int, ScalarField]
    source: ScalarField = 0.0
    outward_fluxes: Mapping[int, ScalarField] = field(default_factory=dict)
    convection: Mapping[int, Convection] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.material, ThermalMaterial):
            object.__setattr__(self, "material", tuple(self.material))
        for name in ("dirichlet", "outward_fluxes", "convection"):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))

    def refined(self, refinement: Refinement) -> HeatProblem:
        materials = self.material
        if not isinstance(materials, ThermalMaterial):
            materials = tuple(materials[int(p)] for p in refinement.cell_parents)

        def inherit(mapping: Mapping[int, ScalarField]) -> dict[int, ScalarField]:
            return {
                i: mapping[int(p)]
                for i, p in enumerate(refinement.facet_parents)
                if int(p) in mapping
            }

        return HeatProblem(
            materials,
            inherit(self.dirichlet),
            self.source,
            inherit(self.outward_fluxes),
            {
                i: self.convection[int(p)]
                for i, p in enumerate(refinement.facet_parents)
                if int(p) in self.convection
            },
        )


def mesh_record(mesh: Mesh) -> dict[str, object]:
    """Canonical topology/selection record; hash includes external identity."""
    return {
        name: (
            dict(
                (key, value.tolist())
                for key, value in sorted(getattr(mesh, name).items())
            )
            if name.endswith("sets")
            else None
            if getattr(mesh, name) is None
            else getattr(mesh, name).tolist()
        )
        for name in (
            "points",
            "tetra4",
            "triangle3",
            "cell_tags",
            "facet_tags",
            "volume_sets",
            "facet_sets",
            "point_ids",
            "tetra4_ids",
            "triangle3_ids",
        )
    }


def mesh_hash(mesh: Mesh) -> str:
    return sha256(
        json.dumps(
            mesh_record(mesh), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


@dataclass(frozen=True)
class AdaptiveStage:
    """Solved immutable mesh/result plus decision and parent evidence."""

    mesh: Mesh
    result: HeatResult
    indicator: HeatIndicator
    marked: tuple[int, ...]
    cell_parents: tuple[int, ...] = ()
    facet_parents: tuple[int, ...] = ()

    def to_dict(self) -> dict[str, object]:
        return {
            "mesh": mesh_record(self.mesh),
            "mesh_sha256": mesh_hash(self.mesh),
            "dofs": len(self.result.node_ids),
            "cells": self.mesh.n_tetra4,
            "indicator": self.indicator.total,
            "cell_indicator": self.indicator.cell.tolist(),
            "marked": list(self.marked),
            "cell_parents": list(self.cell_parents),
            "facet_parents": list(self.facet_parents),
            "temperature": self.result.temperature.tolist(),
            "reactions": self.result.reactions.tolist(),
            "residual": self.result.residual.tolist(),
            "node_ids": self.result.node_ids.tolist(),
            "gradient": self.result.gradient.tolist(),
            "heat_flux": self.result.heat_flux.tolist(),
            "power_balance": self.result.power_balance,
            "source_power": self.result.source_power,
            "dirichlet_inflow": self.result.dirichlet_inflow,
            "prescribed_outflow": self.result.prescribed_outflow,
            "convection_outflow": self.result.convection_outflow,
            "relative_residual": self.result.relative_residual,
            "provenance": self.result.provenance.to_dict(),
        }


@dataclass(frozen=True)
class AdaptiveHeatResult:
    stages: tuple[AdaptiveStage, ...]
    status: Status
    diagnostic: str
    options: AdaptiveOptions

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": 1,
            "kind": "fem-toolbox-adaptive-heat",
            "status": self.status,
            "diagnostic": self.diagnostic,
            "options": dict(vars(self.options)),
            "stages": [stage.to_dict() for stage in self.stages],
        }

    def export(self, path: str | Path) -> Path:
        """Write deterministic strict JSON evidence (meshes and raw SI results).

        This records discrete results and solver input fingerprints. Rerunning
        callable data requires the original problem script, as for solve_heat.
        """
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_dict(), sort_keys=True, indent=2, allow_nan=False)
            + "\n",
            encoding="utf-8",
        )
        return target


def solve_adaptive_heat(
    mesh: Mesh,
    problem: HeatProblem,
    options: AdaptiveOptions | None = None,
    *,
    remesher: Callable[..., Refinement] = refine_tet4,
) -> AdaptiveHeatResult:
    """Solve-estimate-mark-refine until tolerance, budget, stagnation or error.

    No partial or overbudget mesh is accepted. Validation/solver/adapter failures
    return error with previous solved stages intact; programming errors propagate.
    The remesher hook must return exact parent inheritance and conforming Tet4
    geometry, obey max_points, and preserve domain/groups (see refinement ADR).
    """
    control = options or AdaptiveOptions()
    stages: list[AdaptiveStage] = []
    parents: Refinement | None = None

    def finish(status: Status, message: str) -> AdaptiveHeatResult:
        return AdaptiveHeatResult(tuple(stages), status, message, control)

    try:
        validate_tet4(mesh)
        while True:
            if thermal_dofs(mesh).size > control.max_dofs:
                return finish("dof_budget", "active thermal DOF budget reached")
            for mapping in (
                problem.dirichlet,
                problem.outward_fluxes,
                problem.convection,
            ):
                indices(list(mapping), mesh.n_triangle3)
            conditions = []
            dofs = thermal_dofs(mesh)
            for row, value in sorted(problem.dirichlet.items()):
                nodes = mesh.triangle3[row]
                values = [
                    value(point) if callable(value) else value
                    for point in mesh.points[nodes]
                ]
                conditions.append(Dirichlet(dofs.indices("T", nodes), values))
            result = solve_heat(
                mesh,
                problem.material,
                conditions,
                source=problem.source,
                outward_fluxes=problem.outward_fluxes,
                convection=problem.convection,
            )
            full = np.zeros(mesh.n_points)
            full[result.node_ids] = result.temperature
            estimate = heat_indicator(
                mesh,
                full,
                problem.material,
                source=problem.source,
                dirichlet_facets=problem.dirichlet,
                outward_fluxes=problem.outward_fluxes,
                convection=problem.convection,
            )
            marked = tuple(
                int(i) for i in mark_bulk(estimate.cell, fraction=control.fraction)
            )
            stages.append(
                AdaptiveStage(
                    mesh,
                    result,
                    estimate,
                    marked,
                    ()
                    if parents is None
                    else tuple(int(i) for i in parents.cell_parents),
                    ()
                    if parents is None
                    else tuple(int(i) for i in parents.facet_parents),
                )
            )
            if estimate.total <= control.tolerance:
                return finish("tolerance", "heuristic indicator tolerance reached")
            if len(stages) >= control.max_steps:
                return finish("step_budget", "maximum solved mesh count reached")
            window = control.stagnation_window
            if len(stages) > window and all(
                stages[i].indicator.total
                >= (1 - control.min_reduction) * stages[i - 1].indicator.total
                for i in range(len(stages) - window, len(stages))
            ):
                return finish(
                    "stagnation", "indicator reduction below configured minimum"
                )
            try:
                parents = remesher(mesh, marked, max_points=control.max_dofs)
            except ValidationError as error:
                if error.field == "budget":
                    return finish("dof_budget", str(error))
                raise
            validate_refinement(mesh, parents)
            if parents.mesh.n_points == mesh.n_points or mesh_hash(
                parents.mesh
            ) == mesh_hash(mesh):
                return finish("stagnation", "remesher made no progress")
            if thermal_dofs(parents.mesh).size > control.max_dofs:
                return finish(
                    "dof_budget", "remesher exceeds active thermal DOF budget"
                )
            problem = problem.refined(parents)
            mesh = parents.mesh
    except (ValidationError, ValueError, RuntimeError) as error:
        return finish("error", f"{type(error).__name__}: {error}")


__all__ = (
    "AdaptiveHeatResult",
    "AdaptiveOptions",
    "AdaptiveStage",
    "HeatProblem",
    "solve_adaptive_heat",
)
