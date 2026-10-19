"""Thresholded Tet4 geometry and independent full-material elastic reanalysis."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.sparse.csgraph import dijkstra  # type: ignore[import-untyped]

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.io import write_vtu
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.solvers.elasticity import (
    ElasticResult,
    mechanical_dofs,
    solve_elasticity,
)
from fem_toolbox.topology import DensityEvaluation, TopologyProblem


@dataclass(frozen=True)
class GeometryResult:
    """Exact selected-cell geometry with unchanged dead loads and full material."""

    mesh: Mesh
    selected_cells: NDArray[np.int64]
    threshold: float
    objective: float
    compliance_ratio: float
    volume_fraction: float
    states: tuple[ElasticResult, ...]
    von_mises_Pa: NDArray[np.float64]
    added_cells: tuple[int, ...] = ()
    connection_repair: bool = False

    @property
    def maximum_von_mises_Pa(self) -> float:
        """Largest raw cell stress among all cases, not an allowable stress."""
        return float(self.von_mises_Pa.max())

    def to_record(self) -> dict[str, Any]:
        """Predeclared compliance/volume gate and diagnostic cell indices."""
        return {
            "threshold": self.threshold,
            "connection_repair": self.connection_repair,
            "added_original_cells": list(self.added_cells),
            "status": "accepted"
            if self.compliance_ratio <= 2 and self.volume_fraction <= 0.9
            else "gate_failed",
            "cells": self.mesh.n_cells,
            "boundary_triangles": self.mesh.n_facets,
            "selected_original_cells": self.selected_cells.tolist(),
            "compliance_J": self.objective,
            "density_compliance_ratio": self.compliance_ratio,
            "volume_fraction": self.volume_fraction,
            "maximum_von_mises_Pa": self.maximum_von_mises_Pa,
            "critical_original_cells": self.selected_cells[
                np.argsort(self.von_mises_Pa)[-5:][::-1]
            ].tolist(),
            "maximum_relative_residual": max(s.relative_residual for s in self.states),
            "maximum_force_balance_N": max(
                float(np.linalg.norm(s.force_balance)) for s in self.states
            ),
            "maximum_moment_balance_N_m": max(
                float(np.linalg.norm(s.moment_balance)) for s in self.states
            ),
            "limits": "raw Tet4 stress; no stress convergence or strength/manufacturing qualification",
        }


def _connect_cells(
    problem: TopologyProblem, selected: NDArray[np.int64], rho: NDArray[np.float64]
) -> NDArray[np.int64]:
    """Explicit optional repair: add minimum-cost face paths, never fixed voids."""
    owners: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for i, cell in enumerate(problem.mesh.tetra4):
        for face_nodes in combinations(cell, 3):
            owners[tuple(sorted(int(n) for n in face_nodes))].append(i)
    forbidden = set(int(i) for i in problem.void)
    neighbors: dict[int, list[int]] = defaultdict(list)
    rows = []
    columns = []
    weights = []
    for group in owners.values():
        if len(group) == 2 and not forbidden.intersection(group):
            a, b = group
            for source, target in ((a, b), (b, a)):
                neighbors[source].append(target)
                rows.append(source)
                columns.append(target)
                weights.append(
                    float(problem.volumes[target] * (1 - rho[target] + 1e-3))
                )
    graph = sparse.csr_matrix(
        (weights, (rows, columns)), shape=(problem.mesh.n_cells, problem.mesh.n_cells)
    )
    retained = set(int(i) for i in selected)
    while True:
        primary = {min(retained)}
        queue = list(primary)
        while queue:
            current = queue.pop()
            for adjacent in neighbors[current]:
                if adjacent in retained and adjacent not in primary:
                    primary.add(adjacent)
                    queue.append(adjacent)
        targets = sorted(retained - primary)
        if not targets:
            return np.asarray(sorted(retained), dtype=np.int64)
        distances, predecessors, _ = dijkstra(
            graph, indices=sorted(primary), min_only=True, return_predecessors=True
        )
        target = min(targets, key=lambda i: (distances[i], i))
        if not np.isfinite(distances[target]):
            raise ValidationError(
                "fixed void prevents connecting geometry", field="geometry"
            )
        current = target
        while current not in primary:
            retained.add(int(current))
            current = int(predecessors[current])
            if current < 0:
                raise ValidationError(
                    "invalid geometry connection path", field="geometry"
                )


def reanalyze_geometry(
    problem: TopologyProblem,
    field: DensityEvaluation,
    threshold: float,
    *,
    connect: bool = False,
) -> GeometryResult:
    """Keep rho>=threshold, reject disconnected solids or removed loaded faces.

    Newly exposed surfaces are traction-free. Prescribed support components are
    remapped onto retained active nodes; the original solver checks rigid modes.
    Body forces on removed volume are rejected rather than silently changed.
    """
    if not np.isfinite(threshold) or not 0 < threshold < 1:
        raise ValidationError("geometry threshold must lie in (0,1)", field="threshold")
    rho = problem.physical_density(field.density)
    if not np.allclose(rho, field.physical_density, rtol=0, atol=1e-13):
        raise ValidationError("density field does not match problem", field="geometry")
    selected = np.flatnonzero(rho >= threshold)
    if len(selected) == 0:
        raise ValidationError("empty threshold geometry", field="geometry")
    initial_selected = selected.copy()
    if connect:
        selected = _connect_cells(problem, selected, rho)
    cells = problem.mesh.tetra4[selected]
    owners: dict[tuple[int, ...], list[int]] = defaultdict(list)
    parent = list(range(len(cells)))

    def find(i: int) -> int:
        while parent[i] != i:
            i = parent[i]
        return i

    for i, cell in enumerate(cells):
        for face_nodes in combinations(cell, 3):
            face = tuple(sorted(int(n) for n in face_nodes))
            owners[face].append(i)
    for group in owners.values():
        if len(group) > 2:
            raise ValidationError("nonmanifold threshold geometry", field="geometry")
        if len(group) == 2:
            a, b = map(find, group)
            parent[b] = a
    if len({find(i) for i in range(len(cells))}) != 1:
        raise ValidationError("disconnected threshold geometry", field="geometry")
    faces = sorted(face for face, group in owners.items() if len(group) == 1)
    mesh = Mesh(
        points=problem.mesh.points,
        tetra4=cells,
        triangle3=np.asarray(faces, dtype=np.int64),
        tetra4_ids=selected,
        volume_sets={"retained": np.arange(len(cells))},
    )
    inspect_boundary_topology(mesh, require_complete=True)
    lookup = {face: i for i, face in enumerate(faces)}
    new_nodes = np.unique(cells)
    original_nodes = np.unique(problem.mesh.tetra4)
    dofs = mechanical_dofs(mesh)
    states = []
    for case in problem.cases:
        if len(selected) < problem.mesh.n_cells and np.any(
            np.asarray(case.body_force) != 0
        ):
            raise ValidationError(
                "cannot preserve body load on removed volume", field="geometry"
            )
        bcs = []
        for bc in case.conditions:
            old = np.asarray(bc.dofs, dtype=np.int64)
            nodes = original_nodes[old // 3]
            keep = np.isin(nodes, new_nodes)
            remapped = []
            for node, component in zip(nodes[keep], (old % 3)[keep], strict=True):
                remapped.append(int(dofs.indices("u", [node])[int(component)]))
            if remapped:
                bcs.append(Dirichlet(remapped, np.zeros(len(remapped))))
        tractions = {}
        for old_face, vector in (case.tractions or {}).items():
            face = tuple(sorted(int(n) for n in problem.mesh.triangle3[old_face]))
            if face not in lookup:
                raise ValidationError(
                    "loaded face removed by threshold", field="geometry"
                )
            tractions[lookup[face]] = vector
        point_loads = {}
        for node, vector in (case.point_loads or {}).items():
            if node not in new_nodes and np.any(np.asarray(vector) != 0):
                raise ValidationError(
                    "loaded node removed by threshold", field="geometry"
                )
            if node in new_nodes:
                point_loads[node] = vector
        states.append(
            solve_elasticity(
                mesh,
                [problem.materials[int(i)] for i in selected],
                tuple(bcs),
                body_force=case.body_force,
                tractions=tractions,
                point_loads=point_loads,
            )
        )
    objective = float(
        sum(
            w * np.sum(s.applied_loads * s.displacement)
            for w, s in zip(problem.weights, states, strict=True)
        )
    )
    for state in states:
        load_scale = max(
            float(np.linalg.norm(state.applied_loads, axis=1).sum()),
            np.finfo(float).tiny,
        )
        if np.linalg.norm(state.force_balance) > 1e-8 * load_scale or abs(
            2 * state.energy - float(np.sum(state.applied_loads * state.displacement))
        ) > 1e-8 * max(2 * state.energy, np.finfo(float).tiny):
            raise ValidationError(
                "geometry work/energy or force balance failed", field="geometry"
            )
    mises = []
    for state in states:
        s = state.stress
        mises.append(
            np.sqrt(
                0.5
                * (
                    (s[:, 0] - s[:, 1]) ** 2
                    + (s[:, 1] - s[:, 2]) ** 2
                    + (s[:, 2] - s[:, 0]) ** 2
                )
                + 3 * np.sum(s[:, 3:] ** 2, axis=1)
            )
        )
    return GeometryResult(
        mesh,
        selected,
        float(threshold),
        objective,
        objective / field.objective,
        float(problem.volumes[selected].sum() / problem.volumes.sum()),
        tuple(states),
        np.max(mises, axis=0),
        tuple(int(i) for i in np.setdiff1d(selected, initial_selected)),
        connect,
    )


def export_reanalyses(
    problem: TopologyProblem,
    field: DensityEvaluation,
    output: str | Path,
    thresholds: tuple[float, ...] = (0.25, 0.35, 0.45),
) -> list[dict[str, Any]]:
    """Export every valid geometry; retain explicit No-go records for others."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    records = []
    files = {}
    for i, (threshold, connect) in enumerate(
        (t, c) for c in (False, True) for t in thresholds
    ):
        try:
            geometry = reanalyze_geometry(problem, field, threshold, connect=connect)
            record = geometry.to_record()
            name = f"geometry-{i}.vtu"
            write_vtu(
                output / name,
                geometry.mesh,
                cell_data={
                    "von_mises_Pa": {
                        "tetra4": geometry.von_mises_Pa,
                        "triangle3": np.zeros(geometry.mesh.n_facets),
                    }
                },
            )
            record["artifact"] = name
            files[name] = hashlib.sha256((output / name).read_bytes()).hexdigest()
        except ValidationError as error:
            record = {
                "threshold": threshold,
                "connection_repair": connect,
                "status": "invalid_geometry",
                "reason": str(error),
            }
        records.append(record)
    (output / "reanalysis.json").write_text(
        json.dumps(
            {
                "schema": 1,
                "thresholds": records,
                "artifact_sha256": files,
                "gate_passed": any(r["status"] == "accepted" for r in records),
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    return records
