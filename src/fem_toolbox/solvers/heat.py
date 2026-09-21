"""Stationary thermal reference solver with explicit outward-power convention."""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.sparse.csgraph import connected_components  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.assembly import MatrixBlock, VectorBlock, assemble, indices
from fem_toolbox.boundary import Dirichlet, eliminate_dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.heat import conductance_matrix
from fem_toolbox.errors import ThermalNullspaceError, ValidationError
from fem_toolbox.fields import DofMap, Field
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.heat import (
    Convection,
    ScalarField,
    convection_terms,
    flux_load,
    source_load,
)
from fem_toolbox.provenance import RunProvenance


def thermal_dofs(mesh: Mesh) -> DofMap:
    """One T DOF per sorted active Tet4 mesh node, in K; no unused-node DOFs."""
    return DofMap(
        (Field("T", ("scalar",), tuple(int(n) for n in np.unique(mesh.tetra4)), "K"),)
    )


@dataclass(frozen=True)
class HeatResult:
    """Immutable raw thermal results with explicit SI units.

    node_ids (nactive,) map temperature (nactive,) K back to the mesh.
    gradient/heat_flux (ncells,3) are K/m and W/m² in input cell order.
    reactions/residual (nactive,) are W; reaction is zero on free DOFs and
    positive for heat supplied by a prescribed-temperature reservoir.
    source_power, prescribed_outflow, convection_outflow, dirichlet_inflow and
    power_balance are W. Balance = source + Dirichlet - prescribed - convection.
    relative_residual is the free residual norm / || |K| |T| + |f| ||.
    provenance records the effective discrete model, not callable source code.
    """

    node_ids: NDArray[np.int64]
    temperature: NDArray[np.float64]
    gradient: NDArray[np.float64]
    heat_flux: NDArray[np.float64]
    reactions: NDArray[np.float64]
    residual: NDArray[np.float64]
    source_power: float
    prescribed_outflow: float
    convection_outflow: float
    dirichlet_inflow: float
    power_balance: float
    relative_residual: float
    provenance: RunProvenance

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            value = getattr(self, name)
            if isinstance(value, np.ndarray):
                object.__setattr__(
                    self,
                    name,
                    np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(
                        value.shape
                    ),
                )


def _check_nullspace(
    mesh: Mesh,
    active: NDArray[np.int64],
    anchored: NDArray[np.int64],
    rhs: NDArray[np.float64],
) -> None:
    rows, columns = [], []
    for cell in mesh.tetra4:
        for a, b in combinations(np.searchsorted(active, cell), 2):
            rows.extend([a, b])
            columns.extend([b, a])
    graph = sparse.coo_matrix(
        (np.ones(len(rows)), (rows, columns)), shape=(len(active), len(active))
    ).tocsr()
    count, labels = connected_components(graph)
    for component in range(count):
        nodes = np.flatnonzero(labels == component)
        if not np.intersect1d(nodes, anchored).size:
            net = float(rhs[nodes].sum())
            tolerance = 1e-10 * float(np.abs(rhs[nodes]).sum())
            raise ThermalNullspaceError(
                component=component,
                compatible=abs(net) <= tolerance,
                net_power=net,
                tolerance=tolerance,
            )


def solve_heat(
    mesh: Mesh,
    material: ThermalMaterial | Sequence[ThermalMaterial],
    conditions: Iterable[Dirichlet],
    *,
    source: ScalarField = 0.0,
    outward_fluxes: Mapping[int, ScalarField] | None = None,
    convection: Mapping[int, Convection] | None = None,
) -> HeatResult:
    """Solve -div(k grad T)=Q on all Tet4 cells with SI data.

    Material is homogeneous or one per input cell. Dirichlet indices come from
    thermal_dofs(mesh). Surface mappings address exterior triangle row indices.
    Positive outward flux removes heat. Robin law is q_out=h*(T-T_inf); surface
    contributions may add (including prescribed flux on a convecting surface).
    Source/flux/ambient callables use degree-two quadrature. Temperatures share
    the same reference, normally absolute kelvin; no reference is invented.

    ValidationError rejects invalid input, singularity or failed residual/balance
    checks. ThermalNullspaceError diagnoses each first unanchored component as
    compatible or incompatible before sparse factorization. Pure Neumann systems
    are not automatically pinned, even when compatible. No transient coupling.
    """
    if mesh.n_tetra4 == 0:
        raise ValidationError("heat solver requires Tet4 cells", field="mesh")
    materials = (
        [material] * mesh.n_tetra4
        if isinstance(material, ThermalMaterial)
        else list(material)
    )
    if len(materials) != mesh.n_tetra4 or any(
        not isinstance(m, ThermalMaterial) for m in materials
    ):
        raise ValidationError("one ThermalMaterial required per cell", field="material")
    dofs = thermal_dofs(mesh)
    active = np.unique(mesh.tetra4)
    geometries = [
        Tet4Geometry(
            mesh.points[cell],
            entity_id=int(mesh.tetra4_ids[i]) if mesh.tetra4_ids is not None else i,
        )
        for i, cell in enumerate(mesh.tetra4)
    ]
    matrices: list[MatrixBlock] = []
    source_vector = np.zeros(dofs.size)
    flux_vector = np.zeros(dofs.size)
    ambient_vector = np.zeros(dofs.size)
    robin: list[tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64]]] = []
    convective_nodes: list[int] = []
    for cell, g, m in zip(mesh.tetra4, geometries, materials, strict=True):
        ids = dofs.indices("T", cell)
        matrices.append(MatrixBlock(ids, ids, conductance_matrix(g, m)))
        np.add.at(source_vector, ids, source_load(g, source))
    faces = Counter(
        tuple(sorted(int(n) for n in face))
        for cell in mesh.tetra4
        for face in combinations(cell, 3)
    )

    def face_nodes(face: int) -> NDArray[np.int64]:
        face_id = int(indices([face], mesh.n_triangle3)[0])
        nodes = mesh.triangle3[face_id]
        if faces[tuple(sorted(int(n) for n in nodes))] != 1:
            raise ValidationError(
                "thermal boundary must be an exterior Tet4 face", field="surface"
            )
        return np.asarray(nodes, dtype=np.int64)

    for face, flux in (outward_fluxes or {}).items():
        nodes = face_nodes(face)
        np.add.at(
            flux_vector, dofs.indices("T", nodes), flux_load(mesh.points[nodes], flux)
        )
    for face, condition in (convection or {}).items():
        nodes = face_nodes(face)
        if not isinstance(condition, Convection):
            raise ValidationError(
                "expected Convection boundary data", field="convection"
            )
        ids = dofs.indices("T", nodes)
        h, load = convection_terms(mesh.points[nodes], condition)
        matrices.append(MatrixBlock(ids, ids, h))
        np.add.at(ambient_vector, ids, load)
        robin.append((ids, h, load))
        if condition.coefficient > 0:
            convective_nodes.extend(ids.tolist())
    original = assemble(
        dofs.size,
        matrices,
        [
            VectorBlock(
                np.arange(dofs.size), source_vector + flux_vector + ambient_vector
            )
        ],
    )
    bcs = tuple(conditions)
    constrained = eliminate_dirichlet(original, bcs)
    fixed = (
        np.unique(np.concatenate([indices(bc.dofs, dofs.size) for bc in bcs]))
        if bcs
        else np.empty(0, dtype=np.int64)
    )
    anchored = np.union1d(fixed, np.array(convective_nodes, dtype=np.int64))
    _check_nullspace(mesh, active, anchored, original.rhs)
    free = np.setdiff1d(np.arange(dofs.size), fixed)
    temperature = np.zeros(dofs.size)
    temperature[fixed] = constrained.rhs[fixed]
    if len(free):
        reduced = constrained.matrix.tocsr()[free][:, free].tocsc()
        diagonal = reduced.diagonal()
        if np.any(diagonal <= 0):
            raise ValidationError(
                "nonpositive thermal conductance diagonal", field="solver"
            )
        scale = 1 / np.sqrt(diagonal)
        scaling = sparse.diags(scale)
        try:
            factor = splu((scaling @ reduced @ scaling).tocsc())
            if np.min(np.abs(factor.U.diagonal())) <= 1e-12:
                raise ValidationError(
                    "singular or ill-conditioned thermal system", field="solver"
                )
            temperature[free] = scale * factor.solve(scale * constrained.rhs[free])
        except RuntimeError as error:
            raise ValidationError("singular thermal system", field="solver") from error
    if not np.isfinite(temperature).all():
        raise ValidationError("nonfinite temperature solution", field="solver")
    residual = original.residual(temperature)
    roundoff_scale = abs(original.matrix) @ np.abs(temperature) + np.abs(original.rhs)
    norm = float(np.linalg.norm(roundoff_scale))
    relative = float(np.linalg.norm(residual[free])) / norm if norm > 0 else 0.0
    if relative > 1e-9:
        raise ValidationError("free thermal residual exceeds 1e-9", field="solver")
    reaction = np.zeros(dofs.size)
    reaction[fixed] = residual[fixed]
    gradient = np.array(
        [
            g.gradient_matrix
            @ (
                temperature[dofs.indices("T", cell)]
                - temperature[dofs.indices("T", cell)[0]]
            )
            for cell, g in zip(mesh.tetra4, geometries, strict=True)
        ]
    )
    heat_flux = np.array(
        [-m.conductivity @ grad for m, grad in zip(materials, gradient, strict=True)]
    )
    source_power = float(source_vector.sum())
    prescribed_outflow = -float(flux_vector.sum())
    convective_power = [h @ temperature[ids] - load for ids, h, load in robin]
    convection_outflow = sum(float(power.sum()) for power in convective_power)
    dirichlet_inflow = float(reaction.sum())
    balance = source_power + dirichlet_inflow - prescribed_outflow - convection_outflow
    power_scale = float(
        np.abs(source_vector).sum() + np.abs(flux_vector).sum() + np.abs(reaction).sum()
    ) + sum(float(np.abs(power).sum()) for power in convective_power)
    power_tolerance = 1e-9 * power_scale + 100 * np.finfo(float).eps * float(
        roundoff_scale.sum()
    )
    if abs(balance) > power_tolerance:
        raise ValidationError("thermal power balance exceeds tolerance", field="solver")
    provenance = RunProvenance.from_input_bytes(
        json.dumps(
            {
                "points": mesh.points.tolist(),
                "tetra4": mesh.tetra4.tolist(),
                "materials": [m.conductivity.tolist() for m in materials],
                "source_load": source_vector.tolist(),
                "flux_load": flux_vector.tolist(),
                "robin": sorted(
                    [
                        (ids.tolist(), h.tolist(), load.tolist())
                        for ids, h, load in robin
                    ]
                ),
                "fixed_dofs": fixed.tolist(),
                "fixed_values": temperature[fixed].tolist(),
            },
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8"),
        solver_name="fem_toolbox.solvers.heat.solve_heat",
        solver_parameters={
            "input_schema": "heat-discrete-v1",
            "backend": "scipy.sparse.linalg.splu",
            "scaling": "diagonal",
            "pivot_tolerance": 1e-12,
            "neumann_relative_tolerance": 1e-10,
            "relative_residual_tolerance": 1e-9,
            "relative_power_tolerance": 1e-9,
            "power_roundoff_factor": 100,
            "load_quadrature_degree": 2,
            "flux_convention": "positive_outward",
        },
    )
    return HeatResult(
        active,
        temperature,
        gradient,
        heat_flux,
        reaction,
        residual,
        source_power,
        prescribed_outflow,
        convection_outflow,
        dirichlet_inflow,
        balance,
        relative,
        provenance,
    )
