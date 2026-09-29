"""Electrostatic Tet4/Tet10 reference solver, preserving the negative electric block."""

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

from fem_toolbox.assembly import MatrixBlock, VectorBlock, assemble, indices
from fem_toolbox.boundary import Dirichlet, eliminate_dirichlet
from fem_toolbox.elements import Tet4Geometry, Tet10Geometry
from fem_toolbox.errors import ElectricGaugeError, ValidationError
from fem_toolbox.fields import DofMap, Field
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.physics.electrostatics import ChargeDensity
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.electric_elements import (
    dielectric_matrix,
    electric_field,
    stored_energy,
    surface_charge_load,
    volume_charge_load,
)
from fem_toolbox.solvers.linear_backend import LinearBackend, solve_linear


def electric_dofs(mesh: Mesh) -> DofMap:
    """Allocate phi (V) on sorted active Tet4/Tet10 nodes, without unused DOFs."""
    return DofMap(
        (Field("phi", ("scalar",), tuple(int(n) for n in np.unique(mesh.cells)), "V"),)
    )


@dataclass(frozen=True)
class ElectrostaticResult:
    """Immutable nodal and raw element results in SI.

    node_ids (nactive,) map potential (nactive,) V to mesh nodes. electric_field
    and electric_displacement (ncells,3) are E=-grad(phi) V/m and D=kappa E C/m².
    Tet4 fields are elementwise constant (field_points=None); Tet10 fields are
    raw centroid samples located by field_points (ncells,3) in metres.
    element_energy (ncells,) integrates positive stored energy (J) over each
    entire cell, not electric enthalpy or energy evaluated at its centroid.
    reactions and residual (nactive,) are C for the ORIGINAL negative electric
    equation. reactions are unknown outward D-flux contributions at fixed DOFs
    (zero on free DOFs). Adjacent conductor charges are their NEGATIVE.
    volume_charge, prescribed_surface_charge and dirichlet_outward_charge are C.
    charge_balance = prescribed surface + Dirichlet outward - volume charge.
    relative_residual = ||r_free|| / || |A| |phi| + |f| ||, A=-K_phi.
    """

    node_ids: NDArray[np.int64]
    potential: NDArray[np.float64]
    electric_field: NDArray[np.float64]
    electric_displacement: NDArray[np.float64]
    element_energy: NDArray[np.float64]
    reactions: NDArray[np.float64]
    electrode_charges: NDArray[np.float64]
    residual: NDArray[np.float64]
    volume_charge: float
    prescribed_surface_charge: float
    dirichlet_outward_charge: float
    charge_balance: float
    relative_residual: float
    provenance: RunProvenance
    field_points: NDArray[np.float64] | None = None

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

    @property
    def energy(self) -> float:
        """Total stored electrostatic field energy in J."""
        return float(self.element_energy.sum())


def _check_gauge(
    mesh: Mesh,
    active: NDArray[np.int64],
    fixed: NDArray[np.int64],
    rhs: NDArray[np.float64],
) -> None:
    rows, columns = [], []
    for cell in mesh.cells:
        for a, b in combinations(np.searchsorted(active, cell), 2):
            rows.extend([a, b])
            columns.extend([b, a])
    graph = sparse.coo_matrix(
        (np.ones(len(rows)), (rows, columns)), shape=(len(active), len(active))
    ).tocsr()
    count, labels = connected_components(graph)
    for component in range(count):
        nodes = np.flatnonzero(labels == component)
        if not np.intersect1d(nodes, fixed).size:
            net = -float(rhs[nodes].sum())
            tolerance = 1e-10 * float(np.abs(rhs[nodes]).sum())
            raise ElectricGaugeError(
                component=component,
                compatible=abs(net) <= tolerance,
                net_charge=net,
                tolerance=tolerance,
            )


def solve_electrostatics(
    mesh: Mesh,
    material: DielectricMaterial | Sequence[DielectricMaterial],
    conditions: Iterable[Dirichlet],
    *,
    backend: LinearBackend | None = None,
    volume_charge: ChargeDensity = 0.0,
    surface_charges: Mapping[int, ChargeDensity] | None = None,
) -> ElectrostaticResult:
    """Solve div(D)=rho, D=-kappa grad(phi), with the negative electric block.

    ``backend`` is an experimental linear-backend hook; None selects SciPy.
    Reduction, physical residuals and result units remain this solver's contract.

    Permittivity is homogeneous or one material per input cell (absolute F/m).
    Dirichlet DOFs use electric_dofs(mesh). volume_charge is rho C/m³; surface
    mappings address exterior triangle rows and prescribe omega=D.n_out C/m²,
    not adjacent conductor charge. Callables receive x (3,) m. Tet4 load rules
    have degree two; Tet10 volume/surface rules have degree five/four. Only
    straight Tet10 is supported. Unspecified boundaries have D.n=0.

    Each connected component needs a supplied potential reference. Missing gauges
    raise ElectricGaugeError before factorization; no automatic pinning occurs.
    Invalid input, conditioning or failed residual/charge balance raises
    ValidationError. No piezoelectric coupling or floating electrodes.
    """
    if mesh.n_cells == 0:
        raise ValidationError(
            "electrostatics requires Tet4 or Tet10 cells", field="mesh"
        )
    materials = (
        [material] * mesh.n_cells
        if isinstance(material, DielectricMaterial)
        else list(material)
    )
    if len(materials) != mesh.n_cells or any(
        not isinstance(m, DielectricMaterial) for m in materials
    ):
        raise ValidationError(
            "one DielectricMaterial required per cell", field="material"
        )
    if mesh.n_tetra10:
        inspect_boundary_topology(mesh)
    dofs = electric_dofs(mesh)
    active = np.unique(mesh.cells)
    geometries = [
        (Tet10Geometry if mesh.n_tetra10 else Tet4Geometry)(
            mesh.points[cell],
            entity_id=int(mesh.cell_ids[i]) if mesh.cell_ids is not None else i,
        )
        for i, cell in enumerate(mesh.cells)
    ]
    matrices = []
    volume_vector = np.zeros(dofs.size)
    surface_vector = np.zeros(dofs.size)
    for cell, g, m in zip(mesh.cells, geometries, materials, strict=True):
        ids = dofs.indices("phi", cell)
        matrices.append(MatrixBlock(ids, ids, -dielectric_matrix(g, m)))
        np.add.at(volume_vector, ids, volume_charge_load(g, volume_charge))
    if surface_charges:
        faces = Counter(
            tuple(sorted(int(n) for n in face))
            for cell in mesh.cells
            for face in combinations(cell[:4], 3)
        )
        for face, charge in surface_charges.items():
            face_id = int(indices([face], mesh.n_facets)[0])
            nodes = mesh.facets[face_id]
            if faces[tuple(sorted(int(n) for n in nodes[:3]))] != 1:
                raise ValidationError(
                    "electrical boundary must be an exterior tetrahedral face",
                    field="surface",
                )
            np.add.at(
                surface_vector,
                dofs.indices("phi", nodes),
                surface_charge_load(mesh.points[nodes], charge),
            )
    original = assemble(
        dofs.size,
        matrices,
        [VectorBlock(np.arange(dofs.size), volume_vector + surface_vector)],
    )
    bcs = tuple(conditions)
    constrained = eliminate_dirichlet(original, bcs)
    fixed = (
        np.unique(np.concatenate([indices(bc.dofs, dofs.size) for bc in bcs]))
        if bcs
        else np.empty(0, dtype=np.int64)
    )
    _check_gauge(mesh, active, fixed, original.rhs)
    free = np.setdiff1d(np.arange(dofs.size), fixed)
    potential = np.zeros(dofs.size)
    potential[fixed] = constrained.rhs[fixed]
    if len(free):
        # Flip BOTH free equation sides solely for positive diagonal scaling.
        reduced = -constrained.matrix.tocsr()[free][:, free].tocsc()
        diagonal = reduced.diagonal()
        if np.any(diagonal <= 0):
            raise ValidationError("nonpositive dielectric diagonal", field="solver")
        scale = 1 / np.sqrt(diagonal)
        potential[free] = solve_linear(
            reduced,
            -constrained.rhs[free],
            fields={"phi": np.arange(len(free))},
            scale=scale,
            pivot_tolerance=1e-12,
            backend=backend,
        ).values
    if not np.isfinite(potential).all():
        raise ValidationError("nonfinite potential", field="solver")
    residual = original.residual(potential)
    roundoff_scale = abs(original.matrix) @ np.abs(potential) + np.abs(original.rhs)
    norm = float(np.linalg.norm(roundoff_scale))
    relative = float(np.linalg.norm(residual[free])) / norm if norm > 0 else 0.0
    if relative > 1e-9:
        raise ValidationError("free electric residual exceeds 1e-9", field="solver")
    reaction = np.zeros(dofs.size)
    reaction[fixed] = residual[fixed]
    field = np.array(
        [
            electric_field(g, potential[dofs.indices("phi", cell)])
            for cell, g in zip(mesh.cells, geometries, strict=True)
        ]
    )
    displacement = np.array(
        [m.permittivity @ e for m, e in zip(materials, field, strict=True)]
    )
    energy = np.array(
        [
            stored_energy(g, m, potential[dofs.indices("phi", cell)])
            for cell, g, m in zip(mesh.cells, geometries, materials, strict=True)
        ]
    )
    volume_total = -float(volume_vector.sum())
    surface_total = float(surface_vector.sum())
    outward = float(reaction.sum())
    balance = surface_total + outward - volume_total
    charge_scale = float(
        np.abs(volume_vector).sum()
        + np.abs(surface_vector).sum()
        + np.abs(reaction).sum()
    )
    tolerance = 1e-9 * charge_scale + 100 * np.finfo(float).eps * float(
        roundoff_scale.sum()
    )
    if abs(balance) > tolerance:
        raise ValidationError(
            "electric charge balance exceeds tolerance", field="solver"
        )
    provenance = RunProvenance.from_input_bytes(
        json.dumps(
            {
                "points": mesh.points.tolist(),
                mesh.cell_type: mesh.cells.tolist(),
                "materials": [m.permittivity.tolist() for m in materials],
                "volume_load": volume_vector.tolist(),
                "surface_load": surface_vector.tolist(),
                "fixed_dofs": fixed.tolist(),
                "fixed_values": potential[fixed].tolist(),
            },
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8"),
        solver_name="fem_toolbox.solvers.electrostatics.solve_electrostatics",
        solver_parameters={
            "input_schema": "electrostatics-discrete-v1",
            "backend": "scipy.sparse.linalg.splu" if backend is None else backend.name,
            "scaling": "diagonal",
            "pivot_tolerance": 1e-12,
            "gauge_relative_tolerance": 1e-10,
            "relative_residual_tolerance": 1e-9,
            "relative_charge_tolerance": 1e-9,
            "charge_roundoff_factor": 100,
            "load_quadrature_degree": 5 if mesh.n_tetra10 else 2,
            "surface_quadrature_degree": 4 if mesh.n_tetra10 else 2,
            "field_location": "cell_centroid" if mesh.n_tetra10 else "element_constant",
            "energy_evaluation": "cell_integral",
            "equation": "-K_phi phi = surface(D.n) - volume(rho)",
        },
    )
    return ElectrostaticResult(
        active,
        potential,
        field,
        displacement,
        energy,
        reaction,
        -reaction,
        residual,
        volume_total,
        surface_total,
        outward,
        balance,
        relative,
        provenance,
        mesh.points[mesh.cells[:, :4]].mean(axis=1) if mesh.n_tetra10 else None,
    )
