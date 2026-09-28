"""Small-model Tet4/Tet10 piezo modal reference with algebraic electric DOFs.

The dense Schur reference is intentionally limited to 2048 active [u, phi]
DOFs. No electric mass, artificial grounding, or dynamic thermal feedback.
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigh  # type: ignore[import-untyped]

from fem_toolbox.assembly import MatrixBlock, assemble, indices, real_array
from fem_toolbox.boundary import Dirichlet, Electrode, eliminate_dirichlet
from fem_toolbox.elements import Tet4Geometry, Tet10Geometry, quadratic
from fem_toolbox.elements.storage import mass_matrix
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.elasticity import _check_supports
from fem_toolbox.solvers.electrostatics import _check_gauge
from fem_toolbox.solvers.linear_backend import LinearBackend, solve_linear
from fem_toolbox.solvers.piezoelectric import (
    _electrode_reduction,
    _materials,
    piezoelectric_dofs,
)
from fem_toolbox.solvers.piezoelectric_elements import piezoelectric_matrix

MAX_DENSE_DOFS = 2048


def _readonly(value: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.frombuffer(value.tobytes(), dtype=np.float64).reshape(value.shape)


@dataclass(frozen=True, eq=False)
class PiezoelectricModes:
    """Mass-normalized real eigenvectors in sorted active-node order.

    Arrays: displacement (nmodes,nactive,3), potential (nmodes,nactive),
    angular_frequencies (nmodes,) rad/s, frequencies_hz (nmodes,) Hz.
    Modal amplitudes satisfy u.T M u=1 (SI mass metric); they are shapes,
    not a physical displacement/voltage excitation. The largest mechanical
    component is positive. Degenerate eigenspaces have no unique basis.
    Electrical residuals are evaluated after tying conducting electrodes;
    conductor charge is minus the sum of original electrical residuals.
    All arrays are owned immutable snapshots. No experimental validation.
    """

    node_ids: NDArray[np.int64]
    displacement: NDArray[np.float64]
    potential: NDArray[np.float64]
    angular_frequencies: NDArray[np.float64]
    frequencies_hz: NDArray[np.float64]
    relative_residual_u: NDArray[np.float64]
    relative_residual_phi: NDArray[np.float64]
    electrode_names: tuple[str, ...]
    electrode_potentials: NDArray[np.float64]
    electrode_charges: NDArray[np.float64]
    provenance: RunProvenance

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "node_ids", np.frombuffer(self.node_ids.tobytes(), dtype=np.int64)
        )
        for name in (
            "displacement",
            "potential",
            "angular_frequencies",
            "frequencies_hz",
            "relative_residual_u",
            "relative_residual_phi",
            "electrode_potentials",
            "electrode_charges",
        ):
            object.__setattr__(self, name, _readonly(getattr(self, name)))


@dataclass
class _DynamicSystem:
    active: NDArray[np.int64]
    matrix: NDArray[np.float64]
    mass: NDArray[np.float64]
    reduced: NDArray[np.float64]
    index: NDArray[np.int64]
    fixed: NDArray[np.int64]
    free_u: NDArray[np.int64]
    free_phi: NDArray[np.int64]
    electric_map: NDArray[np.float64]
    stiffness: NDArray[np.float64]
    free_mass: NDArray[np.float64]
    electrodes: list[tuple[str, NDArray[np.int64], float | None, float | None]]
    canonical: bytes


def _build_dynamic_system(
    mesh: Mesh,
    material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial],
    density: ArrayLike,
    conditions: Iterable[Dirichlet],
    electrodes: Sequence[Electrode] | None,
    backend: LinearBackend | None = None,
) -> _DynamicSystem:
    if mesh.n_cells == 0:
        raise ValidationError(
            "piezo dynamics requires Tet4 or Tet10 cells", field="mesh"
        )
    dofs = piezoelectric_dofs(mesh)
    if dofs.size > MAX_DENSE_DOFS:
        raise ValidationError(
            f"dense reference supports at most {MAX_DENSE_DOFS} active DOFs",
            field="solver",
        )
    if mesh.n_tetra10:
        inspect_boundary_topology(mesh)
    active = np.unique(mesh.cells)
    nu = 3 * len(active)
    materials = _materials(material, mesh.n_cells)
    rho = real_array(density, np.asarray(density).shape)
    if rho.ndim == 0:
        rho = np.full(mesh.n_cells, float(rho))
    rho = real_array(rho, (mesh.n_cells,))
    if np.any(rho <= 0):
        raise ValidationError("density must be positive in kg/m³", field="density")
    bcs = tuple(conditions)
    for bc in bcs:
        ids = indices(bc.dofs, dofs.size)
        if np.any(real_array(bc.values, (len(ids),)) != 0):
            raise ValidationError(
                "dynamic topology requires homogeneous conditions", field="conditions"
            )
    for electrode in electrodes or ():
        if not isinstance(electrode, Electrode):
            raise ValidationError("expected Electrode", field="electrodes")
        for value in (electrode.potential, electrode.charge):
            if value is not None and float(real_array(value, ())) != 0:
                raise ValidationError(
                    "dynamic topology requires homogeneous electrodes",
                    field="electrodes",
                )
    blocks: list[MatrixBlock] = []
    masses: list[MatrixBlock] = []
    for cell, mat, r in zip(mesh.cells, materials, rho, strict=True):
        geometry = (Tet10Geometry if mesh.n_tetra10 else Tet4Geometry)(
            mesh.points[cell]
        )
        u_ids = dofs.indices("u", cell)
        ids = np.concatenate((u_ids, dofs.indices("phi", cell)))
        blocks.append(MatrixBlock(ids, ids, piezoelectric_matrix(geometry, mat)))
        local_mass = (
            quadratic.mass_matrix(geometry, float(r))
            if isinstance(geometry, Tet10Geometry)
            else mass_matrix(geometry, float(r))
        )
        masses.append(MatrixBlock(u_ids, u_ids, local_mass))
    original = assemble(dofs.size, blocks)
    reduction = _electrode_reduction(mesh, dofs, original, bcs, electrodes)
    # Also validates duplicate/conflicting Dirichlet declarations.
    eliminate_dirichlet(reduction.system, reduction.conditions)
    fixed = (
        np.unique(
            np.concatenate(
                [indices(bc.dofs, dofs.size) for bc in reduction.original_conditions]
            )
        )
        if reduction.original_conditions
        else np.empty(0, dtype=np.int64)
    )
    _check_supports(mesh, dofs, fixed)
    _check_gauge(mesh, active, fixed[fixed >= nu] - nu, original.rhs[nu:])
    reduced_fixed = np.unique(reduction.index[fixed])
    free = np.setdiff1d(np.arange(reduction.system.size), reduced_fixed)
    free_u, free_phi = free[free < nu], free[free >= nu]
    if not len(free_u):
        raise ValidationError("no free mechanical DOFs", field="conditions")
    matrix = np.asarray(original.matrix.toarray(), dtype=np.float64)
    mass = np.asarray(assemble(nu, masses).matrix.toarray(), dtype=np.float64)
    reduced = np.asarray(reduction.system.matrix.toarray(), dtype=np.float64)
    coupling = reduced[np.ix_(free_phi, free_u)]
    if len(free_phi):
        dielectric = -reduced[np.ix_(free_phi, free_phi)]
        if mesh.n_tetra10:
            # Exact bilinear form is symmetric. Quadrature cancellation near
            # zero can leave unequal roundoff in transposed entries. Physical
            # residuals below still use the original assembled operator.
            dielectric = (dielectric + dielectric.T) / 2
        electric_map = np.asarray(
            solve_linear(
                dielectric,
                coupling,
                fields={"phi": np.arange(len(free_phi))},
                scale=1 / np.sqrt(np.diag(dielectric)),
                kind="positive_definite",
                backend=backend,
            ).values,
            dtype=np.float64,
        )
    else:
        electric_map = np.empty((0, len(free_u)))
    stiffness = reduced[np.ix_(free_u, free_u)] + coupling.T @ electric_map
    canonical = json.dumps(
        {
            "points": mesh.points.tolist(),
            mesh.cell_type: mesh.cells.tolist(),
            "materials": [mat.to_provenance_dict() for mat in materials],
            "density": rho.tolist(),
            "fixed": fixed.tolist(),
            "index": reduction.index.tolist(),
            "electrodes": [
                (name, ids.tolist(), v, q) for name, ids, v, q in reduction.data
            ],
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return _DynamicSystem(
        active,
        matrix,
        mass,
        reduced,
        reduction.index,
        reduced_fixed,
        free_u,
        free_phi,
        electric_map,
        stiffness,
        mass[np.ix_(free_u, free_u)],
        reduction.data,
        canonical,
    )


def solve_piezoelectric_modes(
    mesh: Mesh,
    material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial],
    density: ArrayLike,
    conditions: Iterable[Dirichlet],
    *,
    backend: LinearBackend | None = None,
    count: int = 6,
    electrodes: Sequence[Electrode] | None = None,
) -> PiezoelectricModes:
    """Solve the lowest ``count`` supported, undamped natural modes.

    ``backend`` is an experimental linear-backend hook; None selects SciPy.
    Reduction, physical residuals and result units remain this solver's contract.

    Constant density: positive scalar or (ncells,) kg/m³. Dirichlet indices
    use piezoelectric_dofs(mesh); all prescribed values/charges must be zero
    (homogeneous perturbations). Electrode(potential=0) is shorted to the
    reference, Electrode(charge=0) is a floating equipotential open circuit.
    Every connected component needs an explicit potential reference and
    supports eliminating rigid-body modes. Uncovered electric boundaries
    have D.n=0; these are not automatically conducting electrodes.

    Only mechanical DOFs carry consistent Tet4/Tet10 mass. Eliminate free electric
    DOFs via K_eff=Kuu+B D^-1 B.T, then solve K_eff u=omega² M u using a
    symmetric definite eigensolver. Reconstruct phi=D^-1 B.T u. The dense
    reference rejects more than 2048 active total DOFs before allocation.
    Invalid inputs, nonpositive eigenvalues or failed field residuals raise
    ValidationError; a missing gauge raises ElectricGaugeError.
    """
    system = _build_dynamic_system(
        mesh, material, density, conditions, electrodes, backend
    )
    if (
        isinstance(count, bool)
        or not isinstance(count, int)
        or not 1 <= count <= len(system.free_u)
    ):
        raise ValidationError(
            "count must be between 1 and the free mechanical DOF count", field="count"
        )
    kscale = float(np.max(np.diag(system.stiffness)))
    mscale = float(np.max(np.diag(system.free_mass)))
    try:
        eigenvalues, vectors = eigh(
            system.stiffness / kscale,
            system.free_mass / mscale,
            subset_by_index=[0, count - 1],
        )
    except np.linalg.LinAlgError as error:
        raise ValidationError(
            "invalid mechanical eigenproblem", field="solver"
        ) from error
    if not np.isfinite(eigenvalues).all() or np.any(eigenvalues <= 0):
        raise ValidationError("nonpositive or nonfinite eigenvalues", field="solver")
    omega = np.sqrt(eigenvalues * kscale / mscale)
    if not np.isfinite(omega).all():
        raise ValidationError("unrepresentable modal spectrum", field="solver")
    vectors /= np.sqrt(mscale)
    for j in range(count):
        if vectors[np.argmax(np.abs(vectors[:, j])), j] < 0:
            vectors[:, j] *= -1
    reduced_modes = np.zeros((system.reduced.shape[0], count))
    reduced_modes[system.free_u] = vectors
    reduced_modes[system.free_phi] = system.electric_map @ vectors
    modes = reduced_modes[system.index]
    nu = 3 * len(system.active)
    inertia = system.mass @ modes[:nu] * omega**2
    residual = system.reduced @ reduced_modes
    residual[:nu] -= inertia
    denominator = np.abs(system.reduced) @ np.abs(reduced_modes)
    denominator[:nu] += np.abs(system.mass) @ np.abs(modes[:nu]) * omega**2
    relative = []
    for free in (system.free_u, system.free_phi):
        num = np.linalg.norm(residual[free], axis=0)
        den = np.linalg.norm(denominator[free], axis=0)
        relative.append(
            np.divide(num, den, out=np.where(num == 0, 0.0, np.inf), where=den != 0)
        )
    if not np.isfinite(modes).all() or any(
        not np.isfinite(r).all() or np.any(r > 1e-9) for r in relative
    ):
        raise ValidationError("modal field residual exceeds 1e-9", field="solver")
    electric_residual = system.matrix[nu:] @ modes
    electrode_potentials = np.empty((count, len(system.electrodes)))
    electrode_charges = np.empty_like(electrode_potentials)
    for i, (_, ids, _, _) in enumerate(system.electrodes):
        electrode_potentials[:, i] = modes[ids[0]]
        electrode_charges[:, i] = -electric_residual[ids - nu].sum(axis=0)
    return PiezoelectricModes(
        system.active,
        modes[:nu].T.reshape(count, -1, 3),
        modes[nu:].T,
        omega,
        omega / (2 * np.pi),
        relative[0],
        relative[1],
        tuple(e[0] for e in system.electrodes),
        electrode_potentials,
        electrode_charges,
        RunProvenance.from_input_bytes(
            system.canonical,
            solver_name="fem_toolbox.solvers.piezo_dynamics.solve_piezoelectric_modes",
            solver_parameters={
                "backend": "scipy.linalg.eigh",
                "linear_backend": "scipy" if backend is None else backend.name,
                "count": count,
                "electrical_mass": 0,
                "mass": "consistent",
                "normalization": "u.T M u = 1",
                "dense_dof_limit": MAX_DENSE_DOFS,
                "field_residual_tolerance": 1e-9,
            },
        ),
    )
