"""Monolithic quasistatic piezoelectric Tet4 reference solver."""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.assembly import (
    LinearSystem,
    MatrixBlock,
    VectorBlock,
    assemble,
    indices,
    real_array,
)
from fem_toolbox.boundary import (
    Dirichlet,
    Electrode,
    ElectrodeState,
    eliminate_dirichlet,
)
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.piezoelectric import (
    piezoelectric_matrix,
    piezoelectric_response,
    piezoelectric_thermal_load,
)
from fem_toolbox.errors import ValidationError
from fem_toolbox.fields import DofMap, Field
from fem_toolbox.materials import PiezoelectricMaterial, ThermalExpansion
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.physics.elasticity import VectorField, body_load, traction_load
from fem_toolbox.physics.electrostatics import (
    ChargeDensity,
    surface_charge_load,
    volume_charge_load,
)
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.elasticity import _check_supports
from fem_toolbox.solvers.electrostatics import _check_gauge


def piezoelectric_dofs(mesh: Mesh) -> DofMap:
    """Return field-major [u (m), phi (V)] DOFs on sorted active Tet4 nodes."""
    active = tuple(int(node) for node in np.unique(mesh.tetra4))
    return DofMap(
        (
            Field("u", ("x", "y", "z"), active, "m"),
            Field("phi", ("scalar",), active, "V"),
        )
    )


@dataclass(frozen=True, eq=False)
class PiezoelectricResult:
    """Raw monolithic solution and physical checks in active-node/cell order.

    Nodal displacement is (nactive,3) m and potential (nactive,) V; cell fields
    follow stress-charge SI units. Reactions and residuals belong to the
    original signed equations, with mechanical rows in N and electric rows in
    C. ``electrode_charges`` is the negative electric residual on fixed and
    floating electrode nodes; ``electrode_states`` sums these by conductor.
    Fieldwise relative residuals use independent backward-error denominators.
    ``mechanical_scale`` and ``electric_scale`` are positive characteristic
    diagonal block values in N/m and C/V; they equal 1 when scaling is off.
    """

    node_ids: NDArray[np.int64]
    displacement: NDArray[np.float64]
    potential: NDArray[np.float64]
    strain: NDArray[np.float64]
    stress: NDArray[np.float64]
    electric_field: NDArray[np.float64]
    electric_displacement: NDArray[np.float64]
    mechanical_energy: NDArray[np.float64]
    coupling_energy: NDArray[np.float64]
    dielectric_enthalpy: NDArray[np.float64]
    residual_u: NDArray[np.float64]
    residual_phi: NDArray[np.float64]
    reactions_u: NDArray[np.float64]
    reactions_phi: NDArray[np.float64]
    electrode_charges: NDArray[np.float64]
    applied_force: NDArray[np.float64]
    applied_electric_load: NDArray[np.float64]
    force_balance: NDArray[np.float64]
    moment_balance: NDArray[np.float64]
    charge_balance: float
    relative_residual_u: float
    relative_residual_phi: float
    mechanical_scale: float
    electric_scale: float
    electrode_states: tuple[ElectrodeState, ...]
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

    @property
    def total_enthalpy(self) -> float:
        """Integrated electric enthalpy in J, not positive stored energy."""
        return float(
            self.mechanical_energy.sum()
            + self.coupling_energy.sum()
            + self.dielectric_enthalpy.sum()
        )


def _materials(
    material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial], count: int
) -> list[PiezoelectricMaterial]:
    if isinstance(material, PiezoelectricMaterial):
        return [material] * count
    if not isinstance(material, Sequence):
        raise ValidationError(
            "one PiezoelectricMaterial required per cell", field="material"
        )
    materials = list(material)
    if len(materials) != count or any(
        not isinstance(item, PiezoelectricMaterial) for item in materials
    ):
        raise ValidationError(
            "one PiezoelectricMaterial required per cell", field="material"
        )
    return materials


def _thermal_data(
    mesh: Mesh,
    expansion: ThermalExpansion | Sequence[ThermalExpansion] | None,
    temperature: ArrayLike | None,
    electric_thermal_strain: bool,
) -> tuple[
    NDArray[np.float64], list[ThermalExpansion] | None, NDArray[np.float64] | None
]:
    if (expansion is None) != (temperature is None):
        raise ValidationError(
            "expansion and temperature must be supplied together", field="temperature"
        )
    if electric_thermal_strain and expansion is None:
        raise ValidationError(
            "electric thermal strain requires expansion and temperature",
            field="thermal",
        )
    strains = np.zeros((mesh.n_tetra4, 6))
    if expansion is None or temperature is None:
        return strains, None, None
    expansions = (
        [expansion] * mesh.n_tetra4
        if isinstance(expansion, ThermalExpansion)
        else list(expansion)
    )
    if len(expansions) != mesh.n_tetra4 or any(
        not isinstance(item, ThermalExpansion) for item in expansions
    ):
        raise ValidationError(
            "one ThermalExpansion required per cell", field="expansion"
        )
    temperatures = real_array(temperature, (mesh.n_points,))
    if np.any(temperatures < 0):
        raise ValidationError(
            "absolute temperature must be nonnegative K", field="temperature"
        )
    for i, (cell, item) in enumerate(zip(mesh.tetra4, expansions, strict=True)):
        strains[i] = item.alpha * (
            float(np.mean(temperatures[cell])) - item.reference_temperature
        )
    return strains, expansions, temperatures


def _block_scales(
    matrix: sparse.spmatrix,
    free: NDArray[np.int64],
    phi_offset: int,
    enabled: bool,
) -> tuple[NDArray[np.float64], float, float]:
    diagonal = np.asarray(matrix.diagonal(), dtype=np.float64)
    mechanical = diagonal[free < phi_offset]
    electric = -diagonal[free >= phi_offset]
    if np.any(mechanical <= 0) or np.any(electric <= 0):
        raise ValidationError(
            "nonpositive piezoelectric block diagonal", field="solver"
        )
    if not enabled:
        return np.ones(len(free)), 1.0, 1.0
    mechanical_scale = float(np.median(mechanical)) if len(mechanical) else 1.0
    electric_scale = float(np.median(electric)) if len(electric) else 1.0
    scale = np.where(
        free < phi_offset,
        1 / np.sqrt(mechanical_scale),
        1 / np.sqrt(electric_scale),
    )
    return scale, mechanical_scale, electric_scale


def _relative_residual(
    residual: NDArray[np.float64],
    matrix: sparse.spmatrix,
    solution: NDArray[np.float64],
    rhs: NDArray[np.float64],
    free: NDArray[np.int64],
) -> float:
    if not len(free):
        return 0.0
    scale = np.asarray(abs(matrix) @ np.abs(solution) + np.abs(rhs), dtype=np.float64)
    denominator = float(np.linalg.norm(scale[free]))
    numerator = float(np.linalg.norm(residual[free]))
    return (
        numerator / denominator if denominator else (0.0 if numerator == 0 else np.inf)
    )


def solve_piezoelectric(
    mesh: Mesh,
    material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial],
    conditions: Iterable[Dirichlet],
    *,
    body_force: VectorField = (0.0, 0.0, 0.0),
    tractions: Mapping[int, VectorField] | None = None,
    point_loads: Mapping[int, ArrayLike] | None = None,
    volume_charge: ChargeDensity = 0.0,
    surface_charges: Mapping[int, ChargeDensity] | None = None,
    electrodes: Sequence[Electrode] | None = None,
    expansion: ThermalExpansion | Sequence[ThermalExpansion] | None = None,
    temperature: ArrayLike | None = None,
    electric_thermal_strain: bool = False,
    block_scaling: bool = True,
    method: str = "direct",
) -> PiezoelectricResult:
    """Solve the symmetric-indefinite static u-phi Tet4 system with sparse LU.

    Dirichlet DOFs use ``piezoelectric_dofs(mesh)`` and retain their original
    signed reactions. Mechanical loads are N/m³, N/m², and nodal N;
    ``volume_charge`` is C/m³ and ``surface_charges`` prescribe outward D.n
    in C/m² on exterior triangle rows. A potential reference and sufficient
    mechanical supports are required on every disconnected component.
    ``expansion`` and nodal absolute temperature must be supplied together.
    Fixed electrodes share a prescribed potential; floating electrodes share
    an unknown potential and prescribe their total adjacent conductor charge.
    The optional thermal-electric contribution to D is explicitly selected.
    The only supported method is direct LU: CG is invalid for this saddle
    system; MINRES awaits a verified block preconditioner. Invalid inputs,
    missing constraints, singularity, or excessive field residuals raise
    ValidationError (ElectricGaugeError for a missing potential reference).
    """
    if mesh.n_tetra4 == 0:
        raise ValidationError("piezoelectric solve requires Tet4 cells", field="mesh")
    if method != "direct":
        raise ValidationError(
            "only direct LU is supported for the piezoelectric saddle system",
            field="method",
        )
    if not isinstance(block_scaling, bool) or not isinstance(
        electric_thermal_strain, bool
    ):
        raise ValidationError("expected boolean solver switches", field="solver")
    materials = _materials(material, mesh.n_tetra4)
    eigenstrains, expansions, _ = _thermal_data(
        mesh, expansion, temperature, electric_thermal_strain
    )
    dofs = piezoelectric_dofs(mesh)
    active = np.unique(mesh.tetra4)
    phi_offset = dofs.block("phi").start
    geometries = [
        Tet4Geometry(
            mesh.points[cell],
            entity_id=int(mesh.tetra4_ids[i]) if mesh.tetra4_ids is not None else i,
        )
        for i, cell in enumerate(mesh.tetra4)
    ]
    matrices: list[MatrixBlock] = []
    loads: list[VectorBlock] = []
    volume_vector = np.zeros(dofs.size)
    surface_vector = np.zeros(dofs.size)
    for i, (cell, geometry, mat) in enumerate(
        zip(mesh.tetra4, geometries, materials, strict=True)
    ):
        u_ids = dofs.indices("u", cell)
        phi_ids = dofs.indices("phi", cell)
        ids = np.concatenate([u_ids, phi_ids])
        matrices.append(MatrixBlock(ids, ids, piezoelectric_matrix(geometry, mat)))
        loads.append(VectorBlock(u_ids, body_load(geometry, body_force)))
        np.add.at(volume_vector, phi_ids, volume_charge_load(geometry, volume_charge))
        if expansions is not None:
            loads.append(
                VectorBlock(
                    ids,
                    piezoelectric_thermal_load(
                        geometry,
                        mat,
                        eigenstrains[i],
                        electric_thermal_strain=electric_thermal_strain,
                    ),
                )
            )
    if tractions or surface_charges:
        faces = Counter(
            tuple(sorted(int(node) for node in face))
            for cell in mesh.tetra4
            for face in combinations(cell, 3)
        )

        def exterior_face(face_id: int, field: str) -> NDArray[np.int64]:
            index = int(indices([face_id], mesh.n_triangle3)[0])
            nodes = mesh.triangle3[index]
            if faces[tuple(sorted(int(node) for node in nodes))] != 1:
                raise ValidationError(
                    "load face must be an exterior Tet4 face", field=field
                )
            return np.asarray(nodes, dtype=np.int64)

        if tractions:
            for face_id, value in tractions.items():
                nodes = exterior_face(face_id, "traction")
                loads.append(
                    VectorBlock(
                        dofs.indices("u", nodes),
                        traction_load(mesh.points[nodes], value),
                    )
                )
        if surface_charges:
            for face_id, value in surface_charges.items():
                nodes = exterior_face(face_id, "surface")
                np.add.at(
                    surface_vector,
                    dofs.indices("phi", nodes),
                    surface_charge_load(mesh.points[nodes], value),
                )
    if point_loads:
        for node, force in point_loads.items():
            loads.append(
                VectorBlock(dofs.indices("u", [node]), real_array(force, (3,)))
            )
    loads.append(VectorBlock(np.arange(dofs.size), volume_vector + surface_vector))
    original = assemble(dofs.size, matrices, loads)
    user_bcs = tuple(conditions)
    bcs = list(user_bcs)
    electrode_data: list[tuple[str, NDArray[np.int64], float | None, float | None]] = []
    occupied: set[int] = set()
    names: set[str] = set()
    representative = np.arange(dofs.size)
    if electrodes:
        inspect_boundary_topology(mesh)
    user_fixed = (
        np.concatenate([indices(bc.dofs, dofs.size) for bc in user_bcs])
        if user_bcs
        else np.empty(0, dtype=np.int64)
    )
    for electrode in electrodes or ():
        if not isinstance(electrode, Electrode):
            raise ValidationError("expected Electrode", field="electrodes")
        if (
            not isinstance(electrode.name, str)
            or not electrode.name
            or electrode.name.strip() != electrode.name
            or electrode.name in names
        ):
            raise ValidationError(
                "electrode names must be unique and nonempty", field="electrodes"
            )
        names.add(electrode.name)
        if (electrode.potential is None) == (electrode.charge is None):
            raise ValidationError(
                "electrode requires either potential or charge", field="electrodes"
            )
        if isinstance(electrode.facets, str):
            if electrode.facets not in mesh.facet_sets:
                raise ValidationError(
                    f"unknown electrode facet set {electrode.facets!r}",
                    field="electrodes",
                )
            facets = mesh.facet_sets[electrode.facets]
        else:
            facets = indices(electrode.facets, mesh.n_triangle3)
        if not len(facets) or len(np.unique(facets)) != len(facets):
            raise ValidationError(
                "electrode facets must be nonempty and unique", field="electrodes"
            )
        if surface_charges and set(int(face) for face in facets).intersection(
            surface_charges
        ):
            raise ValidationError(
                f"electrode {electrode.name!r} conflicts with prescribed surface charge",
                field="electrodes",
            )
        nodes = np.unique(mesh.triangle3[facets])
        if occupied.intersection(int(node) for node in nodes):
            raise ValidationError("electrode node overlap", field="electrodes")
        occupied.update(int(node) for node in nodes)
        phi_dofs = dofs.indices("phi", nodes)
        potential_value = (
            None
            if electrode.potential is None
            else float(real_array(electrode.potential, ()))
        )
        charge_value = (
            None
            if electrode.charge is None
            else float(real_array(electrode.charge, ()))
        )
        if charge_value is not None:
            if np.intersect1d(phi_dofs, user_fixed).size:
                raise ValidationError(
                    "floating electrode conflicts with Dirichlet potential",
                    field="electrodes",
                )
            representative[phi_dofs] = int(phi_dofs[0])
        else:
            bcs.append(Dirichlet(phi_dofs, np.full(len(phi_dofs), potential_value)))
        electrode_data.append((electrode.name, phi_dofs, potential_value, charge_value))

    unique_representatives = np.unique(representative)
    reduced_index = np.searchsorted(unique_representatives, representative)
    has_floating = any(charge is not None for _, _, _, charge in electrode_data)
    if has_floating:
        projection = sparse.coo_matrix(
            (np.ones(dofs.size), (np.arange(dofs.size), reduced_index)),
            shape=(dofs.size, len(unique_representatives)),
        ).tocsr()
        reduced_rhs = np.asarray(projection.T @ original.rhs, dtype=np.float64)
        for _, electrode_dofs, _, charge in electrode_data:
            if charge is not None:
                reduced_rhs[reduced_index[electrode_dofs[0]]] -= charge
        reduced_system = LinearSystem(
            projection.T @ original.matrix @ projection, reduced_rhs
        )
    else:
        reduced_system = original
    reduced_bcs = tuple(
        Dirichlet(reduced_index[indices(bc.dofs, dofs.size)], bc.values) for bc in bcs
    )
    constrained = eliminate_dirichlet(reduced_system, reduced_bcs)
    fixed_reduced = (
        np.unique(
            np.concatenate(
                [indices(bc.dofs, reduced_system.size) for bc in reduced_bcs]
            )
        )
        if reduced_bcs
        else np.empty(0, dtype=np.int64)
    )
    fixed = (
        np.unique(np.concatenate([indices(bc.dofs, dofs.size) for bc in bcs]))
        if bcs
        else np.empty(0, dtype=np.int64)
    )
    _check_supports(mesh, dofs, fixed)
    fixed_phi = fixed[fixed >= phi_offset] - phi_offset
    _check_gauge(mesh, active, fixed_phi, original.rhs[phi_offset:])
    free_reduced = np.setdiff1d(np.arange(reduced_system.size), fixed_reduced)
    reduced_solution = np.zeros(reduced_system.size)
    reduced_solution[fixed_reduced] = constrained.rhs[fixed_reduced]
    mechanical_scale, electric_scale = 1.0, 1.0
    if len(free_reduced):
        reduced = constrained.matrix.tocsr()[free_reduced][:, free_reduced].tocsc()
        scale, mechanical_scale, electric_scale = _block_scales(
            reduced,
            free_reduced,
            int(np.count_nonzero(unique_representatives < phi_offset)),
            block_scaling,
        )
        scaling = sparse.diags(scale)
        try:
            factor = splu((scaling @ reduced @ scaling).tocsc())
            pivots = np.abs(factor.U.diagonal())
            if not np.isfinite(pivots).all() or np.min(pivots) <= 1e-12 * np.max(
                pivots
            ):
                raise ValidationError(
                    "singular or ill-conditioned piezoelectric system", field="solver"
                )
            reduced_solution[free_reduced] = scale * factor.solve(
                scale * constrained.rhs[free_reduced]
            )
        except RuntimeError as error:
            raise ValidationError(
                "singular piezoelectric system", field="solver"
            ) from error
    solution = reduced_solution[reduced_index]
    if not np.isfinite(solution).all():
        raise ValidationError("nonfinite piezoelectric solution", field="solver")
    residual = original.residual(solution)
    rhs = original.rhs
    matrix = original.matrix
    floating_dofs = (
        np.concatenate(
            [
                electrode_dofs
                for _, electrode_dofs, _, charge in electrode_data
                if charge is not None
            ]
        )
        if has_floating
        else np.empty(0, dtype=np.int64)
    )
    free = np.setdiff1d(np.arange(dofs.size), np.union1d(fixed, floating_dofs))
    free_u = free[free < phi_offset]
    free_phi = free[free >= phi_offset]
    relative_u = _relative_residual(residual, matrix, solution, rhs, free_u)
    relative_phi = _relative_residual(residual, matrix, solution, rhs, free_phi)
    if relative_u > 1e-9 or relative_phi > 1e-9:
        raise ValidationError(
            "mechanical or electric free residual exceeds 1e-9", field="solver"
        )
    reaction = np.zeros(dofs.size)
    reaction[fixed] = residual[fixed]
    electrode_states = tuple(
        ElectrodeState(
            name=name,
            potential=float(solution[electrode_dofs[0]]),
            charge=-float(residual[electrode_dofs].sum()),
            prescribed_charge=charge,
        )
        for name, electrode_dofs, _, charge in electrode_data
    )
    for state in electrode_states:
        if state.prescribed_charge is not None:
            tolerance = (
                1e-9 * max(abs(state.prescribed_charge), abs(state.charge)) + 1e-12
            )
            if abs(state.charge - state.prescribed_charge) > tolerance:
                raise ValidationError(
                    f"electrode {state.name!r} total charge exceeds tolerance",
                    field="electrodes",
                )
    displacement = solution[:phi_offset].reshape(-1, 3)
    potential = solution[phi_offset:]
    responses = [
        piezoelectric_response(
            geometry,
            mat,
            displacement[np.searchsorted(active, cell)],
            potential[np.searchsorted(active, cell)],
            thermal_strain=eigenstrains[i] if expansions is not None else None,
            electric_thermal_strain=electric_thermal_strain,
        )
        for i, (cell, geometry, mat) in enumerate(
            zip(mesh.tetra4, geometries, materials, strict=True)
        )
    ]
    applied_force = rhs[:phi_offset].reshape(-1, 3)
    reactions_u = reaction[:phi_offset].reshape(-1, 3)
    reactions_phi = reaction[phi_offset:]
    force_balance = (applied_force + reactions_u).sum(axis=0)
    moment_balance = np.cross(mesh.points[active], applied_force + reactions_u).sum(
        axis=0
    )
    charge_balance = float(
        surface_vector.sum()
        + reactions_phi.sum()
        + volume_vector.sum()
        + residual[floating_dofs].sum()
    )
    charge_scale = float(
        np.abs(surface_vector).sum()
        + np.abs(reactions_phi).sum()
        + np.abs(volume_vector).sum()
        + np.abs(residual[floating_dofs]).sum()
    )
    # In open circuit, net charge and fixed reactions can vanish by physical
    # cancellation. Roundoff still scales with the un-cancelled electric row
    # products, even when the prescribed electric RHS is exactly zero.
    electric_action = np.asarray(
        abs(matrix[phi_offset:]) @ np.abs(solution) + np.abs(rhs[phi_offset:]),
        dtype=np.float64,
    )
    charge_tolerance = 1e-9 * charge_scale + 100 * np.finfo(float).eps * float(
        electric_action.sum()
    )
    if abs(charge_balance) > charge_tolerance:
        raise ValidationError(
            "electric charge balance exceeds tolerance", field="solver"
        )
    provenance_input = {
        "points": mesh.points.tolist(),
        "tetra4": mesh.tetra4.tolist(),
        "materials": [item.to_provenance_dict() for item in materials],
        "thermal_eigenstrain": eigenstrains.tolist(),
        "electric_thermal_strain": electric_thermal_strain,
        "loads": rhs.tolist(),
        "fixed_dofs": fixed.tolist(),
        "fixed_values": solution[fixed].tolist(),
    }
    if electrode_data:
        provenance_input["electrodes"] = [
            (name, electrode_dofs.tolist(), potential, charge)
            for name, electrode_dofs, potential, charge in electrode_data
        ]
    provenance = RunProvenance.from_input_bytes(
        json.dumps(
            provenance_input,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8"),
        solver_name="fem_toolbox.solvers.piezoelectric.solve_piezoelectric",
        solver_parameters={
            "input_schema": (
                "piezoelectric-discrete-v2"
                if electrode_data
                else "piezoelectric-discrete-v1"
            ),
            "backend": "scipy.sparse.linalg.splu",
            "block_scaling": block_scaling,
            "mechanical_scale": mechanical_scale,
            "electric_scale": electric_scale,
            "pivot_relative_tolerance": 1e-12,
            "field_relative_residual_tolerance": 1e-9,
            "equation": "[Kuu Kuphi; Kuphi.T -Kphi] [u; phi] = [force; charge]",
        },
    )
    return PiezoelectricResult(
        node_ids=active,
        displacement=displacement,
        potential=potential,
        strain=np.array([item.strain for item in responses]),
        stress=np.array([item.stress for item in responses]),
        electric_field=np.array([item.electric_field for item in responses]),
        electric_displacement=np.array(
            [item.electric_displacement for item in responses]
        ),
        mechanical_energy=np.array([item.mechanical_energy for item in responses]),
        coupling_energy=np.array([item.coupling_energy for item in responses]),
        dielectric_enthalpy=np.array([item.dielectric_enthalpy for item in responses]),
        residual_u=residual[:phi_offset].reshape(-1, 3),
        residual_phi=residual[phi_offset:],
        reactions_u=reactions_u,
        reactions_phi=reactions_phi,
        electrode_charges=-reactions_phi
        - np.where(
            np.isin(np.arange(phi_offset, dofs.size), floating_dofs),
            residual[phi_offset:],
            0.0,
        ),
        applied_force=applied_force,
        applied_electric_load=rhs[phi_offset:],
        force_balance=force_balance,
        moment_balance=moment_balance,
        charge_balance=charge_balance,
        relative_residual_u=relative_u,
        relative_residual_phi=relative_phi,
        mechanical_scale=mechanical_scale,
        electric_scale=electric_scale,
        electrode_states=electrode_states,
        provenance=provenance,
    )
