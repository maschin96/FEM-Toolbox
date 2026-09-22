"""Static small-strain elasticity with a SciPy sparse direct reference solve."""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.sparse.csgraph import connected_components  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.assembly import MatrixBlock, VectorBlock, assemble, indices, real_array
from fem_toolbox.boundary import Dirichlet, eliminate_dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.elasticity import stiffness_matrix, strain_matrix
from fem_toolbox.errors import ValidationError
from fem_toolbox.fields import DofMap, Field
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.elasticity import VectorField, body_load, traction_load
from fem_toolbox.provenance import RunProvenance


def mechanical_dofs(mesh: Mesh) -> DofMap:
    """Allocate u=(x,y,z) on sorted active Tet4 node indices, in metres.

    Mesh points unused by Tet4 cells receive no displacement DOFs.
    """
    return DofMap(
        (
            Field(
                "u", ("x", "y", "z"), tuple(int(n) for n in np.unique(mesh.tetra4)), "m"
            ),
        )
    )


@dataclass(frozen=True)
class ElasticResult:
    """Raw, immutable results in active-node / input-cell order.

    node_ids: (nactive,) mesh indices; displacement: (nactive,3) m;
    strain/stress: (ncells,6), total engineering Voigt strain and Pa stress;
    element_energy: (ncells,) mechanical-strain energy J;
    reactions/applied_loads/residual: (nactive,3) N.
    Reactions are zero at free DOFs; residual retains their equilibrium error.
    force_balance (3,) N and moment_balance (3,) N m include external loads and
    reactions. relative_residual is the free residual norm divided by
    || |K| |u| + |f| || (zero for an exactly unloaded state). This backward-error
    scale remains meaningful for rigid motions with cancelling internal forces.
    """

    node_ids: NDArray[np.int64]
    displacement: NDArray[np.float64]
    strain: NDArray[np.float64]
    stress: NDArray[np.float64]
    element_energy: NDArray[np.float64]
    reactions: NDArray[np.float64]
    applied_loads: NDArray[np.float64]
    residual: NDArray[np.float64]
    force_balance: NDArray[np.float64]
    moment_balance: NDArray[np.float64]
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

    @property
    def energy(self) -> float:
        """Total strain energy in joules."""
        return float(self.element_energy.sum())


def _check_supports(mesh: Mesh, dofs: DofMap, fixed: NDArray[np.int64]) -> None:
    """Check six rigid modes on every node-connected component before solving."""
    active = np.unique(mesh.tetra4)
    position = {int(node): i for i, node in enumerate(active)}
    rows, columns = [], []
    for cell in mesh.tetra4:
        for a, b in combinations(cell, 2):
            rows.extend([position[int(a)], position[int(b)]])
            columns.extend([position[int(b)], position[int(a)]])
    graph = sparse.coo_matrix(
        (np.ones(len(rows)), (rows, columns)), shape=(len(active), len(active))
    ).tocsr()
    count, labels = connected_components(graph)
    for component in range(count):
        nodes = active[labels == component]
        xyz = mesh.points[nodes]
        xyz = xyz - xyz.mean(axis=0)
        xyz /= np.max(np.linalg.norm(xyz, axis=1))
        rigid = np.zeros((len(nodes) * 3, 6))
        for i, (x, y, z) in enumerate(xyz):
            rigid[3 * i : 3 * i + 3, :3] = np.eye(3)
            rigid[3 * i : 3 * i + 3, 3:] = [[0, z, -y], [-z, 0, x], [y, -x, 0]]
        constrained = np.isin(dofs.indices("u", nodes), fixed)
        if (
            np.count_nonzero(constrained) < 6
            or np.linalg.matrix_rank(rigid[constrained], tol=1e-10) < 6
        ):
            raise ValidationError(
                f"unrestrained rigid modes in component {component}", field="supports"
            )


def solve_elasticity(
    mesh: Mesh,
    material: ElasticMaterial | Sequence[ElasticMaterial],
    conditions: Iterable[Dirichlet],
    *,
    body_force: VectorField = (0.0, 0.0, 0.0),
    tractions: Mapping[int, VectorField] | None = None,
    point_loads: Mapping[int, ArrayLike] | None = None,
    expansion: ThermalExpansion | Sequence[ThermalExpansion] | None = None,
    temperature: ArrayLike | None = None,
) -> ElasticResult:
    """Solve static linear elasticity on all Tet4 cells, in SI units.

    material is homogeneous or one material per input cell. conditions use
    mechanical_dofs(mesh); nonzero displacement values are allowed. tractions
    maps boundary triangle row indices to global traction vectors/callables.
    point_loads maps active mesh node indices to forces (3,) N. Contributions
    add, including loads on prescribed DOFs. No pressure or follower loads.

    Invalid materials, loads, conflicting conditions, missing supports and
    numerically singular systems raise ValidationError. Optional temperature is
    nodal absolute K; expansion may be homogeneous or per input cell. The thermal
    eigenstrain is alpha times the cell-average temperature difference.
    Diagonally equilibrated LU rejects relative pivots <= 1e-12; very poorly
    conditioned models require reformulation rather than a silent solution.
    """
    if mesh.n_tetra4 == 0:
        raise ValidationError("elasticity requires Tet4 cells", field="mesh")
    materials = (
        [material] * mesh.n_tetra4
        if isinstance(material, ElasticMaterial)
        else list(material)
    )
    if len(materials) != mesh.n_tetra4 or any(
        not isinstance(m, ElasticMaterial) for m in materials
    ):
        raise ValidationError("one ElasticMaterial required per cell", field="material")
    if (expansion is None) != (temperature is None):
        raise ValidationError(
            "expansion and temperature must be supplied together", field="temperature"
        )
    expansions: list[ThermalExpansion] | None = None
    temperatures: NDArray[np.float64] | None = None
    if expansion is not None and temperature is not None:
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
    dofs = mechanical_dofs(mesh)
    active = np.unique(mesh.tetra4)
    geometries = [
        Tet4Geometry(
            mesh.points[cell],
            entity_id=int(mesh.tetra4_ids[i]) if mesh.tetra4_ids is not None else i,
        )
        for i, cell in enumerate(mesh.tetra4)
    ]
    matrices, loads = [], []
    eigenstrains = np.zeros((mesh.n_tetra4, 6))
    for i, (cell, geometry, mat) in enumerate(
        zip(mesh.tetra4, geometries, materials, strict=True)
    ):
        ids = dofs.indices("u", cell)
        matrices.append(MatrixBlock(ids, ids, stiffness_matrix(geometry, mat)))
        loads.append(VectorBlock(ids, body_load(geometry, body_force)))
        if expansions is not None and temperatures is not None:
            eigenstrains[i] = expansions[i].alpha * (
                float(np.mean(temperatures[cell])) - expansions[i].reference_temperature
            )
            loads.append(
                VectorBlock(
                    ids,
                    geometry.volume
                    * strain_matrix(geometry).T
                    @ mat.stiffness
                    @ eigenstrains[i],
                )
            )
    if tractions:
        faces = Counter(
            tuple(sorted(int(n) for n in face))
            for cell in mesh.tetra4
            for face in combinations(cell, 3)
        )
        for face, traction in tractions.items():
            face_id = int(indices([face], mesh.n_triangle3)[0])
            nodes = mesh.triangle3[face_id]
            if faces[tuple(sorted(int(n) for n in nodes))] != 1:
                raise ValidationError(
                    "traction face must be an exterior Tet4 face", field="traction"
                )
            loads.append(
                VectorBlock(
                    dofs.indices("u", nodes),
                    traction_load(mesh.points[nodes], traction),
                )
            )
    if point_loads:
        for node, force in point_loads.items():
            loads.append(
                VectorBlock(dofs.indices("u", [node]), real_array(force, (3,)))
            )
    original = assemble(dofs.size, matrices, loads)
    bcs = tuple(conditions)
    constrained = eliminate_dirichlet(original, bcs)
    fixed = (
        np.unique(np.concatenate([indices(bc.dofs, dofs.size) for bc in bcs]))
        if bcs
        else np.empty(0, dtype=np.int64)
    )
    _check_supports(mesh, dofs, fixed)
    free = np.setdiff1d(np.arange(dofs.size), fixed)
    solution = np.zeros(dofs.size)
    solution[fixed] = constrained.rhs[fixed]
    if len(free):
        reduced = constrained.matrix.tocsr()[free][:, free].tocsc()
        diagonal = reduced.diagonal()
        if np.any(diagonal <= 0):
            raise ValidationError(
                "nonpositive elastic stiffness diagonal", field="solver"
            )
        scale = 1 / np.sqrt(diagonal)
        scaling = sparse.diags(scale)
        try:
            factor = splu((scaling @ reduced @ scaling).tocsc())
            if np.min(np.abs(factor.U.diagonal())) <= 1e-12:
                raise ValidationError(
                    "singular or ill-conditioned elastic system", field="solver"
                )
            solution[free] = scale * factor.solve(scale * constrained.rhs[free])
        except RuntimeError as error:
            raise ValidationError("singular elastic system", field="solver") from error
    if not np.isfinite(solution).all():
        raise ValidationError("nonfinite elastic solution", field="solver")
    residual = original.residual(solution)
    applied = original.rhs
    norm = float(
        np.linalg.norm(abs(original.matrix) @ np.abs(solution) + np.abs(applied))
    )
    relative = float(np.linalg.norm(residual[free])) / norm if norm > 0 else 0.0
    if relative > 1e-9:
        raise ValidationError("free mechanical residual exceeds 1e-9", field="solver")
    reaction = np.zeros(dofs.size)
    reaction[fixed] = residual[fixed]
    strain = np.array(
        [
            strain_matrix(g) @ solution[dofs.indices("u", cell)]
            for cell, g in zip(mesh.tetra4, geometries, strict=True)
        ]
    )
    mechanical_strain = strain - eigenstrains
    stress = np.array(
        [m.stiffness @ e for m, e in zip(materials, mechanical_strain, strict=True)]
    )
    energy = np.array(
        [
            0.5 * g.volume * e @ s
            for g, e, s in zip(geometries, mechanical_strain, stress, strict=True)
        ]
    )
    external = (applied + reaction).reshape(-1, 3)
    provenance = RunProvenance.from_input_bytes(
        json.dumps(
            {
                "points": mesh.points.tolist(),
                "tetra4": mesh.tetra4.tolist(),
                "materials": [m.stiffness.tolist() for m in materials],
                "thermal_eigenstrain": eigenstrains.tolist(),
                "temperature": temperatures.tolist()
                if temperatures is not None
                else None,
                "expansion": [
                    (e.alpha.tolist(), e.reference_temperature) for e in expansions
                ]
                if expansions is not None
                else None,
                "loads": applied.tolist(),
                "fixed_dofs": fixed.tolist(),
                "fixed_values": solution[fixed].tolist(),
            },
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8"),
        solver_name="fem_toolbox.solvers.elasticity.solve_elasticity",
        solver_parameters={
            "input_schema": "elasticity-discrete-v1",
            "backend": "scipy.sparse.linalg.splu",
            "scaling": "diagonal",
            "pivot_tolerance": 1e-12,
            "rigid_rank_tolerance": 1e-10,
            "relative_residual_tolerance": 1e-9,
            "load_quadrature_degree": 2,
        },
    )
    return ElasticResult(
        active,
        solution.reshape(-1, 3),
        strain,
        stress,
        energy,
        reaction.reshape(-1, 3),
        applied.reshape(-1, 3),
        residual.reshape(-1, 3),
        external.sum(axis=0),
        np.cross(mesh.points[active], external).sum(axis=0),
        relative,
        provenance,
    )
