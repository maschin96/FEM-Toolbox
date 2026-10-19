"""Filtered SIMP design for stationary Tet4 elasticity; SI and fixed dead loads.

The reference implementation delegates every state to the verified elasticity
solver. It is intentionally limited to small, deterministic local studies.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.spatial import cKDTree  # type: ignore[import-untyped]

from fem_toolbox.assembly import indices, real_array
from fem_toolbox.boundary import Dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ElasticMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.elasticity import ElasticResult, solve_elasticity

if TYPE_CHECKING:
    from fem_toolbox.model import Model


def _immutable(value: ArrayLike) -> NDArray[np.float64]:
    array = np.asarray(value, dtype=float)
    return np.frombuffer(array.tobytes(), dtype=float).reshape(array.shape)


def _vectors(value: Mapping[int, ArrayLike]) -> Mapping[int, tuple[float, ...]]:
    result = {}
    for key, vector in value.items():
        if isinstance(key, bool) or not isinstance(key, (int, np.integer)) or key < 0:
            raise ValidationError(
                "expected nonnegative integer load index", field="loads"
            )
        result[int(key)] = tuple(float(x) for x in real_array(vector, (3,)))
    return MappingProxyType(result)


@dataclass(frozen=True)
class LoadCase:
    """Named homogeneous supports and density-independent SI forces.

    conditions use mechanical_dofs(mesh). body_force is N/m³, tractions map
    Triangle3 rows to N/m², point_loads map mesh nodes to N. No callable loads.
    """

    name: str
    conditions: tuple[Dirichlet, ...]
    body_force: ArrayLike = (0.0, 0.0, 0.0)
    tractions: Mapping[int, ArrayLike] | None = None
    point_loads: Mapping[int, ArrayLike] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValidationError("load case needs a name", field="cases")
        conditions = []
        for bc in self.conditions:
            dofs = indices(bc.dofs)
            values = real_array(bc.values, (len(dofs),))
            if np.any(values != 0):
                raise ValidationError(
                    "compliance requires homogeneous supports", field="supports"
                )
            frozen_ids = np.frombuffer(dofs.tobytes(), dtype=np.int64)
            conditions.append(Dirichlet(frozen_ids, _immutable(values)))
        object.__setattr__(self, "conditions", tuple(conditions))
        object.__setattr__(
            self, "body_force", _immutable(real_array(self.body_force, (3,)))
        )
        object.__setattr__(self, "tractions", _vectors(self.tractions or {}))
        object.__setattr__(self, "point_loads", _vectors(self.point_loads or {}))


@dataclass(frozen=True)
class DensityEvaluation:
    """Densities/gradients in input-cell order; compliance J, volume fraction 1."""

    density: NDArray[np.float64]
    physical_density: NDArray[np.float64]
    objective: float
    volume_fraction: float
    objective_gradient: NDArray[np.float64]
    volume_gradient: NDArray[np.float64]
    states: tuple[ElasticResult, ...]

    def __post_init__(self) -> None:
        for key in (
            "density",
            "physical_density",
            "objective_gradient",
            "volume_gradient",
        ):
            object.__setattr__(self, key, _immutable(getattr(self, key)))


class TopologyProblem:
    """Volume-weighted linear density filter and SIMP on a fixed Tet4 mesh.

    Variables have shape (ncells,), fixed values included. Radius is metres;
    weights are positive, not implicitly normalized. Materials are homogeneous
    or one immutable positive tensor per cell. No projection or continuation.
    """

    def __init__(
        self,
        mesh: Mesh,
        material: ElasticMaterial | Sequence[ElasticMaterial],
        cases: Sequence[LoadCase],
        *,
        filter_radius: float,
        weights: ArrayLike | None = None,
        fixed_solid: ArrayLike = (),
        fixed_void: ArrayLike = (),
        penalty: float = 3.0,
        minimum_stiffness: float = 1e-3,
    ) -> None:
        if mesh.n_cells == 0 or mesh.n_tetra10 or mesh.n_triangle6:
            raise ValidationError("topology requires nonempty Tet4 mesh", field="mesh")
        values = np.asarray([filter_radius, penalty, minimum_stiffness], dtype=float)
        if (
            not np.isfinite(values).all()
            or filter_radius <= 0
            or penalty < 1
            or not 0 < minimum_stiffness < 1
        ):
            raise ValidationError(
                "invalid radius, SIMP penalty or minimum stiffness", field="topology"
            )
        self.mesh = mesh
        self.materials = tuple(
            [material] * mesh.n_cells
            if isinstance(material, ElasticMaterial)
            else material
        )
        if len(self.materials) != mesh.n_cells or not all(
            isinstance(m, ElasticMaterial) for m in self.materials
        ):
            raise ValidationError(
                "one elastic material required per cell", field="material"
            )
        self.cases = tuple(cases)
        if (
            not self.cases
            or not all(isinstance(c, LoadCase) for c in self.cases)
            or len({c.name for c in self.cases}) != len(self.cases)
        ):
            raise ValidationError("unique nonempty load cases required", field="cases")
        self.weights = _immutable(
            np.ones(len(cases))
            if weights is None
            else real_array(weights, (len(cases),))
        )
        if np.any(self.weights <= 0):
            raise ValidationError("load weights must be positive", field="weights")
        self.filter_radius = float(filter_radius)
        self.penalty = float(penalty)
        self.minimum_stiffness = float(minimum_stiffness)
        self.solid = indices(fixed_solid, mesh.n_cells)
        self.void = indices(fixed_void, mesh.n_cells)
        if len(np.intersect1d(self.solid, self.void)):
            raise ValidationError("solid and void regions overlap", field="design")
        self.free = np.setdiff1d(
            np.arange(mesh.n_cells), np.union1d(self.solid, self.void)
        )
        geometry = [Tet4Geometry(mesh.points[cell]) for cell in mesh.tetra4]
        self.volumes = _immutable([g.volume for g in geometry])
        centers = mesh.points[mesh.tetra4].mean(axis=1)
        tree = cKDTree(centers)
        rows, columns, coefficients = [], [], []
        for i, neighbors in enumerate(tree.query_ball_point(centers, filter_radius)):
            for j in sorted(neighbors):
                coefficient = (
                    max(
                        0.0,
                        filter_radius - float(np.linalg.norm(centers[i] - centers[j])),
                    )
                    * self.volumes[j]
                )
                if coefficient > 0:
                    rows.append(i)
                    columns.append(j)
                    coefficients.append(coefficient)
        h = sparse.csr_matrix(
            (coefficients, (rows, columns)), shape=(mesh.n_cells, mesh.n_cells)
        )
        self._filter = sparse.diags(1 / np.asarray(h.sum(axis=1)).ravel()) @ h
        mask = np.zeros(mesh.n_cells)
        mask[self.free] = 1
        self._derivative = sparse.diags(mask) @ self._filter @ sparse.diags(mask)
        self.volume_gradient = _immutable(
            self._derivative.T @ self.volumes / self.volumes.sum()
        )

    @classmethod
    def from_models(
        cls,
        models: Sequence[Model],
        *,
        weights: ArrayLike | None = None,
        **options: Any,
    ) -> TopologyProblem:
        """Share the validated Model parser; reject different meshes/materials."""
        if not models:
            raise ValidationError("at least one Model required", field="models")
        snapshots = [model.elasticity_inputs() for model in models]
        first = snapshots[0]
        cases = []
        for i, item in enumerate(snapshots):
            for key in ("points", "tetra4", "triangle3"):
                if not np.array_equal(
                    getattr(first["mesh"], key), getattr(item["mesh"], key)
                ):
                    raise ValidationError("load cases must share mesh", field="models")
            if len(first["material"]) != len(item["material"]) or any(
                not np.array_equal(a.stiffness, b.stiffness)
                for a, b in zip(first["material"], item["material"], strict=True)
            ):
                raise ValidationError("load cases must share material", field="models")
            cases.append(
                LoadCase(
                    f"case-{i}",
                    item["conditions"],
                    item["body_force"],
                    item["tractions"],
                )
            )
        return cls(first["mesh"], first["material"], cases, weights=weights, **options)

    def physical_density(self, density: ArrayLike) -> NDArray[np.float64]:
        """Map bounded design variables to physical densities, no threshold."""
        x = real_array(density, (self.mesh.n_cells,))
        if (
            np.any(x < 0)
            or np.any(x > 1)
            or np.any(x[self.solid] != 1)
            or np.any(x[self.void] != 0)
        ):
            raise ValidationError(
                "density bounds or fixed values violated", field="density"
            )
        rho = np.asarray(self._filter @ x).ravel()
        rho[self.solid] = 1
        rho[self.void] = 0
        return np.asarray(np.clip(rho, 0, 1), dtype=np.float64)

    def evaluate(self, density: ArrayLike) -> DensityEvaluation:
        """Solve independent load states and exact filtered compliance gradients."""
        rho = self.physical_density(density)
        a = self.minimum_stiffness + (1 - self.minimum_stiffness) * rho**self.penalty
        materials = tuple(
            ElasticMaterial(factor * m.stiffness)
            for factor, m in zip(a, self.materials, strict=True)
        )
        states = tuple(
            solve_elasticity(
                self.mesh,
                materials,
                case.conditions,
                body_force=case.body_force,
                tractions=case.tractions,
                point_loads=case.point_loads,
            )
            for case in self.cases
        )
        objective = 0.0
        energy = np.zeros(self.mesh.n_cells)
        for weight, state in zip(self.weights, states, strict=True):
            work = float(np.sum(state.applied_loads * state.displacement))
            scale = max(abs(work), 2 * state.energy, np.finfo(float).tiny)
            force_scale = max(
                float(np.linalg.norm(state.applied_loads, axis=1).sum()),
                np.finfo(float).tiny,
            )
            if (
                abs(work - 2 * state.energy) > 1e-8 * scale
                or np.linalg.norm(state.force_balance) > 1e-8 * force_scale
            ):
                raise ValidationError(
                    "state work/energy or force balance failed", field="solver"
                )
            objective += float(weight) * work
            energy += float(weight) * state.element_energy
        if not np.isfinite(objective) or objective <= 0:
            raise ValidationError(
                "positive compliance requires nonzero free loading", field="loads"
            )
        derivative = (
            -2
            * energy
            / a
            * self.penalty
            * (1 - self.minimum_stiffness)
            * rho ** (self.penalty - 1)
        )
        return DensityEvaluation(
            np.asarray(density),
            rho,
            objective,
            float(self.volumes @ rho / self.volumes.sum()),
            np.asarray(self._derivative.T @ derivative).ravel(),
            self.volume_gradient,
            states,
        )

    def initial_density(self, volume_fraction: float) -> NDArray[np.float64]:
        """Find feasible uniform free variables including fixed/filter effects."""
        if not np.isfinite(volume_fraction) or not 0 < volume_fraction <= 1:
            raise ValidationError("volume fraction must lie in (0,1]", field="volume")
        lower = np.zeros(self.mesh.n_cells)
        lower[self.solid] = 1
        upper = lower.copy()
        upper[self.free] = 1

        def volume(x: NDArray[np.float64]) -> float:
            return float(self.volumes @ self.physical_density(x) / self.volumes.sum())

        if volume(lower) > volume_fraction + 1e-12:
            raise ValidationError(
                "volume budget infeasible for fixed regions/filter", field="volume"
            )
        if volume(upper) <= volume_fraction:
            return upper
        lo, hi = 0.0, 1.0
        for _ in range(60):
            mid = (lo + hi) / 2
            x = lower.copy()
            x[self.free] = mid
            if volume(x) > volume_fraction:
                hi = mid
            else:
                lo = mid
        lower[self.free] = lo
        return lower
