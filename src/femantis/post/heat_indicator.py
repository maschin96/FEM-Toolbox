"""Residual indicator for stationary Tet4 heat; not a guaranteed error bound."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.assembly import indices, real_array
from femantis.elements import Tet4Geometry
from femantis.elements.quadrature import tetrahedron_quadrature, triangle_quadrature
from femantis.errors import ValidationError
from femantis.materials import ThermalMaterial
from femantis.mesh import Mesh, inspect_boundary_topology
from femantis.physics.heat import Convection, ScalarField


@dataclass(frozen=True)
class HeatIndicator:
    """Cell-row indicators sqrt(W K), squared contributions W K, global l2 norm.

    Arrays own immutable storage. Interior-face contributions are divided equally
    between their two owners. This is a heuristic energy-error indicator, without
    reliability/efficiency constants or robustness guarantees for high contrast.
    """

    volume_squared: NDArray[np.float64]
    interior_squared: NDArray[np.float64]
    boundary_squared: NDArray[np.float64]

    def __post_init__(self) -> None:
        shape = np.asarray(self.volume_squared).shape
        if len(shape) != 1:
            raise ValidationError("expected one value per cell", field="indicator")
        for name in ("volume_squared", "interior_squared", "boundary_squared"):
            values = real_array(getattr(self, name), shape)
            if np.any(values < 0):
                raise ValidationError("negative squared indicator", field="indicator")
            object.__setattr__(self, name, np.frombuffer(values.tobytes(), dtype=float))

    @property
    def cell(self) -> NDArray[np.float64]:
        """Local indicator eta_K in sqrt(W K), ordered as mesh.tetra4."""
        values = np.sqrt(
            self.volume_squared + self.interior_squared + self.boundary_squared
        )
        return np.frombuffer(values.tobytes(), dtype=float)

    @property
    def total(self) -> float:
        """Global indicator sqrt(sum_K eta_K²) in sqrt(W K)."""
        return float(np.linalg.norm(self.cell))


def mark_bulk(indicators: ArrayLike, *, fraction: float = 0.5) -> NDArray[np.int64]:
    """Minimal-cardinality Dörfler set capturing fraction of sum eta_K².

    Nonnegative finite eta_K; 0 < fraction <= 1. Ties favor lower cell rows;
    returned rows are sorted ascending. Zero total gives an empty set. Rescaling
    before squaring avoids overflow/underflow under a common unit-scale change.
    No numerical noise cutoff is inferred from dimensional values.
    """
    eta = real_array(indicators, np.asarray(indicators).shape)
    theta = float(real_array(fraction, ()))
    if eta.ndim != 1 or np.any(eta < 0) or not 0 < theta <= 1:
        raise ValidationError("invalid indicators or bulk fraction", field="marking")
    positive = np.flatnonzero(eta > 0)
    if not len(positive):
        return np.empty(0, dtype=np.int64)
    order = positive[np.argsort(-eta[positive], kind="stable")]
    if theta == 1:
        return np.sort(order)
    squared = (eta[order] / eta[order[0]]) ** 2
    cumulative = np.cumsum(squared)
    count = int(np.searchsorted(cumulative, theta * cumulative[-1])) + 1
    return np.sort(order[:count])


def _value(field: ScalarField, point: NDArray[np.float64]) -> float:
    return float(real_array(field(point) if callable(field) else field, ()))


def _diameter(points: NDArray[np.float64]) -> float:
    return max(float(np.linalg.norm(a - b)) for a, b in combinations(points, 2))


def heat_indicator(
    mesh: Mesh,
    temperature: ArrayLike,
    material: ThermalMaterial | Sequence[ThermalMaterial],
    *,
    source: ScalarField = 0.0,
    dirichlet_facets: Iterable[int] = (),
    outward_fluxes: Mapping[int, ScalarField] | None = None,
    convection: Mapping[int, Convection] | None = None,
) -> HeatIndicator:
    """Estimate stationary Tet4 heat energy error from full-mesh nodal T (K).

    Temperature has shape (mesh.n_points,). Conductivity is cellwise constant,
    symmetric positive definite, in W/(m K). Source Q is W/m³, prescribed outward
    flux W/m². Boundary arguments address triangle3 rows; full exterior coverage
    is mandatory. Explicit Dirichlet facets have no residual term; other facets
    default to insulation. Flux and Robin terms add, as in solve_heat. Dirichlet
    cannot overlap either. Caller must supply the actual problem's boundary data.

    eta_K² = h_K²/lambda_K * integral_K Q²
      + 1/2 sum_inner h_F/lambda_F * integral_F jump(q.n)²
      + sum_N/R h_F/lambda_K * integral_F (q.n-q_pres-h*(T-T_inf))².
    h denotes diameter, lambda_K=min eig(k_K), lambda_F=min of adjacent values.
    q=-k grad(T_h); div(q)=0 inside affine Tet4 cells. Volume degree-5 and surface
    degree-4 quadrature integrate squared residuals; nonpolynomial data require
    separate quadrature assessment. Norm targets sqrt(integral grad(e).k.grad(e))
    with units sqrt(W K); no certified error bound, data-oscillation correction,
    transient, coupled, Tet10 or automatic remeshing support is implied.
    """
    if mesh.n_tetra4 == 0 or mesh.n_tetra10:
        raise ValidationError("indicator requires a nonempty Tet4 mesh", field="mesh")
    temperatures = real_array(temperature, (mesh.n_points,))
    materials = (
        [material] * mesh.n_tetra4
        if isinstance(material, ThermalMaterial)
        else list(material)
    )
    if len(materials) != mesh.n_tetra4 or any(
        not isinstance(m, ThermalMaterial) for m in materials
    ):
        raise ValidationError("one ThermalMaterial required per cell", field="material")
    topology = inspect_boundary_topology(mesh, require_complete=True)
    fixed = set(indices(list(dirichlet_facets), mesh.n_triangle3).tolist())
    fluxes = dict(outward_fluxes or {})
    robin = dict(convection or {})
    indices(list(fluxes), mesh.n_triangle3)
    indices(list(robin), mesh.n_triangle3)
    if fixed.intersection(fluxes) or fixed.intersection(robin):
        raise ValidationError(
            "Dirichlet and flux/Robin facets overlap", field="surface"
        )
    if any(not isinstance(value, Convection) for value in robin.values()):
        raise ValidationError("expected Convection boundary data", field="convection")
    geometries = [Tet4Geometry(mesh.points[cell]) for cell in mesh.tetra4]
    lambdas = np.array([np.linalg.eigvalsh(m.conductivity)[0] for m in materials])
    flux = np.array(
        [
            -m.conductivity
            @ g.gradient_matrix
            @ (temperatures[cell] - temperatures[cell[0]])
            for cell, g, m in zip(mesh.tetra4, geometries, materials, strict=True)
        ]
    )
    volume = np.zeros(mesh.n_tetra4)
    interior = np.zeros(mesh.n_tetra4)
    boundary = np.zeros(mesh.n_tetra4)
    rule = tetrahedron_quadrature()
    bary = np.column_stack((1 - rule.points.sum(axis=1), rule.points))
    # Store internal rows rather than external IDs, which need not be contiguous.
    faces: dict[tuple[int, ...], list[tuple[int, NDArray[np.float64], float]]] = {}
    for i, (cell, g) in enumerate(zip(mesh.tetra4, geometries, strict=True)):
        values = np.array([_value(source, p) for p in bary @ mesh.points[cell]])
        volume[i] = (
            _diameter(mesh.points[cell]) ** 2
            / lambdas[i]
            * 6
            * g.volume
            * float(rule.weights @ values**2)
        )
        for face in g.faces:
            key = tuple(sorted(int(cell[n]) for n in face.local_nodes))
            faces.setdefault(key, []).append((i, face.unit_normal, face.area))
    for nodes, owners in faces.items():
        if len(owners) == 2:
            (a, normal_a, area), (b, normal_b, _) = owners
            jump = float(flux[a] @ normal_a + flux[b] @ normal_b)
            term = _diameter(mesh.points[list(nodes)]) / min(lambdas[a], lambdas[b])
            term *= area * jump**2 / 2
            interior[a] += term
            interior[b] += term
    surface_rule = triangle_quadrature()
    surface_bary = np.column_stack(
        (1 - surface_rule.points.sum(axis=1), surface_rule.points)
    )
    for row, boundary_face in enumerate(topology.triangle_faces):
        if row in fixed:
            continue
        surface_nodes = list(boundary_face.nodes)
        points = mesh.points[surface_nodes]
        i = boundary_face.cell_index
        normal_flux = float(flux[i] @ boundary_face.unit_normal)
        residuals = []
        for n, p in zip(surface_bary, surface_bary @ points, strict=True):
            prescribed = _value(fluxes.get(row, 0.0), p)
            if row in robin:
                condition = robin[row]
                prescribed += condition.coefficient * (
                    float(n @ temperatures[surface_nodes])
                    - _value(condition.ambient, p)
                )
            residuals.append(normal_flux - prescribed)
        boundary[i] += (
            _diameter(points)
            / lambdas[i]
            * 2
            * boundary_face.area
            * float(surface_rule.weights @ np.square(residuals))
        )
    return HeatIndicator(volume, interior, boundary)
