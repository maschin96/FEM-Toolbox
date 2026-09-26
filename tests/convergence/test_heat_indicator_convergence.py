"""A-priori criteria: first-order trend, bounded effectivity, localization."""

from itertools import product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.post import heat_indicator, mark_bulk
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs


def _energy_errors(mesh, result, conductivity, exact_gradient, *, gauss_order=None):
    # Independent symmetric four-point rule, exact for quadratic gradient error.
    bary = np.full((4, 4), (5 - np.sqrt(5)) / 20)
    np.fill_diagonal(bary, (5 + 3 * np.sqrt(5)) / 20)
    weights = np.full(4, 0.25)
    if gauss_order is not None:
        roots, gauss_weights = np.polynomial.legendre.leggauss(gauss_order)
        roots, gauss_weights = (roots + 1) / 2, gauss_weights / 2
        points, weights = [], []
        for a, b, c in product(range(gauss_order), repeat=3):
            r, t, u = roots[[a, b, c]]
            xyz = np.array([r, (1 - r) * t, (1 - r) * (1 - t) * u])
            points.append([1 - xyz.sum(), *xyz])
            weights.append(
                6
                * gauss_weights[a]
                * gauss_weights[b]
                * gauss_weights[c]
                * (1 - r) ** 2
                * (1 - t)
            )
        bary, weights = np.array(points), np.array(weights)
    values = []
    for i, cell in enumerate(mesh.tetra4):
        geometry = Tet4Geometry(mesh.points[cell])
        error = result.gradient[i] - exact_gradient(bary @ mesh.points[cell])
        values.append(
            geometry.volume
            * np.dot(weights, np.einsum("ni,ij,nj->n", error, conductivity, error))
        )
    return np.array(values)


def _solve(mesh, material, exact, source):
    nodes = np.unique(mesh.triangle3)
    result = solve_heat(
        mesh,
        material,
        [Dirichlet(thermal_dofs(mesh).indices("T", nodes), exact(mesh.points[nodes]))],
        source=source,
    )
    estimate = heat_indicator(
        mesh,
        result.temperature,
        material,
        source=source,
        dirichlet_facets=np.arange(mesh.n_triangle3),
    )
    return result, estimate


def test_indicator_four_mesh_energy_trend_and_effectivity(cube_mesh):
    material = ThermalMaterial([[2, 0.3, 0.2], [0.3, 3, 0.1], [0.2, 0.1, 4]])
    weights = np.array([1.0, 0.5, 0.25])

    def exact(x):
        return 300 + (x * (1 - x)) @ weights

    def gradient(x):
        return (1 - 2 * x) * weights

    levels = np.array([2, 4, 6, 8])
    errors, indicators = [], []
    for n in levels:
        mesh = cube_mesh(n)
        result, estimate = _solve(mesh, material, exact, 9.0)
        errors.append(
            np.sqrt(_energy_errors(mesh, result, material.conductivity, gradient).sum())
        )
        indicators.append(estimate.total)
        assert abs(result.power_balance) < 1e-9
    errors, indicators = np.array(errors), np.array(indicators)
    effectivity = indicators / errors
    rates = np.log(indicators[:-1] / indicators[1:]) / np.log(levels[1:] / levels[:-1])
    print(
        "levels, energy errors, eta, effectivity, rates:",
        levels,
        errors,
        indicators,
        effectivity,
        rates,
    )
    # Broad absolute bound: a heuristic need not tend to effectivity one.
    # A stable nonzero effectivity and order-one decay are required in advance.
    assert np.all((effectivity > 0.1) & (effectivity < 50))
    assert effectivity.max() / effectivity.min() < 1.3
    assert np.all((rates[-2:] > 0.8) & (rates[-2:] < 1.2))
    assert np.all(np.diff(indicators) < 0)


def test_localized_gradient_marks_the_error_region(cube_mesh):
    mesh = cube_mesh(8)
    material = ThermalMaterial.isotropic(1.0)

    def exact(x):
        return np.exp(-80 * (x[..., 0] - 0.5) ** 2)

    def source(x):
        return float((160 - 25600 * (x[0] - 0.5) ** 2) * exact(x))

    def gradient(x):
        result = np.zeros_like(x)
        result[:, 0] = -160 * (x[:, 0] - 0.5) * exact(x)
        return result

    result, estimate = _solve(mesh, material, exact, source)
    error = _energy_errors(mesh, result, material.conductivity, gradient, gauss_order=8)
    lower_order = _energy_errors(
        mesh, result, material.conductivity, gradient, gauss_order=7
    )
    np.testing.assert_allclose(error.sum(), lower_order.sum(), rtol=1e-5, atol=0)
    band = abs(mesh.points[mesh.tetra4].mean(axis=1)[:, 0] - 0.5) < 0.25
    marked = mark_bulk(estimate.cell, fraction=0.5)
    print(
        "band eta/error shares, marked fraction:",
        np.sum(estimate.cell[band] ** 2) / estimate.total**2,
        error[band].sum() / error.sum(),
        band[marked].mean(),
    )
    assert np.sum(estimate.cell[band] ** 2) / estimate.total**2 > 0.65
    assert error[band].sum() / error.sum() > 0.65
    assert band[marked].mean() > 0.8
    assert np.sum(estimate.cell[marked] ** 2) >= 0.5 * estimate.total**2
