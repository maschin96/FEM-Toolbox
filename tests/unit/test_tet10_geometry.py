"""Independent polynomial identities for straight quadratic simplex geometry."""

from itertools import permutations, product
from math import factorial

import numpy as np
import pytest
from numpy.testing import assert_allclose

EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
CORNERS = np.array(
    [[0.2, -0.3, 0.1], [2, 0.1, 0.3], [0.4, 1.2, -0.2], [-0.1, 0.2, 1.4]]
)


def points(corners=CORNERS):
    return np.vstack([corners, [(corners[a] + corners[b]) / 2 for a, b in EDGES]])


def test_nodal_basis_partition_and_quadratic_reproduction():
    from fem_toolbox.elements import Tet10Geometry

    xyz = points()
    geometry = Tet10Geometry(xyz)
    assert_allclose(
        [geometry.shape_values(p) for p in xyz], np.eye(10), atol=2e-14, rtol=0
    )
    p = np.array([0.2, 0.3, 0.1, 0.4]) @ CORNERS
    n, g = geometry.shape_values(p), geometry.shape_gradients(p)
    assert_allclose(n.sum(), 1, atol=2e-14, rtol=0)
    assert_allclose(g.sum(axis=0), 0, atol=2e-14, rtol=0)
    for powers in product(range(3), repeat=3):
        if sum(powers) > 2:
            continue
        values = np.prod(xyz**powers, axis=1)
        expected = np.prod(p**powers)
        derivative = np.zeros(3)
        for axis, power in enumerate(powers):
            if power:
                reduced = list(powers)
                reduced[axis] -= 1
                derivative[axis] = power * np.prod(p**reduced)
        assert_allclose(n @ values, expected, atol=2e-14, rtol=2e-14)
        assert_allclose(g.T @ values, derivative, atol=2e-14, rtol=2e-14)
    assert_allclose(geometry.gradient_matrix(p), g.T, atol=0, rtol=0)


@pytest.mark.parametrize("dimension", [2, 3])
def test_quadrature_integrates_all_reference_monomials_through_advertised_degree(
    dimension,
):
    from fem_toolbox.elements import tetrahedron_quadrature, triangle_quadrature

    rule = (triangle_quadrature if dimension == 2 else tetrahedron_quadrature)()
    assert rule.degree >= 4
    assert np.all(rule.weights > 0)
    assert np.all(rule.points > 0) and np.all(rule.points.sum(axis=1) < 1)
    for powers in product(range(rule.degree + 1), repeat=dimension):
        if sum(powers) > rule.degree:
            continue
        exact = np.prod([factorial(i) for i in powers]) / factorial(
            sum(powers) + dimension
        )
        actual = rule.weights @ np.prod(rule.points**powers, axis=1)
        assert_allclose(actual, exact, rtol=2e-14, atol=1e-16)
    with pytest.raises(ValueError):
        rule.weights[0] = 0


@pytest.mark.parametrize("order", list(permutations(range(4))))
def test_orientation_all_corner_permutations_and_surface_trace(order):
    from fem_toolbox.elements import Tet10Geometry, Triangle6Geometry

    xyz = points(CORNERS[list(order)])
    geometry = Tet10Geometry(xyz)
    positive = Tet10Geometry(geometry.oriented_coordinates)
    assert positive.signed_volume > 0
    ids = np.arange(10) + 100
    assert_allclose(
        geometry.oriented_connectivity(ids),
        ids[list(geometry.oriented_node_order)],
        atol=0,
        rtol=0,
    )
    assert_allclose(positive.volume, geometry.volume, rtol=2e-14)
    volume_rule = geometry.quadrature()
    assert_allclose(volume_rule.weights.sum(), geometry.volume, rtol=2e-14)
    flux = np.zeros(3)
    for face in geometry.faces:
        nodes = list(face.local_nodes)
        surface = Triangle6Geometry(xyz[nodes])
        rule = surface.quadrature()
        assert_allclose(rule.weights.sum(), face.area, rtol=2e-14)
        assert_allclose(surface.unit_normal, face.unit_normal, atol=2e-14, rtol=0)
        p = np.array([0.2, 0.3, 0.5]) @ xyz[nodes[:3]]
        trace = np.zeros(10)
        trace[nodes] = surface.shape_values(p)
        assert_allclose(geometry.shape_values(p), trace, atol=2e-14, rtol=0)
        assert_allclose(surface.shape_gradients(p).sum(axis=0), 0, atol=2e-14, rtol=0)
        flux += face.area * face.unit_normal
    assert_allclose(flux, 0, atol=3e-14, rtol=0)


def test_curved_degenerate_and_invalid_inputs_rejected():
    from fem_toolbox.elements import Tet10Geometry
    from fem_toolbox.errors import ValidationError

    for bad in (
        points()[:9],
        points().astype(complex),
        np.full((10, 3), np.nan),
        np.zeros((10, 3)),
    ):
        with pytest.raises(ValidationError):
            Tet10Geometry(bad)
    curved = points()
    curved[4, 2] += 0.01
    with pytest.raises(ValidationError, match="curved"):
        Tet10Geometry(curved)


def test_triangle_nodal_basis_tangential_derivatives_and_input_validation():
    from fem_toolbox.elements import Triangle6Geometry
    from fem_toolbox.errors import ValidationError

    xyz = points()[[0, 1, 2, 4, 5, 6]]
    surface = Triangle6Geometry(xyz)
    assert_allclose(
        [surface.shape_values(p) for p in xyz], np.eye(6), atol=2e-14, rtol=0
    )
    p = np.array([0.2, 0.3, 0.5]) @ xyz[:3]
    projector = np.eye(3) - np.outer(surface.unit_normal, surface.unit_normal)
    for powers in product(range(3), repeat=3):
        if sum(powers) > 2:
            continue
        nodal = np.prod(xyz**powers, axis=1)
        exact_gradient = np.zeros(3)
        for axis, power in enumerate(powers):
            if power:
                reduced = list(powers)
                reduced[axis] -= 1
                exact_gradient[axis] = power * np.prod(p**reduced)
        assert_allclose(
            surface.shape_values(p) @ nodal, np.prod(p**powers), atol=2e-14, rtol=2e-14
        )
        assert_allclose(
            surface.shape_gradients(p).T @ nodal,
            projector @ exact_gradient,
            atol=2e-14,
            rtol=2e-14,
        )
    with pytest.raises(ValidationError, match="plane"):
        surface.shape_values(p + 0.1 * surface.unit_normal)
    for bad in (xyz[:5], np.zeros((6, 3)), np.full((6, 3), np.inf)):
        with pytest.raises(ValidationError):
            Triangle6Geometry(bad)
    curved = xyz.copy()
    curved[3] += 0.01 * surface.unit_normal
    with pytest.raises(ValidationError, match="curved"):
        Triangle6Geometry(curved)


@pytest.mark.parametrize("scale", [1e-6, 1.0, 1e6])
def test_units_immutable_storage_and_gradient_scaling(scale):
    from fem_toolbox.elements import Tet10Geometry
    from fem_toolbox.errors import ValidationError

    original = points()
    translated = scale * (original + np.array([3, -2, 1]))
    geometry = Tet10Geometry(translated)
    translated[:] = 0
    p = CORNERS.mean(axis=0)
    baseline = Tet10Geometry(original)
    assert_allclose(geometry.volume, baseline.volume * scale**3, rtol=2e-14)
    assert_allclose(
        geometry.shape_values(scale * (p + np.array([3, -2, 1]))),
        baseline.shape_values(p),
        atol=3e-14,
        rtol=0,
    )
    assert_allclose(
        geometry.shape_gradients(scale * (p + np.array([3, -2, 1]))) * scale,
        baseline.shape_gradients(p),
        atol=3e-14,
        rtol=0,
    )
    with pytest.raises(ValueError):
        geometry.coordinates.setflags(write=True)
    for bad in ([0, 0, 0], [np.nan, 0, 0], np.zeros(10), np.arange(10, dtype=float)):
        with pytest.raises(ValidationError):
            geometry.oriented_connectivity(bad)
    for bad in ([0, 0], [np.inf, 0, 0], [1j, 0, 0]):
        with pytest.raises(ValidationError):
            geometry.shape_gradients(bad)
