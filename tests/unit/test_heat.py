"""Thermal material, operator and load references derived by integration."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.heat import conductance_matrix
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.physics.heat import (
    Convection,
    convection_terms,
    flux_load,
    source_load,
)


def test_conductivity_and_constant_null_mode():
    data = np.array([[2.0, 0.5, 0], [0.5, 3, 0], [0, 0, 4]])
    mat = ThermalMaterial(data)
    data[:] = 0
    g = Tet4Geometry([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    expected = (
        np.array(
            [[10, -2.5, -3.5, -4], [-2.5, 2, 0.5, 0], [-3.5, 0.5, 3, 0], [-4, 0, 0, 4]]
        )
        / 6
    )
    k = conductance_matrix(g, mat)
    assert_allclose(k, expected, atol=2e-15, rtol=0)
    assert_allclose(k @ np.ones(4), 0, atol=2e-15)
    assert np.count_nonzero(np.linalg.eigvalsh(k) > 1e-12) == 3
    assert_allclose(ThermalMaterial.isotropic(5).conductivity, np.eye(3) * 5, atol=0)
    with pytest.raises(ValueError):
        mat.conductivity[0, 0] = 1
    order = [0, 2, 1, 3]
    assert_allclose(
        conductance_matrix(Tet4Geometry(g.coordinates[order]), mat),
        k[np.ix_(order, order)],
        atol=2e-15,
    )


@pytest.mark.parametrize(
    "data",
    [
        np.eye(2),
        np.zeros((3, 3)),
        np.diag([1, 1, -1]),
        np.eye(3) * np.nan,
        np.eye(3, dtype=complex),
        [[1, 1, 0], [0, 1, 0], [0, 0, 1]],
    ],
)
def test_invalid_conductivity(data):
    with pytest.raises(ValidationError):
        ThermalMaterial(data)


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
def test_invalid_isotropic_conductivity(value):
    with pytest.raises(ValidationError):
        ThermalMaterial.isotropic(value)


def test_volume_surface_and_robin_integrals():
    g = Tet4Geometry([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert_allclose(source_load(g, 24), np.ones(4), atol=2e-15)
    assert_allclose(
        source_load(g, lambda x: x[0]), np.array([1, 2, 1, 1]) / 120, atol=2e-16
    )
    p = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    # Positive outward flux removes heat from the RHS.
    assert_allclose(flux_load(p, 6), -np.ones(3), atol=2e-15)
    assert_allclose(flux_load(p, lambda x: x[0]), -np.array([1, 2, 1]) / 24, atol=2e-16)
    h, f = convection_terms(p, Convection(12, 300))
    assert_allclose(h, (np.ones((3, 3)) + np.eye(3)) * 0.5, atol=2e-15)
    assert_allclose(f, np.full(3, 600), atol=2e-13)
    assert_allclose(h @ np.full(3, 300) - f, 0, atol=2e-13)
    h2, f2 = convection_terms(p, Convection(12, lambda x: 300 + x[0]))
    assert_allclose(h2, h, atol=0)
    assert_allclose(f2, h @ [300, 301, 300], atol=2e-13)
    with pytest.raises(ValidationError):
        Convection(-1, 300)
    with pytest.raises(ValidationError):
        flux_load(np.zeros((3, 3)), 1)
    with pytest.raises(ValidationError):
        source_load(g, lambda x: [1, 2])


def test_rotated_anisotropic_operator_and_surface_permutation():
    points = np.array([[0.1, 0.2, 0.3], [1.2, 0, 0.1], [0, 1.5, 0.2], [0.2, 0.3, 1.7]])
    rotation = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
    conductivity = np.array([[2.0, 0.3, 0.2], [0.3, 3, 0.1], [0.2, 0.1, 4]])
    original = conductance_matrix(Tet4Geometry(points), ThermalMaterial(conductivity))
    rotated = conductance_matrix(
        Tet4Geometry(points @ rotation.T + np.array([2, 3, 4])),
        ThermalMaterial(rotation @ conductivity @ rotation.T),
    )
    assert_allclose(rotated, original, atol=2e-14, rtol=2e-13)
    p = points[:3]
    order = [0, 2, 1]
    assert_allclose(
        flux_load(p[order], lambda x: x[0]),
        flux_load(p, lambda x: x[0])[order],
        atol=2e-15,
        rtol=0,
    )
