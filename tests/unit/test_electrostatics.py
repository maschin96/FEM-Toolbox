"""Electrical operators and charges, independent hand-integrated references."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.electrostatics import dielectric_matrix
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.physics.electrostatics import surface_charge_load, volume_charge_load


def test_dielectric_matrix_nullmode_and_orientation():
    data = 1e-11 * np.array([[2, 0.5, 0], [0.5, 3, 0], [0, 0, 4]])
    m = DielectricMaterial(data)
    data[:] = 0
    g = Tet4Geometry([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    expected = (
        1e-11
        * np.array(
            [[10, -2.5, -3.5, -4], [-2.5, 2, 0.5, 0], [-3.5, 0.5, 3, 0], [-4, 0, 0, 4]]
        )
        / 6
    )
    k = dielectric_matrix(g, m)
    assert_allclose(k, expected, atol=2e-26, rtol=0)
    assert_allclose(k @ np.ones(4), 0, atol=2e-26, rtol=0)
    assert np.count_nonzero(np.linalg.eigvalsh(k) > 1e-23) == 3
    order = [0, 2, 1, 3]
    assert_allclose(
        dielectric_matrix(Tet4Geometry(g.coordinates[order]), m),
        k[np.ix_(order, order)],
        atol=2e-26,
        rtol=0,
    )
    with pytest.raises(ValueError):
        m.permittivity[0, 0] = 1
    assert_allclose(
        DielectricMaterial.isotropic(2e-11).permittivity, np.eye(3) * 2e-11, atol=0
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
def test_invalid_permittivity(data):
    with pytest.raises(ValidationError):
        DielectricMaterial(data)


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
def test_invalid_isotropic_permittivity(value):
    with pytest.raises(ValidationError):
        DielectricMaterial.isotropic(value)


def test_signed_charge_integrals():
    g = Tet4Geometry([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert_allclose(volume_charge_load(g, 24e-9), np.full(4, -1e-9), atol=1e-24, rtol=0)
    assert_allclose(
        volume_charge_load(g, lambda x: x[0]),
        -np.array([1, 2, 1, 1]) / 120,
        atol=2e-16,
        rtol=0,
    )
    p = g.coordinates[:3]
    assert_allclose(surface_charge_load(p, 6e-9), np.full(3, 1e-9), atol=1e-24, rtol=0)
    assert_allclose(
        surface_charge_load(p, lambda x: x[0]),
        np.array([1, 2, 1]) / 24,
        atol=2e-16,
        rtol=0,
    )
    with pytest.raises(ValidationError):
        surface_charge_load(np.zeros((3, 3)), 1)
    with pytest.raises(ValidationError):
        volume_charge_load(g, lambda x: [1, 2])
