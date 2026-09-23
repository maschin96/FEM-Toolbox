"""Analytical simplex integrals and rigid-translation storage invariants."""

import numpy as np
import pytest

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.storage import capacity_matrix, mass_matrix
from fem_toolbox.errors import ValidationError


@pytest.fixture
def geometry():
    # Volume = 1 m³, independently from base area 3 and height 1.
    return Tet4Geometry(
        np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]])
    )


@pytest.mark.parametrize("lumped", [False, True])
def test_storage_integrals_and_translation(geometry, lumped):
    mass = mass_matrix(geometry, 12.0, lumped=lumped)
    capacity = capacity_matrix(geometry, 12.0, 5.0, lumped=lumped)
    np.testing.assert_allclose(capacity.sum(axis=1), 15.0, rtol=1e-14)
    np.testing.assert_allclose(capacity.sum(), 60.0, rtol=1e-14)
    velocity = np.tile([2.0, -3.0, 4.0], 4)
    assert velocity @ mass @ velocity / 2 == pytest.approx(174.0, rel=1e-14)
    np.testing.assert_allclose(mass, mass.T, atol=0)
    assert np.linalg.eigvalsh(mass).min() > 0
    if lumped:
        np.testing.assert_allclose(capacity, np.eye(4) * 15.0, rtol=1e-14)
    else:
        # Integral Ni² = V/10; integral Ni Nj = V/20 for i != j.
        assert capacity[0, 0] == pytest.approx(6.0, rel=1e-14)
        assert capacity[0, 1] == pytest.approx(3.0, rel=1e-14)
        # Integral x² = V * 2²/10 for this tetrahedron.
        x = np.array([0.0, 2.0, 0.0, 0.0])
        assert x @ capacity @ x == pytest.approx(24.0, rel=1e-14)


@pytest.mark.parametrize("bad", [0.0, -1.0, np.nan, np.inf, True, 1j, [1.0]])
def test_invalid_storage_parameters(geometry, bad):
    with pytest.raises(ValidationError):
        mass_matrix(geometry, bad)
    with pytest.raises(ValidationError):
        capacity_matrix(geometry, 1.0, bad)


def test_lumping_is_explicit_boolean(geometry):
    with pytest.raises(ValidationError):
        mass_matrix(geometry, 1.0, lumped="false")


def test_unrepresentable_storage_is_rejected(geometry):
    with pytest.raises(ValidationError):
        capacity_matrix(geometry, 1e300, 1e300)


def test_node_permutation_and_volume_scaling():
    points = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 1.0]]
    )
    order = np.array([2, 0, 3, 1])
    original = mass_matrix(Tet4Geometry(points), 7.0)
    permuted = mass_matrix(Tet4Geometry(points[order]), 7.0)
    dofs = (3 * order[:, None] + np.arange(3)).ravel()
    np.testing.assert_allclose(permuted, original[np.ix_(dofs, dofs)], rtol=1e-14)
    scaled = mass_matrix(Tet4Geometry(2 * points), 7.0)
    np.testing.assert_allclose(scaled, 8 * original, rtol=1e-14)
