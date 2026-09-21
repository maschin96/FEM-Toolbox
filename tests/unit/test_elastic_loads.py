import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError
from fem_toolbox.physics.elasticity import body_load, traction_load


def test_constant_and_linear_body_force():
    g = Tet4Geometry([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert_allclose(
        body_load(g, [4, 8, 12]).reshape(4, 3),
        np.tile([1, 2, 3], (4, 1)) / 6,
        atol=2e-16,
    )
    # integral Ni*x = V*(1+delta_i1)/20, since x=N1.
    expected = np.zeros((4, 3))
    expected[:, 0] = np.array([1, 2, 1, 1]) / 120
    assert_allclose(
        body_load(g, lambda x: [x[0], 0, 0]).reshape(4, 3), expected, atol=2e-16
    )


def test_surface_force_orientation_and_linear_quadrature():
    p = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0]])
    expected = np.tile([1.0, 2, 3], (3, 1)) / 6
    assert_allclose(traction_load(p, [1, 2, 3]).reshape(3, 3), expected, atol=2e-16)
    assert_allclose(
        traction_load(p[[0, 2, 1]], [1, 2, 3]).reshape(3, 3), expected, atol=2e-16
    )
    expected = np.zeros((3, 3))
    expected[:, 0] = np.array([1, 2, 1]) / 24
    assert_allclose(
        traction_load(p, lambda x: [x[0], 0, 0]).reshape(3, 3), expected, atol=2e-16
    )
    with pytest.raises(ValidationError):
        traction_load(np.zeros((3, 3)), [1, 0, 0])
    with pytest.raises(ValidationError):
        traction_load(p, [np.nan, 0, 0])
