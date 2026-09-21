"""Independent constitutive and element references, SI units."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.elasticity import stiffness_matrix, strain_matrix
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ElasticMaterial

POINTS = np.array([[0.0, 0, 0], [2, 0, 0], [0.2, 1, 0], [0.1, 0.2, 1.5]])


def test_isotropic_stiffness_and_pure_shear():
    material = ElasticMaterial.isotropic(120.0, 0.2)
    # lambda=100/3, mu=50; engineering shear stress = mu*gamma.
    expected = np.diag([400 / 3] * 3 + [50.0] * 3)
    expected[:3, :3] += (np.ones((3, 3)) - np.eye(3)) * (100 / 3)
    assert_allclose(material.stiffness, expected, rtol=2e-15, atol=2e-14)
    assert_allclose(
        material.stiffness @ [0, 0, 0, 0.1, 0, 0], [0, 0, 0, 5, 0, 0], atol=1e-14
    )
    with pytest.raises(ValueError):
        material.stiffness[0, 0] = 1


@pytest.mark.parametrize(
    "young,poisson",
    [(0, 0.2), (-1, 0.2), (1, -1), (1, 0.5), (np.inf, 0.2), (1, np.nan), (1, 0.499)],
)
def test_invalid_isotropic_material(young, poisson):
    with pytest.raises(ValidationError):
        ElasticMaterial.isotropic(young, poisson)


@pytest.mark.parametrize(
    "matrix",
    [
        np.eye(3),
        np.zeros((6, 6)),
        -np.eye(6),
        np.eye(6) * np.nan,
        np.eye(6, dtype=complex),
        np.eye(6) + np.triu(np.ones((6, 6)), 1),
    ],
)
def test_invalid_anisotropic_material(matrix):
    with pytest.raises(ValidationError):
        ElasticMaterial(matrix)


def test_anisotropic_energy_and_copy():
    c = np.diag([10.0, 20, 30, 4, 5, 6])
    c[0, 1] = c[1, 0] = 2
    material = ElasticMaterial(c)
    c[:] = 0
    e = np.array([0.1, 0.2, 0, 0, 0.3, 0])
    assert_allclose(0.5 * e @ material.stiffness @ e, 0.715, atol=2e-15)


def test_affine_strain_and_six_rigid_modes():
    geometry = Tet4Geometry(POINTS)
    b = strain_matrix(geometry)
    a = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]])
    u = POINTS @ a.T + [1, 2, 3]
    assert_allclose(b @ u.ravel(), [0.1, 0.5, 0.9, 1.4, 1.0, 0.6], atol=2e-15)
    k = stiffness_matrix(geometry, ElasticMaterial.isotropic(120, 0.2))
    assert_allclose(k, k.T, atol=2e-14)
    assert np.count_nonzero(np.linalg.eigvalsh(k) > 1e-10) == 6
    for direction in np.eye(3):
        for motion in [np.tile(direction, (4, 1)), np.cross(direction, POINTS)]:
            assert_allclose(b @ motion.ravel(), 0, atol=2e-15)
            assert_allclose(k @ motion.ravel(), 0, atol=5e-14)
    permutation = [0, 2, 1, 3]
    kp = stiffness_matrix(
        Tet4Geometry(POINTS[permutation]), ElasticMaterial.isotropic(120, 0.2)
    )
    ids = np.array([[3 * n, 3 * n + 1, 3 * n + 2] for n in permutation]).ravel()
    assert_allclose(kp, k[np.ix_(ids, ids)], atol=5e-14)
