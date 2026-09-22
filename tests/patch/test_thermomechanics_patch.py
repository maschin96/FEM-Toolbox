"""Independent SI and material-axis benchmarks for free/restrained expansion."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion
from fem_toolbox.solvers.elasticity import solve_elasticity


def test_rotated_expansion_in_si_free_and_restrained(cube_mesh):
    mesh = cube_mesh(1)
    # Principal alpha_x=20 microstrain/K, rotated 90 degrees about z.
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    expansion = ThermalExpansion.from_principal([20e-6, 0, 0], 293.15, rotation)
    material = ElasticMaterial.isotropic(100e9, 0.25)
    temperature = np.full(mesh.n_points, 318.15)
    expected_u = mesh.points * [0, 500e-6, 0]

    origin = int(np.flatnonzero(np.all(mesh.points == [0, 0, 0], axis=1))[0])
    x_node = int(np.flatnonzero(np.all(mesh.points == [1, 0, 0], axis=1))[0])
    y_node = int(np.flatnonzero(np.all(mesh.points == [0, 1, 0], axis=1))[0])
    fixed = [
        3 * origin,
        3 * origin + 1,
        3 * origin + 2,
        3 * x_node + 1,
        3 * x_node + 2,
        3 * y_node + 2,
    ]
    free = solve_elasticity(
        mesh,
        material,
        [Dirichlet(fixed, expected_u.ravel()[fixed])],
        expansion=expansion,
        temperature=temperature,
    )
    assert_allclose(free.displacement, expected_u, rtol=0, atol=1e-18)
    assert_allclose(free.stress, 0, rtol=0, atol=1e-5)
    assert free.energy < 1e-18

    restrained = solve_elasticity(
        mesh,
        material,
        [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
        expansion=expansion,
        temperature=temperature,
    )
    # With E=100 GPa, nu=1/4, C11=120 GPa and C12=40 GPa.
    # sigma=-C*epsilon_th, where epsilon_yy=500 microstrain.
    assert_allclose(
        restrained.stress,
        np.tile([-20e6, -60e6, -20e6, 0, 0, 0], (mesh.n_tetra4, 1)),
        rtol=0,
        atol=1e-5,
    )
    assert_allclose(restrained.energy, 15_000, rtol=0, atol=1e-8)
    top = mesh.points[restrained.node_ids, 1] == 1
    assert_allclose(restrained.reactions[top, 1].sum(), -60e6, rtol=0, atol=1e-5)
    assert_allclose(restrained.reactions[~top, 1].sum(), 60e6, rtol=0, atol=1e-5)
    assert np.linalg.norm(restrained.force_balance) < 1e-5
