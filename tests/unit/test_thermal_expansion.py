import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.solvers.heat import thermal_dofs
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics


def test_free_anisotropic_expansion(cube_mesh):
    mesh = cube_mesh(1)
    alpha = ThermalExpansion([0.01, 0.02, 0.03, 0, 0, 0], 290)
    expected = mesh.points * [0.1, 0.2, 0.3]
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
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(1000, 0.25),
        [Dirichlet(fixed, expected.ravel()[fixed])],
        expansion=alpha,
        temperature=np.full(mesh.n_points, 300),
    )
    assert_allclose(result.displacement, expected, atol=2e-15)
    assert_allclose(result.stress, 0, atol=2e-12)
    assert_allclose(result.energy, 0, atol=1e-27)


def test_restrained_expansion_reaction(cube_mesh):
    mesh = cube_mesh(1)
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(1000, 0.25),
        [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
        expansion=ThermalExpansion([0.01, 0.01, 0.01, 0, 0, 0], 290),
        temperature=np.full(mesh.n_points, 300),
    )
    assert_allclose(
        result.stress,
        -200 * np.ones((mesh.n_tetra4, 6)) * [1, 1, 1, 0, 0, 0],
        atol=2e-12,
    )
    assert_allclose(result.force_balance, 0, atol=2e-12)
    assert_allclose(result.energy, 30, atol=2e-12)
    right = mesh.points[result.node_ids, 0] == 1
    assert_allclose(result.reactions[right, 0].sum(), -200, atol=2e-12)
    assert_allclose(result.reactions[~right, 0].sum(), 200, atol=2e-12)


def test_invalid_temperature(cube_mesh):
    mesh = cube_mesh(1)
    with pytest.raises(ValidationError):
        solve_elasticity(
            mesh,
            ElasticMaterial.isotropic(100, 0.25),
            [],
            expansion=ThermalExpansion([0.01, 0, 0, 0, 0, 0], 290),
            temperature=[300],
        )


def test_sequential_heat_to_mechanics(cube_mesh):
    mesh = cube_mesh(1)
    heat_dofs = thermal_dofs(mesh)
    mechanical = mechanical_dofs(mesh)
    all_nodes = np.arange(mesh.n_points)
    expected = mesh.points * [0.1, 0.2, 0.3]
    result = solve_thermomechanics(
        mesh,
        ThermalMaterial.isotropic(1),
        ElasticMaterial.isotropic(1000, 0.25),
        [ThermalExpansion([0.01, 0.02, 0.03, 0, 0, 0], 290)] * mesh.n_tetra4,
        [Dirichlet(heat_dofs.indices("T", all_nodes), np.full(mesh.n_points, 300))],
        [Dirichlet(mechanical.indices("u", all_nodes), expected.ravel())],
    )
    assert_allclose(result.heat.temperature, 300)
    assert_allclose(result.elasticity.stress, 0, atol=2e-12)


def test_regionwise_reference_temperature(cube_mesh):
    mesh = cube_mesh(1)
    expansions = [
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 290 if i == 0 else 300)
        for i in range(mesh.n_tetra4)
    ]
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(1000, 0.25),
        [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
        expansion=expansions,
        temperature=np.full(mesh.n_points, 300),
    )
    assert np.linalg.norm(result.stress[0]) > 0
    assert_allclose(result.stress[1:], 0, atol=0)


def test_principal_expansion_rotates_with_material_axes():
    angle = np.pi / 4
    axes = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    expansion = ThermalExpansion.from_principal([0.02, 0, 0], 300, axes)
    assert_allclose(expansion.alpha, [0.01, 0.01, 0, 0, 0, 0.02], atol=1e-15)
    with pytest.raises(ValidationError):
        ThermalExpansion.from_principal([0.02, 0, 0], 300, np.eye(3) * 2)
