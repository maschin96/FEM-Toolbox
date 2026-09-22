"""Global affine piezoelectric patch on a distorted Tet4 mesh."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric


def test_affine_coupled_patch_with_free_interior_node(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    points = mesh.points
    stiffness = np.diag([10.0, 20, 30, 4, 5, 6])
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 2.0
    coupling[1, 1] = 1.0
    material = PiezoelectricMaterial(stiffness, coupling, np.diag([3.0, 4, 5]))
    exact_u = np.column_stack(
        [
            0.2 * points[:, 0] + 0.03 * points[:, 1],
            -0.1 * points[:, 1],
            0.05 * points[:, 2],
        ]
    )
    exact_phi = 0.3 * points[:, 0] - 0.2 * points[:, 1] + 0.1 * points[:, 2] + 2
    boundary = np.unique(mesh.triangle3)
    assert len(boundary) < mesh.n_points
    dofs = piezoelectric_dofs(mesh)
    result = solve_piezoelectric(
        mesh,
        material,
        [
            Dirichlet(dofs.indices("u", boundary), exact_u[boundary].ravel()),
            Dirichlet(dofs.indices("phi", boundary), exact_phi[boundary]),
        ],
    )
    assert_allclose(result.displacement, exact_u, atol=4e-13, rtol=0)
    assert_allclose(result.potential, exact_phi, atol=4e-13, rtol=0)
    assert_allclose(
        result.stress,
        [[2.6, -2.2, 1.5, 0, 0, 0.18]] * mesh.n_tetra4,
        atol=3e-12,
        rtol=0,
    )
    assert_allclose(
        result.electric_displacement,
        [[-0.5, 0.7, -0.5]] * mesh.n_tetra4,
        atol=3e-12,
        rtol=0,
    )
    assert result.relative_residual_u < 1e-13
    assert result.relative_residual_phi < 1e-13
    assert_allclose(result.force_balance, 0, atol=2e-13, rtol=0)
    assert_allclose(result.moment_balance, 0, atol=2e-13, rtol=0)
    assert_allclose(result.charge_balance, 0, atol=2e-13, rtol=0)


def test_affine_patch_with_tractions_and_prescribed_electric_flux(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    points = mesh.points
    stiffness = np.diag([10.0, 20, 30, 4, 5, 6])
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 2.0
    coupling[1, 1] = 1.0
    material = PiezoelectricMaterial(stiffness, coupling, np.diag([3.0, 4, 5]))
    exact_u = np.column_stack(
        [
            0.2 * points[:, 0] + 0.03 * points[:, 1],
            -0.1 * points[:, 1],
            0.05 * points[:, 2],
        ]
    )
    exact_phi = 0.3 * points[:, 0] - 0.2 * points[:, 1] + 0.1 * points[:, 2] + 2
    stress = np.array([[2.6, 0.18, 0], [0.18, -2.2, 0], [0, 0, 1.5]])
    electric_displacement = np.array([-0.5, 0.7, -0.5])
    tractions = {}
    surface_charges = {}
    for i, face in enumerate(mesh.triangle3):
        vertices = points[face]
        area_normal = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
        if np.dot(area_normal, vertices.mean(axis=0) - [0.5, 0.5, 0.5]) < 0:
            area_normal *= -1
        normal = area_normal / np.linalg.norm(area_normal)
        tractions[i] = stress @ normal
        surface_charges[i] = float(electric_displacement @ normal)
    supports = np.flatnonzero(points[:, 0] == 0)
    dofs = piezoelectric_dofs(mesh)
    result = solve_piezoelectric(
        mesh,
        material,
        [
            Dirichlet(dofs.indices("u", supports), exact_u[supports].ravel()),
            Dirichlet(dofs.indices("phi", [0]), exact_phi[[0]]),
        ],
        tractions=tractions,
        surface_charges=surface_charges,
    )
    assert_allclose(result.displacement, exact_u, atol=2e-12, rtol=0)
    assert_allclose(result.potential, exact_phi, atol=2e-12, rtol=0)
    assert_allclose(result.force_balance, 0, atol=2e-12, rtol=0)
    assert_allclose(result.charge_balance, 0, atol=2e-12, rtol=0)
