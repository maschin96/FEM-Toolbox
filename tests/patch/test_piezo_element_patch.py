"""Coupled affine patch with independent boundary traction and charge integrals."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.piezoelectric import (
    piezoelectric_matrix,
    piezoelectric_response,
)
from fem_toolbox.materials import PiezoelectricMaterial


def test_affine_piezo_patch_on_distorted_mixed_orientation_mesh(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    c = np.diag([10.0, 20, 30, 4, 5, 6])
    e = np.zeros((3, 6))
    e[0, 0] = 2
    e[1, 1] = 1
    material = PiezoelectricMaterial(c, e, np.diag([3.0, 4, 5]))
    points = mesh.points
    displacement = np.column_stack(
        [
            0.2 * points[:, 0] + 0.03 * points[:, 1],
            -0.1 * points[:, 1],
            0.05 * points[:, 2],
        ]
    )
    potential = 0.3 * points[:, 0] - 0.2 * points[:, 1] + 0.1 * points[:, 2] + 2
    sigma = np.array([[2.6, 0.18, 0], [0.18, -2.2, 0], [0, 0, 1.5]])
    electric_displacement = np.array([-0.5, 0.7, -0.5])

    mechanical_residual = np.zeros((mesh.n_points, 3))
    electric_residual = np.zeros(mesh.n_points)
    for index, original_cell in enumerate(mesh.tetra4):
        cell = original_cell if index % 2 == 0 else original_cell[[0, 2, 1, 3]]
        geometry = Tet4Geometry(points[cell])
        local = np.concatenate([displacement[cell].ravel(), potential[cell]])
        residual = piezoelectric_matrix(geometry, material) @ local
        np.add.at(mechanical_residual, cell, residual[:12].reshape(4, 3))
        np.add.at(electric_residual, cell, residual[12:])
        response = piezoelectric_response(
            geometry, material, displacement[cell], potential[cell]
        )
        assert_allclose(
            response.stress,
            [2.6, -2.2, 1.5, 0, 0, 0.18],
            atol=3e-14,
            rtol=0,
        )
        assert_allclose(
            response.electric_displacement,
            electric_displacement,
            atol=3e-14,
            rtol=0,
        )

    expected_force = np.zeros_like(mechanical_residual)
    expected_charge = np.zeros_like(electric_residual)
    for face in mesh.triangle3:
        vertices = points[face]
        area_normal = 0.5 * np.cross(
            vertices[1] - vertices[0], vertices[2] - vertices[0]
        )
        if np.dot(area_normal, vertices.mean(axis=0) - [0.5, 0.5, 0.5]) < 0:
            area_normal *= -1
        expected_force[face] += (sigma @ area_normal) / 3
        expected_charge[face] += (electric_displacement @ area_normal) / 3
    assert_allclose(mechanical_residual, expected_force, atol=3e-14, rtol=0)
    assert_allclose(electric_residual, expected_charge, atol=3e-14, rtol=0)
    assert_allclose(mechanical_residual.sum(axis=0), 0, atol=3e-14, rtol=0)
    assert_allclose(electric_residual.sum(), 0, atol=3e-14, rtol=0)
