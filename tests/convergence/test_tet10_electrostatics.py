"""Quartic manufactured potential; independent P2 interpolation and Duffy norms."""

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.solvers.electrostatics import solve_electrostatics
from tests.reference.quadratic_mesh import elevate


def test_four_mesh_electric_l2_h1_and_integrated_energy(cube_mesh):
    roots, weights = np.polynomial.legendre.leggauss(6)
    roots, weights = (roots + 1) / 2, weights / 2
    bary, quadrature = [], []
    for i, j, k in product(range(6), repeat=3):
        a, b, c = roots[[i, j, k]]
        x, y, z = a, (1 - a) * b, (1 - a) * (1 - b) * c
        bary.append([1 - x - y - z, x, y, z])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - a) ** 2 * (1 - b))
    bary, quadrature = np.array(bary), np.array(quadrature)
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    shapes = np.column_stack(
        [bary * (2 * bary - 1), *[4 * bary[:, i] * bary[:, j] for i, j in edges]]
    )
    coefficients = np.array([1, 0.5, 0.25])
    levels = (1, 2, 3, 4)
    errors = []
    for n in levels:
        mesh = elevate(cube_mesh(n))
        fixed = np.unique(mesh.facets)
        result = solve_electrostatics(
            mesh,
            DielectricMaterial.isotropic(2e-11),
            [Dirichlet(fixed, mesh.points[fixed] ** 4 @ coefficients)],
            volume_charge=lambda p: -24e-11 * (p**2 @ coefficients),
        )
        l2, h1 = 0.0, 0.0
        energies = []
        for cell in mesh.cells:
            corners = mesh.points[cell[:4]]
            jacobian = (corners[1:] - corners[0]).T
            determinant = abs(np.linalg.det(jacobian))
            xyz = bary @ corners
            bary_gradients = np.vstack([-np.ones(3), np.eye(3)]) @ np.linalg.inv(
                jacobian
            )
            gradients = np.concatenate(
                [
                    (4 * bary - 1)[:, :, None] * bary_gradients[None, :, :],
                    np.stack(
                        [
                            4
                            * (
                                bary[:, i, None] * bary_gradients[j]
                                + bary[:, j, None] * bary_gradients[i]
                            )
                            for i, j in edges
                        ],
                        axis=1,
                    ),
                ],
                axis=1,
            )
            values = result.potential[cell]
            gradient = np.einsum("qni,n->qi", gradients, values)
            l2 += determinant * (
                quadrature @ (shapes @ values - xyz**4 @ coefficients) ** 2
            )
            h1 += determinant * (
                quadrature @ np.sum((gradient - 4 * xyz**3 * coefficients) ** 2, axis=1)
            )
            energies.append(
                1e-11 * determinant * (quadrature @ np.sum(gradient**2, axis=1))
            )
        assert_allclose(result.element_energy, energies, rtol=2e-12, atol=1e-25)
        assert abs(result.charge_balance) < 2e-23
        assert abs(result.volume_charge + 14e-11) < 2e-23
        errors.append(np.sqrt([l2, h1]))
    rates = (
        np.log(np.array(errors[:-1]) / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("Tet10 electrostatics L2/H1 errors:", errors, "rates:", rates)
    assert np.all((rates[-2:, 0] > 2.7) & (rates[-2:, 0] < 3.3))
    assert np.all((rates[-2:, 1] > 1.8) & (rates[-2:, 1] < 2.2))
