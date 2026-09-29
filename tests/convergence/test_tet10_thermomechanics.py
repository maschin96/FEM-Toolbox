"""Coupled quartic manufactured temperature/displacement; independent P2 interpolation and Duffy norms."""

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs
from fem_toolbox.solvers.heat import thermal_dofs
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics
from tests.reference.quadratic_mesh import elevate


def test_four_mesh_thermomechanical_fields_and_energy(cube_mesh):
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
        prescribed = np.zeros((len(fixed), 3))
        prescribed[:, 0] = 0.01 * (mesh.points[fixed] ** 4 @ coefficients)
        coupled = solve_thermomechanics(
            mesh,
            ThermalMaterial.isotropic(2),
            ElasticMaterial.isotropic(100, 0.25),
            ThermalExpansion([0.002, 0, 0, 0, 0, 0], 1),
            [
                Dirichlet(
                    thermal_dofs(mesh).indices("T", fixed),
                    1 + mesh.points[fixed] ** 4 @ coefficients,
                )
            ],
            [Dirichlet(mechanical_dofs(mesh).indices("u", fixed), prescribed.ravel())],
            source=lambda p: -24 * (p**2 @ coefficients),
            body_force=lambda p: [
                -0.01 * (1440 * p[0] ** 2 + 240 * p[1] ** 2 + 120 * p[2] ** 2)
                + 0.96 * p[0] ** 3,
                0.16 * p[1] ** 3,
                0.08 * p[2] ** 3,
            ],
        )
        result = coupled.elasticity
        l2, h1, t_l2, t_h1, stress_l2 = 0.0, 0.0, 0.0, 0.0, 0.0
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
            values = result.displacement[cell]
            gradient = np.einsum("qni,nj->qji", gradients, values)
            exact = np.zeros_like(xyz)
            exact[:, 0] = 0.01 * (xyz**4 @ coefficients)
            exact_gradient = np.zeros_like(gradient)
            exact_gradient[:, 0, :] = 0.04 * xyz**3 * coefficients
            l2 += determinant * (
                quadrature @ np.sum((shapes @ values - exact) ** 2, axis=1)
            )
            h1 += determinant * (
                quadrature @ np.sum((gradient - exact_gradient) ** 2, axis=(1, 2))
            )
            symmetric = (gradient + gradient.transpose(0, 2, 1)) / 2
            temperature = shapes @ coupled.heat.temperature[cell]
            thermal_gradient = np.einsum(
                "qni,n->qi", gradients, coupled.heat.temperature[cell]
            )
            exact_temperature = 1 + xyz**4 @ coefficients
            t_l2 += determinant * (quadrature @ (temperature - exact_temperature) ** 2)
            t_h1 += determinant * (
                quadrature
                @ np.sum((thermal_gradient - 4 * xyz**3 * coefficients) ** 2, axis=1)
            )
            mechanical = symmetric.copy()
            mechanical[:, 0, 0] -= 0.002 * (temperature - 1)
            exact_mechanical = (exact_gradient + exact_gradient.transpose(0, 2, 1)) / 2
            exact_mechanical[:, 0, 0] -= 0.002 * (exact_temperature - 1)
            error = mechanical - exact_mechanical
            stress_error = (
                40 * np.trace(error, axis1=1, axis2=2)[:, None, None] * np.eye(3)
                + 80 * error
            )
            stress_l2 += determinant * (
                quadrature @ np.sum(stress_error**2, axis=(1, 2))
            )
            density = 20 * np.trace(mechanical, axis1=1, axis2=2) ** 2 + 40 * np.sum(
                mechanical**2, axis=(1, 2)
            )
            energies.append(determinant * (quadrature @ density))
        assert_allclose(result.element_energy, energies, rtol=2e-12, atol=1e-16)
        assert_allclose(result.force_balance, 0, rtol=0, atol=2e-11)
        assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-11)
        assert abs(coupled.heat.power_balance) < 2e-11
        errors.append(np.sqrt([t_l2, t_h1, l2, h1, stress_l2]))
    rates = (
        np.log(np.array(errors[:-1]) / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print(
        "Tet10 coupled T_L2/T_H1/u_L2/u_H1/stress_L2 errors:", errors, "rates:", rates
    )
    assert np.all((rates[-2:, 0] > 2.7) & (rates[-2:, 0] < 3.3))
    assert np.all((rates[-2:, 1] > 1.8) & (rates[-2:, 1] < 2.2))
    assert np.all((rates[-2:, 2] > 2.7) & (rates[-2:, 2] < 3.3))
    assert np.all((rates[-2:, 3:] > 1.8) & (rates[-2:, 3:] < 2.2))
