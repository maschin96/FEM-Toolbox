"""Quartic manufactured displacement; independent P2 interpolation and Duffy norms."""

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import ElasticMaterial, PiezoelectricMaterial
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from tests.reference.quadratic_mesh import elevate


def test_four_mesh_coupled_l2_h1_and_enthalpy(cube_mesh):
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
        coupling = np.zeros((3, 6))
        coupling[0, 0] = 2
        material = PiezoelectricMaterial(
            ElasticMaterial.isotropic(100, 0.25).stiffness, coupling, 3 * np.eye(3)
        )
        dofs = piezoelectric_dofs(mesh)
        # P=x^4+y^4/2+z^4/4, u=(.01P,0,0), phi=P.
        # sigma=C eps+e.T grad(phi), D=e eps-3 grad(phi).
        result = solve_piezoelectric(
            mesh,
            material,
            [
                Dirichlet(dofs.indices("u", fixed), prescribed.ravel()),
                Dirichlet(
                    dofs.indices("phi", fixed), mesh.points[fixed] ** 4 @ coefficients
                ),
            ],
            body_force=lambda p: [
                -38.4 * p[0] ** 2 - 2.4 * p[1] ** 2 - 1.2 * p[2] ** 2,
                0,
                0,
            ],
            volume_charge=lambda p: -35.76 * p[0] ** 2 - 18 * p[1] ** 2 - 9 * p[2] ** 2,
        )
        l2, h1, phi_l2, phi_h1 = 0.0, 0.0, 0.0, 0.0
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
            density = 20 * np.trace(gradient, axis1=1, axis2=2) ** 2 + 40 * np.sum(
                symmetric**2, axis=(1, 2)
            )
            phi = result.potential[cell]
            phi_gradient = np.einsum("qni,n->qi", gradients, phi)
            phi_l2 += determinant * (
                quadrature @ (shapes @ phi - xyz**4 @ coefficients) ** 2
            )
            phi_h1 += determinant * (
                quadrature
                @ np.sum((phi_gradient - 4 * xyz**3 * coefficients) ** 2, axis=1)
            )
            coupling_density = 2 * gradient[:, 0, 0] * phi_gradient[:, 0]
            electric_density = -1.5 * np.sum(phi_gradient**2, axis=1)
            energies.append(
                determinant
                * np.array(
                    [
                        quadrature @ density,
                        quadrature @ coupling_density,
                        quadrature @ electric_density,
                    ]
                )
            )
        assert_allclose(
            np.column_stack(
                [
                    result.mechanical_energy,
                    result.coupling_energy,
                    result.dielectric_enthalpy,
                ]
            ),
            energies,
            rtol=2e-12,
            atol=1e-16,
        )
        assert abs(result.charge_balance) < 2e-11
        assert_allclose(result.force_balance, 0, rtol=0, atol=2e-11)
        assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-11)
        errors.append(np.sqrt([l2, h1, phi_l2, phi_h1]))
    rates = (
        np.log(np.array(errors[:-1]) / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("Tet10 piezo u/phi L2/H1 errors:", errors, "rates:", rates)
    assert np.all((rates[-2:, [0, 2]] > 2.7) & (rates[-2:, [0, 2]] < 3.3))
    assert np.all((rates[-2:, [1, 3]] > 1.8) & (rates[-2:, [1, 3]] < 2.2))
