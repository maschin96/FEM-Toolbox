"""Independent quadratic manufactured electrostatic solution on four meshes."""

from itertools import product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics


def test_electric_l2_and_energy_convergence(cube_mesh):
    material = DielectricMaterial(
        1e-11 * np.array([[2, 0.3, 0.2], [0.3, 3, 0.1], [0.2, 0.1, 4]])
    )
    weights_xyz = np.array([1.0, 0.5, 0.25])

    # phi=x(1-x)+y(1-y)/2+z(1-z)/4; rho=9e-11 C/m³.
    def exact(x):
        return (x * (1 - x)) @ weights_xyz

    roots, weights = np.polynomial.legendre.leggauss(4)
    roots, weights = (roots + 1) / 2, weights / 2
    bary, quadrature = [], []
    for i, j, k in product(range(4), repeat=3):
        r, s, t = roots[[i, j, k]]
        x, y, z = r, (1 - r) * s, (1 - r) * (1 - s) * t
        bary.append([1 - x - y - z, x, y, z])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - r) ** 2 * (1 - s))
    bary, quadrature = np.array(bary), np.array(quadrature)
    errors = []
    levels = (2, 4, 6, 8)
    for n in levels:
        mesh = cube_mesh(n)
        boundary = np.flatnonzero(
            np.any((mesh.points == 0) | (mesh.points == 1), axis=1)
        )
        result = solve_electrostatics(
            mesh,
            material,
            [
                Dirichlet(
                    electric_dofs(mesh).indices("phi", boundary),
                    exact(mesh.points[boundary]),
                )
            ],
            volume_charge=9e-11,
        )
        l2, energy = 0.0, 0.0
        for cell_id, nodes in enumerate(mesh.tetra4):
            g = Tet4Geometry(mesh.points[nodes])
            xyz = bary @ mesh.points[nodes]
            error = bary @ result.potential[nodes] - exact(xyz)
            gradient_error = (
                -result.electric_field[cell_id] - (1 - 2 * xyz) * weights_xyz
            )
            l2 += 6 * g.volume * np.dot(quadrature, error**2)
            energy += (
                6
                * g.volume
                * np.dot(
                    quadrature,
                    np.einsum(
                        "ni,ij,nj->n",
                        gradient_error,
                        material.permittivity,
                        gradient_error,
                    ),
                )
            )
        errors.append([np.sqrt(l2), np.sqrt(energy)])
        assert abs(result.charge_balance) < 1e-20
        assert abs(result.dirichlet_outward_charge - 9e-11) < 1e-20
        assert result.relative_residual < 1e-12
    errors = np.array(errors)
    rates = (
        np.log(errors[:-1] / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("n, L2 error, energy error:", list(zip(levels, errors.tolist(), strict=True)))
    print("rates:", rates.tolist())
    assert np.all((rates[-2:, 0] > 1.8) & (rates[-2:, 0] < 2.2))
    assert np.all((rates[-2:, 1] > 0.85) & (rates[-2:, 1] < 1.15))
