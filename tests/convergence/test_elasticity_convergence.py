"""Manufactured quadratic solution on four meshes; no sampled golden data."""

from itertools import product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ElasticMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity


def test_quadratic_solution_l2_and_energy_convergence(cube_mesh):
    material = ElasticMaterial.isotropic(100.0, 0.25)
    # u=(x²,0,0), epsilon=(2x,0,0,0,0,0), -div sigma=(-240,0,0).
    # Independent tensor-product Gauss rule under the Duffy map: degree-four
    # displacement error squared, including its Jacobian, is integrated exactly.
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
        prescribed = np.zeros((len(boundary), 3))
        prescribed[:, 0] = mesh.points[boundary, 0] ** 2
        dofs = mechanical_dofs(mesh)
        result = solve_elasticity(
            mesh,
            material,
            [Dirichlet(dofs.indices("u", boundary), prescribed.ravel())],
            body_force=[-240, 0, 0],
        )
        l2, energy = 0.0, 0.0
        for cell_id, nodes in enumerate(mesh.tetra4):
            g = Tet4Geometry(mesh.points[nodes])
            xyz = bary @ mesh.points[nodes]
            exact = np.zeros_like(xyz)
            exact[:, 0] = xyz[:, 0] ** 2
            difference = bary @ result.displacement[nodes] - exact
            exact_strain = np.zeros((len(bary), 6))
            exact_strain[:, 0] = 2 * xyz[:, 0]
            strain_error = result.strain[cell_id] - exact_strain
            l2 += 6 * g.volume * np.dot(quadrature, np.sum(difference**2, axis=1))
            energy += (
                6
                * g.volume
                * np.dot(
                    quadrature,
                    np.einsum(
                        "ni,ij,nj->n", strain_error, material.stiffness, strain_error
                    ),
                )
            )
        errors.append([np.sqrt(l2), np.sqrt(energy)])
        assert result.relative_residual < 1e-12
        assert np.linalg.norm(result.force_balance) < 1e-9
        assert np.linalg.norm(result.moment_balance) < 1e-9
    errors = np.array(errors)
    rates = (
        np.log(errors[:-1] / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("n, L2 error, energy error:", list(zip(levels, errors.tolist(), strict=True)))
    print("successive rates:", rates.tolist())
    # Asymptotic windows: include truncation effects but exclude wrong P1 order.
    assert np.all((rates[-2:, 0] > 1.8) & (rates[-2:, 0] < 2.2))
    assert np.all((rates[-2:, 1] > 0.85) & (rates[-2:, 1] < 1.15))
