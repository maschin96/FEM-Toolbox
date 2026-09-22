"""Independent 1D FE cross-check and four-level coupled solution convergence."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from tests.reference.piezo_bar_1d import solve_bar


def material(c, e, k):
    stiffness = np.diag([20.0, 20.0, c, 4.0, 4.0, 4.0])
    coupling = np.zeros((3, 6))
    coupling[2, 2] = e
    return PiezoelectricMaterial(stiffness, coupling, np.diag([5.0, 5.0, k]))


def conditions(mesh, *, top_u=None, top_phi=None):
    nodes = np.arange(mesh.n_points)
    z = mesh.points[:, 2]
    bottom = np.flatnonzero(z == 0)
    top = np.flatnonzero(z == 1)
    dofs = piezoelectric_dofs(mesh)
    bcs = [
        Dirichlet(dofs.indices("u", nodes, ["x", "y"]), np.zeros(2 * len(nodes))),
        Dirichlet(dofs.indices("u", bottom, ["z"]), np.zeros(len(bottom))),
        Dirichlet(dofs.indices("phi", bottom), np.zeros(len(bottom))),
    ]
    if top_u is not None:
        bcs.append(Dirichlet(dofs.indices("u", top, ["z"]), np.full(len(top), top_u)))
    if top_phi is not None:
        bcs.append(Dirichlet(dofs.indices("phi", top), np.full(len(top), top_phi)))
    return bcs, top


def test_heterogeneous_bar_matches_separate_dense_1d_fe(cube_mesh):
    n = 4
    mesh = cube_mesh(n)
    lower = (12.0, 1.5, 2.5)
    upper = (18.0, 2.0, 3.0)
    layers = [lower] * 2 + [upper] * 2
    materials = [
        material(*(lower if mesh.points[cell, 2].mean() < 0.5 else upper))
        for cell in mesh.tetra4
    ]
    bcs, top = conditions(mesh, top_phi=0.8)
    traction_faces = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 2] == 1, axis=1))
    result = solve_piezoelectric(
        mesh,
        materials,
        bcs,
        tractions={int(face): [0, 0, 0.6] for face in traction_faces},
    )
    expected_u, expected_phi, expected_charge, expected_h = solve_bar(
        layers, force=0.6, voltage=0.8
    )
    z_level = np.rint(n * mesh.points[:, 2]).astype(int)
    assert_allclose(result.displacement[:, 2], expected_u[z_level], atol=2e-13, rtol=0)
    assert_allclose(result.potential, expected_phi[z_level], atol=2e-13, rtol=0)
    assert_allclose(
        result.electrode_charges[top].sum(), expected_charge, atol=2e-12, rtol=0
    )
    assert_allclose(result.total_enthalpy, expected_h, atol=2e-12, rtol=0)


def test_quadratic_coupled_solution_converges_on_four_meshes(cube_mesh):
    # u_z=a z(1-z), phi=p z(1-z), with b_z=2(C a+e p)
    # and rho=2(k p-e a). Side traction and side D-flux vanish.
    c, e, k = 10.0, 2.0, 3.0
    a, p = 0.2, 0.3
    levels = (2, 4, 6, 8)
    errors = []
    enthalpy_errors = []
    exact_enthalpy = (0.5 * c * a**2 + e * a * p - 0.5 * k * p**2) / 3
    roots, weights = np.polynomial.legendre.leggauss(3)
    for n in levels:
        mesh = cube_mesh(n)
        bcs, _ = conditions(mesh, top_u=0.0, top_phi=0.0)
        result = solve_piezoelectric(
            mesh,
            material(c, e, k),
            bcs,
            body_force=[0, 0, 2 * (c * a + e * p)],
            volume_charge=2 * (k * p - e * a),
        )
        z = mesh.points[:, 2]
        z_level = np.rint(n * z).astype(int)
        u = np.array(
            [result.displacement[z_level == i, 2].mean() for i in range(n + 1)]
        )
        phi = np.array([result.potential[z_level == i].mean() for i in range(n + 1)])
        assert_allclose(result.displacement[:, 2], u[z_level], atol=2e-12, rtol=0)
        assert_allclose(result.potential, phi[z_level], atol=2e-12, rtol=0)
        l2_sq = h1_sq = 0.0
        for i in range(n):
            h = 1 / n
            zs = (i + (roots + 1) / 2) * h
            shape = (roots + 1) / 2
            uh = (1 - shape) * u[i] + shape * u[i + 1]
            phih = (1 - shape) * phi[i] + shape * phi[i + 1]
            l2_sq += (
                h
                / 2
                * np.dot(
                    weights,
                    (uh - a * zs * (1 - zs)) ** 2 + (phih - p * zs * (1 - zs)) ** 2,
                )
            )
            du = (u[i + 1] - u[i]) / h - a * (1 - 2 * zs)
            dp = (phi[i + 1] - phi[i]) / h - p * (1 - 2 * zs)
            h1_sq += h / 2 * np.dot(weights, c * du**2 + k * dp**2)
        errors.append([np.sqrt(l2_sq), np.sqrt(h1_sq)])
        enthalpy_errors.append(abs(result.total_enthalpy - exact_enthalpy))
        assert result.relative_residual_u < 1e-12
        assert result.relative_residual_phi < 1e-12
        assert_allclose(result.force_balance, 0, atol=2e-11)
        assert abs(result.charge_balance) < 2e-12
    errors = np.asarray(errors)
    rates = (
        np.log(errors[:-1] / errors[1:])
        / np.log(np.asarray(levels[1:]) / levels[:-1])[:, None]
    )
    assert np.all((rates[-2:, 0] > 1.8) & (rates[-2:, 0] < 2.2))
    assert np.all((rates[-2:, 1] > 0.85) & (rates[-2:, 1] < 1.15))
    assert np.all(np.diff(enthalpy_errors) < 0)
