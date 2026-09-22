"""Manufactured T-to-u solution with independently differentiated source terms."""

from itertools import product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs
from fem_toolbox.solvers.heat import thermal_dofs
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics


def test_coupled_temperature_displacement_stress_convergence(cube_mesh):
    # T = 300 + x(1-x) K, u = (0.01 x², 0, 0) m.
    # k=2 W/(m K) gives Q=4 W/m³. For E=100 Pa, nu=1/4,
    # C11=120 Pa; alpha_x=0.002/K gives
    # f_x = -C11 [0.02 - 0.002(1-2x)] N/m³.
    conductivity = ThermalMaterial.isotropic(2)
    elasticity = ElasticMaterial.isotropic(100, 0.25)
    expansion = ThermalExpansion([0.002, 0, 0, 0, 0, 0], 300)

    def exact_temperature(x):
        return 300 + x[..., 0] * (1 - x[..., 0])

    def exact_displacement(x):
        u = np.zeros_like(x)
        u[..., 0] = 0.01 * x[..., 0] ** 2
        return u

    def exact_stress(x):
        mechanical_xx = 0.02 * x[..., 0] - 0.002 * x[..., 0] * (1 - x[..., 0])
        sigma = np.zeros((*x.shape[:-1], 6))
        sigma[..., :3] = mechanical_xx[..., None] * [120, 40, 40]
        return sigma

    def body_force(x):
        return [-120 * (0.02 - 0.002 * (1 - 2 * x[0])), 0, 0]

    # Fourth-order tensor-product Gauss rule under a Duffy map integrates the
    # squared polynomial errors, including the map Jacobian, exactly.
    roots, weights = np.polynomial.legendre.leggauss(4)
    roots, weights = (roots + 1) / 2, weights / 2
    barycentric, quadrature = [], []
    for i, j, k in product(range(4), repeat=3):
        r, s, t = roots[[i, j, k]]
        a, b, c = r, (1 - r) * s, (1 - r) * (1 - s) * t
        barycentric.append([1 - a - b - c, a, b, c])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - r) ** 2 * (1 - s))
    barycentric, quadrature = np.array(barycentric), np.array(quadrature)

    levels = (2, 4, 6, 8)
    errors = []
    for n in levels:
        mesh = cube_mesh(n)
        boundary = np.flatnonzero(
            np.any((mesh.points == 0) | (mesh.points == 1), axis=1)
        )
        heat_bcs = [
            Dirichlet(
                thermal_dofs(mesh).indices("T", boundary),
                exact_temperature(mesh.points[boundary]),
            )
        ]
        elastic_bcs = [
            Dirichlet(
                mechanical_dofs(mesh).indices("u", boundary),
                exact_displacement(mesh.points[boundary]).ravel(),
            )
        ]
        result = solve_thermomechanics(
            mesh,
            conductivity,
            elasticity,
            expansion,
            heat_bcs,
            elastic_bcs,
            source=4,
            body_force=body_force,
        )
        t_l2, u_l2, stress_l2 = 0.0, 0.0, 0.0
        for cell_id, nodes in enumerate(mesh.tetra4):
            geometry = Tet4Geometry(mesh.points[nodes])
            xyz = barycentric @ mesh.points[nodes]
            dt = barycentric @ result.heat.temperature[nodes] - exact_temperature(xyz)
            du = barycentric @ result.elasticity.displacement[
                nodes
            ] - exact_displacement(xyz)
            ds = result.elasticity.stress[cell_id] - exact_stress(xyz)
            measure = 6 * geometry.volume * quadrature
            t_l2 += np.dot(measure, dt**2)
            u_l2 += np.dot(measure, np.sum(du**2, axis=1))
            stress_l2 += np.dot(measure, np.sum(ds**2, axis=1))
        errors.append(np.sqrt([t_l2, u_l2, stress_l2]))
        assert abs(result.heat.power_balance) < 1e-9
        assert abs(result.heat.dirichlet_inflow + 4) < 1e-9
        assert np.linalg.norm(result.elasticity.force_balance) < 1e-9
        assert np.linalg.norm(result.elasticity.moment_balance) < 1e-9
        assert result.heat.relative_residual < 1e-12
        assert result.elasticity.relative_residual < 1e-12
    errors = np.array(errors)
    rates = (
        np.log(errors[:-1] / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("levels, L2(T,u,sigma):", list(zip(levels, errors.tolist(), strict=True)))
    print("rates:", rates.tolist())
    assert np.all((rates[-2:, :2] > 1.8) & (rates[-2:, :2] < 2.2))
    assert np.all((rates[-2:, 2] > 0.85) & (rates[-2:, 2] < 1.15))
