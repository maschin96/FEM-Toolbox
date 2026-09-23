"""Independent cooling and spatial manufactured-solution convergence."""

from itertools import combinations, product

import numpy as np
import pytest

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs
from fem_toolbox.solvers.transient_heat import (
    TemperatureBC,
    TransientConvection,
    solve_transient_heat,
)
from fem_toolbox.time import Ramp


def test_uniform_cooling_first_order_and_steady_limit():
    # Regular tetrahedron: uniform temperature is an exact FE eigenmode.
    points = np.array(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0.5, np.sqrt(3) / 2, 0],
            [0.5, np.sqrt(3) / 6, np.sqrt(2 / 3)],
        ]
    )
    mesh = Mesh(
        points=points, tetra4=[[0, 1, 2, 3]], triangle3=list(combinations(range(4), 3))
    )
    volume = np.sqrt(2) / 12
    area = np.sqrt(3) / 4
    h = volume / (4 * area)  # rho*c=1, exact decay rate h*A_total/(rho*c*V)=1/s
    mat = ThermalMaterial.isotropic(1, density=1, specific_heat=1)
    robin = {i: TransientConvection(h, 300.0) for i in range(4)}
    errors = []
    for steps in [10, 20, 40, 80]:
        result = solve_transient_heat(
            mesh,
            mat,
            np.linspace(0, 1, steps + 1),
            310.0,
            convection=robin,
            output_every=steps,
        )
        errors.append(abs(result.frames[-1].temperature[0] - (300 + 10 * np.exp(-1))))
        assert max(abs(b.balance) for b in result.balances) < 1e-10
    rates = np.log2(np.array(errors[:-1]) / errors[1:])
    print("cooling errors/rates:", errors, rates)
    assert np.all((rates > 0.9) & (rates < 1.1))
    final = solve_transient_heat(
        mesh, mat, np.linspace(0, 40, 81), 310.0, convection=robin, output_every=80
    )
    np.testing.assert_allclose(final.frames[-1].temperature, 300.0, rtol=0, atol=2e-10)


def test_affine_time_boundary_has_correct_capacity_cross_terms(cube_mesh):
    mesh = cube_mesh(2)
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    bcs = [
        TemperatureBC(
            [int(i)], Ramp(0, 1, 300 + mesh.points[i, 0], 302 + mesh.points[i, 0])
        )
        for i in boundary
    ]
    result = solve_transient_heat(
        mesh,
        ThermalMaterial.isotropic(2, density=3, specific_heat=4),
        [0, 0.2, 0.7, 1],
        300 + mesh.points[:, 0],
        conditions=bcs,
        source=24.0,
    )
    np.testing.assert_allclose(
        result.frames[-1].temperature, 302 + mesh.points[:, 0], atol=2e-11, rtol=0
    )
    for b in result.balances:
        assert b.storage_power == pytest.approx(24.0, abs=2e-9)
        assert b.dirichlet_inflow == pytest.approx(0.0, abs=2e-9)


def test_spatial_convergence_with_time_exact_linear_manufactured_solution(cube_mesh):
    # T=300+t+x(1-x)+.5y(1-y)+.25z(1-z), k=I, rho*c=1.
    # Q=1+3.5=4.5. Linear time dependence has zero backward-Euler truncation error.
    # Initial discrete elliptic projection avoids contaminating spatial rates
    # with a startup layer; reference norms are independent analytical integrals.
    def exact(x, t=0):
        return 300 + t + (x * (1 - x)) @ np.array([1.0, 0.5, 0.25])

    roots, weights = np.polynomial.legendre.leggauss(3)
    roots, weights = (roots + 1) / 2, weights / 2
    bary, quadrature = [], []
    for i, j, k in product(range(3), repeat=3):
        r, s, t = roots[[i, j, k]]
        x, y, z = r, (1 - r) * s, (1 - r) * (1 - s) * t
        bary.append([1 - x - y - z, x, y, z])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - r) ** 2 * (1 - s))
    bary, quadrature = np.array(bary), np.array(quadrature)
    errors = []
    levels = [2, 4, 6, 8]
    mat = ThermalMaterial.isotropic(1, density=1, specific_heat=1)
    for n in levels:
        mesh = cube_mesh(n)
        boundary = np.flatnonzero(
            np.any((mesh.points == 0) | (mesh.points == 1), axis=1)
        )
        initial = solve_heat(
            mesh,
            mat,
            [
                Dirichlet(
                    thermal_dofs(mesh).indices("T", boundary),
                    exact(mesh.points[boundary]),
                )
            ],
            source=3.5,
        ).temperature
        bcs = [
            TemperatureBC(
                [int(i)],
                Ramp(
                    0, 1, float(exact(mesh.points[i])), float(exact(mesh.points[i], 1))
                ),
            )
            for i in boundary
        ]
        result = solve_transient_heat(
            mesh,
            mat,
            [0, 0.1, 0.2],
            initial,
            conditions=bcs,
            source=4.5,
            output_every=2,
        )
        frame = result.frames[-1]
        l2, h1 = 0.0, 0.0
        for ci, nodes in enumerate(mesh.tetra4):
            g = Tet4Geometry(mesh.points[nodes])
            xyz = bary @ mesh.points[nodes]
            error = bary @ frame.temperature[nodes] - exact(xyz, 0.2)
            grad_error = frame.gradient[ci] - (1 - 2 * xyz) * np.array([1.0, 0.5, 0.25])
            l2 += 6 * g.volume * np.dot(quadrature, error**2)
            h1 += 6 * g.volume * np.dot(quadrature, np.sum(grad_error**2, axis=1))
        errors.append([np.sqrt(l2), np.sqrt(h1)])
    rates = (
        np.log(np.array(errors[:-1]) / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("spatial errors/rates:", errors, rates)
    assert np.all((rates[-2:, 0] > 1.8) & (rates[-2:, 0] < 2.2))
    assert np.all((rates[-2:, 1] > 0.85) & (rates[-2:, 1] < 1.15))
