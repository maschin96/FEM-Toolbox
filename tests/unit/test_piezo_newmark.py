"""Newmark verification against an independent scalar oscillator."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezo_transient import (
    KinematicValues,
    solve_transient_piezoelectric,
)


def problem():
    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], tetra4=[[0, 1, 2, 3]]
    )
    e = np.zeros((3, 6))
    e[0, 0] = 2
    mat = PiezoelectricMaterial(np.diag([10, 20, 30, 4, 5, 6]), e, np.diag([3, 4, 5]))
    fixed = np.setdiff1d(np.arange(16), [3, 13])
    return mesh, mat, [Dirichlet(fixed, np.zeros(len(fixed)))]


def test_free_oscillator_energy_and_second_order():
    mesh, mat, bc = problem()
    omega = np.sqrt((10 + 4 / 3) / 0.6)
    u0 = np.zeros((4, 3))
    u0[1, 0] = 0.01
    errors = []
    for n in [40, 80, 160]:
        result = solve_transient_piezoelectric(
            mesh, mat, 6.0, bc, np.linspace(0, 1, n + 1), initial_displacement=u0
        )
        errors.append(abs(result.displacement[-1, 1, 0] - 0.01 * np.cos(omega)))
        assert_allclose(result.total_energy, result.total_energy[0], rtol=3e-12)
        assert_allclose(
            result.potential[:, 1], 2 / 3 * result.displacement[:, 1, 0], atol=1e-14
        )
        assert_allclose(
            result.acceleration[:, 1, 0],
            -(omega**2) * result.displacement[:, 1, 0],
            atol=2e-13,
        )
    assert 3.8 < errors[0] / errors[1] < 4.2
    assert 3.8 < errors[1] / errors[2] < 4.2


def test_force_step_no_early_impulse():
    mesh, mat, bc = problem()

    def load(t):
        f = np.zeros(16)
        f[3] = float(t >= 0.5)
        return f

    result = solve_transient_piezoelectric(
        mesh, mat, 6.0, bc, [0, 1], loads=load, events=[0.5]
    )
    assert_allclose(result.times, [0, 0.5, 1])
    assert_allclose(result.displacement[1], 0, atol=1e-15)
    assert_allclose(result.velocity[1], 0, atol=1e-15)
    assert_allclose(result.acceleration[1, 1, 0], 10.0, rtol=1e-14)
    assert result.displacement[2, 1, 0] > 0


def test_moving_support_inertia_and_consistent_initial_state():
    mesh, mat, bc = problem()

    def prescribed(t):
        u = np.zeros((4, 3))
        v = u.copy()
        a = u.copy()
        u[0, 0] = 0.003 * t * t
        v[0, 0] = 0.006 * t
        a[0, 0] = 0.006
        return KinematicValues(u, v, a, np.zeros(4))

    result = solve_transient_piezoelectric(
        mesh, mat, 6.0, bc, np.linspace(0, 0.3, 61), prescribed=prescribed
    )
    # M10=rho*V/20=.05; M11=.1. At t0 all displacements vanish.
    assert_allclose(result.acceleration[0, 1, 0], -0.05 / 0.1 * 0.006, atol=1e-14)
    assert_allclose(result.displacement[:, 0, 0], 0.003 * result.times**2, atol=1e-14)
    assert_allclose(result.velocity[:, 0, 0], 0.006 * result.times, atol=1e-14)
    assert max(result.relative_residual_u) < 1e-10
    assert max(result.relative_residual_phi) < 1e-10


def test_inconsistent_initial_and_discontinuous_support_rejected():
    mesh, mat, bc = problem()
    bad = np.ones((4, 3))
    with pytest.raises(ValidationError, match="initial"):
        solve_transient_piezoelectric(
            mesh, mat, 6.0, bc, [0, 1], initial_displacement=bad
        )

    def prescribed(t):
        u = np.zeros((4, 3))
        u[0, 0] = float(t >= 0.5)
        return KinematicValues(u, np.zeros_like(u), np.zeros_like(u), np.zeros(4))

    with pytest.raises(ValidationError, match="discontinuous"):
        solve_transient_piezoelectric(
            mesh, mat, 6.0, bc, [0, 1], prescribed=prescribed, events=[0.5]
        )


def test_constant_force_equilibrium_and_signed_charge():
    mesh, mat, bc = problem()
    f = np.zeros(16)
    f[3] = 1
    f[13] = 0.1
    # [10/6, 2/6; 2/6, -3/6] [u,phi]=[1,.1].
    exact = np.linalg.solve([[10 / 6, 2 / 6], [2 / 6, -3 / 6]], [1, 0.1])
    u = np.zeros((4, 3))
    u[1, 0] = exact[0]
    result = solve_transient_piezoelectric(
        mesh, mat, 6.0, bc, [0, 0.03, 0.2, 1], loads=lambda t: f, initial_displacement=u
    )
    assert_allclose(result.displacement[:, 1, 0], exact[0], atol=5e-14)
    assert_allclose(result.potential[:, 1], exact[1], atol=5e-14)
    assert_allclose(result.velocity, 0, atol=1e-13)


def test_prescribed_voltage_step_keeps_mechanics_continuous():
    mesh, mat, _ = problem()
    fixed = np.setdiff1d(np.arange(16), [3])

    def prescribed(t):
        phi = np.zeros(4)
        phi[1] = float(t >= 0.5)
        return KinematicValues(
            np.zeros((4, 3)), np.zeros((4, 3)), np.zeros((4, 3)), phi
        )

    result = solve_transient_piezoelectric(
        mesh,
        mat,
        6.0,
        [Dirichlet(fixed, np.zeros(len(fixed)))],
        [0, 0.5, 1],
        prescribed=prescribed,
    )
    assert_allclose(result.displacement[1], 0, atol=1e-14)
    assert_allclose(result.velocity[1], 0, atol=1e-14)
    assert_allclose(result.potential[:, 1], [0, 1, 1], atol=1e-14)
    assert_allclose(result.acceleration[1, 1, 0], -(2 / 6) / 0.1, atol=1e-14)


@pytest.mark.parametrize("times", [[0, 0], [1, 0], [0, np.nan], [0]])
def test_invalid_grid(times):
    mesh, mat, bc = problem()
    with pytest.raises(ValidationError):
        solve_transient_piezoelectric(mesh, mat, 6.0, bc, times)


def test_quadratic_rigid_translation_manufactured_solution():
    mesh, mat, bc = problem()

    def prescribed(t):
        u = np.zeros((4, 3))
        v = u.copy()
        a = u.copy()
        u[[0, 2, 3], 0] = 0.003 * t * t
        v[[0, 2, 3], 0] = 0.006 * t
        a[[0, 2, 3], 0] = 0.006
        return KinematicValues(u, v, a, np.zeros(4))

    force = np.zeros(16)
    # Uniform body acceleration: rho*V/4 per vertex, all four x components.
    force[[0, 3, 6, 9]] = 0.25 * 0.006
    result = solve_transient_piezoelectric(
        mesh,
        mat,
        6.0,
        bc,
        np.linspace(0, 1, 41),
        loads=lambda t: force,
        prescribed=prescribed,
    )
    assert_allclose(result.displacement[:, 1, 0], 0.003 * result.times**2, atol=3e-15)
    assert_allclose(result.velocity[:, 1, 0], 0.006 * result.times, atol=3e-15)
    assert_allclose(result.acceleration[:, 1, 0], 0.006, atol=3e-14)
    assert_allclose(result.potential, 0, atol=3e-15)
