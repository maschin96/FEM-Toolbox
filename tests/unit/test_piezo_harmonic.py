"""Complex oscillator, static limit, resonance and time-domain equivalence."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
from fem_toolbox.solvers.piezo_transient import solve_transient_piezoelectric


def problem():
    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], tetra4=[[0, 1, 2, 3]]
    )
    e = np.zeros((3, 6))
    e[0, 0] = 2
    mat = PiezoelectricMaterial(np.diag([10, 20, 30, 4, 5, 6]), e, np.diag([3, 4, 5]))
    fixed = np.setdiff1d(np.arange(16), [3, 13])
    return mesh, mat, [Dirichlet(fixed, np.zeros(len(fixed)))]


@pytest.mark.parametrize("hz", [0.0, 0.2, 1.0])
def test_complex_force_and_static_low_frequency_limit(hz):
    mesh, mat, bc = problem()
    force = np.zeros(16, dtype=complex)
    force[3] = 1 + 2j
    result = solve_harmonic_piezoelectric(mesh, mat, 6.0, bc, hz, loads=force)
    omega = 2 * np.pi * hz
    expected = (1 + 2j) / ((10 + 4 / 3) / 6 - 0.1 * omega**2)
    assert_allclose(result.displacement[1, 0], expected, rtol=3e-14)
    assert_allclose(result.potential[1], 2 / 3 * expected, rtol=3e-14)
    assert_allclose(result.velocity, 1j * omega * result.displacement, atol=1e-14)
    assert_allclose(result.acceleration, -(omega**2) * result.displacement, atol=1e-14)
    assert result.relative_residual_u < 1e-12
    assert result.relative_residual_phi < 1e-12
    assert result.provenance.solver_parameters["phasor"] == "Re(x_hat exp(+i omega t))"


def test_complex_voltage_and_moving_support_inertia():
    mesh, mat, _ = problem()
    fixed = np.setdiff1d(np.arange(16), [3])
    g = np.zeros(16, dtype=complex)
    g[13] = 2 - 3j
    # Moving y-axis vertex in x adds shear but no x normal strain contribution
    # to u1x stiffness (C is diagonal). Its consistent mass coupling is .05 kg.
    g[6] = 1 + 0.5j
    omega = 1.7
    result = solve_harmonic_piezoelectric(
        mesh,
        mat,
        6.0,
        [Dirichlet(fixed, np.zeros(len(fixed)))],
        omega / (2 * np.pi),
        prescribed=g,
    )
    expected = (-(2 / 6) * g[13] + omega**2 * 0.05 * g[6]) / (10 / 6 - 0.1 * omega**2)
    assert_allclose(result.displacement[1, 0], expected, rtol=3e-14)
    assert_allclose(result.potential[1], g[13], rtol=3e-14)


def test_resonance_rejected_even_with_zero_force():
    mesh, mat, bc = problem()
    hz = np.sqrt((10 + 4 / 3) / 0.6) / (2 * np.pi)
    with pytest.raises(ValidationError, match="resonan"):
        solve_harmonic_piezoelectric(mesh, mat, 6.0, bc, hz)


@pytest.mark.parametrize("hz", [-1, np.nan, 1j, True])
def test_invalid_frequency(hz):
    mesh, mat, bc = problem()
    with pytest.raises(ValidationError):
        solve_harmonic_piezoelectric(mesh, mat, 6.0, bc, hz)


def test_harmonic_and_time_domain_agree_off_resonance():
    mesh, mat, bc = problem()
    omega = 1.3
    force = np.zeros(16, dtype=complex)
    force[3] = 0.01 + 0.02j
    harmonic = solve_harmonic_piezoelectric(
        mesh, mat, 6.0, bc, omega / (2 * np.pi), loads=force
    )
    errors = []
    for n in [100, 200]:
        transient = solve_transient_piezoelectric(
            mesh,
            mat,
            6.0,
            bc,
            np.linspace(0, 1, n + 1),
            loads=lambda t: np.real(force * np.exp(1j * omega * t)),
            initial_displacement=harmonic.displacement.real,
            initial_velocity=harmonic.velocity.real,
        )
        expected = np.real(
            harmonic.displacement[1, 0] * np.exp(1j * omega * transient.times)
        )
        errors.append(np.max(np.abs(transient.displacement[:, 1, 0] - expected)))
    assert errors[1] < 1e-6
    assert 3.8 < errors[0] / errors[1] < 4.2


def test_unrepresentable_frequency_is_validation_error():
    mesh, mat, bc = problem()
    with pytest.raises(ValidationError, match="representable"):
        solve_harmonic_piezoelectric(mesh, mat, 6.0, bc, 1e200)


def test_unrepresentable_residual_metric_is_rejected():
    mesh, mat, bc = problem()
    force = np.zeros(16, dtype=complex)
    force[3] = 1e200 + 2e200j
    with np.errstate(over="ignore", invalid="ignore"), pytest.raises(ValidationError):
        solve_harmonic_piezoelectric(mesh, mat, 6.0, bc, 0.2, loads=force)
