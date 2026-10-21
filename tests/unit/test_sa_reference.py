"""The independent reference must solve its own wall/centre balance and scale."""

import numpy as np
import pytest

from tests.reference.sa_channel_1d import solve_channel_reference


def test_reference_branch_residual_wall_stress_and_initial_guess_independence():
    first = solve_channel_reference()
    other = solve_channel_reference(initial_amplitude=0.7, tolerance=2e-9)
    points = np.linspace(0, 2, 101)
    u, a = first.evaluate(points)
    u_other, a_other = other.evaluate(points)
    assert first.wall_shear_Pa == 1
    assert first.re_tau == pytest.approx(587.19)
    assert first.solution.rms_residuals.max() <= 1e-8
    assert u[0] == u[-1] == a[0] == a[-1] == 0
    assert np.min(a) >= 0
    np.testing.assert_allclose(u, u[::-1], atol=1e-12)
    np.testing.assert_allclose(a, a[::-1], atol=1e-12)
    np.testing.assert_allclose(u, u_other, atol=1e-7, rtol=0)
    np.testing.assert_allclose(a, a_other, atol=1e-9, rtol=0)
    assert first.bulk_velocity_m_s == pytest.approx(other.bulk_velocity_m_s, abs=1e-7)


def test_reference_dimensional_similarity_without_changing_model_constants():
    base = solve_channel_reference()
    # delta x0.05 and u_tau x0.01 keep Re_tau fixed when nu x0.0005.
    scaled = solve_channel_reference(
        half_height_m=0.05,
        pressure_drive_N_m3=0.002,
        kinematic_viscosity_m2_s=0.0005 / 587.19,
    )
    u, a = base.evaluate([0.1, 0.5, 1.0])
    us, as_ = scaled.evaluate(np.array([0.1, 0.5, 1.0]) * 0.05)
    np.testing.assert_allclose(us, 0.01 * u, rtol=1e-9, atol=1e-10)
    np.testing.assert_allclose(as_, 0.0005 * a, rtol=1e-9, atol=1e-12)
    assert scaled.re_tau == pytest.approx(base.re_tau)
    assert scaled.flow_per_span_m2_s == pytest.approx(
        0.0005 * base.flow_per_span_m2_s, rel=1e-9
    )


def test_fixed_flow_height_change_keeps_reynolds_and_scalar_profile():
    base = solve_channel_reference()
    factor = 1.01
    changed = solve_channel_reference(
        half_height_m=factor, pressure_drive_N_m3=1 / factor**3
    )
    y = np.linspace(0, 2, 101)
    u, a = base.evaluate(y)
    uc, ac = changed.evaluate(factor * y)
    np.testing.assert_allclose(uc, u / factor, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(ac, a, rtol=1e-8, atol=1e-10)
    assert changed.flow_per_span_m2_s == pytest.approx(
        base.flow_per_span_m2_s, rel=1e-8
    )
    assert changed.re_tau == pytest.approx(base.re_tau)
