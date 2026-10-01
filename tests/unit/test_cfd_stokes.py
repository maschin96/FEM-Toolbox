"""Pressure/velocity coupling checked independently of assembled matrices."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError


def solve(**kwargs):
    from fem_toolbox.cfd import solve_stokes

    return solve_stokes(**kwargs)


def test_hydrostatic_pressure_and_gauge():
    r = solve(
        cells=(6, 5),
        size_m=(2, 1),
        viscosity_Pa_s=3,
        force_x_N_m3=np.ones((5, 5)),
        force_y_N_m3=np.full((4, 6), 2),
        pressure_reference_Pa=7,
    )
    y, x = np.meshgrid(
        (np.arange(5) + 0.5) / 5, (np.arange(6) + 0.5) * 2 / 6, indexing="ij"
    )
    expected = 7 + x - x[0, 0] + 2 * (y - y[0, 0])
    np.testing.assert_allclose(r.pressure_Pa, expected, rtol=0, atol=1e-11)
    assert np.max(np.abs(r.u_m_s)) < 1e-12
    assert np.max(np.abs(r.v_m_s)) < 1e-12
    assert np.max(np.abs(r.divergence_s_inv)) < 1e-12
    assert r.converged
    with pytest.raises(ValueError):
        r.pressure_Pa.setflags(write=True)


def test_moving_lid_conservation_and_gauge_invariance():
    options = dict(
        cells=(8, 6),
        size_m=(1, 2),
        viscosity_Pa_s=2,
        tangential_wall_velocities_m_s=(0, 1, 0, 0),
    )
    r = solve(**options)
    shifted = solve(**options, pressure_reference_Pa=13)
    assert np.max(np.abs(r.u_m_s)) > 0.05
    assert np.max(np.abs(r.v_m_s)) > 0.01
    np.testing.assert_allclose(r.u_m_s, shifted.u_m_s, rtol=0, atol=1e-12)
    np.testing.assert_allclose(
        shifted.pressure_Pa - r.pressure_Pa, 13, rtol=0, atol=1e-10
    )
    assert np.max(np.abs(r.divergence_s_inv)) < 1e-11
    assert np.max(np.abs(r.momentum_x_residual_N_m3)) < 1e-9
    assert np.max(np.abs(r.momentum_y_residual_N_m3)) < 1e-9
    assert abs(r.power_balance_W_m) < 1e-10
    assert r.dissipation_W_m > 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"cells": (1, 3)},
        {"cells": (True, 3)},
        {"cells": (3.0, 3)},
        {"size_m": (0, 1)},
        {"size_m": (1, float("inf"))},
        {"force_x_N_m3": np.zeros((3, 3))},
        {"force_y_N_m3": float("nan")},
        {"viscosity_Pa_s": 0},
        {"pressure_reference_Pa": float("nan")},
        {"tangential_wall_velocities_m_s": (0, 1)},
        {"max_unknowns": 2},
        {"relative_tolerance": 0},
    ],
)
def test_invalid_input(kwargs):
    options = dict(cells=(3, 3), size_m=(1, 1), viscosity_Pa_s=1)
    options.update(kwargs)
    with pytest.raises(ValidationError):
        solve(**options)


def test_manufactured_four_meshes():
    errors = []
    pressure_errors = []
    pi = np.pi
    for n in (8, 16, 32, 64):
        yu, xu = np.meshgrid(
            (np.arange(n) + 0.5) / n, np.arange(1, n) / n, indexing="ij"
        )
        yv, xv = np.meshgrid(
            np.arange(1, n) / n, (np.arange(n) + 0.5) / n, indexing="ij"
        )
        fx = -4 * pi**3 * (2 * np.cos(2 * pi * xu) - 1) * np.sin(
            2 * pi * yu
        ) + 2 * pi * np.cos(2 * pi * xu) * np.sin(2 * pi * yu)
        fy = -4 * pi**3 * np.sin(2 * pi * xv) * (
            1 - 2 * np.cos(2 * pi * yv)
        ) + 2 * pi * np.sin(2 * pi * xv) * np.cos(2 * pi * yv)
        r = solve(
            cells=(n, n),
            size_m=(1, 1),
            viscosity_Pa_s=1,
            force_x_N_m3=fx,
            force_y_N_m3=fy,
        )
        u = pi * (1 - np.cos(2 * pi * xu)) * np.sin(2 * pi * yu)
        v = -pi * np.sin(2 * pi * xv) * (1 - np.cos(2 * pi * yv))
        errors.append(
            np.sqrt(
                np.mean((r.u_m_s[:, 1:-1] - u) ** 2)
                + np.mean((r.v_m_s[1:-1, :] - v) ** 2)
            )
        )
        yp, xp = np.meshgrid(
            (np.arange(n) + 0.5) / n, (np.arange(n) + 0.5) / n, indexing="ij"
        )
        p = np.sin(2 * pi * xp) * np.sin(2 * pi * yp)
        p -= p[0, 0]
        pressure_errors.append(np.sqrt(np.mean((r.pressure_Pa - p) ** 2)))
        assert np.max(np.abs(r.divergence_s_inv)) < 1e-9
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 1.8)
    assert np.all(np.log2(np.array(pressure_errors[:-1]) / pressure_errors[1:]) > 1.5)
    assert errors[-1] < 0.01
    assert pressure_errors[-1] < 0.01


def test_unresolvable_area_is_validation_error():
    with pytest.raises(ValidationError):
        solve(cells=(2, 2), size_m=(1e-200, 1e-200), viscosity_Pa_s=5e-324)


def test_zero_forcing_smallest_grid():
    r = solve(cells=(2, 2), size_m=(1, 1), viscosity_Pa_s=1, pressure_reference_Pa=-3)
    np.testing.assert_array_equal(r.pressure_Pa, np.full((2, 2), -3))
    np.testing.assert_array_equal(r.divergence_s_inv, np.zeros((2, 2)))
    assert r.power_balance_W_m == 0
