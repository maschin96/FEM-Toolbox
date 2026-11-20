"""Four-mesh independent accuracy gates for issue #86."""

import numpy as np
import pytest

from fem_toolbox.cfd import (
    CartesianFluidGrid,
    ChannelBoundary,
    ClosedWalls,
    IncompressibleFluid,
    solve_stationary_flow,
)


def test_channel_four_meshes():
    errors = []
    gradient_errors = []
    for n in (8, 16, 32, 64):
        y = (np.arange(n) + 0.5) / n
        r = solve_stationary_flow(
            CartesianFluidGrid((2 * n, n), (2, 1)),
            IncompressibleFluid(1, 1),
            ChannelBoundary(y * (1 - y)),
        )
        errors.append(np.max(np.abs(r.u_m_s - y[:, None] * (1 - y[:, None]))))
        gradient = np.diff(np.mean(r.pressure_Pa, axis=0)) / (1 / n)
        gradient_errors.append(np.sqrt(np.mean((gradient + 2) ** 2)))
        assert abs(r.global_mass_balance_kg_s) < 1e-9 * r.inlet_mass_flow_kg_s
        assert np.max(np.abs(r.divergence_s_inv)) < 1e-9
    assert errors[-1] < 0.001
    assert gradient_errors[-1] < 0.02
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 0.8)


@pytest.mark.parametrize(
    "scheme,density,error_limit,order_limit",
    [
        ("central", 0.1, 0.01, 1.8),
        ("linear-upwind", 1.0, 0.1, 0.7),
    ],
)
def test_nonlinear_manufactured_four_meshes(scheme, density, error_limit, order_limit):
    pi = np.pi
    errors = []
    pressure_errors = []
    for n in (16, 32, 64, 128):
        yu, xu = np.meshgrid(
            (np.arange(n) + 0.5) / n, np.arange(1, n) / n, indexing="ij"
        )
        yv, xv = np.meshgrid(
            np.arange(1, n) / n, (np.arange(n) + 0.5) / n, indexing="ij"
        )
        u = pi * (1 - np.cos(2 * pi * xu)) * np.sin(2 * pi * yu)
        v_at_u = -pi * np.sin(2 * pi * xu) * (1 - np.cos(2 * pi * yu))
        v = -pi * np.sin(2 * pi * xv) * (1 - np.cos(2 * pi * yv))
        u_at_v = pi * (1 - np.cos(2 * pi * xv)) * np.sin(2 * pi * yv)
        convx = u * 2 * pi**2 * np.sin(2 * pi * xu) * np.sin(
            2 * pi * yu
        ) + v_at_u * 2 * pi**2 * (1 - np.cos(2 * pi * xu)) * np.cos(2 * pi * yu)
        convy = u_at_v * (
            -2 * pi**2 * np.cos(2 * pi * xv) * (1 - np.cos(2 * pi * yv))
        ) + v * (-2 * pi**2 * np.sin(2 * pi * xv) * np.sin(2 * pi * yv))
        fx = (
            -4 * pi**3 * (2 * np.cos(2 * pi * xu) - 1) * np.sin(2 * pi * yu)
            + 2 * pi * np.cos(2 * pi * xu) * np.sin(2 * pi * yu)
            + density * convx
        )
        fy = (
            -4 * pi**3 * np.sin(2 * pi * xv) * (1 - 2 * np.cos(2 * pi * yv))
            + 2 * pi * np.sin(2 * pi * xv) * np.cos(2 * pi * yv)
            + density * convy
        )
        r = solve_stationary_flow(
            CartesianFluidGrid((n, n), (1, 1)),
            IncompressibleFluid(density, 1),
            ClosedWalls(),
            force_x_N_m3=fx,
            force_y_N_m3=fy,
            convection_scheme=scheme,
        )
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
        assert np.max(np.abs(r.momentum_x_residual_N_m3)) < 1e-6
    assert errors[-1] < error_limit
    assert pressure_errors[-1] < error_limit
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > order_limit)
    assert np.all(
        np.log2(np.array(pressure_errors[:-1]) / pressure_errors[1:]) > order_limit
    )
