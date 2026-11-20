"""Stationary nonlinear flow and explicit inlet/outlet acceptance tests."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError


def api():
    from fem_toolbox.cfd import (
        CartesianFluidGrid,
        ChannelBoundary,
        ClosedWalls,
        IncompressibleFluid,
        solve_stationary_flow,
    )

    return (
        CartesianFluidGrid,
        ChannelBoundary,
        ClosedWalls,
        IncompressibleFluid,
        solve_stationary_flow,
    )


def test_channel_inlet_outlet_and_mass():
    Grid, Channel, _, Fluid, solve = api()
    n = 16
    y = (np.arange(n) + 0.5) / n
    r = solve(Grid((32, n), (2, 1)), Fluid(1, 1), Channel(y * (1 - y), 5))
    assert r.converged and r.iterations > 0
    assert r.outlet_mass_flow_kg_s == pytest.approx(
        r.inlet_mass_flow_kg_s, rel=1e-9, abs=1e-12
    )
    assert np.max(np.abs(r.divergence_s_inv)) < 1e-9
    assert np.max(np.abs(r.u_m_s - y[:, None] * (1 - y[:, None]))) < 0.003
    assert r.pressure_Pa[0, 0] > r.pressure_Pa[0, -1] > 5
    assert abs(r.power_balance_W_m) < 1e-8
    assert np.max(np.abs(r.momentum_x_residual_N_m3)) < 1e-7
    with pytest.raises(ValueError):
        r.u_m_s.setflags(write=True)


def test_nonzero_convection_and_iteration_history():
    Grid, _, Walls, Fluid, solve = api()
    r = solve(Grid((10, 10), (1, 1)), Fluid(1, 1), Walls((0, 1, 0, 0)))
    assert r.iterations > 1
    assert len(r.residual_history) == r.iterations
    assert r.residual_history[-1] < r.residual_history[0]
    assert abs(r.convective_power_W_m) > 1e-5
    assert abs(r.power_balance_W_m) < 1e-8
    assert np.max(np.abs(r.divergence_s_inv)) < 1e-9
    assert max(abs(x) for x in r.global_momentum_balance_N_m) < 1e-8
    assert r.provenance.to_dict()["solver"]["parameters"]["convection"] == "donor-cell"


def test_gauge_and_zero_flow():
    Grid, _, Walls, Fluid, solve = api()
    r = solve(Grid((3, 3), (1, 1)), Fluid(1, 1), Walls(), pressure_reference_Pa=9)
    np.testing.assert_array_equal(r.pressure_Pa, np.full((3, 3), 9))
    np.testing.assert_array_equal(r.u_m_s, np.zeros((3, 4)))
    assert r.iterations == 1


def test_nonconvergence_reports_last_residual():
    Grid, _, Walls, Fluid, solve = api()
    with pytest.raises(ValidationError, match="converg") as error:
        solve(Grid((8, 8), (1, 1)), Fluid(1, 1), Walls((0, 1, 0, 0)), max_iterations=1)
    assert error.value.field == "iterations"


@pytest.mark.parametrize(
    "kind", ["grid", "fluid", "inlet", "gauge", "relaxation", "iterations", "unknowns"]
)
def test_invalid_configuration(kind):
    Grid, Channel, Walls, Fluid, solve = api()
    with pytest.raises(ValidationError):
        if kind == "grid":
            Grid((1, 3), (1, 1))
        elif kind == "fluid":
            Fluid(-1, 1)
        elif kind == "inlet":
            solve(Grid((3, 3), (1, 1)), Fluid(1, 1), Channel([1, 2]))
        elif kind == "gauge":
            solve(
                Grid((3, 3), (1, 1)),
                Fluid(1, 1),
                Channel([1, 1, 1]),
                pressure_reference_Pa=4,
            )
        elif kind == "relaxation":
            solve(Grid((3, 3), (1, 1)), Fluid(1, 1), Walls(), relaxation=0)
        elif kind == "iterations":
            solve(Grid((3, 3), (1, 1)), Fluid(1, 1), Walls(), max_iterations=True)
        else:
            solve(Grid((3, 3), (1, 1)), Fluid(1, 1), Walls(), max_unknowns=2)


def test_outlet_backflow_rejected():
    Grid, Channel, _, Fluid, solve = api()
    y = (np.arange(8) + 0.5) / 8
    fx = np.broadcast_to(100 * (y[:, None] - 0.5), (8, 8)).copy()
    with pytest.raises(ValidationError, match="backflow"):
        solve(
            Grid((8, 8), (1, 1)),
            Fluid(1, 1),
            Channel(np.ones(8) * 0.1),
            force_x_N_m3=fx,
        )


def test_closed_gauge_invariance_and_detached_provenance():
    Grid, _, Walls, Fluid, solve = api()
    options = (Grid((6, 6), (1, 1)), Fluid(1, 1), Walls((0, 1, 0, 0)))
    a = solve(*options)
    b = solve(*options, pressure_reference_Pa=100)
    np.testing.assert_allclose(a.u_m_s, b.u_m_s, rtol=0, atol=1e-12)
    np.testing.assert_allclose(b.pressure_Pa - a.pressure_Pa, 100, rtol=0, atol=1e-11)
    assert a.provenance.input_sha256 != b.provenance.input_sha256


@pytest.mark.parametrize(
    "options",
    [
        {"convection_scheme": "unknown"},
        {"convection_scheme": "central"},
        {"momentum_tolerance_N_m3": -1},
        {"force_x_N_m3": float("nan")},
    ],
)
def test_numerical_contract_errors(options):
    Grid, _, Walls, Fluid, solve = api()
    with pytest.raises(ValidationError):
        solve(Grid((4, 4), (1, 1)), Fluid(100, 1), Walls((0, 1, 0, 0)), **options)


def test_stokes_limit_matches_independent_solver():
    Grid, _, Walls, Fluid, solve = api()
    from fem_toolbox.cfd import solve_stokes

    flow = solve(Grid((8, 6), (1, 2)), Fluid(1e-9, 2), Walls((0, 1, 0, 0)))
    stokes = solve_stokes(
        cells=(8, 6),
        size_m=(1, 2),
        viscosity_Pa_s=2,
        tangential_wall_velocities_m_s=(0, 1, 0, 0),
    )
    np.testing.assert_allclose(flow.u_m_s, stokes.u_m_s, rtol=0, atol=1e-8)
    np.testing.assert_allclose(flow.pressure_Pa, stokes.pressure_Pa, rtol=0, atol=1e-7)


def test_pressure_outlet_shift_and_owned_inlet():
    Grid, Channel, _, Fluid, solve = api()
    profile = np.array([0.1, 0.2, 0.2, 0.1])
    boundary = Channel(profile, 3)
    profile[:] = 99
    a = solve(Grid((8, 4), (2, 1)), Fluid(1, 1), boundary)
    b = solve(Grid((8, 4), (2, 1)), Fluid(1, 1), Channel([0.1, 0.2, 0.2, 0.1], 8))
    np.testing.assert_array_equal(a.u_m_s[:, 0], [0.1, 0.2, 0.2, 0.1])
    np.testing.assert_allclose(a.u_m_s, b.u_m_s, rtol=0, atol=1e-12)
    np.testing.assert_allclose(b.pressure_Pa - a.pressure_Pa, 5, rtol=0, atol=1e-11)
