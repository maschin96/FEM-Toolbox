"""RC/RL/RLC dynamics, state consistency and discrete energy identity (#117)."""

import numpy as np
import pytest

from fem_toolbox.circuits import Circuit
from fem_toolbox.errors import ValidationError
from fem_toolbox.time import Step, TimeSeries


def solve(text, times, **kwargs):
    from fem_toolbox.circuits import solve_transient

    return solve_transient(Circuit.from_netlist(text), times, **kwargs)


def test_rc_step_analytic_charge_and_energy():
    times = np.linspace(0, 2, 201)
    r = solve(
        "V1 in 0 0\nR1 in out 2\nC1 out 0 .5", times, sources={"V1": Step(0, 0, 1)}
    )
    # tau=R*C=1 s. BE is O(dt): beforehand fixed 0.002 V error at dt=.01 s.
    assert np.max(abs(r.node_potentials_V["out"] - (1 - np.exp(-times)))) < 0.002
    np.testing.assert_allclose(
        np.diff(r.capacitor_charge_C["C1"]),
        r.branch_currents_A["C1"][1:] * np.diff(times),
        rtol=1e-12,
        atol=1e-14,
    )
    assert np.max(abs(r.energy_balance_J)) < 1e-12
    assert np.max(abs(r.power_balance_W)) < 1e-12
    assert np.max(abs(r.kcl_residual_A["out"])) < 1e-12
    assert np.all(r.numerical_dissipation_J >= 0)
    with pytest.raises(ValueError):
        r.times[0] = 99
    with pytest.raises(ValueError):
        r.node_potentials_V["out"][0] = 99


def test_rl_step_analytic_flux():
    times = np.linspace(0, 2, 201)
    r = solve(
        "V1 in 0 0\nR1 in out 2\nL1 out 0 2", times, sources={"V1": Step(0, 0, 1)}
    )
    assert np.max(abs(r.branch_currents_A["L1"] - 0.5 * (1 - np.exp(-times)))) < 0.001
    np.testing.assert_allclose(
        np.diff(r.inductor_flux_Wb["L1"]),
        r.branch_voltages_V["L1"][1:] * np.diff(times),
        rtol=1e-12,
        atol=1e-14,
    )
    assert np.max(abs(r.energy_balance_J)) < 1e-12


def test_explicit_capacitor_and_inductor_initial_states():
    times = np.linspace(0, 1, 101)
    c = solve("R1 a 0 1\nC1 a 0 1", times, capacitor_voltages_V={"c1": 2})
    inductive = solve("R1 a 0 1\nL1 a 0 1", times, inductor_currents_A={"L1": 1})
    assert np.max(abs(c.node_potentials_V["a"] - 2 * np.exp(-times))) < 0.004
    assert np.max(abs(inductive.branch_currents_A["L1"] - np.exp(-times))) < 0.002
    assert c.stored_energy_J[0] == pytest.approx(2, rel=0, abs=1e-12)
    assert inductive.stored_energy_J[0] == pytest.approx(0.5, rel=0, abs=1e-12)
    assert inductive.node_potentials_V["a"][0] == pytest.approx(-1, rel=0, abs=1e-12)


def test_passive_rlc_energy_and_first_order_convergence():
    errors = []
    for steps in (100, 200, 400):
        times = np.linspace(0, 3, steps + 1)
        r = solve("C1 a 0 1\nR1 a b 1\nL1 b 0 1", times, capacitor_voltages_V={"C1": 1})
        omega = np.sqrt(3) / 2
        reference = np.exp(-times / 2) * (
            np.cos(omega * times) + np.sin(omega * times) / (2 * omega)
        )
        errors.append(np.max(abs(r.node_potentials_V["a"] - reference)))
        assert np.all(np.diff(r.stored_energy_J) <= 1e-12)
        assert np.max(abs(r.energy_balance_J)) < 1e-12
        np.testing.assert_allclose(
            np.diff(r.stored_energy_J),
            -r.resistor_loss_J[1:] - r.numerical_dissipation_J[1:],
            rtol=1e-11,
            atol=1e-12,
        )
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 0.9)
    assert errors[-1] < 0.004


def test_interior_event_is_inserted_without_pre_event_forcing():
    r = solve(
        "V1 in 0 0\nR1 in out 1\nC1 out 0 1",
        [0, 0.2, 0.4, 0.6],
        sources={"V1": Step(0.35, 0, 1)},
    )
    np.testing.assert_allclose(r.times, [0, 0.2, 0.35, 0.4, 0.6], rtol=0, atol=1e-15)
    assert r.node_potentials_V["out"][2] == pytest.approx(0, rel=0, abs=1e-14)
    assert r.node_potentials_V["out"][3] == pytest.approx(0.05 / 1.05, rel=0, abs=1e-12)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"capacitor_voltages_V": {"R1": 1}},
        {"capacitor_voltages_V": {"C1": float("nan")}},
        {"inductor_currents_A": {"missing": 0}},
        {"sources": {"R1": 1}},
        {"sources": {"V1": 1j}},
        {"sources": {"V1": TimeSeries([0, 1], [0, 1], unit="K")}},
        {"max_steps": 1},
        {"max_history_bytes": 1},
    ],
)
def test_invalid_initial_values_profiles_and_limits(kwargs):
    with pytest.raises(ValidationError):
        solve("V1 in 0 0\nR1 in out 1\nC1 out 0 1", [0, 0.5, 1], **kwargs)


def test_unsupported_initial_constraint_loop_is_diagnosed():
    with pytest.raises(ValidationError, match="initial"):
        solve("V1 a 0 1\nC1 a 0 1", [0, 1], capacitor_voltages_V={"C1": 0})


@pytest.mark.parametrize("times", [[0], [0, 0], [1, 0], [0, float("nan")]])
def test_invalid_times(times):
    with pytest.raises(ValidationError):
        solve("R1 a 0 1", times)


def test_deterministic_provenance_tracks_initial_state_and_profiles():
    text = "R1 a 0 1\nC1 a 0 1"
    a = solve(text, [0, 0.1], capacitor_voltages_V={"C1": 1})
    b = solve(text, [0, 0.1], capacitor_voltages_V={"c1": 1})
    c = solve(text, [0, 0.1], capacitor_voltages_V={"C1": 2})
    assert a.provenance.to_dict() == b.provenance.to_dict()
    assert a.provenance.input_sha256 != c.provenance.input_sha256
