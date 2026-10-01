"""Complex peak phasors checked against series impedances and power (#118)."""

import numpy as np
import pytest

from fem_toolbox.circuits import Circuit
from fem_toolbox.errors import ValidationError


def solve(text, frequency, **kwargs):
    from fem_toolbox.circuits import solve_ac

    return solve_ac(Circuit.from_netlist(text), frequency, **kwargs)


@pytest.mark.parametrize("frequency", [0.1, 1, 100])
@pytest.mark.parametrize("kind,value", [("C", 0.01), ("L", 0.2)])
def test_series_rc_rl_impedance_and_phase(frequency, kind, value):
    r = solve(f"V1 in 0 2\nR1 in out 10\n{kind}1 out 0 {value}", frequency)
    omega = 2 * np.pi * frequency
    impedance = 1 / (1j * omega * value) if kind == "C" else 1j * omega * value
    current = 2 / (10 + impedance)
    assert r.branch_currents_A["R1"] == pytest.approx(current, rel=1e-12, abs=1e-14)
    assert r.node_potentials_V["out"] == pytest.approx(
        current * impedance, rel=1e-12, abs=1e-12
    )
    assert r.branch_powers_VA["R1"] == pytest.approx(
        0.5 * abs(current) ** 2 * 10, rel=1e-12
    )
    assert (r.branch_powers_VA[kind + "1"].imag > 0) == (kind == "L")
    assert abs(r.power_balance_VA) < 1e-12
    assert max(abs(v) for v in r.kcl_residual_A.values()) < 1e-12


def test_complex_sources_and_average_stored_energy():
    r = solve("V1 a 0 1\nC1 a 0 .02\nR1 a 0 10", 2, sources={"v1": 2 + 3j})
    assert r.node_potentials_V["a"] == pytest.approx(2 + 3j, rel=0, abs=1e-12)
    assert r.capacitor_charge_C["C1"] == pytest.approx(
        0.02 * (2 + 3j), rel=0, abs=1e-12
    )
    assert r.stored_energy_J == pytest.approx(0.25 * 0.02 * 13, rel=0, abs=1e-12)
    assert r.branch_powers_VA["C1"].real == pytest.approx(0, rel=0, abs=1e-12)
    with pytest.raises(TypeError):
        r.branch_currents_A["R1"] = 100


def test_series_rlc_multiple_frequencies_and_damped_resonance():
    # R=2, L=.5, C=.02, omega0=10 rad/s. At resonance impedance is 2 Ohm.
    for omega in (1, 9, 10, 11, 100):
        r = solve("V1 in 0 4\nR1 in a 2\nL1 a b .5\nC1 b 0 .02", omega / (2 * np.pi))
        expected = 4 / (2 + 1j * (omega * 0.5 - 1 / (omega * 0.02)))
        assert r.branch_currents_A["R1"] == pytest.approx(
            expected, rel=1e-12, abs=1e-13
        )
        assert r.inductor_flux_Wb["L1"] == pytest.approx(0.5 * expected, rel=1e-12)
        assert abs(r.power_balance_VA) < 1e-11
        expected_energy = (
            0.25 * 0.02 * abs(r.branch_voltages_V["C1"]) ** 2
            + 0.25 * 0.5 * abs(expected) ** 2
        )
        assert r.stored_energy_J == pytest.approx(expected_energy, rel=1e-12, abs=1e-12)
        assert r.branch_currents_A["C1"] == pytest.approx(
            1j * omega * r.capacitor_charge_C["C1"], rel=1e-12, abs=1e-12
        )


def test_zero_frequency_dc_values_and_power_convention():
    r = solve("V1 a 0 2\nR1 a 0 10\nC1 a 0 .1", 0)
    assert r.branch_currents_A["R1"] == pytest.approx(0.2, rel=0, abs=1e-12)
    assert r.branch_powers_VA["R1"] == pytest.approx(0.4, rel=0, abs=1e-12)
    assert r.stored_energy_J == pytest.approx(0.2, rel=0, abs=1e-12)
    with pytest.raises(ValidationError):
        solve("V1 a 0 1", 0, sources={"V1": 1j})


def test_lossless_resonance_even_with_zero_forcing():
    from fem_toolbox.circuits import CircuitTopologyError

    for amplitude in (0, 1):
        with pytest.raises(CircuitTopologyError, match="resonance"):
            solve(f"I1 0 a {amplitude}\nL1 a 0 1\nC1 a 0 1", 1 / (2 * np.pi))


@pytest.mark.parametrize("frequency", [-1, True, float("inf"), float("nan"), 1e308])
def test_bad_frequency(frequency):
    with pytest.raises(ValidationError):
        solve("R1 a 0 1", frequency)


@pytest.mark.parametrize(
    "sources",
    [
        {"R1": 1},
        {"unknown": 1},
        {"V1": complex(float("nan"))},
        {"V1": True},
        {"V1": 1, "v1": 2},
    ],
)
def test_invalid_source_override(sources):
    with pytest.raises(ValidationError):
        solve("V1 a 0 1\nR1 a 0 1", 1, sources=sources)


def test_ac_capacitive_ground_path_without_dc_path():
    r = solve("I1 0 a 1\nC1 a 0 1", 1 / (2 * np.pi))
    assert r.node_potentials_V["a"] == pytest.approx(-1j, rel=0, abs=1e-12)


def test_source_scalar_overflow_is_specific_validation_error():
    with pytest.raises(ValidationError):
        solve("V1 a 0 1", 1, sources={"V1": 10**400})
