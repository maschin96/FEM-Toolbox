"""Independent RC/RL step response convergence with predeclared budgets (#119)."""

import numpy as np
import pytest

from fem_toolbox.circuits import Circuit, solve_transient
from fem_toolbox.time import Step


@pytest.mark.parametrize("kind", ["C", "L"])
def test_rc_rl_step_first_order_and_energy(kind):
    circuit = Circuit.from_netlist(f"V1 in 0 0\nR1 in out 1\n{kind}1 out 0 1")
    errors = []
    for steps in (40, 80, 160):
        times = np.linspace(0, 2, steps + 1)
        result = solve_transient(circuit, times, sources={"V1": Step(0, 0, 1)})
        actual = (
            result.node_potentials_V["out"]
            if kind == "C"
            else result.branch_currents_A["L1"]
        )
        errors.append(np.max(abs(actual - (1 - np.exp(-times)))))
        assert np.max(abs(result.energy_balance_J)) < 1e-12
        expected_energy = 0.5 * actual**2
        np.testing.assert_allclose(
            result.stored_energy_J, expected_energy, rtol=0, atol=1e-12
        )
    # A priori BE expectation: first order, with room for finite-step effects.
    orders = np.log2(np.array(errors[:-1]) / errors[1:])
    assert np.all((orders > 0.9) & (orders < 1.1))
    assert errors[-1] < 0.003  # V for RC, A for RL; dt=.0125 s, tau=1 s.
