"""Run the experimental DC divider and check independent SI reference values."""

from math import isclose
from pathlib import Path

from fem_toolbox.circuits import Circuit, solve_dc


def main() -> None:
    circuit = Circuit.from_netlist(
        Path(__file__).with_name("divider.cir").read_text(encoding="utf-8")
    )
    result = solve_dc(circuit)
    # O(10 V), mA and mW: tolerance leaves ample float64 roundoff allowance
    # while remaining many orders below the analytic reference signals.
    assert isclose(result.node_potentials_V["out"], 8, rel_tol=0, abs_tol=1e-12)
    assert isclose(result.branch_currents_A["R1"], 0.004, rel_tol=0, abs_tol=1e-14)
    assert isclose(result.branch_powers_W["V1"], -0.048, rel_tol=0, abs_tol=1e-13)
    assert max(abs(v) for v in result.kcl_residual_A.values()) < 1e-12
    assert abs(result.power_balance_W) < 1e-12
    print(f"U_out = {result.node_potentials_V['out']:.6f} V")
    print(f"I_R1 = {result.branch_currents_A['R1']:.6f} A")
    print(f"P_V1 = {result.branch_powers_W['V1']:.6f} W")
    print(f"Power balance = {result.power_balance_W:.3e} W")


if __name__ == "__main__":
    main()
