"""Analytic driven RLC transient/AC references with CSV and provenance output."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from fem_toolbox.circuits import Circuit, solve_ac, solve_transient
from fem_toolbox.time import Step


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("circuit_analysis_results"))
    output = parser.parse_args().output
    circuit = Circuit.from_file(Path(__file__).with_name("circuit.toml"))
    times = np.linspace(0, 3, 301)
    transient = solve_transient(circuit, times, sources={"V1": Step(0, 0, 1)})
    omega = np.sqrt(3) / 2
    expected = 1 - np.exp(-times / 2) * (
        np.cos(omega * times) + np.sin(omega * times) / (2 * omega)
    )
    error = float(np.max(abs(transient.node_potentials_V["out"] - expected)))
    # Predeclared BE budget: dt=.01 s, tau/order scales near 1 s, 0.004 V.
    assert error < 0.004
    assert np.max(abs(transient.energy_balance_J)) < 1e-12
    harmonic = solve_ac(circuit, 1 / (2 * np.pi), sources={"V1": 1})
    # At omega=1, X_L+X_C=0, Z=1 Ohm, I=1 A, V_C=-i V.
    assert abs(harmonic.branch_currents_A["R1"] - 1) < 1e-12
    assert abs(harmonic.node_potentials_V["out"] + 1j) < 1e-12
    assert abs(harmonic.power_balance_VA) < 1e-12
    output.mkdir(parents=True, exist_ok=False)
    with (output / "transient.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "time_s",
                "out_V",
                "R1_A",
                "stored_energy_J",
                "source_work_J",
                "resistor_loss_J",
                "numerical_dissipation_J",
                "energy_balance_J",
            ]
        )
        writer.writerows(
            zip(
                times,
                transient.node_potentials_V["out"],
                transient.branch_currents_A["R1"],
                transient.stored_energy_J,
                transient.source_work_J,
                transient.resistor_loss_J,
                transient.numerical_dissipation_J,
                transient.energy_balance_J,
                strict=True,
            )
        )
    report = {
        "example_schema_version": 1,
        "transient_max_error_V": error,
        "ac_current_R1_A": [
            harmonic.branch_currents_A["R1"].real,
            harmonic.branch_currents_A["R1"].imag,
        ],
        "ac_absorbed_power_R1_W": harmonic.branch_powers_VA["R1"].real,
        "transient_provenance": transient.provenance.to_dict(),
        "ac_provenance": harmonic.provenance.to_dict(),
    }
    (output / "report.json").write_text(
        json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"Transient error: {error:.6g} V; AC current: {harmonic.branch_currents_A['R1']:.6g} A"
    )


if __name__ == "__main__":
    main()
