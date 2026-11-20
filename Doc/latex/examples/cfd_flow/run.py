"""Run a complete declared channel problem and export auditable SI fields."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from fem_toolbox.cfd import FlowProblem

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=Path("cfd_flow_results"))
args = parser.parse_args()
problem = FlowProblem.from_file(Path(__file__).with_name("flow.toml"))
result = problem.solve()
nx, ny = problem.grid.cells
lx, ly = problem.grid.size_m
hx, hy = lx / nx, ly / ny
y = (np.arange(ny) + 0.5) * hy
error = float(np.max(np.abs(result.u_m_s - y[:, None] * (1 - y[:, None]))))
if error >= 0.003 or abs(result.global_mass_balance_kg_s) >= 1e-9:
    raise RuntimeError("channel reference gate failed")
args.output.mkdir(parents=True, exist_ok=False)
with (args.output / "pressure.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(
        ["x_m", "y_m", "pressure_Pa", "divergence_s_inv", "cell_mass_balance_kg_s"]
    )
    for j in range(ny):
        for i in range(nx):
            writer.writerow(
                [
                    (i + 0.5) * hx,
                    y[j],
                    result.pressure_Pa[j, i],
                    result.divergence_s_inv[j, i],
                    result.cell_mass_balance_kg_s[j, i],
                ]
            )
with (args.output / "velocity.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(["component", "x_m", "y_m", "velocity_m_s"])
    for j in range(ny):
        for i in range(nx + 1):
            writer.writerow(["u", i * hx, y[j], result.u_m_s[j, i]])
    for j in range(ny + 1):
        for i in range(nx):
            writer.writerow(["v", (i + 0.5) * hx, j * hy, result.v_m_s[j, i]])
report = dict(
    converged=result.converged,
    iterations=result.iterations,
    max_velocity_error_m_s=error,
    mass_balance_kg_s=result.global_mass_balance_kg_s,
    inlet_mass_flow_kg_s=result.inlet_mass_flow_kg_s,
    outlet_mass_flow_kg_s=result.outlet_mass_flow_kg_s,
    momentum_balance_N_m=result.global_momentum_balance_N_m,
    power_balance_W_m=result.power_balance_W_m,
    residual_history=result.residual_history,
    provenance=result.provenance.to_dict(),
)
(args.output / "report.json").write_text(
    json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
)
print(
    f"converged in {result.iterations} iterations; max velocity error {error:.6g} m/s"
)
