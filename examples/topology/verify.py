"""Run the predeclared three-mesh/two-radius M11 verification matrix."""

import argparse
import hashlib
import json
import platform
import runpy
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy

from fem_toolbox import __version__
from fem_toolbox.topology_study import TopologyStudy, read_topology_result

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument(
    "--quick", action="store_true", help="smoke only; not the three-mesh gate"
)
args = parser.parse_args()
output = args.output.resolve()
output.mkdir(parents=True, exist_ok=False)
prepare = runpy.run_path(str(Path(__file__).parent / "prepare.py"))["make_inputs"]
rows = []
for n in [2] if args.quick else [2, 3, 4]:
    for radius in (0.3, 0.6):
        name = f"n{n}-r{radius}"
        config = prepare(output / (name + "-input"), n, radius)
        start = time.perf_counter()
        process = subprocess.run(
            [
                sys.executable,
                "-I",
                str(Path(__file__).parent / "run.py"),
                "--config",
                str(config),
                "--output",
                str(output / name),
            ],
            capture_output=True,
            text=True,
        )
        (output / (name + ".log")).write_text(process.stdout + process.stderr)
        record = read_topology_result(output / name)
        residual = max(s["relative_residual"] for s in record["states"])
        force = max(
            float(np.linalg.norm(s["force_balance_N"])) for s in record["states"]
        )
        work = max(
            abs(s["compliance_J"] - 2 * s["strain_energy_J"]) / s["compliance_J"]
            for s in record["states"]
        )
        passed = (
            process.returncode == 0
            and record["status"] == "converged"
            and record["compliance_J"] <= 0.99 * record["initial_compliance_J"]
            and record["volume_fraction"] <= 0.6 + 1e-8
            and residual <= 1e-9
            and force <= 1e-8
            and work <= 1e-8
        )
        row = {
            "name": name,
            "divisions": [2 * n, n, n],
            "filter_radius_m": radius,
            "cells": 12 * n**3,
            "seconds": time.perf_counter() - start,
            "compliance_J": record["compliance_J"],
            "initial_compliance_J": record["initial_compliance_J"],
            "volume_fraction": record["volume_fraction"],
            "gray_measure": record["gray_measure"],
            "intermediate_cell_volume_fraction": record[
                "intermediate_cell_volume_fraction"
            ],
            "maximum_relative_residual": residual,
            "maximum_force_balance_N": force,
            "maximum_relative_work_energy_error": work,
            "evaluations": record["evaluations"],
            "status": record["status"],
            "gate": "passed" if passed else "failed",
            "result_sha256": hashlib.sha256(
                (output / name / "result.json").read_bytes()
            ).hexdigest(),
        }
        rows.append(row)
        print(json.dumps(row), flush=True)
        (output / "verification.json").write_text(
            json.dumps(
                {
                    "schema": 1,
                    "complete_matrix": not args.quick,
                    "toolbox_version": __version__,
                    "python": platform.python_version(),
                    "numpy": np.__version__,
                    "scipy": scipy.__version__,
                    "rows": rows,
                },
                indent=2,
            )
            + "\n"
        )
# Repeat a complete study, using identical input bytes and effective SI data.
study = TopologyStudy.from_file(output / "n2-r0.3-input/study.toml")
repeat = study.run(output / "repeat")
first = read_topology_result(output / "n2-r0.3")
repeat_error = (
    abs(repeat.final.objective - first["compliance_J"]) / first["compliance_J"]
)
# Directional finite differences independently re-solve every state.
x = study.problem.initial_density(0.6)
exact = study.problem.evaluate(x)
direction = np.random.default_rng(81).normal(size=len(x))
direction[study.problem.solid] = 0
# Avoid stepping outside boxes at the prescribed cells.
direction[study.problem.void] = 0
direction /= np.linalg.norm(direction)
steps = []
for h in (1e-3, 1e-4, 1e-5):
    a, b = (
        study.problem.evaluate(x + h * direction),
        study.problem.evaluate(x - h * direction),
    )
    fd = (a.objective - b.objective) / (2 * h)
    analytical = float(exact.objective_gradient @ direction)
    steps.append(
        {
            "step": h,
            "central_difference_J": fd,
            "adjoint_J": analytical,
            "relative_error": abs(fd - analytical) / max(abs(analytical), 1e-15),
            "absolute_volume_derivative_error": abs(
                (a.volume_fraction - b.volume_fraction) / (2 * h)
                - exact.volume_gradient @ direction
            ),
        }
    )
report = json.loads((output / "verification.json").read_text())
report.update(
    repeat_relative_error=repeat_error,
    gradient_steps=steps,
    gates_passed=all(r["gate"] == "passed" for r in rows)
    and repeat_error <= 1e-10
    and steps[-1]["relative_error"] <= 2e-4
    and steps[-1]["absolute_volume_derivative_error"] <= 2e-9,
)
(output / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
raise SystemExit(0 if report["gates_passed"] else 2)
