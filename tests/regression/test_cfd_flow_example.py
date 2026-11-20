"""The shipped TOML reference must produce auditable result artifacts."""

import json
import os
import subprocess
import sys
from pathlib import Path


def test_cfd_example_runs_from_clean_directory(tmp_path):
    root = Path(__file__).resolve().parents[2]
    env = dict(os.environ, PYTHONPATH=str(root / "src"))
    output = tmp_path / "results"
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "examples/cfd_flow/run.py"),
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    report = json.loads((output / "report.json").read_text())
    assert report["converged"] is True
    assert abs(report["mass_balance_kg_s"]) < 1e-9
    assert report["max_velocity_error_m_s"] < 0.003
    assert len(report["provenance"]["input_sha256"]) == 64
    assert (output / "pressure.csv").read_text().startswith("x_m,y_m,pressure_Pa")
