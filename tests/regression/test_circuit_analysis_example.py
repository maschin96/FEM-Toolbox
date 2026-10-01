"""Reproduce the complete TOML/CSV/JSON circuit recipe in a fresh output dir."""

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


def test_circuit_analysis_recipe(tmp_path):
    root = Path(__file__).parents[2]
    output = tmp_path / "results"
    subprocess.run(
        [
            sys.executable,
            str(root / "examples/circuit_analyses/run.py"),
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(root / "src")},
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads((output / "report.json").read_text())
    assert report["transient_max_error_V"] < 0.004
    assert report["ac_absorbed_power_R1_W"] == pytest.approx(0.5, rel=0, abs=1e-12)
    with (output / "transient.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 301
    assert max(abs(float(r["energy_balance_J"])) for r in rows) < 1e-12
