"""Prepared file-only Quad-FV solve, balances and overwrite protection."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]


def test_prepared_quad_channel_solve(tmp_path):
    command = [
        sys.executable,
        str(ROOT / "examples/quad_flow/run.py"),
        "--output",
        str(tmp_path / "result"),
    ]
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    subprocess.run(command, cwd=tmp_path, env=env, check=True, capture_output=True)
    report = json.loads((tmp_path / "result/report.json").read_text())
    assert report["converged"]
    assert report["inlet_mass_flow_kg_s"] == pytest.approx(0.06875, abs=1e-14)
    assert abs(report["global_mass_balance_kg_s"]) < 2e-11
    assert max(map(abs, report["global_momentum_balance_N"])) < 2e-9
    assert abs(report["power_balance_W"]) < 2e-7
    assert len(report["provenance"]["solver"]["parameters"]["geometry_sha256"]) == 64
    fields = json.loads((tmp_path / "result/fields.json").read_text())
    assert len(fields["velocity_m_s"]) == 16
    assert len(fields["mass_residual_kg_m3_s"]) == 16
    assert len(fields["wall_traction_Pa"]) == 8
    history = (tmp_path / "result/residuals.csv").read_text().splitlines()
    assert len(history) == report["iterations"] + 1
    original = (tmp_path / "result/report.json").read_bytes()
    again = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True)
    assert again.returncode != 0
    assert (tmp_path / "result/report.json").read_bytes() == original
