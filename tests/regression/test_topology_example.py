"""Installed entry point and exported result remain usable from another cwd."""

import json
import os
import subprocess
import sys
from pathlib import Path


def test_topology_model_recipe_and_geometry_from_installed_package(tmp_path):
    root = Path(__file__).parents[2]
    result = subprocess.run(
        [
            sys.executable,
            "-I",
            str(root / "examples/topology/run.py"),
            "--output",
            str(tmp_path / "run"),
        ],
        cwd=tmp_path,
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1"},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads((tmp_path / "run/result.json").read_text())
    geometry = json.loads((tmp_path / "run/geometry/reanalysis.json").read_text())
    assert data["status"] == "converged"
    assert data["volume_fraction"] <= 0.6 + 1e-8
    assert geometry["gate_passed"]
    assert len(geometry["thresholds"]) == 6
