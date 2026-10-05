"""Prepared file-only fluid import; independent O(1) SI preflight values."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]


def test_prepared_channel_import(tmp_path):
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    output = tmp_path / "preflight.json"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples/fluid_mesh/preflight.py"),
            "--output",
            str(output),
        ],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    data = json.loads(output.read_text())
    assert data["cells"] == 16
    assert data["boundary_faces"] == 16
    # Four midpoint samples of 4y(1-y), rho=1, b=0.1, H=1.
    assert data["inlet_mass_flow_kg_s"] == pytest.approx(0.06875, rel=0, abs=1e-14)
    assert data["minimum_wall_distance_m"] == pytest.approx(0.125, rel=0, abs=1e-14)
    assert data["maximum_layer_growth"] == pytest.approx(1, rel=0, abs=1e-14)
    assert data["wall_resolution"] == "not_checked_without_wall_shear"
    assert data["solver_status"] == "not_run"
