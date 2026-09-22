"""Installed-package style regression for the declarative reference example."""

import json
import os
import subprocess
import sys
from pathlib import Path

import meshio
import numpy as np
from numpy.testing import assert_allclose

ROOT = Path(__file__).parents[2]
SCRIPT = ROOT / "examples" / "thermopiezoelectric" / "run.py"
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


def test_declarative_example_exports_fields_and_provenance(tmp_path: Path) -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(tmp_path)],
        cwd=tmp_path,
        env=ENV,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    report = json.loads((tmp_path / "report.json").read_text(encoding="utf-8"))
    mesh = meshio.read(tmp_path / "fields.vtu")
    assert report["model"]["physics"] == "thermopiezoelectric"
    assert len(report["limitations"]) == 3
    assert report["provenance"]["schema_version"] == 1
    assert len(report["provenance"]["input_sha256"]) == 64
    assert_allclose(
        mesh.point_data["temperature_K"], [300, 310, 300, 300], atol=1e-12, rtol=0
    )
    # The signed electric row uses -rho: the free node has D_x=-rho/4=-1 C/m².
    assert_allclose(
        mesh.point_data["displacement_m"][1, 0], 0.025 - 1 / 22, atol=1e-12, rtol=0
    )
    assert_allclose(mesh.point_data["potential_V"][1], 1 / 4.4, atol=1e-12, rtol=0)
    assert_allclose(
        mesh.cell_data_dict["heat_flux_W_m2"]["tetra"][0],
        [-20, 0, 0],
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        mesh.cell_data_dict["stress_Pa"]["tetra"][0, 0], 0, atol=1e-12, rtol=0
    )
    assert_allclose(
        mesh.cell_data_dict["electric_displacement_C_m2"]["tetra"][0, 0],
        -1,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        report["metrics"]["free_displacement_x_m"],
        mesh.point_data["displacement_m"][1, 0],
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        report["metrics"]["free_potential_V"],
        mesh.point_data["potential_V"][1],
        atol=1e-12,
        rtol=0,
    )
    assert abs(report["balances"]["heat_power_W"]) < 1e-10
    assert np.linalg.norm(report["balances"]["force_N"]) < 1e-10
    assert np.linalg.norm(report["balances"]["moment_N_m"]) < 1e-10
    assert abs(report["balances"]["charge_C"]) < 1e-10
    assert max(report["relative_residuals"].values()) < 1e-10

    second = tmp_path / "repeated"
    repeated = subprocess.run(
        [sys.executable, str(SCRIPT), "--output", str(second)],
        cwd=tmp_path,
        env=ENV,
        capture_output=True,
        text=True,
        check=False,
    )
    assert repeated.returncode == 0, repeated.stderr
    assert (second / "report.json").read_bytes() == (
        tmp_path / "report.json"
    ).read_bytes()
