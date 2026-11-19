"""Prepared FSI lid examples: file-only runs, exports, repetition and inputs (#93).

Reference values are the coarsest #92 mesh pairs (Doc/93). The coupling stops
at max|r_d| <= 1e-8 Lref = 1e-9 m, so an independent repetition, also via
the Gmsh/JSON round trip, must reproduce the lid displacement within 2e-9 m.
Regenerated JSON inputs agree to 1e-12 (the SA BVP is libm dependent).
Repeating the identical files on the same platform must reproduce the
history exactly, because every solver step is a deterministic direct solve.
"""

import csv
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from tests.reference.fsi_lid import write_example

ROOT = Path(__file__).parents[2]
EXAMPLE = ROOT / "examples/fsi_lid"
REFERENCE_LID_M = {"laminar": 5.71964e-4, "sa": 9.64275e-4}


def run(case, output, cwd):
    command = [sys.executable, str(EXAMPLE / "run.py"), case, "--output", str(output)]
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    return subprocess.run(command, cwd=cwd, env=env, capture_output=True, text=True)


def lid_displacement(output):
    fields = json.loads((output / "structure_fields.json").read_text())
    points = np.asarray(fields["reference_points_m"])
    displacement = np.asarray(fields["structure_displacement_m"])
    return displacement[np.isclose(points[:, 1], 0.1, atol=1e-15), 1]


@pytest.mark.parametrize("case", ["laminar", "sa"])
def test_prepared_lid_runs_from_files_and_exports_units(case, tmp_path):
    completed = run(case, tmp_path / "result", tmp_path)
    assert completed.returncode == 0, completed.stderr
    output = tmp_path / "result"
    report = json.loads((output / "report.json").read_text())
    assert report["status"] == "converged"
    assert report["fsi_schema_version"] == 1
    assert set(report["input_files_sha256"]) >= {
        "config",
        "fluid.domain",
        "structure.model",
        "structure.mesh",
    }
    for name, digest in report["input_files_sha256"].items():
        assert len(digest) == 64, name
    assert set(report["versions"]) == {"fem_toolbox", "numpy", "scipy", "python"}
    assert report["settings"]["fluid_model"] == (
        "laminar" if case == "laminar" else "sa-noft2"
    )
    assert (
        report["coupling"]["displacement_residual_m"]
        <= report["settings"]["displacement_gate_m"]
    )
    assert report["coupling"]["load_change_N"] <= report["settings"]["load_gate_N"]
    assert report["coupling"]["transfer_relative_error"] <= 1e-10
    assert report["structure"]["relative_residual"] <= 1e-10
    assert len(report["fluid"]["provenance"]["input_sha256"]) == 64
    lid = lid_displacement(output)
    assert np.ptp(lid) <= 1e-10
    assert lid.mean() == pytest.approx(REFERENCE_LID_M[case], abs=2e-9)
    with (output / "history.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == report["iterations"]
    assert "displacement_residual_m" in rows[0] and "load_change_N" in rows[0]
    fluid = json.loads((output / "fluid_fields.json").read_text())
    cells = len(fluid["quad4"])
    assert np.asarray(fluid["velocity_m_s"]).shape == (cells, 2)
    assert np.asarray(fluid["pressure_Pa"]).shape == (cells,)
    if case == "sa":
        assert np.asarray(fluid["nu_tilde_m2_s"]).shape == (cells,)
        assert max(fluid["wall_y_plus"]) <= 1
        assert report["sa_noft2_constants"]["kappa"] == 0.41
        assert report["fluid"]["maximum_y_plus"] <= 1
    # Existing output is refused before any solve and stays untouched.
    original = (output / "report.json").read_bytes()
    again = run(case, output, tmp_path)
    assert again.returncode != 0 and "already exists" in again.stderr
    assert (output / "report.json").read_bytes() == original


def test_identical_repetition_reproduces_history_exactly(tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    assert run("laminar", first, tmp_path).returncode == 0
    assert run("laminar", second, tmp_path).returncode == 0
    assert (first / "history.csv").read_text() == (second / "history.csv").read_text()
    a = json.loads((first / "report.json").read_text())
    b = json.loads((second / "report.json").read_text())
    assert a["input_sha256"] == b["input_sha256"]
    assert a["mesh"] == b["mesh"]


@pytest.mark.parametrize("case", ["laminar", "sa"])
def test_versioned_inputs_are_reproducible_from_the_reference_generator(case, tmp_path):
    generated = write_example(tmp_path / case, case).parent
    names = sorted(p.name for p in (EXAMPLE / case).iterdir())
    assert names == sorted(p.name for p in generated.iterdir())
    for name in names:
        produced, versioned = generated / name, EXAMPLE / case / name
        if name.endswith(".json"):
            # The SA inlet profile and start field come from a BVP solve whose
            # last bit depends on the platform libm; compare numbers to 1e-12.
            assert_same_json(
                json.loads(produced.read_text()), json.loads(versioned.read_text())
            )
        else:
            assert produced.read_bytes() == versioned.read_bytes(), name


def assert_same_json(produced, versioned, path="$"):
    if isinstance(versioned, dict):
        assert isinstance(produced, dict) and produced.keys() == versioned.keys(), path
        for key in versioned:
            assert_same_json(produced[key], versioned[key], f"{path}.{key}")
    elif isinstance(versioned, list):
        assert isinstance(produced, list) and len(produced) == len(versioned), path
        for index, (a, b) in enumerate(zip(produced, versioned, strict=True)):
            assert_same_json(a, b, f"{path}[{index}]")
    elif isinstance(versioned, float):
        assert isinstance(produced, float), path
        assert produced == pytest.approx(versioned, rel=1e-12, abs=1e-300), path
    else:
        assert type(produced) is type(versioned) and produced == versioned, path
