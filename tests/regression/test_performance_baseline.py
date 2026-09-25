"""Benchmark reports must carry real solves and independently checked fields."""

import json
import os
import subprocess
import sys
from pathlib import Path

import meshio
import numpy as np
import pytest

SCRIPT = Path(__file__).parents[2] / "benchmarks" / "profile_baseline.py"


@pytest.mark.parametrize("case", ["heat", "electrode", "transient"])
def test_baseline_worker(case, tmp_path):
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--worker",
            case,
            "--size",
            "2",
            "--steps",
            "3",
            "--output",
            str(tmp_path / "run"),
        ],
        check=True,
        env={**os.environ, "PYTHONPATH": "src"},
    )
    report = json.loads((tmp_path / "run" / "report.json").read_text())
    assert report["cells"] == 48
    assert report["nodes"] == 27
    assert report["dofs"] == (108 if case == "electrode" else 27)
    assert report["matrix_solves"] == (3 if case == "transient" else 1)
    assert report["factorized_systems"][0]["nnz"] > 0
    assert report["max_absolute_error"] < 1e-9
    assert report["instrumented_parity_error"] == 0
    assert report["verification"]["balance_error"] < 1e-9
    if case == "transient":
        manifest = json.loads(
            (tmp_path / "run" / "series" / "manifest.json").read_text()
        )
        assert len(manifest["frames"]) == 4
        assert manifest["status"] == "complete"
    assert report["traced_peak_bytes"] > 0
    assert report["process_peak_rss_bytes"] > 0
    for key in ["geometry", "assembly", "boundary", "factorization", "solve", "export"]:
        assert report["phases_s"][key] > 0
    assert report["environment"]["versions"]["scipy"]
    exported = meshio.read(tmp_path / "run" / "solution.vtu")
    x = exported.points[:, 0]
    field = exported.point_data["field"]
    expected = {
        "heat": 300 + x,
        "electrode": x / 20.25,
        "transient": np.full(len(x), 301.0),
    }[case]
    np.testing.assert_allclose(field, expected, rtol=0, atol=1e-9)


def test_invalid_benchmark_size(tmp_path):
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--size", "0", "--output", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "positive integer" in result.stderr
