"""Actual FEM -> quantity table -> surrogate artifact end-to-end contract."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.surrogates import PolynomialSurrogate

ROOT = Path(__file__).parents[2]


def test_heat_surrogate_example(tmp_path):
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples/surrogate_heat/run.py"),
            "--output",
            str(tmp_path / "out"),
        ],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
    report = json.loads((tmp_path / "out/report.json").read_text())
    assert len(report["split"]["train_ids"]) == 6
    assert len(report["split"]["test_ids"]) == 2
    assert max(report["test"]["max_abs"]) < 1e-10
    model = PolynomialSurrogate.load(tmp_path / "out/surrogate.json")
    assert_allclose(
        model.predict([[3], [7]]).values,
        [[300 + 3 / 8, 3 / 6], [300 + 7 / 8, 7 / 6]],
        atol=1e-10,
        rtol=0,
    )
    assert np.isfinite(model.predict([5]).values).all()
