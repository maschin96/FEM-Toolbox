"""Independent norm and selection checks for the M7 accuracy/cost benchmark."""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from benchmarks.accuracy_cost import comparison, make_mesh, relative_errors


@pytest.mark.parametrize("order", [1, 2])
def test_error_norms_with_independent_affine_reference(order):
    mesh = make_mesh(2, order)

    def exact(p):
        return p[:, 0, None]

    def gradient(p):
        return np.broadcast_to([[[1.0, 0, 0]]], (len(p), 1, 3))

    values = mesh.points[:, 0, None] + 1
    errors = relative_errors(mesh, values, exact, gradient)
    assert_allclose(errors["l2_relative"], np.sqrt(3), rtol=2e-14, atol=0)
    assert errors["h1_relative"] < 1e-13


def test_budget_selection_reports_unattained_target():
    rows = [
        {
            "case": "heat",
            "backend": "scipy",
            "order": order,
            "size": n,
            "wall_s": time,
            "rss_bytes": 100,
            "h1_relative": error,
        }
        for order, n, time, error in [
            (1, 2, 1.0, 0.4),
            (1, 4, 4.0, 0.2),
            (2, 2, 2.0, 0.1),
        ]
        for _ in range(3)
    ]
    report = comparison(rows, [0.25, 0.05], [1.5, 3.0], 200)
    assert report["accuracy"][0]["selected"]["size"] == 4
    assert report["accuracy"][1]["selected"] is None
    assert report["cost"][0]["selected"]["size"] == 2
    assert report["cost"][1]["selected"]["size"] == 2
    assert report["groups"][0]["wall_s"]["stdev"] == 0


@pytest.mark.parametrize("case", ["heat", "bending"])
@pytest.mark.parametrize("order", [1, 2])
def test_accuracy_worker(case, order, tmp_path):
    script = Path(__file__).parents[2] / "benchmarks/accuracy_cost.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--worker",
            case,
            "--order",
            str(order),
            "--size",
            "2",
            "--backend",
            "scipy",
            "--output",
            str(tmp_path / "sample"),
        ],
        check=True,
        env={**os.environ, "PYTHONPATH": "src"},
    )
    report = json.loads((tmp_path / "sample/report.json").read_text())
    assert 0 < report["h1_relative"] < 1
    assert report["residual"] < 1e-9
    assert report["systems"][0]["rows"] > 0
    assert report["systems"][0]["nnz"] > 0
    assert report["wall_s"] > 0
    assert (tmp_path / "sample/field.npy").is_file()
