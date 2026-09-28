"""Native execution of the accuracy/cost worker, without timing thresholds."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from benchmarks.accuracy_cost import worker


@pytest.mark.parametrize("case", ["heat", "bending"])
@pytest.mark.parametrize("order", [1, 2])
def test_native_accuracy_worker(case, order, petsc_backend, tmp_path):
    fields = []
    for backend in ("scipy", "petsc"):
        output = tmp_path / backend
        worker(
            SimpleNamespace(
                worker=case,
                size=2,
                order=order,
                backend=backend,
                petsc_scalar=None,
                output=output,
            )
        )
        report = json.loads((output / "report.json").read_text())
        assert report["residual"] < 1e-9
        assert 0 < report["h1_relative"] < 1
        if backend == "petsc":
            assert report["backend_name"] == petsc_backend.name
        fields.append(np.load(output / "field.npy"))
    assert_allclose(fields[0], fields[1], rtol=1e-8, atol=1e-10)
