"""Independent barycentric integration for the handbook's local extension."""

import runpy
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

ROOT = Path(__file__).resolve().parents[2]


def test_developer_source_integrates_linear_density_independently():
    module = runpy.run_path(str(ROOT / "Doc/latex/examples/developer_extension.py"))
    from fem_toolbox.elements import Tet4Geometry
    from fem_toolbox.physics.heat import source_load

    geometry = Tet4Geometry(np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    source = module["LinearHeatSource"](2, (3, 0, 0))
    # Integral lambda_i*lambda_j = V/20 for i != j, V/10 for i == j.
    assert_allclose(
        source_load(geometry, source),
        [13 / 120, 2 / 15, 13 / 120, 13 / 120],
        atol=2e-15,
        rtol=0,
    )
    with pytest.raises(ValueError):
        module["LinearHeatSource"](2, (np.nan, 0, 0))
