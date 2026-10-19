"""Persisted design inputs and results are hash bound and reproducible."""

import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ValidationError
from fem_toolbox.topology_study import TopologyStudy, read_topology_result
from tests.unit.test_topology import problem


def test_export_roundtrip_and_drift_detection(cube_mesh, tmp_path):
    study = TopologyStudy(problem(cube_mesh(1)), 0.6)
    result = study.run(tmp_path / "result")
    saved = read_topology_result(tmp_path / "result")
    assert saved["schema"] == 1
    assert saved["status"] == result.status
    assert_allclose(saved["density"], result.final.density, atol=0, rtol=0)
    assert saved["units"]["compliance"] == "J"
    assert saved["provenance"]["input_sha256"]
    assert (tmp_path / "result/design.vtu").is_file()
    source = tmp_path / "result/input.json"
    source.write_text(source.read_text() + " ")
    with pytest.raises(ValidationError, match="hash"):
        read_topology_result(tmp_path / "result")


def test_repeat_has_same_input_hash_and_export(cube_mesh, tmp_path):
    study = TopologyStudy(problem(cube_mesh(1)), 0.6)
    study.run(tmp_path / "a")
    study.run(tmp_path / "b")
    assert (tmp_path / "a/result.json").read_bytes() == (
        tmp_path / "b/result.json"
    ).read_bytes()
    with pytest.raises(FileExistsError):
        study.run(tmp_path / "a")


def test_unknown_or_escaping_configuration_is_rejected(tmp_path):
    path = tmp_path / "study.toml"
    path.write_text(
        'schema_version=1\nmodels=["../elsewhere.toml"]\n[design]\nvolume_fraction=0.6\nfilter_radius=0.3\n'
    )
    with pytest.raises((ValueError, ValidationError)):
        TopologyStudy.from_file(path)
    path.write_text("schema_version=1\nsurprise=true\nmodels=[]\n")
    with pytest.raises(ValidationError, match="unknown"):
        TopologyStudy.from_file(path)


def test_model_snapshot_and_declared_recipe_roundtrip():
    from pathlib import Path
    from fem_toolbox.model import Model

    root = Path(__file__).parents[2] / "examples/topology/reference"
    original = Model.from_file(root / "y.toml")
    snapshot = original.elasticity_inputs()
    snapshot["body_force"][0] = 100
    assert original.elasticity_inputs()["body_force"][0] == 0
    study = TopologyStudy.from_file(root / "study.toml")
    assert len(study.problem.cases) == 2
    assert len(study.problem.solid) > 0
    state = study.problem.evaluate(study.problem.initial_density(0.6))
    assert state.volume_fraction <= 0.6 + 1e-10
