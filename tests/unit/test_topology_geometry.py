"""A SIMP field is not a discrete solid until separately reanalysed."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ValidationError
from fem_toolbox.topology_geometry import reanalyze_geometry
from tests.unit.test_topology import problem


def test_full_solid_geometry_matches_independent_elastic_state(cube_mesh):
    design = problem(cube_mesh(1))
    field = design.evaluate(np.ones(6))
    result = reanalyze_geometry(design, field, 0.5)
    assert result.mesh.n_cells == 6
    assert result.volume_fraction == pytest.approx(1)
    assert_allclose(result.objective, field.objective, rtol=2e-10)
    assert result.maximum_von_mises_Pa > 0
    assert_allclose(result.states[0].force_balance, 0, atol=1e-10, rtol=0)


def test_empty_or_invalid_threshold_is_diagnosed(cube_mesh):
    design = problem(cube_mesh(1))
    field = design.evaluate(np.full(6, 0.2))
    with pytest.raises(ValidationError, match="empty"):
        reanalyze_geometry(design, field, 0.5)
    for threshold in (0, 1, float("nan")):
        with pytest.raises(ValidationError, match="threshold"):
            reanalyze_geometry(design, field, threshold)


def test_face_disconnected_remnants_are_diagnosed(cube_mesh):
    from fem_toolbox.topology import TopologyProblem

    base = problem(cube_mesh(1))
    design = TopologyProblem(base.mesh, base.materials, base.cases, filter_radius=1e-4)
    x = np.zeros(6)
    x[[0, 3]] = 1
    field = design.evaluate(x)
    with pytest.raises(ValidationError, match="disconnected"):
        reanalyze_geometry(design, field, 0.5)


def test_removed_point_load_is_never_silently_dropped(cube_mesh):
    from fem_toolbox.topology import TopologyProblem

    base = problem(cube_mesh(1))
    design = TopologyProblem(base.mesh, base.materials, base.cases, filter_radius=1e-4)
    x = np.zeros(6)
    x[0] = 1
    field = design.evaluate(x)
    with pytest.raises(ValidationError, match="loaded node"):
        reanalyze_geometry(design, field, 0.5)


def test_explicit_connection_repair_is_reported_on_reference(tmp_path):
    from pathlib import Path

    from fem_toolbox.topology_study import TopologyStudy

    study = TopologyStudy.from_file(
        Path(__file__).parents[2] / "examples/topology/reference/study.toml"
    )
    result = study.run(tmp_path / "density")
    with pytest.raises(ValidationError, match="disconnected"):
        reanalyze_geometry(study.problem, result.final, 0.45)
    repaired = reanalyze_geometry(study.problem, result.final, 0.45, connect=True)
    record = repaired.to_record()
    assert record["connection_repair"] is True
    assert len(record["added_original_cells"]) > 0
    assert record["status"] == "accepted"
    assert repaired.volume_fraction <= 0.9
    assert repaired.compliance_ratio <= 2
