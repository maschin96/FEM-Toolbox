"""M8 joint sampling, physical heat reference, recovery and provenance."""

import json
import shutil
from pathlib import Path

import pytest
from numpy.testing import assert_allclose

from fem_toolbox.study_reports import Quantity
from fem_toolbox.uncertainty import (
    Normal,
    SamplingPlan,
    UncertainParameter,
    UncertaintyStudy,
    Uniform,
)

ROOT = Path(__file__).resolve().parents[2]


def fixture(tmp_path):
    inputs = tmp_path / "inputs"
    shutil.copytree(ROOT / "examples/surrogate_heat", inputs)
    plan = SamplingPlan(
        (
            UncertainParameter(
                "source",
                "loads.source.value",
                "W/m^3",
                Uniform(4, 12),
                "synthetic teaching assumption",
            ),
            UncertainParameter(
                "temperature",
                "boundary.dirichlet.0.value",
                "K",
                Normal(300, 1),
                "synthetic teaching assumption",
            ),
        ),
        size=16,
        seed=64,
    )
    q = Quantity(
        "Tmax",
        "point",
        "temperature_K",
        "K",
        "maximum nodal temperature",
        reduction="max",
    )
    return inputs, plan, q


def test_joint_rows_heat_reference_resume_and_sensitivity(tmp_path):
    inputs, plan, q = fixture(tmp_path)
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    output = tmp_path / "output"
    report = study.run(
        output,
        quantities=[q],
        numerical_error_bounds={"Tmax": 1e-10},
        sensitivity_center=[8, 300],
        sensitivity_steps=[0.1, 0.01],
    )
    assert report["complete"]
    assert len(report["samples"]) == 16  # Joint rows, never 16**2 Cartesian cases.
    x = plan.sample()
    actual = [r["values"]["Tmax"] for r in report["samples"]]
    assert_allclose(actual, x[:, 1] + x[:, 0] / 8, atol=1e-11, rtol=0)
    assert_allclose(
        report["sensitivity"]["Tmax"]["gradient"], [1 / 8, 1], atol=2e-10, rtol=0
    )
    assert report["summaries"]["Tmax"]["numerical_error_bound"] == 1e-10
    assert len(report["sample_sha256"]) == 64
    first = json.loads((output / "uncertainty.json").read_text())
    resumed = study.run(
        output,
        quantities=[q],
        numerical_error_bounds={"Tmax": 1e-10},
        sensitivity_center=[8, 300],
        sensitivity_steps=[0.1, 0.01],
    )
    assert resumed == first
    for sample in resumed["samples"]:
        catalog = json.loads((output / sample["catalog"]).read_text())
        assert len(catalog["cases"]) == 1
        assert catalog["cases"][0]["attempts"] == 1
        assert sample["row"]["provenance"]["input_sha256"]


def test_failed_physical_case_remains_visible_and_no_conditional_statistics(tmp_path):
    inputs, _, q = fixture(tmp_path)
    plan = SamplingPlan(
        (
            UncertainParameter(
                "k",
                "materials.T.thermal.conductivity",
                "W/(m K)",
                Uniform(-2, -1),
                "deliberately invalid negative material",
            ),
        ),
        size=3,
        seed=1,
    )
    report = UncertaintyStudy(plan, model="model.toml", root=inputs).run(
        tmp_path / "output", quantities=[q], numerical_error_bounds={"Tmax": None}
    )
    assert not report["complete"]
    assert len(report["samples"]) == 3
    assert all(s["row"]["status"] in {"failed", "invalid"} for s in report["samples"])
    assert report["summaries"] == {}


def test_new_plan_cannot_reuse_existing_output(tmp_path):
    inputs, plan, q = fixture(tmp_path)
    output = tmp_path / "output"
    UncertaintyStudy(plan, model="model.toml", root=inputs).run(
        output, quantities=[q], numerical_error_bounds={"Tmax": None}
    )
    other = SamplingPlan(plan.parameters, size=plan.size, seed=65)
    with pytest.raises(ValueError, match="sampling"):
        UncertaintyStudy(other, model="model.toml", root=inputs).run(
            output, quantities=[q], numerical_error_bounds={"Tmax": None}
        )


def test_modified_model_invalidates_results(tmp_path):
    inputs, plan, q = fixture(tmp_path)
    output = tmp_path / "output"
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    first = study.run(output, quantities=[q], numerical_error_bounds={"Tmax": None})
    path = inputs / "model.toml"
    path.write_text(
        path.read_text().replace("conductivity = 2.0", "conductivity = 4.0")
    )
    second = study.run(output, quantities=[q], numerical_error_bounds={"Tmax": None})
    assert second["complete"]
    assert (
        second["samples"][0]["catalog_sha256"] != first["samples"][0]["catalog_sha256"]
    )
    x = plan.sample()
    assert_allclose(
        [s["values"]["Tmax"] for s in second["samples"]],
        x[:, 1] + x[:, 0] / 16,
        atol=1e-11,
        rtol=0,
    )


def test_missing_budget_and_invalid_sensitivity_fail_before_solving(tmp_path):
    inputs, plan, q = fixture(tmp_path)
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    with pytest.raises(ValueError, match="bounds"):
        study.run(tmp_path / "output", quantities=[q], numerical_error_bounds={})
    assert not (tmp_path / "output").exists()
    with pytest.raises(ValueError, match="steps"):
        study.run(
            tmp_path / "output",
            quantities=[q],
            numerical_error_bounds={"Tmax": None},
            sensitivity_center=[8, 300],
            sensitivity_steps=[0, 0.1],
        )
    assert not (tmp_path / "output").exists()


def test_missing_output_remains_visible_and_lock_protects_writer(tmp_path):
    inputs, plan, _ = fixture(tmp_path)
    q = Quantity(
        "absent", "point", "missing_field", "K", "negative fixture", reduction="max"
    )
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    output = tmp_path / "output"
    report = study.run(output, quantities=[q], numerical_error_bounds={"absent": None})
    assert not report["complete"]
    assert report["summaries"] == {}
    assert all(s["row"]["missing"]["absent"] for s in report["samples"])
    previous = (output / "uncertainty.json").read_bytes()
    (output / ".uncertainty.lock").write_text("other writer")
    with pytest.raises(FileExistsError):
        study.run(output, quantities=[q], numerical_error_bounds={"absent": None})
    assert (output / "uncertainty.json").read_bytes() == previous
    assert (output / ".uncertainty.lock").exists()


def test_sensitivity_step_refinement_on_independent_heat_solution(tmp_path):
    inputs, plan, q = fixture(tmp_path)
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    first = study.run(
        tmp_path / "output",
        quantities=[q],
        numerical_error_bounds={"Tmax": None},
        sensitivity_center=[8, 300],
        sensitivity_steps=[0.2, 0.02],
    )
    second = study.run(
        tmp_path / "output",
        quantities=[q],
        numerical_error_bounds={"Tmax": None},
        sensitivity_center=[8, 300],
        sensitivity_steps=[0.1, 0.01],
    )
    assert_allclose(
        first["sensitivity"]["Tmax"]["gradient"], [1 / 8, 1], atol=2e-10, rtol=0
    )
    assert_allclose(
        second["sensitivity"]["Tmax"]["gradient"], [1 / 8, 1], atol=2e-10, rtol=0
    )


def test_table_receipt_must_match_current_catalog(tmp_path):
    inputs, plan, q = fixture(tmp_path)
    output = tmp_path / "output"
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    report = study.run(output, quantities=[q], numerical_error_bounds={"Tmax": None})
    path = output / report["samples"][0]["table"] / "metadata.json"
    data = json.loads(path.read_text())
    data["catalog_sha256"] = "0" * 64
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="catalog"):
        study.run(output, quantities=[q], numerical_error_bounds={"Tmax": None})


def test_failed_rebuild_does_not_leave_prior_success_report(tmp_path, monkeypatch):
    inputs, plan, q = fixture(tmp_path)
    output = tmp_path / "output"
    study = UncertaintyStudy(plan, model="model.toml", root=inputs)
    study.run(output, quantities=[q], numerical_error_bounds={"Tmax": None})
    from fem_toolbox.studies import Study

    def fail(*args, **kwargs):
        raise OSError("interrupted filesystem")

    monkeypatch.setattr(Study, "run", fail)
    with pytest.raises(OSError, match="interrupted"):
        study.run(output, quantities=[q], numerical_error_bounds={"Tmax": None})
    assert not (output / "uncertainty.json").exists()
    assert not (output / ".uncertainty.lock").exists()
