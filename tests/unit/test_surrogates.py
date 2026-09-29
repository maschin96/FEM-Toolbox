"""Independent polynomial functions and strict data/split contracts."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.surrogates import (
    DataSplit,
    PolynomialSurrogate,
    SurrogateData,
    Variable,
    fit_surrogate,
    split_cases,
)


def linear_data():
    x = np.array([[-2.0], [-1.0], [0.0], [1.0], [2.0], [0.5], [1.5]])
    y = np.column_stack([3 + 2 * x[:, 0], np.full(len(x), 7.0)])
    return SurrogateData(
        x,
        y,
        [str(i) for i in range(len(x))],
        [Variable("x", "m", "coordinate")],
        [
            Variable("temperature", "K", "linear reference"),
            Variable("constant", "W", "constant output"),
        ],
        provenance={"source": "independent analytical function"},
    )


def test_linear_multioutput_and_report():
    data = linear_data()
    split = DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    model = fit_surrogate(data, split, degree=1, ridge=0)
    prediction = model.predict([[0.5], [1.5]])
    assert_allclose(prediction.values, [[4, 7], [6, 7]], atol=1e-12, rtol=0)
    assert prediction.extrapolated.tolist() == [False, False]
    assert prediction.targets[0].unit == "K"
    assert model.report["test"]["max_abs"][0] < 1e-12
    assert model.report["split"]["train_ids"] == list(split.train_ids)
    assert len(model.report["dataset_sha256"]) == 64


def test_quadratic_interaction_unseen_points():
    x = np.array(
        [(a, b) for a in [-1.0, 0.0, 1.0] for b in [-1.0, 0.0, 1.0]]
        + [(0.2, 0.3), (-0.4, 0.7)]
    )
    y = (
        2
        + 3 * x[:, 0]
        - x[:, 1]
        + 4 * x[:, 0] ** 2
        + 5 * x[:, 0] * x[:, 1]
        - 2 * x[:, 1] ** 2
    )[:, None]
    data = SurrogateData(
        x,
        y,
        [str(i) for i in range(len(x))],
        [Variable("a", "1", "a"), Variable("b", "1", "b")],
        [Variable("y", "1", "analytic quadratic")],
        provenance={"source": "formula"},
    )
    model = fit_surrogate(
        data, DataSplit(tuple(str(i) for i in range(9)), ("9", "10")), degree=2
    )
    assert_allclose(model.predict(x[9:]).values, y[9:], atol=2e-13, rtol=0)


@pytest.mark.parametrize(
    "split",
    [
        DataSplit(("0", "1"), ("1", "2")),
        DataSplit(("0", "1"), ()),
        DataSplit(("missing",), ("1",)),
    ],
)
def test_bad_splits(split):
    with pytest.raises(ValueError):
        fit_surrogate(linear_data(), split)


def test_extrapolation_is_explicit():
    model = fit_surrogate(
        linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    )
    with pytest.raises(ValueError, match="outside"):
        model.predict([3.0])
    prediction = model.predict([[-2.0], [3.0]], allow_extrapolation=True)
    assert prediction.extrapolated.tolist() == [False, True]
    assert_allclose(prediction.values, [[-1, 7], [9, 7]], atol=1e-12, rtol=0)


def test_ridge_preserves_intercept_and_reports_mismatch():
    x = np.array([[-2.0], [-1.0], [0.0], [1.0], [2.0], [0.5]])
    data = SurrogateData(
        x,
        np.column_stack([x[:, 0] ** 2, np.full(6, 7.0)]),
        [str(i) for i in range(6)],
        [Variable("x", "1", "x")],
        [Variable("square", "1", "square"), Variable("constant", "1", "constant")],
        provenance={"source": "formula"},
    )
    model = fit_surrogate(
        data, DataSplit(("0", "1", "2", "3", "4"), ("5",)), degree=1, ridge=1000
    )
    assert_allclose(model.predict([0.5]).values, [[2.0, 7.0]], atol=1e-12, rtol=0)
    # Symmetric training gives intercept 2; test is .25 -> absolute error 1.75.
    for metric in ["mae", "rmse", "max_abs"]:
        assert model.report["test"][metric][0] == pytest.approx(1.75, abs=1e-12)


def test_scaling_never_uses_test_data():
    original = linear_data()
    split = DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    a = fit_surrogate(original, split)
    x = original.X.copy()
    y = original.Y.copy()
    x[-1] = 100
    y[-1] = [999, 999]
    data = SurrogateData(
        x,
        y,
        original.case_ids,
        original.parameters,
        original.targets,
        provenance={"source": "changed test"},
    )
    b = fit_surrogate(data, split)
    for name in ["center", "scale", "minimum", "maximum", "coefficients"]:
        assert a.to_document()[name] == b.to_document()[name]
    assert b.report["test_extrapolated"] == [False, True]


def test_leakage_duplicate_points_and_seeded_groups():
    d = linear_data()
    x = d.X.copy()
    x[-1] = x[0]
    data = SurrogateData(
        x, d.Y, d.case_ids, d.parameters, d.targets, provenance={"source": "duplicate"}
    )
    with pytest.raises(ValueError, match="identical"):
        fit_surrogate(data, DataSplit(("0", "1", "2", "3", "4"), ("5", "6")))
    a = split_cases(data, test_fraction=0.3, seed=11)
    assert a == split_cases(data, test_fraction=0.3, seed=11)
    assert ("0" in a.train_ids) == ("6" in a.train_ids)
    assert a.seed == 11


@pytest.mark.parametrize("kind", ["constant", "under", "rank"])
def test_degenerate_training(kind):
    if kind == "constant":
        x = np.array([[1.0], [1.0], [1.0], [2.0]])
        degree = 1
    elif kind == "under":
        x = np.array([[0.0], [1.0], [2.0], [3.0]])
        degree = 2
    else:
        x = np.array([[0.0, 0.0], [1.0, 1.0], [2.0, 2.0], [3.0, 3.0]])
        degree = 1
    data = SurrogateData(
        x,
        np.arange(4.0)[:, None],
        ["a", "b", "c", "d"],
        [Variable(f"x{i}", "1", "parameter") for i in range(x.shape[1])],
        [Variable("y", "1", "target")],
        provenance={"source": "fixture"},
    )
    split = (
        DataSplit(("a", "b"), ("c", "d"))
        if kind == "under"
        else DataSplit(("a", "b", "c"), ("d",))
    )
    with pytest.raises(
        ValueError,
        match={
            "constant": "constant",
            "under": "underdetermined",
            "rank": "rank-deficient",
        }[kind],
    ):
        fit_surrogate(data, split, degree=degree)
    if kind == "rank":
        model = fit_surrogate(data, split, ridge=1)
        assert np.isfinite(model.predict([0.5, 0.5]).values).all()


@pytest.mark.parametrize("bad", [complex(1, 1), float("nan"), float("inf"), True])
def test_invalid_data(bad):
    with pytest.raises(ValueError):
        SurrogateData(
            [[bad]],
            [[1]],
            ["id"],
            [Variable("x", "1", "x")],
            [Variable("y", "1", "y")],
            provenance={"source": "invalid"},
        )


def test_json_roundtrip_overwrite_and_corruption(tmp_path):
    import json

    model = fit_surrogate(
        linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    )
    path = tmp_path / "model.json"
    model.save(path)
    loaded = PolynomialSurrogate.load(path)
    assert loaded.report == model.report
    assert_allclose(
        loaded.predict([[0.5]]).values, model.predict([[0.5]]).values, atol=0, rtol=0
    )
    with pytest.raises(FileExistsError):
        model.save(path)
    model.save(path, overwrite=True)
    doc = json.loads(path.read_text())
    doc["model"]["coefficients"][0][0] += 1
    path.write_text(json.dumps(doc))
    with pytest.raises(ValueError, match="checksum"):
        PolynomialSurrogate.load(path)


def test_artifact_structural_validation():
    model = fit_surrogate(
        linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    )
    for key, value in [
        ("surrogate_schema_version", 99),
        ("coefficients", [[1.0]]),
        ("scale", [0]),
        ("basis", [[0], [2]]),
    ]:
        doc = model.to_document()
        doc[key] = value
        with pytest.raises(ValueError):
            PolynomialSurrogate(doc)


def test_real_table_adapter_and_explicit_exclusion(tmp_path):
    import json
    from pathlib import Path
    from shutil import copytree

    from fem_toolbox.studies import Study
    from fem_toolbox.study_reports import Quantity, export_study_table

    copytree(Path(__file__).parents[2] / "examples/parameter_study", tmp_path / "input")
    source = tmp_path / "input/study.toml"
    source.write_text(source.read_text().replace("[2.0, 4.0]", "[2.0, 4.0, -1.0]"))
    Study.from_file(source).run(tmp_path / "runs")
    export_study_table(
        tmp_path / "runs/catalog.json",
        tmp_path / "table",
        [
            Quantity(
                "temperature",
                "point",
                "temperature_K",
                "K",
                "maximum temperature",
                reduction="max",
            )
        ],
    )
    with pytest.raises(ValueError, match=r"excluded|invalid"):
        SurrogateData.from_table(
            tmp_path / "table",
            parameter_columns=["parameter:conductivity"],
            target_columns=["temperature"],
        )
    data = SurrogateData.from_table(
        tmp_path / "table",
        parameter_columns=["parameter:conductivity"],
        target_columns=["temperature"],
        exclude_invalid=True,
    )
    assert data.X.shape == (2, 1)
    assert_allclose(data.Y[:, 0], [301, 300.5], rtol=0, atol=1e-11)
    assert data.parameters[0].unit == "W/(m K)"
    assert len(data.provenance["exclusions"]) == 1
    assert data.provenance["exclusions"][0]["reasons"]
    table = json.loads((tmp_path / "table/metadata.json").read_text())
    assert data.case_ids == tuple(r["case_id"] for r in table["rows"][:2])
    with pytest.raises(ValueError, match="role"):
        SurrogateData.from_table(
            tmp_path / "table",
            parameter_columns=["temperature"],
            target_columns=["parameter:conductivity"],
        )


def test_extreme_ridge_rejects_numerical_rank_loss():
    with pytest.raises(ValueError, match="rank"):
        fit_surrogate(
            linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6")), ridge=1e100
        )


def test_hand_computed_metrics_multiple_test_points():
    d = linear_data()
    data = SurrogateData(
        d.X,
        d.X**2,
        d.case_ids,
        d.parameters,
        [Variable("square", "m^2", "coordinate squared")],
        provenance={"source": "formula"},
    )
    model = fit_surrogate(data, DataSplit(("0", "1", "2", "3", "4"), ("5", "6")))
    assert model.report["test"]["mae"] == pytest.approx([1.0], abs=1e-12)
    assert model.report["test"]["rmse"] == pytest.approx([1.25], abs=1e-12)
    assert model.report["test"]["max_abs"] == pytest.approx([1.75], abs=1e-12)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"degree": True},
        {"degree": 3},
        {"ridge": -1},
        {"ridge": float("nan")},
        {"ridge": True},
    ],
)
def test_invalid_fit_options(kwargs):
    with pytest.raises(ValueError):
        fit_surrogate(
            linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6")), **kwargs
        )


def test_owned_inputs_and_metadata():
    x = np.array([[0.0], [1.0], [2.0]])
    y = np.array([[3.0], [4.0], [5.0]])
    provenance = {"source": {"dataset": "original"}}
    d = SurrogateData(
        x,
        y,
        ["a", "b", "c"],
        [Variable("x", "m", "position")],
        [Variable("y", "K", "temperature")],
        provenance=provenance,
    )
    fingerprint = d.fingerprint
    x[:] = 100
    y[:] = 200
    provenance["source"]["dataset"] = "changed"
    assert d.fingerprint == fingerprint
    with pytest.raises(ValueError):
        d.X.setflags(write=True)
    bad = d.provenance
    bad["source"] = "changed"
    assert d.provenance["source"]["dataset"] == "original"
    with pytest.raises(ValueError):
        SurrogateData(
            [[1], [2]],
            [[3]],
            ["a", "b"],
            d.parameters,
            d.targets,
            provenance={"source": "bad"},
        )
    with pytest.raises(ValueError):
        SurrogateData(
            [[1], [2]],
            [[3], [4]],
            ["a", "a"],
            d.parameters,
            d.targets,
            provenance={"source": "bad"},
        )


def test_loaded_prediction_has_no_solver_or_data_dependency(tmp_path, monkeypatch):
    from fem_toolbox.model import Model

    model = fit_surrogate(
        linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    )
    path = tmp_path / "model.json"
    model.save(path)

    def forbidden(*args):
        raise AssertionError("FEM must not execute")

    monkeypatch.setattr(Model, "solve", forbidden)
    loaded = PolynomialSurrogate.load(path)
    assert_allclose(loaded.predict([0.5]).values, [[4, 7]], atol=1e-12, rtol=0)
    doc = loaded.to_document()
    doc["coefficients"][0][0] = -999
    report = loaded.report
    report["test"]["mae"][0] = 999
    assert loaded.report["test"]["mae"][0] < 1e-12
    assert_allclose(loaded.predict([0.5]).values, [[4, 7]], atol=1e-12, rtol=0)
    for bad in ([float("nan")], [[1, 2]], [], [1j]):
        with pytest.raises(ValueError):
            loaded.predict(bad)


def test_artifact_missing_version_and_bad_json(tmp_path):
    model = fit_surrogate(
        linear_data(), DataSplit(("0", "1", "2", "3", "4"), ("5", "6"))
    )
    doc = model.to_document()
    del doc["report"]["versions"]["numpy"]
    with pytest.raises(ValueError, match="version"):
        PolynomialSurrogate(doc)
    path = tmp_path / "broken.json"
    path.write_text("{not json")
    with pytest.raises(ValueError):
        PolynomialSurrogate.load(path)
