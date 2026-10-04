"""Real study exports, independent reductions, status and provenance."""

import csv
import json
from pathlib import Path
from shutil import copytree

import pytest

from fem_toolbox.studies import Study
from fem_toolbox.study_reports import Quantity, export_study_table, read_study_table

ROOT = Path(__file__).parents[2]


def study(tmp_path):
    copytree(ROOT / "examples/parameter_study", tmp_path / "input")
    return Study.from_file(tmp_path / "input/study.toml")


def test_real_export_without_solve_and_missing_cases(tmp_path, monkeypatch):
    from fem_toolbox.model import Model

    runner = study(tmp_path)
    runner.run(tmp_path / "runs", max_cases=1)

    def forbidden(*args):
        raise AssertionError("export must not solve")

    monkeypatch.setattr(Model, "solve", forbidden)
    quantities = [
        Quantity(
            "temperature",
            "point",
            "temperature_K",
            "K",
            "maximum active-node temperature",
            reduction="max",
        ),
        Quantity(
            "balance",
            "metric",
            "heat.power_balance_W",
            "W",
            "signed heat power imbalance",
            absolute_tolerance=1e-11,
        ),
    ]
    export_study_table(tmp_path / "runs/catalog.json", tmp_path / "table", quantities)
    table = read_study_table(tmp_path / "table")
    assert table["rows"][0]["values"]["parameter:conductivity"] == 2
    assert table["rows"][0]["values"]["temperature"] == pytest.approx(301, abs=1e-11)
    assert table["rows"][0]["balance_status"] == "passed"
    assert table["rows"][0]["mesh_convergence"] == "not_checked"
    assert table["rows"][1]["values"]["temperature"] is None
    assert table["rows"][1]["missing"]["temperature"]
    rows = list(csv.DictReader((tmp_path / "table/values.csv").open()))
    assert rows[1]["temperature"] == ""
    export_study_table(tmp_path / "runs/catalog.json", tmp_path / "again", quantities)
    for filename in ["values.csv", "metadata.json", "report.md"]:
        assert (tmp_path / "table" / filename).read_bytes() == (
            tmp_path / "again" / filename
        ).read_bytes()


def test_corruption_unit_and_overwrite(tmp_path):
    runner = study(tmp_path)
    c = runner.run(tmp_path / "runs")
    (tmp_path / "runs" / c["cases"][0]["result"]["fields"]).write_text("corrupt")
    q = Quantity(
        "temperature", "point", "temperature_K", "Pa", "wrong unit", reduction="max"
    )
    export_study_table(tmp_path / "runs/catalog.json", tmp_path / "table", [q])
    t = read_study_table(tmp_path / "table")
    assert t["rows"][0]["status"] == "invalid"
    assert "unit" in t["rows"][1]["missing"]["temperature"]
    with pytest.raises(FileExistsError):
        export_study_table(tmp_path / "runs/catalog.json", tmp_path / "table", [q])
    (tmp_path / "table/values.csv").write_text("corrupt")
    with pytest.raises(ValueError, match="hash"):
        read_study_table(tmp_path / "table")


def test_coupled_quantities_and_failed_balance(tmp_path):
    import meshio
    import numpy as np

    copytree(ROOT / "examples/thermopiezoelectric", tmp_path / "input")
    runner = Study(
        {"study_schema_version": 1, "model": "case.toml", "parameters": []},
        root=tmp_path / "input",
    )
    c = runner.run(tmp_path / "runs")
    bundle = json.loads(
        (tmp_path / "runs" / c["cases"][0]["result"]["manifest"]).read_text()
    )
    fields = meshio.read(tmp_path / "runs" / c["cases"][0]["result"]["fields"])
    quantities = [
        Quantity(
            "u",
            "point",
            "displacement_m",
            "m",
            "node 1 x displacement",
            component="x",
            indices=(1,),
        ),
        Quantity(
            "mechanical",
            "cell",
            "mechanical_energy_J",
            "J",
            "sum of mechanical cell energies",
            reduction="sum",
        ),
        Quantity(
            "coupling",
            "cell",
            "coupling_energy_J",
            "J",
            "sum of coupling cell energies",
            reduction="sum",
        ),
        Quantity(
            "dielectric",
            "cell",
            "dielectric_enthalpy_J",
            "J",
            "sum of dielectric cell enthalpies",
            reduction="sum",
        ),
        Quantity(
            "charge",
            "metric",
            "piezoelectric.charge_balance_C",
            "C",
            "signed global charge imbalance",
        ),
        Quantity(
            "force",
            "metric",
            "piezoelectric.force_balance_N",
            "N",
            "force-balance norm",
            component="norm",
            absolute_tolerance=1e-10,
        ),
        Quantity(
            "power",
            "metric",
            "heat.source_power_W",
            "W",
            "deliberately strict check of source power",
            absolute_tolerance=0,
        ),
    ]
    export_study_table(tmp_path / "runs/catalog.json", tmp_path / "table", quantities)
    row = read_study_table(tmp_path / "table")["rows"][0]
    assert not row["missing"], row["missing"]
    assert row["values"]["u"] == pytest.approx(-9 / 440, abs=1e-12)
    for key, name in [
        ("mechanical", "mechanical_energy_J"),
        ("coupling", "coupling_energy_J"),
        ("dielectric", "dielectric_enthalpy_J"),
    ]:
        assert row["values"][key] == np.sum(fields.cell_data_dict[name]["tetra"])
    assert (
        row["values"]["charge"]
        == bundle["metrics"]["piezoelectric"]["charge_balance_C"]
    )
    assert row["balance_checks"]["force"] == "passed"


def test_failure_missing_metric_and_invalid_selection(tmp_path):
    runner = study(tmp_path)
    data = runner._manifest.copy()
    data["parameters"][0]["values"] = [2.0, -1.0]
    runner = Study(data, root=tmp_path / "input")
    runner.run(tmp_path / "runs")
    quantities = [
        Quantity(
            "bad_node", "point", "temperature_K", "K", "absent node", indices=(100,)
        ),
        Quantity("unknown", "metric", "heat.absent_W", "W", "absent metric"),
        Quantity(
            "source",
            "metric",
            "heat.source_power_W",
            "W",
            "zero source test",
            absolute_tolerance=0,
        ),
    ]
    export_study_table(tmp_path / "runs/catalog.json", tmp_path / "table", quantities)
    rows = read_study_table(tmp_path / "table")["rows"]
    assert rows[0]["balance_status"] == "failed"
    assert rows[0]["missing"]["bad_node"]
    assert rows[0]["missing"]["unknown"]
    assert rows[1]["status"] == "failed"
    assert rows[1]["values"]["source"] is None
    assert rows[1]["balance_status"] == "not_available"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"reduction": "bad"},
        {"indices": (-1,)},
        {"indices": (1, 1)},
        {"absolute_tolerance": float("nan")},
        {"name": "case_id"},
    ],
)
def test_invalid_definitions(kwargs):
    base = dict(
        name="q",
        source="metric",
        field="heat.source_power_W",
        unit="W",
        definition="source",
    )
    base.update(kwargs)
    with pytest.raises(ValueError):
        Quantity(**base)


@pytest.mark.parametrize("relocate", [False, True])
def test_swapped_results_rejected_by_report_and_resume(tmp_path, relocate):
    runner = study(tmp_path)
    root = tmp_path / "runs"
    catalog = runner.run(root)
    first, second = catalog["cases"]
    first["result"], second["result"] = second["result"], first["result"]
    (root / "catalog.json").write_text(json.dumps(catalog))
    if relocate:
        copytree(root, tmp_path / "copied")
        root = tmp_path / "copied"
    q = Quantity(
        "temperature", "point", "temperature_K", "K", "maximum", reduction="max"
    )
    export_study_table(root / "catalog.json", tmp_path / "table", [q])
    rows = read_study_table(tmp_path / "table")["rows"]
    assert all(row["status"] == "invalid" for row in rows)
    assert all(row["values"]["temperature"] is None for row in rows)
    assert all(row["missing"]["temperature"] for row in rows)
    repaired = runner.run(root)
    assert all(case["attempts"] == 2 for case in repaired["cases"])
    assert all(
        case["history"][-1]["result"] == old["result"]
        for case, old in zip(repaired["cases"], catalog["cases"], strict=True)
    )
    export_study_table(root / "catalog.json", tmp_path / "repaired", [q])
    values = [
        row["values"]["temperature"]
        for row in read_study_table(tmp_path / "repaired")["rows"]
    ]
    # Analytic unit-tetrahedron solution; 1e-11 K covers float64 roundoff at 300 K.
    assert values == pytest.approx([301, 300.5], rel=0, abs=1e-11)
    assert runner.run(root) == repaired


@pytest.mark.parametrize("damage", ["foreign_bundle", "document", "legacy"])
def test_result_binding_checked_independently_of_paths(tmp_path, damage):
    from hashlib import sha256

    from fem_toolbox.studies import _valid_result

    runner = study(tmp_path)
    root = tmp_path / "runs"
    catalog = runner.run(root)
    case = catalog["cases"][0]
    result = case["result"]
    manifest = root / result["manifest"]
    if damage == "foreign_bundle":
        foreign = catalog["cases"][1]["result"]
        manifest.write_bytes((root / foreign["manifest"]).read_bytes())
        (root / result["fields"]).write_bytes((root / foreign["fields"]).read_bytes())
        result["fields_sha256"] = foreign["fields_sha256"]
        result["model_input_sha256"] = foreign["model_input_sha256"]
    elif damage == "document":
        case["document"]["materials"]["T"]["thermal"]["conductivity"] = 9.0
    else:
        bundle = json.loads(manifest.read_text())
        del bundle["provenance"]["solver"]["parameters"]["study_case"]
        manifest.write_text(json.dumps(bundle))
    result["manifest_sha256"] = sha256(manifest.read_bytes()).hexdigest()
    assert not _valid_result(root, case)


def test_intact_copied_study_reused_without_solve(tmp_path, monkeypatch):
    from fem_toolbox.model import Model

    runner = study(tmp_path)
    catalog = runner.run(tmp_path / "runs")
    copytree(tmp_path / "runs", tmp_path / "copied")

    def forbidden(*args):
        raise AssertionError("intact relocated results must be reused")

    monkeypatch.setattr(Model, "solve", forbidden)
    assert runner.run(tmp_path / "copied") == catalog
