"""Stationary studies: independent inputs, identity and durable resume."""

import json
from pathlib import Path
from shutil import copyfile

import pytest

from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.studies import Study

FIXTURE = Path(__file__).parents[1] / "fixtures" / "overlapping_groups.msh"


def inputs(root: Path) -> dict:
    copyfile(FIXTURE, root / "mesh.msh")
    (root / "model.toml").write_text("""schema_version = 1
units = "SI"
physics = "heat"
[mesh]
path = "mesh.msh"
length_scale_to_m = 1.0
[materials.T.solid]
conductivity = 2.0
unit = "W/(m K)"
[regions.T]
ceramic = "solid"
[[boundary.dirichlet]]
field = "T"
nodes = [0, 1, 2, 3, 4]
value = [300.0, 310.0, 300.0, 300.0, 310.0]
unit = "K"
""")
    return {
        "study_schema_version": 1,
        "model": "model.toml",
        "parameters": [
            {
                "name": "k",
                "target": "materials.T.solid.conductivity",
                "unit": "W/(m K)",
                "values": [2.0, 4.0],
            },
        ],
    }


def test_resolve_identity_and_input_changes(tmp_path: Path) -> None:
    manifest = inputs(tmp_path)
    study = Study(manifest, root=tmp_path)
    first = study.resolve()
    assert [c["parameters"]["k"]["value"] for c in first] == [2, 4]
    assert first == study.resolve()
    assert len({c["case_id"] for c in first}) == 2
    path = tmp_path / "model.toml"
    path.write_text(path.read_text().replace("300.0", "301.0"))
    changed = study.resolve()
    assert [c["case_id"] for c in changed] == [c["case_id"] for c in first]
    assert changed[0]["fingerprints"]["input"] != first[0]["fingerprints"]["input"]
    # Returned data never mutates the study.
    changed[0]["document"]["units"] = "bad"
    assert study.resolve()[0]["document"]["units"] == "SI"


@pytest.mark.parametrize("change", ["unit", "target", "nan", "duplicate", "schema"])
def test_bad_manifest(tmp_path: Path, change: str) -> None:
    manifest = inputs(tmp_path)
    p = manifest["parameters"][0]
    if change == "unit":
        p["unit"] = "Pa"
    elif change == "target":
        p["target"] = "materials.T.solid.missing"
    elif change == "nan":
        p["values"] = [float("nan")]
    elif change == "duplicate":
        manifest["parameters"].append(dict(p))
    else:
        manifest["study_schema_version"] = True
    with pytest.raises(ModelSchemaError):
        Study(manifest, root=tmp_path).resolve()


def test_toml_manifest(tmp_path: Path) -> None:
    inputs(tmp_path)
    source = tmp_path / "study.toml"
    source.write_text("""study_schema_version = 1
model = "model.toml"
[[parameters]]
name = "k"
target = "materials.T.solid.conductivity"
unit = "W/(m K)"
values = [2.0, 4.0]
""")
    assert len(Study.from_file(source).resolve()) == 2


def test_run_parity_resume_and_corruption(tmp_path: Path, monkeypatch) -> None:
    import meshio
    import numpy as np
    from numpy.testing import assert_allclose

    from fem_toolbox.io import export_result_bundle
    from fem_toolbox.model import Model

    study = Study(inputs(tmp_path), root=tmp_path)
    output = tmp_path / "results"
    catalog = study.run(output)
    assert [c["status"] for c in catalog["cases"]] == ["succeeded", "succeeded"]
    # Build independent direct cases, without using Study.resolve's documents.
    base = Model.from_file(tmp_path / "model.toml")
    for i, conductivity in enumerate([2.0, 4.0]):
        doc = base.to_document()
        doc["materials"]["T"]["solid"]["conductivity"] = conductivity
        direct = Model(mesh=base.mesh, document=doc).solve()
        _, reference = export_result_bundle(tmp_path / f"direct-{i}", base.mesh, direct)
        actual_path = output / catalog["cases"][i]["result"]["manifest"]
        actual = json.loads(actual_path.read_text())
        expected = json.loads(reference.read_text())
        assert actual["metrics"] == expected["metrics"]
        fields = meshio.read(actual_path.with_suffix(".vtu"))
        reference_fields = meshio.read(reference.with_suffix(".vtu"))
        for name in fields.point_data:
            assert_allclose(
                fields.point_data[name],
                reference_fields.point_data[name],
                rtol=1e-12,
                atol=1e-12,
            )
        for name, arrays in fields.cell_data.items():
            for a, b in zip(arrays, reference_fields.cell_data[name], strict=True):
                assert_allclose(a, b, rtol=1e-12, atol=1e-12)
        assert np.isfinite(fields.points).all()
    original_solve = Model.solve
    calls = []

    def counted(model):
        calls.append(model)
        return original_solve(model)

    monkeypatch.setattr(Model, "solve", counted)
    assert study.run(output) == catalog
    assert not calls
    broken = output / catalog["cases"][1]["result"]["manifest"]
    broken.write_text("{}")
    repaired = study.run(output)
    assert len(calls) == 1
    assert repaired["cases"][0] == catalog["cases"][0]
    assert repaired["cases"][1]["attempts"] == 2
    assert repaired["cases"][1]["history"][-1]["status"] == "invalid"


def test_interruption_during_export_and_failure_are_resumable(
    tmp_path: Path, monkeypatch
) -> None:
    import fem_toolbox.studies as module

    manifest = inputs(tmp_path)
    manifest["parameters"][0]["values"] = [2.0, 4.0, -1.0]
    study = Study(manifest, root=tmp_path)
    output = tmp_path / "results"
    export = module.export_result_bundle
    count = 0

    def interrupted(prefix, mesh, run):
        nonlocal count
        count += 1
        if count == 2:
            prefix.with_suffix(".vtu").write_text("partial")
            raise KeyboardInterrupt
        return export(prefix, mesh, run)

    monkeypatch.setattr(module, "export_result_bundle", interrupted)
    with pytest.raises(KeyboardInterrupt):
        study.run(output)
    partial = json.loads((output / "catalog.json").read_text())
    assert [c["status"] for c in partial["cases"]] == [
        "succeeded",
        "interrupted",
        "pending",
    ]
    monkeypatch.setattr(module, "export_result_bundle", export)
    done = study.run(output)
    assert [c["status"] for c in done["cases"]] == ["succeeded", "succeeded", "failed"]
    assert [c["attempts"] for c in done["cases"]] == [1, 2, 1]
    assert done["cases"][2]["error"]
    assert partial["cases"][0] == done["cases"][0]
    assert (
        output / partial["cases"][1]["attempt_directory"] / "result.vtu"
    ).read_text() == "partial"


def test_changed_input_invalidates_and_removed_cases_are_retained(
    tmp_path: Path,
) -> None:
    manifest = inputs(tmp_path)
    study = Study(manifest, root=tmp_path)
    output = tmp_path / "results"
    old = study.run(output, max_cases=1)
    assert [c["status"] for c in old["cases"]] == ["succeeded", "pending"]
    model = tmp_path / "model.toml"
    model.write_text(model.read_text().replace("300.0", "301.0"))
    manifest["parameters"][0]["values"] = [2.0]
    new = Study(manifest, root=tmp_path).run(output)
    assert new["cases"][0]["case_id"] == old["cases"][0]["case_id"]
    assert new["cases"][0]["attempts"] == 2
    assert new["cases"][0]["history"][-1]["status"] == "invalid"
    assert len(new["retired_cases"]) == 1


def test_product_units_and_relocated_study(tmp_path: Path) -> None:
    from shutil import copytree

    original = tmp_path / "original"
    original.mkdir()
    manifest = inputs(original)
    manifest["parameters"].append(
        {
            "name": "temperature",
            "target": "boundary.dirichlet.0.value",
            "unit": "K",
            "values": [
                [300.0, 310.0, 300.0, 300.0, 310.0],
                [301.0, 311.0, 301.0, 301.0, 311.0],
            ],
        }
    )
    first = Study(manifest, root=original).resolve()
    assert [c["parameters"]["k"]["value"] for c in first] == [2.0, 2.0, 4.0, 4.0]
    relocated = tmp_path / "relocated"
    copytree(original, relocated)
    assert Study(manifest, root=relocated).resolve() == first


def test_mesh_content_and_scale_invalidate_fingerprints(tmp_path: Path) -> None:
    study = Study(inputs(tmp_path), root=tmp_path)
    first = study.resolve()
    mesh = tmp_path / "mesh.msh"
    mesh.write_text(mesh.read_text() + "\n")
    changed = study.resolve()
    assert changed[0]["case_id"] == first[0]["case_id"]
    assert changed[0]["fingerprints"]["mesh"] != first[0]["fingerprints"]["mesh"]
    model = tmp_path / "model.toml"
    model.write_text(
        model.read_text().replace(
            "length_scale_to_m = 1.0", "length_scale_to_m = 0.001"
        )
    )
    assert (
        study.resolve()[0]["fingerprints"]["input"]
        != changed[0]["fingerprints"]["input"]
    )


@pytest.mark.parametrize("damage", ["missing", "vtu", "running"])
def test_resume_damaged_or_crashed_case(tmp_path: Path, damage: str) -> None:
    study = Study(inputs(tmp_path), root=tmp_path)
    output = tmp_path / "results"
    first = study.run(output)
    record = first["cases"][0]
    if damage == "missing":
        (output / record["result"]["manifest"]).unlink()
    elif damage == "vtu":
        (output / record["result"]["fields"]).write_text("broken")
    else:
        record["status"] = "running"
        (output / "catalog.json").write_text(json.dumps(first))
    second = study.run(output)
    assert [c["attempts"] for c in second["cases"]] == [2, 1]
    assert all(c["status"] == "succeeded" for c in second["cases"])


def test_catalog_corruption_and_concurrent_writer_fail_without_overwrite(
    tmp_path: Path,
) -> None:
    study = Study(inputs(tmp_path), root=tmp_path)
    output = tmp_path / "results"
    study.run(output, max_cases=0)
    path = output / "catalog.json"
    before = path.read_bytes()
    lock = output / ".study.lock"
    lock.write_text("other writer")
    with pytest.raises(FileExistsError, match="locked"):
        study.run(output)
    assert path.read_bytes() == before
    lock.unlink()
    path.write_text("{bad")
    with pytest.raises(ModelSchemaError, match="catalog"):
        study.run(output)
    assert path.read_text() == "{bad"
    assert not lock.exists()


@pytest.mark.parametrize(
    "target", ["boundary.dirichlet.-1.value", "boundary.dirichlet.00.value"]
)
def test_aliased_array_targets_rejected(tmp_path: Path, target: str) -> None:
    manifest = inputs(tmp_path)
    manifest["parameters"] = [
        {
            "name": "temperature",
            "target": target,
            "unit": "K",
            "values": [[300.0, 310.0, 300.0, 300.0, 310.0]],
        }
    ]
    with pytest.raises(ModelSchemaError):
        Study(manifest, root=tmp_path).resolve()


def test_paths_and_nonstationary_rejected(tmp_path: Path) -> None:
    manifest = inputs(tmp_path)
    manifest["model"] = "../model.toml"
    with pytest.raises(ModelSchemaError, match="escapes"):
        Study(manifest, root=tmp_path).resolve()
    manifest["model"] = "model.toml"
    model = tmp_path / "model.toml"
    model.write_text(
        model.read_text().replace('physics = "heat"', 'physics = "transient_heat"')
    )
    with pytest.raises(ModelSchemaError, match="stationary"):
        Study(manifest, root=tmp_path).resolve()


def test_atomic_catalog_failure_preserves_previous(tmp_path: Path, monkeypatch) -> None:
    import fem_toolbox.studies as module

    study = Study(inputs(tmp_path), root=tmp_path)
    output = tmp_path / "results"
    study.run(output, max_cases=0)
    before = (output / "catalog.json").read_bytes()

    def fail_replace(*args):
        raise OSError("injected disk failure")

    monkeypatch.setattr(module.os, "replace", fail_replace)
    with pytest.raises(OSError, match="disk failure"):
        study.run(output)
    assert (output / "catalog.json").read_bytes() == before
    assert not (output / ".study.lock").exists()


def test_coupled_solver_settings_and_mesh_variants(tmp_path: Path) -> None:
    from shutil import copytree

    from fem_toolbox.io import export_result_bundle
    from fem_toolbox.model import Model

    source = Path(__file__).parents[2] / "examples" / "thermopiezoelectric"
    copytree(source, tmp_path, dirs_exist_ok=True)
    model = Model.from_file(tmp_path / "case.toml")
    mesh_name = model.to_document()["mesh"]["path"]
    copyfile(tmp_path / mesh_name, tmp_path / "variant.msh")
    manifest = {
        "study_schema_version": 1,
        "model": "case.toml",
        "parameters": [
            {
                "name": "scaled",
                "target": "solver.block_scaling",
                "unit": "1",
                "values": [True, False],
            },
            {
                "name": "geometry",
                "target": "mesh.path",
                "unit": "1",
                "values": [mesh_name, "variant.msh"],
            },
        ],
    }
    study = Study(manifest, root=tmp_path)
    cases = study.resolve()
    assert cases[0]["fingerprints"]["solver"] != cases[2]["fingerprints"]["solver"]
    assert cases[0]["fingerprints"] == cases[1]["fingerprints"]
    result = study.run(tmp_path / "out")
    assert all(c["status"] == "succeeded" for c in result["cases"])
    for i, scaled in enumerate([True, False]):
        document = model.to_document()
        document["solver"]["block_scaling"] = scaled
        direct = Model(mesh=model.mesh, document=document).solve()
        _, ref = export_result_bundle(tmp_path / f"ref-{i}", model.mesh, direct)
        actual = tmp_path / "out" / result["cases"][2 * i]["result"]["manifest"]
        assert (
            json.loads(actual.read_text())["metrics"]
            == json.loads(ref.read_text())["metrics"]
        )


def test_tet10_bundle_reused(tmp_path: Path) -> None:
    from shutil import copytree

    source = Path(__file__).parents[2] / "examples" / "tet10_heat"
    copytree(source, tmp_path, dirs_exist_ok=True)
    study = Study(
        {"study_schema_version": 1, "model": "stationary.toml", "parameters": []},
        root=tmp_path,
    )
    first = study.run(tmp_path / "out")
    assert first["cases"][0]["status"] == "succeeded"
    result = json.loads(
        (tmp_path / "out" / first["cases"][0]["result"]["manifest"]).read_text()
    )
    assert result["bundle_schema_version"] == 2
    assert study.run(tmp_path / "out") == first


def test_toolbox_revision_invalidates(tmp_path: Path, monkeypatch) -> None:
    import fem_toolbox.studies as module

    study = Study(inputs(tmp_path), root=tmp_path)
    first = study.resolve()
    monkeypatch.setattr(module, "_STUDY_REVISION", 999)
    second = study.resolve()
    assert first[0]["case_id"] == second[0]["case_id"]
    assert first[0]["fingerprints"]["input"] != second[0]["fingerprints"]["input"]


def test_manifest_detached_and_invalid_limits(tmp_path: Path) -> None:
    manifest = inputs(tmp_path)
    study = Study(manifest, root=tmp_path)
    manifest["parameters"].clear()
    assert len(study.resolve()) == 2
    for limit in (-1, True, 1.5):
        with pytest.raises(ModelSchemaError, match="max_cases"):
            study.run(tmp_path / "out", max_cases=limit)


@pytest.mark.parametrize("bad", ["physics", "mesh", "material_units"])
def test_malformed_base_tables_report_schema_error(tmp_path: Path, bad: str) -> None:
    manifest = inputs(tmp_path)
    model = tmp_path / "model.toml"
    content = model.read_text()
    if bad == "physics":
        content = content.replace('physics = "heat"', 'physics = ["heat"]')
    elif bad == "mesh":
        content = content.replace(
            '[mesh]\npath = "mesh.msh"\nlength_scale_to_m = 1.0', "mesh = []"
        )
    else:
        content = content.replace('unit = "W/(m K)"', 'unit = ["W/(m K)"]')
    model.write_text(content)
    with pytest.raises(ModelSchemaError):
        Study(manifest, root=tmp_path).resolve()
