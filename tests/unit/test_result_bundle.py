"""Stable, raw SI result bundle contract."""

import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.io import export_result_bundle
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
)
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics
from fem_toolbox.solvers.thermopiezoelectric import solve_thermopiezoelectric

CASE = Path(__file__).parents[2] / "examples" / "thermopiezoelectric" / "case.toml"


def test_model_bundle_roundtrip_has_raw_fields_and_signed_metrics(
    tmp_path: Path,
) -> None:
    model = Model.from_file(CASE)
    solved = model.solve()
    fields, manifest = export_result_bundle(tmp_path / "result", model.mesh, solved)
    data = json.loads(manifest.read_text(encoding="utf-8"))
    visual = meshio.read(fields)
    assert data["bundle_schema_version"] == 1
    assert data["toolbox_version"]
    assert data["physics"] == "thermopiezoelectric"
    assert data["provenance"]["input_sha256"] == solved.provenance.input_sha256
    assert data["fields"]["point"]["temperature_K"]["kind"] == "raw"
    assert data["fields"]["tetra4"]["stress_Pa"]["components"] == [
        "xx",
        "yy",
        "zz",
        "yz",
        "xz",
        "xy",
    ]
    assert_array_equal(
        visual.point_data["temperature_K"], solved.field_result.heat.temperature
    )
    assert_array_equal(
        visual.point_data["displacement_m"],
        solved.field_result.piezoelectric.displacement,
    )
    assert_array_equal(
        visual.cell_data_dict["stress_Pa"]["tetra"],
        solved.field_result.piezoelectric.stress,
    )
    assert data["metrics"]["piezoelectric"]["total_enthalpy_J"] == pytest.approx(
        solved.field_result.piezoelectric.total_enthalpy
    )
    assert data["metrics"]["piezoelectric"]["charge_balance_C"] == pytest.approx(
        solved.field_result.piezoelectric.charge_balance
    )
    assert data["limitations"]
    assert_allclose(visual.point_data["potential_V"][1], 1 / 4.4, atol=1e-12, rtol=0)
    assert_allclose(
        visual.cell_data_dict["electric_displacement_C_m2"]["tetra"][0, 0],
        -1,
        atol=1e-12,
        rtol=0,
    )
    assert abs(data["metrics"]["piezoelectric"]["charge_balance_C"]) < 1e-10
    reader = Path(__file__).parents[2] / "scripts" / "read_result_bundle.py"
    inspected = subprocess.run(
        [sys.executable, "-I", str(reader), str(manifest)],
        cwd=tmp_path,
        check=True,
        capture_output=True,
        text=True,
    )
    summary = json.loads(inspected.stdout)
    assert summary["input_sha256"] == solved.provenance.input_sha256
    assert summary["metrics"]["piezoelectric"]["charge_balance_C"] == pytest.approx(
        solved.field_result.piezoelectric.charge_balance
    )
    fields.write_bytes(fields.read_bytes() + b"corrupt")
    corrupted = subprocess.run(
        [sys.executable, "-I", str(reader), str(manifest)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert corrupted.returncode != 0
    assert "VTU checksum mismatch" in corrupted.stderr


def test_bundle_rejects_existing_files_without_modifying_them(tmp_path: Path) -> None:
    model = Model.from_file(CASE)
    target = tmp_path / "result.json"
    target.write_text("existing", encoding="utf-8")
    with pytest.raises(ValueError, match="exists"):
        export_result_bundle(tmp_path / "result", model.mesh, model.solve())
    assert target.read_text(encoding="utf-8") == "existing"


@pytest.mark.parametrize(
    "mode",
    [
        "heat",
        "elasticity",
        "electrostatics",
        "thermomechanics",
        "piezoelectric",
        "thermopiezoelectric",
    ],
)
def test_direct_stationary_results_have_consistent_entity_contract(
    tmp_path: Path, mode: str
) -> None:
    mesh = Mesh(
        points=np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 2, 2]], dtype=float
        ),
        tetra4=np.array([[0, 1, 2, 3]]),
        triangle3=np.array([[1, 2, 3]]),
        volume_sets={"body": np.array([0])},
        facet_sets={"surface": np.array([0])},
    )
    nodes = np.array([0, 1, 2, 3])
    heat_condition = Dirichlet(
        thermal_dofs(mesh).indices("T", nodes), np.full(4, 300.0)
    )
    mechanical_condition = Dirichlet(
        mechanical_dofs(mesh).indices("u", nodes), np.zeros(12)
    )
    electric_condition = Dirichlet(
        electric_dofs(mesh).indices("phi", nodes), np.zeros(4)
    )
    piezo_dofs = piezoelectric_dofs(mesh)
    piezo_conditions = [
        Dirichlet(piezo_dofs.indices("u", nodes), np.zeros(12)),
        Dirichlet(piezo_dofs.indices("phi", nodes), np.zeros(4)),
    ]
    heat = ThermalMaterial.isotropic(2)
    elastic = ElasticMaterial.isotropic(10, 0.25)
    electric = DielectricMaterial.isotropic(4)
    piezo = PiezoelectricMaterial(np.eye(6) * 10, np.zeros((3, 6)), np.eye(3) * 4)
    expansion = ThermalExpansion(np.zeros(6), 300)
    if mode == "heat":
        result = solve_heat(mesh, heat, [heat_condition])
        expected_field = "temperature_K"
    elif mode == "elasticity":
        result = solve_elasticity(mesh, elastic, [mechanical_condition])
        expected_field = "displacement_m"
    elif mode == "electrostatics":
        result = solve_electrostatics(mesh, electric, [electric_condition])
        expected_field = "potential_V"
    elif mode == "thermomechanics":
        result = solve_thermomechanics(
            mesh, heat, elastic, expansion, [heat_condition], [mechanical_condition]
        )
        expected_field = "temperature_K"
    elif mode == "piezoelectric":
        result = solve_piezoelectric(mesh, piezo, piezo_conditions)
        expected_field = "potential_V"
    else:
        result = solve_thermopiezoelectric(
            mesh, heat, piezo, expansion, [heat_condition], piezo_conditions
        )
        expected_field = "temperature_K"
    fields, manifest = export_result_bundle(tmp_path / mode, mesh, result)
    data = json.loads(manifest.read_text())
    visual = meshio.read(fields)
    assert data["physics"] == mode
    assert data["fields"]["point"][expected_field]["active_node_ids"] == nodes.tolist()
    assert visual.point_data[expected_field].shape[0] == 5
    assert visual.point_data[expected_field][-1].sum() == 0
    assert data["fields"]["tetra4"]
    assert "tetra" in visual.cell_data_dict[next(iter(data["fields"]["tetra4"]))]
    assert "triangle" in visual.cell_data_dict[next(iter(data["fields"]["tetra4"]))]


def test_nonfinite_field_fails_before_creating_bundle(tmp_path: Path) -> None:
    model = Model.from_file(CASE)
    solved = model.solve()
    bad_heat = replace(solved.field_result.heat, temperature=np.array([np.nan] * 4))
    bad = replace(solved.field_result, heat=bad_heat)
    with pytest.raises(ValueError, match="temperature_K"):
        export_result_bundle(tmp_path / "bad", model.mesh, bad)
    assert not (tmp_path / "bad.vtu").exists()
    assert not (tmp_path / "bad.json").exists()


def test_incomplete_cell_field_fails_before_creating_bundle(tmp_path: Path) -> None:
    model = Model.from_file(CASE)
    solved = model.solve()
    bad_piezo = replace(solved.field_result.piezoelectric, stress=np.empty((0, 6)))
    bad = replace(solved.field_result, piezoelectric=bad_piezo)
    with pytest.raises(ValueError, match="stress_Pa"):
        export_result_bundle(tmp_path / "bad", model.mesh, bad)
    assert not (tmp_path / "bad.vtu").exists()
    assert not (tmp_path / "bad.json").exists()


def test_vtu_write_error_leaves_no_partial_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import fem_toolbox.io.result_bundle as bundle_module

    model = Model.from_file(CASE)

    def broken_writer(*args: object, **kwargs: object) -> None:
        raise OSError("disk unavailable")

    monkeypatch.setattr(bundle_module, "write_vtu", broken_writer)
    with pytest.raises(OSError, match="disk unavailable"):
        export_result_bundle(tmp_path / "bad", model.mesh, model.solve())
    assert not (tmp_path / "bad.vtu").exists()
    assert not (tmp_path / "bad.json").exists()


def test_manifest_publication_error_rolls_back_vtu(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import fem_toolbox.io.result_bundle as bundle_module

    original = bundle_module.os.link

    def broken_second_link(source: Path, destination: Path) -> None:
        if destination.suffix == ".json":
            raise OSError("manifest unavailable")
        original(source, destination)

    monkeypatch.setattr(bundle_module.os, "link", broken_second_link)
    model = Model.from_file(CASE)
    with pytest.raises(OSError, match="manifest unavailable"):
        export_result_bundle(tmp_path / "bad", model.mesh, model.solve())
    assert not (tmp_path / "bad.vtu").exists()
    assert not (tmp_path / "bad.json").exists()


def test_electrode_charge_and_signed_energy_are_preserved(
    tmp_path: Path, cube_mesh
) -> None:
    mesh = cube_mesh(1)
    left = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 0, axis=1))
    right = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 1, axis=1))
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.5
    material = PiezoelectricMaterial(np.eye(6) * 10, coupling, np.eye(3) * 2)
    displacement = np.column_stack(
        [0.1 * mesh.points[:, 0], np.zeros(mesh.n_points), np.zeros(mesh.n_points)]
    )
    bc = Dirichlet(
        piezoelectric_dofs(mesh).indices("u", np.arange(mesh.n_points)),
        displacement.ravel(),
    )
    solved = solve_piezoelectric(
        mesh,
        material,
        [bc],
        electrodes=[
            Electrode("left", left, potential=0),
            Electrode("right", right, potential=1),
        ],
    )
    _, manifest = export_result_bundle(tmp_path / "electrodes", mesh, solved)
    metrics = json.loads(manifest.read_text())["metrics"]["piezoelectric"]
    by_name = {item["name"]: item for item in metrics["electrodes"]}
    assert by_name["right"]["charge_C"] == pytest.approx(1.95, abs=1e-12)
    assert metrics["total_enthalpy_J"] == pytest.approx(solved.total_enthalpy)
