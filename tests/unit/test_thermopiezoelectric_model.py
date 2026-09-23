"""Declarative end-to-end T -> (u, phi) case and schema errors."""

import json
from pathlib import Path
from shutil import copyfile

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model


def _case() -> tuple[Mesh, dict]:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        volume_sets={"body": np.array([0])},
    )
    e = np.zeros((3, 6))
    e[0, 0] = 2
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "thermopiezoelectric",
        "materials": {
            "T": {"heat": {"conductivity": 2, "unit": "W/(m K)"}},
            "piezo": {
                "solid": {
                    "stiffness": (np.eye(6) * 10).tolist(),
                    "stress_coupling": e.tolist(),
                    "permittivity": (np.eye(3) * 4).tolist(),
                    "units": {
                        "stiffness": "Pa",
                        "stress_coupling": "C/m^2",
                        "permittivity": "F/m",
                    },
                }
            },
            "expansion": {
                "solid": {
                    "alpha": [0.01, 0, 0, 0, 0, 0],
                    "alpha_unit": "1/K",
                    "reference_temperature": 300,
                    "reference_temperature_unit": "K",
                }
            },
        },
        "regions": {
            "T": {"body": "heat"},
            "piezo": {"body": "solid"},
            "expansion": {"body": "solid"},
        },
        "boundary": {
            "dirichlet": [
                {"field": "T", "nodes": [0, 1, 2, 3], "value": 310, "unit": "K"},
                {"field": "u", "nodes": [0, 1, 2, 3], "value": 0, "unit": "m"},
                {"field": "phi", "nodes": [0, 2, 3], "value": 0, "unit": "V"},
            ]
        },
        "solver": {"electric_thermal_strain": True, "block_scaling": True},
    }
    return mesh, document


def test_model_solves_and_reports_three_field_balances() -> None:
    mesh, document = _case()
    model = Model(mesh=mesh, document=document)
    run = model.solve()
    assert run.physics == "thermopiezoelectric"
    assert_allclose(run.field_result.heat.temperature, 310, atol=1e-12, rtol=0)
    assert_allclose(run.field_result.piezoelectric.potential[1], -0.05, atol=1e-12)
    assert abs(run.field_result.heat.power_balance) < 1e-12
    assert np.linalg.norm(run.field_result.piezoelectric.force_balance) < 1e-12
    assert np.linalg.norm(run.field_result.piezoelectric.moment_balance) < 1e-12
    assert abs(run.field_result.piezoelectric.charge_balance) < 1e-12
    replay = Model(mesh=mesh, document=model.to_document()).solve()
    assert replay.provenance.input_sha256 == run.provenance.input_sha256


def test_expansion_unit_and_missing_region_are_rejected() -> None:
    mesh, document = _case()
    document["materials"]["expansion"]["solid"]["alpha_unit"] = "1/C"
    with pytest.raises(ModelSchemaError, match="alpha_unit"):
        Model(mesh=mesh, document=document)
    document["materials"]["expansion"]["solid"]["alpha_unit"] = "1/K"
    del document["regions"]["expansion"]
    with pytest.raises(ModelSchemaError, match="expansion"):
        Model(mesh=mesh, document=document)


def test_principal_expansion_and_thermal_option_validation() -> None:
    mesh, document = _case()
    expansion = document["materials"]["expansion"]["solid"]
    del expansion["alpha"]
    expansion["principal_alpha"] = [0.01, 0, 0]
    expansion["axes"] = np.eye(3).tolist()
    result = Model(mesh=mesh, document=document).solve()
    assert_allclose(result.field_result.piezoelectric.potential[1], -0.05, atol=1e-12)
    document["solver"]["electric_thermal_strain"] = "yes"
    with pytest.raises(ModelSchemaError, match=r"solver\.electric_thermal_strain"):
        Model(mesh=mesh, document=document)


def test_model_rejects_material_gaps_and_conflicts_before_solve() -> None:
    mesh, document = _case()
    document["regions"]["T"] = {}
    with pytest.raises(ModelSchemaError, match=r"regions\.T.*missing"):
        Model(mesh=mesh, document=document)

    mesh = Mesh(
        points=mesh.points,
        tetra4=mesh.tetra4,
        volume_sets={"body": np.array([0]), "overlap": np.array([0])},
    )
    document["regions"]["T"] = {"body": "heat", "overlap": "other"}
    document["materials"]["T"]["other"] = {
        "conductivity": 3,
        "unit": "W/(m K)",
    }
    with pytest.raises(ModelSchemaError, match=r"regions\.T.*conflict"):
        Model(mesh=mesh, document=document)


def test_bundle_reports_material_axes_and_thermal_solver_option(tmp_path: Path) -> None:
    mesh, document = _case()
    expansion = document["materials"]["expansion"]["solid"]
    del expansion["alpha"]
    expansion["principal_alpha"] = [0.01, 0, 0]
    expansion["axes"] = np.eye(3).tolist()
    solved = Model(mesh=mesh, document=document).solve()
    _, manifest = export_result_bundle(tmp_path / "result", mesh, solved)
    report = json.loads(manifest.read_text(encoding="utf-8"))
    parameters = report["provenance"]["solver"]["parameters"]
    assert (
        parameters["effective_materials"]["expansion"]["catalog"]["solid"][
            "axes_local_to_global"
        ]
        == np.eye(3).tolist()
    )
    assert (
        report["solver_provenance"]["piezoelectric"]["solver"]["parameters"][
            "electric_thermal_strain"
        ]
        is True
    )


def test_from_toml_file_transfers_temperature_on_imported_mesh(tmp_path: Path) -> None:
    source = Path(__file__).parents[1] / "fixtures" / "overlapping_groups.msh"
    copyfile(source, tmp_path / "mesh.msh")
    _, document = _case()
    solid = document["materials"]["piezo"]["solid"]
    path = tmp_path / "case.toml"
    path.write_text(
        f"""schema_version = 1
units = "SI"
physics = "thermopiezoelectric"
[mesh]
path = "mesh.msh"
length_scale_to_m = 1.0
[materials.T.heat]
conductivity = 2.0
unit = "W/(m K)"
[materials.piezo.solid]
stiffness = {json.dumps(solid["stiffness"])}
stress_coupling = {json.dumps(solid["stress_coupling"])}
permittivity = {json.dumps(solid["permittivity"])}
[materials.piezo.solid.units]
stiffness = "Pa"
stress_coupling = "C/m^2"
permittivity = "F/m"
[materials.expansion.solid]
alpha = [0.01, 0, 0, 0, 0, 0]
alpha_unit = "1/K"
reference_temperature = 300.0
reference_temperature_unit = "K"
[regions.T]
ceramic = "heat"
[regions.piezo]
ceramic = "solid"
[regions.expansion]
ceramic = "solid"
[[boundary.dirichlet]]
field = "T"
nodes = [0, 1, 2, 3, 4]
value = 310.0
unit = "K"
[[boundary.dirichlet]]
field = "u"
nodes = [0, 1, 2, 3, 4]
value = 0.0
unit = "m"
[[boundary.dirichlet]]
field = "phi"
nodes = [0, 2, 3]
value = 0.0
unit = "V"
[solver]
electric_thermal_strain = true
""",
        encoding="utf-8",
    )
    result = Model.from_file(path).solve().field_result
    assert_allclose(result.heat.temperature, 310, atol=1e-12, rtol=0)
    assert_allclose(
        result.piezoelectric.potential,
        -0.05 * Model.from_file(path).mesh.points[result.piezoelectric.node_ids, 0],
        atol=1e-12,
        rtol=0,
    )
