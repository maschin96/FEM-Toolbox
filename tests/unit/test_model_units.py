"""Explicit datasheet units must produce the same effective SI material."""

import copy
from pathlib import Path
from shutil import copyfile

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.io import read_gmsh
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model


def _case() -> tuple[Mesh, dict]:
    mesh = Mesh(
        points=np.eye(4, 3),
        tetra4=np.array([[0, 1, 2, 3]]),
        volume_sets={"body": np.array([0])},
    )
    s = np.diag([1 / 10e9] * 3 + [1 / 4e9] * 3)
    d = np.zeros((3, 6))
    d[2, 2] = 200e-12
    k = np.eye(3) * 1000 * 8.8541878188e-12
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "piezoelectric",
        "materials": {
            "piezo": {
                "pzt": {
                    "form": "strain_charge",
                    "compliance": (s * 1e9).tolist(),
                    "strain_coupling": (d * 1e12).tolist(),
                    "permittivity_at_stress": (k / 8.8541878188e-12).tolist(),
                    "units": {
                        "compliance": "1/GPa",
                        "strain_coupling": "pC/N",
                        "permittivity_at_stress": "relative",
                    },
                }
            }
        },
        "regions": {"piezo": {"body": "pzt"}},
        "boundary": {
            "dirichlet": [
                {"field": "u", "nodes": [0, 1, 2, 3], "value": 0, "unit": "m"},
                {"field": "phi", "nodes": [0, 1, 2, 3], "value": 0, "unit": "V"},
            ]
        },
    }
    return mesh, document


def test_datasheet_units_match_independent_si_reference() -> None:
    mesh, document = _case()
    converted = Model(mesh=mesh, document=document)
    reference = copy.deepcopy(document)
    material = reference["materials"]["piezo"]["pzt"]
    material["compliance"] = (np.asarray(material["compliance"]) * 1e-9).tolist()
    material["strain_coupling"] = (
        np.asarray(material["strain_coupling"]) * 1e-12
    ).tolist()
    material["permittivity_at_stress"] = (
        np.asarray(material["permittivity_at_stress"]) * 8.8541878188e-12
    ).tolist()
    material["units"] = {
        "compliance": "1/Pa",
        "strain_coupling": "C/N",
        "permittivity_at_stress": "F/m",
    }
    expected = Model(mesh=mesh, document=reference)
    a = converted.assignment.cell_materials[0]
    b = expected.assignment.cell_materials[0]
    assert_allclose(a.stiffness, b.stiffness, rtol=1e-13)
    assert_allclose(a.stress_coupling, b.stress_coupling, rtol=1e-13)
    assert_allclose(a.permittivity, b.permittivity, rtol=1e-13)
    assert_allclose(
        converted.solve().field_result.potential,
        expected.solve().field_result.potential,
        atol=1e-14,
    )
    data = converted.solve().provenance.solver_parameters["material_inputs"]["piezo"][
        "pzt"
    ]
    assert data["units"]["strain_coupling"] == "pC/N"
    assert data["si_factors"]["strain_coupling"] == 1e-12
    effective = converted.solve().provenance.solver_parameters["effective_materials"]
    assert effective["piezo"]["catalog"]["pzt"]["shear_strain"] == "engineering"


@pytest.mark.parametrize(
    "key,unit",
    [
        ("compliance", "Pa"),
        ("strain_coupling", "C/m^2"),
        ("permittivity_at_stress", ""),
    ],
)
def test_incompatible_datasheet_unit_has_field_path(key: str, unit: str) -> None:
    mesh, document = _case()
    document["materials"]["piezo"]["pzt"]["units"][key] = unit
    with pytest.raises(ModelSchemaError, match=rf"materials\.piezo\.pzt\.units\.{key}"):
        Model(mesh=mesh, document=document)


def test_absolute_celsius_and_temperature_difference_are_distinct() -> None:
    mesh, document = _case()
    document["physics"] = "thermopiezoelectric"
    document["materials"]["T"] = {"heat": {"conductivity": 2, "unit": "W/(m K)"}}
    document["materials"]["expansion"] = {
        "thermal": {
            "alpha": [1e-5, 0, 0, 0, 0, 0],
            "alpha_unit": "1/delta_degC",
            "reference_temperature": 20,
            "reference_temperature_unit": "degC",
        }
    }
    document["regions"]["T"] = {"body": "heat"}
    document["regions"]["expansion"] = {"body": "thermal"}
    document["boundary"]["dirichlet"].append(
        {"field": "T", "nodes": [0, 1, 2, 3], "value": 30, "unit": "degC"}
    )
    model = Model(mesh=mesh, document=document)
    assert model.assignments["expansion"].cell_materials[
        0
    ].reference_temperature == pytest.approx(293.15)
    assert_allclose(model.solve().field_result.heat.temperature, 303.15, atol=1e-12)
    document["materials"]["expansion"]["thermal"]["alpha_unit"] = "1/degC"
    with pytest.raises(ModelSchemaError, match="alpha_unit"):
        Model(mesh=mesh, document=document)


def test_explicit_geometry_unit_must_match_scale() -> None:
    mesh, document = _case()
    document["mesh"] = {"length_unit": "mm", "length_scale_to_m": 1}
    with pytest.raises(ModelSchemaError, match=r"mesh\.length_scale_to_m"):
        Model(mesh=mesh, document=document)


def test_local_polarization_axis_rotates_into_global_frame() -> None:
    mesh, document = _case()
    rotation = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
    document["materials"]["piezo"]["pzt"]["axes"] = rotation
    material = Model(mesh=mesh, document=document).assignment.cell_materials[0]
    assert_allclose(material.axes, rotation, atol=0)
    assert material.stress_coupling[0, 0] != 0
    provenance = material.to_provenance_dict()
    assert provenance["voigt_order"] == ["xx", "yy", "zz", "yz", "xz", "xy"]
    assert provenance["shear_strain"] == "engineering"


def test_stress_charge_gpa_and_relative_permittivity_match_si() -> None:
    mesh, document = _case()
    material = document["materials"]["piezo"]["pzt"]
    material.clear()
    material.update(
        {
            "stiffness": (np.eye(6) * 10).tolist(),
            "stress_coupling": np.zeros((3, 6)).tolist(),
            "permittivity": (np.eye(3) * 1000).tolist(),
            "units": {
                "stiffness": "GPa",
                "stress_coupling": "C/m^2",
                "permittivity": "relative",
            },
        }
    )
    converted = Model(mesh=mesh, document=document)
    reference = copy.deepcopy(document)
    si = reference["materials"]["piezo"]["pzt"]
    si["stiffness"] = (np.eye(6) * 10e9).tolist()
    si["permittivity"] = (np.eye(3) * 1000 * 8.8541878188e-12).tolist()
    si["units"]["stiffness"] = "Pa"
    si["units"]["permittivity"] = "F/m"
    expected = Model(mesh=mesh, document=reference)
    assert_allclose(
        converted.assignment.cell_materials[0].stiffness,
        expected.assignment.cell_materials[0].stiffness,
        rtol=1e-15,
    )
    assert_allclose(
        converted.assignment.cell_materials[0].permittivity,
        expected.assignment.cell_materials[0].permittivity,
        rtol=1e-15,
    )


def test_mm_mesh_file_matches_independently_scaled_reference(tmp_path: Path) -> None:
    source = Path(__file__).parents[1] / "fixtures" / "overlapping_groups.msh"
    copyfile(source, tmp_path / "mesh.msh")
    template = """schema_version = 1
units = "SI"
physics = "heat"
[mesh]
path = "mesh.msh"
length_unit = "mm"
length_scale_to_m = 0.001
[materials.T.solid]
conductivity = 2.0
unit = "W/(m K)"
[regions.T]
ceramic = "solid"
[[boundary.dirichlet]]
field = "T"
nodes = [0, 1, 2, 3, 4]
value = 300.0
unit = "K"
"""
    path = tmp_path / "case.toml"
    path.write_text(template)
    model = Model.from_file(path)
    assert_allclose(
        model.mesh.points, read_gmsh(tmp_path / "mesh.msh").points * 0.001, atol=0
    )
    assert_allclose(model.solve().field_result.temperature, 300, atol=1e-12)
