"""Versioned TOML model loading, strict validation and solver dispatch."""

from pathlib import Path
from shutil import copyfile

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox import __version__
from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model

FIXTURE = Path(__file__).parents[1] / "fixtures" / "overlapping_groups.msh"


def _heat_file(tmp_path: Path, *, extra: str = "", unit: str = "SI") -> Path:
    copyfile(FIXTURE, tmp_path / "mesh.msh")
    path = tmp_path / "case.toml"
    path.write_text(
        f'''schema_version = 1
units = "{unit}"
physics = "heat"
{extra}
[mesh]
path = "mesh.msh"
length_scale_to_m = 1.0

[materials.T.ceramic]
conductivity = 2.0
unit = "W/(m K)"

[regions.T]
ceramic = "ceramic"

[[boundary.dirichlet]]
field = "T"
nodes = [0, 1, 2, 3, 4]
value = 300.0
unit = "K"
''',
        encoding="utf-8",
    )
    return path


def test_load_solve_and_roundtrip_document(tmp_path: Path) -> None:
    model = Model.from_file(_heat_file(tmp_path))
    result = model.solve()
    assert_allclose(result.field_result.temperature, 300, atol=1e-12, rtol=0)
    assert result.provenance.solver_parameters["model_schema_version"] == 1
    assert result.provenance.toolbox_version == __version__
    replay = Model(mesh=model.mesh, document=model.to_document()).solve()
    assert replay.provenance.input_sha256 == result.provenance.input_sha256
    assert_allclose(replay.field_result.temperature, result.field_result.temperature)


@pytest.mark.parametrize(
    ("extra", "unit", "expected"),
    [
        ("surprise = 1", "SI", "surprise"),
        ("", "mm", "units"),
    ],
)
def test_unknown_key_and_unit_are_rejected(
    tmp_path: Path, extra: str, unit: str, expected: str
) -> None:
    with pytest.raises(ModelSchemaError, match=expected):
        Model.from_file(_heat_file(tmp_path, extra=extra, unit=unit))


def test_relative_mesh_path_cannot_escape_model_directory(tmp_path: Path) -> None:
    path = _heat_file(tmp_path)
    path.write_text(
        path.read_text().replace('path = "mesh.msh"', 'path = "../mesh.msh"')
    )
    with pytest.raises(ModelSchemaError, match=r"mesh\.path"):
        Model.from_file(path)


def test_python_model_runs_electrostatics_and_detects_bad_field_unit() -> None:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        volume_sets={"dielectric": np.array([0])},
    )
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "electrostatics",
        "materials": {"phi": {"eps": {"permittivity": 2.0, "unit": "F/m"}}},
        "regions": {"phi": {"dielectric": "eps"}},
        "boundary": {
            "dirichlet": [
                {"field": "phi", "nodes": [0, 2, 3], "value": 0.0, "unit": "V"},
                {"field": "phi", "nodes": [1], "value": 1.0, "unit": "V"},
            ]
        },
    }
    result = Model(mesh=mesh, document=document).solve()
    assert_allclose(result.field_result.potential, [0, 1, 0, 0], atol=1e-13)
    document["boundary"]["dirichlet"][0]["unit"] = "mV"
    with pytest.raises(ModelSchemaError, match=r"boundary\.dirichlet.*unit"):
        Model(mesh=mesh, document=document)


def _cube_regions(cube_mesh):
    original = cube_mesh(1)
    left = np.flatnonzero(np.all(original.points[original.triangle3, 0] == 0, axis=1))
    right = np.flatnonzero(np.all(original.points[original.triangle3, 0] == 1, axis=1))
    return Mesh(
        points=original.points,
        tetra4=original.tetra4,
        triangle3=original.triangle3,
        volume_sets={"body": np.arange(original.n_tetra4)},
        facet_sets={"left": left, "right": right},
    )


def test_heat_model_uses_named_flux_and_dirichlet_sets(cube_mesh) -> None:
    mesh = _cube_regions(cube_mesh)
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "heat",
        "materials": {"T": {"k": {"conductivity": 2.0, "unit": "W/(m K)"}}},
        "regions": {"T": {"body": "k"}},
        "boundary": {
            "dirichlet": [{"field": "T", "set": "left", "value": 300.0, "unit": "K"}],
            "flux": [{"set": "right", "value": 2.0, "unit": "W/m^2"}],
        },
    }
    result = Model(mesh=mesh, document=document).solve()
    assert_allclose(
        result.field_result.temperature, 300 - mesh.points[:, 0], atol=2e-12, rtol=0
    )
    assert abs(result.field_result.power_balance) < 1e-11


def test_piezo_model_routes_floating_electrode_and_options(cube_mesh) -> None:
    mesh = _cube_regions(cube_mesh)
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "piezoelectric",
        "materials": {
            "piezo": {
                "p": {
                    "stiffness": (np.eye(6) * 10).tolist(),
                    "stress_coupling": np.zeros((3, 6)).tolist(),
                    "permittivity": (np.eye(3) * 2).tolist(),
                    "units": {
                        "stiffness": "Pa",
                        "stress_coupling": "C/m^2",
                        "permittivity": "F/m",
                    },
                }
            }
        },
        "regions": {"piezo": {"body": "p"}},
        "boundary": {
            "dirichlet": [
                {
                    "field": "u",
                    "nodes": list(range(mesh.n_points)),
                    "value": 0.0,
                    "unit": "m",
                },
                {"field": "phi", "set": "left", "value": 0.0, "unit": "V"},
            ],
            "electrodes": [
                {"name": "right", "set": "right", "charge": 2.0, "charge_unit": "C"}
            ],
        },
        "solver": {"method": "direct", "block_scaling": True},
    }
    result = Model(mesh=mesh, document=document).solve()
    assert_allclose(
        result.field_result.potential[mesh.points[:, 0] == 1], 1, atol=2e-12, rtol=0
    )
    assert_allclose(
        result.field_result.electrode_states[0].charge, 2, atol=2e-12, rtol=0
    )
    replay = Model(
        mesh=mesh, document=Model(mesh=mesh, document=document).to_document()
    ).solve()
    assert replay.provenance.input_sha256 == result.provenance.input_sha256


def test_elasticity_model_routes_material_and_displacements(cube_mesh) -> None:
    mesh = _cube_regions(cube_mesh)
    values = np.column_stack(
        [0.1 * mesh.points[:, 0], np.zeros(mesh.n_points), np.zeros(mesh.n_points)]
    ).ravel()
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "elasticity",
        "materials": {"u": {"solid": {"young": 10.0, "poisson": 0.0, "unit": "Pa"}}},
        "regions": {"u": {"body": "solid"}},
        "boundary": {
            "dirichlet": [
                {
                    "field": "u",
                    "nodes": list(range(mesh.n_points)),
                    "value": values.tolist(),
                    "unit": "m",
                }
            ]
        },
    }
    result = Model(mesh=mesh, document=document).solve()
    assert_allclose(
        result.field_result.displacement.ravel(), values, atol=1e-13, rtol=0
    )


def test_schema_version_and_nested_unknown_key_are_concrete(tmp_path: Path) -> None:
    source = _heat_file(tmp_path)
    for replacement, expected in [
        (("schema_version = 1", "schema_version = 2"), "schema_version"),
        (('unit = "W/(m K)"', 'unit = "W/(m K)"\nsecret = 1'), "materials.T.ceramic"),
    ]:
        text = source.read_text()
        source.write_text(text.replace(*replacement))
        with pytest.raises(ModelSchemaError, match=expected):
            Model.from_file(source)
        source.write_text(text)


def test_invalid_field_and_empty_named_boundary_are_schema_errors(cube_mesh) -> None:
    mesh = _cube_regions(cube_mesh)
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "heat",
        "materials": {"T": {"k": {"conductivity": 2.0, "unit": "W/(m K)"}}},
        "regions": {"T": {"body": "k"}},
        "boundary": {
            "dirichlet": [{"field": [], "set": "left", "value": 300.0, "unit": "K"}]
        },
    }
    with pytest.raises(ModelSchemaError, match=r"boundary\.dirichlet.*field"):
        Model(mesh=mesh, document=document)
    document["boundary"]["dirichlet"][0]["field"] = "T"
    mesh_empty = Mesh(
        points=mesh.points,
        tetra4=mesh.tetra4,
        triangle3=mesh.triangle3,
        volume_sets=mesh.volume_sets,
        facet_sets={"left": np.array([], dtype=int)},
    )
    with pytest.raises(ModelSchemaError, match="empty"):
        Model(mesh=mesh_empty, document=document)
