import copy

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError
from fem_toolbox.materials.catalog import thermal_material
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model


def document():
    return {
        "schema_version": 2,
        "units": "SI",
        "physics": "transient_heat",
        "materials": {"al": {"preset": "aluminium_6060_t6_20c"}},
        "regions": {"body": "al"},
        "time": {"values": [0.0, 1.0], "unit": "s"},
        "initial_temperature": {"value": 300.0, "unit": "K"},
        "loads": {
            "source": {
                "value": {
                    "kind": "step",
                    "time": 0.5,
                    "before": 0.0,
                    "after": 2432700.0,
                },
                "unit": "W/m^3",
            }
        },
        "probes": {"center": {"point": [0.5, 0.5, 0.5]}},
    }


def test_model_preset_and_profile(cube_mesh):
    base = cube_mesh(1)
    mesh = Mesh(
        points=base.points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        volume_sets={"body": np.arange(base.n_tetra4)},
    )
    model = Model(mesh=mesh, document=document())
    result = model.solve()
    assert result.physics == "transient_heat"
    np.testing.assert_allclose(
        result.field_result.frames[-1].temperature, 300.5, atol=1e-10
    )
    assert model.to_document()["materials"]["al"]["preset"] == "aluminium_6060_t6_20c"
    assert thermal_material("aluminium_6060_t6_20c").density == 2700.0
    changed = copy.deepcopy(document())
    changed["schema_version"] = 1
    with pytest.raises(ValidationError):
        Model(mesh=mesh, document=changed)
    changed = copy.deepcopy(document())
    changed["loads"]["source"]["value"]["typo"] = 1
    with pytest.raises(ValidationError):
        Model(mesh=mesh, document=changed)


def test_toml_csv_example_and_portable_document():
    from pathlib import Path

    model = Model.from_file(
        Path(__file__).parents[2] / "examples/transient_heat/model.toml"
    )
    document = model.to_document()
    assert document["loads"]["source"]["value"]["kind"] == "series"
    result = model.solve()
    repeated = Model(mesh=model.mesh, document=document).solve()
    assert result.provenance.to_dict() == repeated.provenance.to_dict()
    assert result.field_result.frames[-1].temperature.min() > 293.15
    assert max(abs(b.balance) for b in result.field_result.balances) < 1e-9


def test_csv_cannot_escape_model_directory(tmp_path):
    from fem_toolbox.transient_model import resolve_csv_profiles

    with pytest.raises(ValidationError, match="escapes"):
        resolve_csv_profiles({"kind": "csv", "path": "../outside.csv"}, tmp_path)


@pytest.mark.parametrize(
    "change",
    [
        {"solver": {"lumped": "false"}},
        {"output": {"every": 0}},
        {"time": {"values": [0, 0], "unit": "s"}},
        {"mesh": {"path": "body.msh", "length_unit": "mm", "length_scale_to_m": 1}},
        {"materials": {"al": {"preset": "aluminium_6060_t6_20c", "density": 1}}},
        {
            "loads": {
                "source": {
                    "value": {
                        "kind": "series",
                        "times": [0, 1],
                        "values": [0, 1],
                        "unit": "W",
                    },
                    "unit": "W/m^3",
                }
            }
        },
    ],
)
def test_model_rejects_ambiguous_or_invalid_data(cube_mesh, change):
    base = cube_mesh(1)
    mesh = Mesh(
        points=base.points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        volume_sets={"body": np.arange(base.n_tetra4)},
    )
    raw = document()
    raw.update(change)
    with pytest.raises(ValidationError):
        Model(mesh=mesh, document=raw)
