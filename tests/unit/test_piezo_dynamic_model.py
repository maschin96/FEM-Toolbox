"""Schema 3 routes existing piezo material/region definitions to dynamic solvers."""

import copy
import tomllib
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.io.meshio_adapter import read_gmsh
from fem_toolbox.model import Model


def case(physics="harmonic_piezoelectric"):
    root = Path(__file__).resolve().parents[2] / "examples/thermopiezoelectric"
    base = tomllib.loads((root / "case.toml").read_text())
    base["physics"] = "piezoelectric"
    base["materials"] = {"piezo": base["materials"]["piezo"]}
    base["regions"] = {"piezo": base["regions"]["piezo"]}
    base["boundary"]["dirichlet"] = [
        b for b in base["boundary"]["dirichlet"] if b["field"] != "T"
    ]
    base["loads"] = {"body_force": {"value": [6, 0, 0], "unit": "N/m^3"}}
    base["solver"] = {"method": "direct", "block_scaling": True}
    base.pop("mesh")
    analysis = {
        "frequencies": {"values": [0, 0.2], "unit": "Hz"},
        "load_factor": {"real": 1, "imag": 2},
    }
    if physics == "piezo_modes":
        base["loads"] = {}
        analysis = {"count": 1}
    elif physics == "transient_piezoelectric":
        analysis = {
            "times": {"values": [0, 0.5, 1], "unit": "s"},
            "load_factor": {"kind": "step", "time": 0.5, "before": 0, "after": 1},
        }
    return read_gmsh(root / "mesh.msh"), {
        "schema_version": 3,
        "units": "SI",
        "physics": physics,
        "base": base,
        "density": {"value": 6, "unit": "kg/m^3"},
        "analysis": analysis,
    }


@pytest.mark.parametrize(
    "physics", ["piezo_modes", "transient_piezoelectric", "harmonic_piezoelectric"]
)
def test_model_routes_and_preserves_document(physics):
    mesh, document = case(physics)
    original = copy.deepcopy(document)
    model = Model(mesh=mesh, document=document)
    result = model.solve()
    assert model.to_document() == original
    assert result.provenance.solver_parameters["model_schema_version"] == 3
    assert result.field_result.mesh is mesh
    if physics == "piezo_modes":
        assert_allclose(
            result.field_result.solution.angular_frequencies,
            [np.sqrt((10 + 4 / 4) / 0.6)],
            rtol=2e-14,
        )
    elif physics == "transient_piezoelectric":
        assert_allclose(result.field_result.solution.displacement[1], 0, atol=1e-14)
    else:
        state = result.field_result.solution[0]
        assert_allclose(
            state.displacement[1, 0], 0.25 * (1 + 2j) / ((10 + 4 / 4) / 6), rtol=3e-14
        )


@pytest.mark.parametrize("bad", ["Hz_missing", "extra", "density"])
def test_model_rejects_invalid_units_or_unknown_fields(bad):
    mesh, document = case()
    if bad == "Hz_missing":
        document["analysis"]["frequencies"]["unit"] = "rad/s"
    elif bad == "extra":
        document["analysis"]["ignored_option"] = True
    else:
        document["density"]["value"] = -1
    with pytest.raises(ModelSchemaError):
        Model(mesh=mesh, document=document)


def test_polynomial_support_model_and_voltage_profile():
    mesh, document = case("transient_piezoelectric")
    document["analysis"]["prescribed"] = [
        {"dofs": [0], "unit": "m", "coefficients": [0, 0, 0.003]}
    ]
    result = Model(mesh=mesh, document=document).solve().field_result.solution
    assert_allclose(result.displacement[:, 0, 0], 0.003 * result.times**2, atol=1e-14)
    assert_allclose(result.velocity[:, 0, 0], 0.006 * result.times, atol=1e-14)
    assert_allclose(result.acceleration[:, 0, 0], 0.006, atol=1e-14)


def test_csv_profile_resolution_in_from_file(tmp_path):
    root = Path(__file__).resolve().parents[2] / "examples/thermopiezoelectric"
    (tmp_path / "mesh.msh").write_bytes((root / "mesh.msh").read_bytes())
    # Reuse the checked example below; CSV samples become part of the detached document.
    example = Path(__file__).resolve().parents[2] / "examples/piezo_dynamics"
    (tmp_path / "model.toml").write_bytes((example / "model.toml").read_bytes())
    (tmp_path / "drive.csv").write_text(
        "time_s,factor\n0,0\n0.2,0\n0.3,1\n0.7,1\n1,0\n"
    )
    model = Model.from_file(tmp_path / "model.toml")
    assert model.to_document()["analysis"]["load_factor"]["kind"] == "series"
    result = model.solve().field_result.solution
    assert len(result.times) > 10
    assert result.displacement[-1, 1, 0] != 0
