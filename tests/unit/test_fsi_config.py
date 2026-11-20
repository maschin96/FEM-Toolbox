"""Strict fsi_schema_version=1 TOML: mandatory values, paths and failures (#93)."""

import json
import shutil
from pathlib import Path

import pytest

from fem_toolbox.coupling import FsiCouplingError
from fem_toolbox.coupling.fsi_config import load_fsi_config, run_fsi_config
from fem_toolbox.errors import ModelSchemaError

EXAMPLE = Path(__file__).parents[2] / "examples/fsi_lid"


@pytest.fixture
def case(tmp_path):
    target = tmp_path / "case"
    shutil.copytree(EXAMPLE / "laminar", target)
    return target / "fsi.toml"


def edit(config, old, new):
    text = config.read_text()
    assert old in text
    config.write_text(text.replace(old, new))


def test_loader_builds_the_problem_and_hashes_every_input(case):
    loaded = load_fsi_config(case)
    problem = loaded.problem
    assert problem.wall_name == "top" and problem.facet_set == "wet"
    assert problem.relaxation == 1.0 and problem.max_coupling_iterations == 200
    assert problem.fluid_model == "laminar"
    assert set(loaded.input_files) == {
        "config",
        "fluid.domain",
        "structure.model",
        "structure.mesh",
    }
    edit(case.with_name("domain.json"), '"width_m": 1.0', '"width_m": 1.5')
    with pytest.raises(FsiCouplingError) as caught:
        # The changed width breaks the interface, and the hash changes as well.
        run_fsi_config(case, case.parent / "out")
    assert caught.value.kind == "interface"
    assert (
        load_fsi_config(case).input_files["fluid.domain"]
        != loaded.input_files["fluid.domain"]
    )


@pytest.mark.parametrize(
    ("old", "new", "path"),
    [
        ("fsi_schema_version = 1", "fsi_schema_version = 2", "fsi_schema_version"),
        ("relaxation = 1.0\n", "", "coupling"),
        ("relaxation = 1.0", "relaxation = true", "coupling.relaxation"),
        (
            "max_coupling_iterations = 200",
            "max_coupling_iterations = 2.0",
            "coupling.max",
        ),
        ("preserve_inlet_flow = true", "preserve_inlet_flow = 1", "coupling.preserve"),
        ('model = "laminar"', 'model = "k-omega"', "fluid.model"),
        ('domain = "domain.json"', 'domain = "/etc/domain.json"', "fluid.domain"),
        ('domain = "domain.json"', 'domain = "../domain.json"', "fluid.domain"),
        ('domain = "domain.json"', 'domain = "missing.json"', "fluid.domain"),
        ('model = "structure.toml"', 'model = "structure.msh"', "structure.model"),
        (
            "[coupling]",
            "[fluid.options]\nviscosity_Pa_s = 2.0\n[coupling]",
            "fluid.options",
        ),
        (
            'model = "laminar"',
            'model = "laminar"\ninitial_nu_tilde = 1e-5',
            "only valid",
        ),
        ('model = "laminar"', 'model = "sa-noft2"', "explicit SA start"),
        ("[coupling]", "extra = 1\n[coupling]", "structure"),
    ],
)
def test_invalid_documents_are_rejected_with_a_concrete_path(case, old, new, path):
    edit(case, old, new)
    with pytest.raises(ModelSchemaError, match=path.replace(".", r"\.")):
        load_fsi_config(case)


def test_existing_output_is_refused_before_loading(case, tmp_path):
    (tmp_path / "out").mkdir()
    with pytest.raises(FileExistsError):
        run_fsi_config(case, tmp_path / "out")


def test_coupling_stop_writes_failure_artifacts_and_reraises(case, tmp_path):
    edit(case, "max_coupling_iterations = 200", "max_coupling_iterations = 1")
    with pytest.raises(FsiCouplingError) as caught:
        run_fsi_config(case, tmp_path / "failed")
    assert caught.value.kind == "budget"
    report = json.loads((tmp_path / "failed/report.json").read_text())
    assert report["status"] == "failed"
    assert report["failure"]["kind"] == "budget"
    assert len(report["last_displacement_m"]) > 0
    history = (tmp_path / "failed/history.csv").read_text().splitlines()
    assert len(history) == 2  # header plus the single completed step
    assert not (tmp_path / "failed/fluid_fields.json").exists()
