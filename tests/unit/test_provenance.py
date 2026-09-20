"""Tests for deterministic and validated simulation provenance."""

from hashlib import sha256
from importlib.metadata import version
from math import inf

import pytest

from fem_toolbox import __version__
from fem_toolbox.errors import ProvenanceError
from fem_toolbox.provenance import RunProvenance


def test_version_comes_from_installed_package_metadata() -> None:
    """Version access does not depend on a Git executable or checkout."""
    assert __version__ == version("fem-toolbox")


def test_provenance_is_deterministic_and_detached_from_mutable_input() -> None:
    parameters: dict[str, object] = {
        "ordering": "COLAMD",
        "tolerances": [1.0e-12, 1.0e-9],
    }

    provenance = RunProvenance.from_input_bytes(
        b"deterministic model input",
        solver_name="scipy.sparse.linalg.spsolve",
        solver_parameters=parameters,
    )
    parameters["ordering"] = "NATURAL"
    parameters["tolerances"] = []

    assert provenance.input_sha256 == sha256(b"deterministic model input").hexdigest()
    assert provenance.toolbox_version == __version__
    assert provenance.solver_parameters["ordering"] == "COLAMD"
    assert provenance.to_dict() == {
        "schema_version": 1,
        "toolbox_version": __version__,
        "input_sha256": sha256(b"deterministic model input").hexdigest(),
        "solver": {
            "name": "scipy.sparse.linalg.spsolve",
            "parameters": {
                "ordering": "COLAMD",
                "tolerances": [1.0e-12, 1.0e-9],
            },
        },
    }

    with pytest.raises(TypeError):
        provenance.solver_parameters["ordering"] = "NATURAL"  # type: ignore[index]


@pytest.mark.parametrize(
    ("kwargs", "field"),
    [
        ({"input_sha256": "not-a-digest"}, "input_sha256"),
        ({"solver_name": ""}, "solver_name"),
        ({"toolbox_version": ""}, "toolbox_version"),
        ({"solver_parameters": {"rtol": inf}}, "solver_parameters"),
        ({"solver_parameters": {1: "invalid key"}}, "solver_parameters"),
        ({"solver_parameters": {"callback": object()}}, "solver_parameters"),
    ],
)
def test_invalid_provenance_raises_specific_error(
    kwargs: dict[str, object], field: str
) -> None:
    valid: dict[str, object] = {
        "input_sha256": "0" * 64,
        "solver_name": "reference",
        "solver_parameters": {},
        "toolbox_version": "0.1.0.dev0",
    }
    valid.update(kwargs)

    with pytest.raises(ProvenanceError) as caught:
        RunProvenance(**valid)  # type: ignore[arg-type]

    assert caught.value.field == field
