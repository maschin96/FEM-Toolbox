"""Strict SI TOML configuration, owned ports and provenance boundaries (#115)."""

from copy import deepcopy

import pytest

from fem_toolbox.circuits import Circuit, solve_dc
from fem_toolbox.errors import ValidationError


def document():
    return {
        "circuit_schema_version": 1,
        "components": [
            {"name": "V1", "positive": "a", "negative": "0", "value": 5, "unit": "V"},
            {
                "name": "R1",
                "positive": "a",
                "negative": "0",
                "value": 1000,
                "unit": "Ohm",
            },
        ],
        "electrical_ports": [{"name": "piezo", "node": "a", "electrode": "top"}],
        "thermal_ports": [
            {"name": "cooler", "region": "contact", "components": ["R1"]}
        ],
    }


def test_config_roundtrip_owned_ports_and_no_lossy_netlist():
    data = document()
    c = Circuit.from_document(data)
    assert c.electrical_ports[0].node == "a"
    assert c.thermal_ports[0].components == ("R1",)
    assert Circuit.from_document(c.to_document()) == c
    data["thermal_ports"][0]["components"].clear()
    assert c.thermal_ports[0].components == ("R1",)
    with pytest.raises(ValidationError, match="ports"):
        c.to_netlist()


def test_file_loader_and_port_hash(tmp_path):
    path = tmp_path / "circuit.toml"
    path.write_text("""circuit_schema_version=1
[[components]]
name="R1"
positive="a"
negative="0"
value=1000
unit="Ohm"
[[electrical_ports]]
name="top"
node="a"
electrode="electrode_top"
""")
    c = Circuit.from_file(path)
    assert c.electrical_ports[0].electrode == "electrode_top"
    without_ports = Circuit(c.components)
    assert (
        solve_dc(c).provenance.input_sha256
        != solve_dc(without_ports).provenance.input_sha256
    )


@pytest.mark.parametrize(
    "change",
    [
        lambda d: d.update(circuit_schema_version=True),
        lambda d: d.update(circuit_schema_version=2),
        lambda d: d.update(unknown=1),
        lambda d: d["components"][1].update(unit="kOhm"),
        lambda d: d["components"][1].update(extra=1),
        lambda d: d["electrical_ports"][0].update(node="missing"),
        lambda d: d["electrical_ports"].append(dict(d["electrical_ports"][0])),
        lambda d: d["thermal_ports"][0].update(components=["V1"]),
        lambda d: d["thermal_ports"][0].update(components=["R1", "R1"]),
        lambda d: d["thermal_ports"].append(
            {"name": "second", "region": "other", "components": ["R1"]}
        ),
        lambda d: d["electrical_ports"][0].update(electrode=""),
        lambda d: d["thermal_ports"][0].update(region=""),
    ],
)
def test_invalid_config_is_specific_error(change):
    data = deepcopy(document())
    change(data)
    with pytest.raises(ValidationError):
        Circuit.from_document(data)


def test_invalid_toml_and_missing_file(tmp_path):
    with pytest.raises(ValidationError):
        Circuit.from_file(tmp_path / "absent.toml")
    path = tmp_path / "bad.toml"
    path.write_text("{bad")
    with pytest.raises(ValidationError):
        Circuit.from_file(path)


def test_thermal_component_list_type_is_checked():
    data = document()
    data["thermal_ports"][0]["components"] = 1
    with pytest.raises(ValidationError):
        Circuit.from_document(data)
