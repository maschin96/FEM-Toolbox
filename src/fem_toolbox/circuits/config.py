"""Strict versioned SI configuration; no executable serialization."""

from __future__ import annotations

import tomllib
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from fem_toolbox.circuits.model import Circuit, Component
from fem_toolbox.circuits.ports import ElectricalPort, ThermalPort
from fem_toolbox.errors import ModelSchemaError, UnitError


def _keys(value: object, keys: set[str], *, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        raise ModelSchemaError(path=path, reason=f"expected keys {sorted(keys)}")
    return value


def _records(value: object, *, path: str) -> list[Any]:
    if not isinstance(value, list):
        raise ModelSchemaError(path=path, reason="expected table array")
    return value


def from_document(document: Mapping[str, object]) -> Circuit:
    """Construct an immutable circuit from exact schema keys and SI units."""
    if (
        not isinstance(document, Mapping)
        or not {"circuit_schema_version", "components"} <= set(document)
        or set(document)
        - {"circuit_schema_version", "components", "electrical_ports", "thermal_ports"}
    ):
        raise ModelSchemaError(path="circuit", reason="missing or unknown keys")
    if (
        type(document["circuit_schema_version"]) is not int
        or document["circuit_schema_version"] != 1
    ):
        raise ModelSchemaError(
            path="circuit_schema_version", reason="expected integer 1"
        )
    components = []
    units = {"R": "Ohm", "C": "F", "L": "H", "V": "V", "I": "A"}
    for row in _records(document["components"], path="components"):
        r = _keys(
            row, {"name", "positive", "negative", "value", "unit"}, path="components"
        )
        c = Component(r["name"], r["positive"], r["negative"], r["value"])
        if r["unit"] != units[c.kind]:
            raise UnitError(field=c.name, expected=units[c.kind], actual=str(r["unit"]))
        components.append(c)
    ep = [
        ElectricalPort(
            **_keys(r, {"name", "node", "electrode"}, path="electrical_ports")
        )
        for r in _records(document.get("electrical_ports", []), path="electrical_ports")
    ]
    tp = [
        ThermalPort(**_keys(r, {"name", "region", "components"}, path="thermal_ports"))
        for r in _records(document.get("thermal_ports", []), path="thermal_ports")
    ]
    return Circuit(components, electrical_ports=ep, thermal_ports=tp)


def from_file(path: str | Path) -> Circuit:
    """Read TOML; report file/syntax errors as ModelSchemaError."""
    try:
        document = tomllib.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise ModelSchemaError(path=str(path), reason=str(error)) from error
    return from_document(document)
