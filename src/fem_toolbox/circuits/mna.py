"""Shared linear storage, source assembly and deterministic provenance helpers."""

from __future__ import annotations

import json
from collections.abc import Mapping
from importlib.metadata import version

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_matrix, csc_matrix  # type: ignore[import-untyped]

from fem_toolbox.circuits.model import Circuit
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


def storage(
    circuit: Circuit, nodes: dict[str, int], branches: dict[str, int]
) -> csc_matrix:
    rows, columns, values = [], [], []
    for c in circuit.components:
        if c.kind == "C":
            incidence = [
                (nodes[n], s)
                for n, s in ((c.positive, 1), (c.negative, -1))
                if n != "0"
            ]
            for row, a in incidence:
                for column, b in incidence:
                    rows.append(row)
                    columns.append(column)
                    values.append(a * b * c.value)
        elif c.kind == "L":
            rows.append(branches[c.name])
            columns.append(branches[c.name])
            values.append(-c.value)
    size = len(nodes) + len(branches)
    return coo_matrix((values, (rows, columns)), shape=(size, size)).tocsc()


def source_names(
    circuit: Circuit, overrides: Mapping[str, object] | None
) -> dict[str, object]:
    defaults: dict[str, object] = {
        c.name: c.value for c in circuit.components if c.kind in "VI"
    }
    seen: set[str] = set()
    for name, value in (overrides or {}).items():
        if (
            not isinstance(name, str)
            or name.upper() not in defaults
            or name.upper() in seen
        ):
            raise ValidationError(
                "unknown or duplicate source override", field="sources"
            )
        seen.add(name.upper())
        defaults[name.upper()] = value
    return defaults


def source_rhs(
    circuit: Circuit,
    nodes: dict[str, int],
    branches: dict[str, int],
    values: Mapping[str, complex],
) -> NDArray[np.complex128]:
    rhs = np.zeros(len(nodes) + len(branches), dtype=np.complex128)
    for c in circuit.components:
        if c.kind == "V":
            rhs[branches[c.name]] = values[c.name]
        elif c.kind == "I":
            for n, sign in ((c.positive, 1), (c.negative, -1)):
                if n != "0":
                    rhs[nodes[n]] -= sign * values[c.name]
    return rhs


def provenance(
    circuit: Circuit,
    analysis: str,
    parameters: Mapping[str, object],
    extra: Mapping[str, object],
) -> RunProvenance:
    canonical = json.dumps(
        {"circuit": circuit.to_document(), **extra},
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return RunProvenance.from_input_bytes(
        canonical,
        solver_name=f"fem_toolbox.circuits.solve_{analysis}",
        solver_parameters={
            "revision": 1,
            "numpy": np.__version__,
            "scipy": version("scipy"),
            **parameters,
        },
    )
