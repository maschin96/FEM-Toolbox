"""Experimental circuit input contract; values have type-defined SI units."""

from dataclasses import FrozenInstanceError

import pytest

from fem_toolbox.errors import ValidationError


def api():
    # Import in the test so Red identifies the absent feature, not conftest.
    from fem_toolbox.circuits import Circuit, Component

    return Circuit, Component


def test_netlist_normalization_roundtrip_and_owned_components():
    Circuit, Component = api()
    circuit = Circuit.from_netlist(
        "* divider\nv1 IN 0 12\nR1 IN OUT 1k\nR2 OUT 0 2k\n.end\n"
    )
    assert circuit.nodes == ("0", "in", "out")
    assert circuit.components == (
        Component("R1", "in", "out", 1000),
        Component("R2", "out", "0", 2000),
        Component("V1", "in", "0", 12),
    )
    assert Circuit.from_netlist(circuit.to_netlist()) == circuit
    incoming = list(circuit.components)
    copied = Circuit(incoming)
    incoming.clear()
    assert copied == circuit
    with pytest.raises(FrozenInstanceError):
        copied.components = ()


@pytest.mark.parametrize(
    ("token", "expected"),
    [("1M", 1e-3), ("1MEG", 1e6), ("2.5e-3k", 2.5), ("2u", 2e-6)],
)
def test_exact_engineering_suffixes(token, expected):
    Circuit, _ = api()
    circuit = Circuit.from_netlist(f"R1 a 0 {token}")
    assert circuit.components[0].value == expected


@pytest.mark.parametrize(
    "text",
    [
        "",
        "R1 a b 10",
        "R1 a 0 0",
        "R1 a 0 -1",
        "C1 a 0 0",
        "L1 a 0 -1",
        "R1 a 0 nan",
        "R1 a 0 inf",
        "R1 a 0 1e999",
        "R1 a 0 1e-999",
        "R1 a 0 1kOhm",
        "R1 a 0 1 V",
        "R1 a 0 {1+2}",
        "D1 a 0 1",
        "R1 a 0 1\nr1 b 0 2",
        "R1 a a 1\nV1 a 0 1",
        ".include other.cir",
        ".param a=1",
        "title\nR1 a 0 1",
        "R1 a 0 1\n.end\nV1 a 0 1",
        "R1 a 0 1\n.end\n.end",
        "R1 a-b 0 1",
        "R1 a 0 1 ; comment",
    ],
)
def test_invalid_netlists_fail_explicitly(text):
    Circuit, _ = api()
    with pytest.raises(ValidationError):
        Circuit.from_netlist(text)


@pytest.mark.parametrize("value", [True, "1", complex(1), float("nan"), float("inf")])
def test_direct_component_rejects_nonreal_or_nonfinite_values(value):
    _, Component = api()
    with pytest.raises(ValidationError):
        Component("R1", "a", "0", value)


def test_parser_error_identifies_line():
    Circuit, _ = api()
    with pytest.raises(ValidationError, match="line 2"):
        Circuit.from_netlist("* comment\nR1 a 0 bad")


def test_suffix_applied_before_float_range_check():
    Circuit, _ = api()
    # The final SI value is representable even when the unsuffixed token is not.
    c = Circuit.from_netlist("V1 a 0 1e309m")
    assert c.components[0].value == pytest.approx(1e306, rel=1e-15)


def test_nonzero_source_underflow_is_not_silently_zero():
    Circuit, _ = api()
    with pytest.raises(ValidationError, match="representable"):
        Circuit.from_netlist("V1 a 0 1e-999")
