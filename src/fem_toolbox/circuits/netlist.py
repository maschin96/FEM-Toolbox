"""Strict text adapter: no expression evaluation or external includes."""

from __future__ import annotations

import re

from fem_toolbox.circuits.model import Circuit, Component
from fem_toolbox.errors import ModelSchemaError, ValidationError

_EXPONENTS = {
    "": 0,
    "t": 12,
    "g": 9,
    "meg": 6,
    "k": 3,
    "m": -3,
    "u": -6,
    "n": -9,
    "p": -12,
    "f": -15,
}
_NUMBER = re.compile(
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))(?:[eE]([+-]?\d+))?(meg|[tgkmunpf]?)",
    re.IGNORECASE | re.ASCII,
)


def parse_netlist(text: str) -> Circuit:
    """Parse four-token R/C/L/V/I lines with SI engineering prefixes.

    Full '*' comments, blank lines and optional final '.end' are supported.
    M is milli and MEG is mega; trailing unit text is rejected. No title line,
    directives or inline comments. Raises ModelSchemaError for invalid lines.
    """
    if not isinstance(text, str):
        raise ModelSchemaError(path="netlist", reason="expected text")
    components = []
    ended = False
    names: set[str] = set()
    for number, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("*"):
            continue
        try:
            if ended:
                raise ValueError("content after .end")
            if line.lower() == ".end":
                ended = True
                continue
            tokens = line.split()
            if len(tokens) != 4:
                raise ValueError("expected: name positive negative value")
            name, positive, negative, token = tokens
            match = _NUMBER.fullmatch(token)
            if match is None:
                raise ValueError("invalid SI number or unsupported suffix")
            exponent = int(match[2] or "0") + _EXPONENTS[match[3].lower()]
            value = float(f"{match[1]}e{exponent}")
            if value == 0 and any(digit in "123456789" for digit in match[1]):
                raise ValueError("nonzero SI value is not representable")
            component = Component(name, positive, negative, value)
            if component.name in names:
                raise ValueError("duplicate component name")
            names.add(component.name)
            components.append(component)
        except (ValueError, OverflowError) as error:
            raise ModelSchemaError(
                path=f"netlist line {number}", reason=str(error)
            ) from error
    try:
        return Circuit(tuple(components))
    except ValidationError as error:
        raise ModelSchemaError(path="netlist", reason=str(error)) from error
