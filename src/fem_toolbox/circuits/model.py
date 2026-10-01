"""Owned, immutable SI components and deterministic circuit inputs."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from numbers import Real

from fem_toolbox.errors import ValidationError

CIRCUIT_SCHEMA_VERSION = 1


def _finite_scalar(value: object, *, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValidationError("value must be a real SI scalar", field=field)
    try:
        result = float(value)
    except OverflowError as error:
        raise ValidationError("value is not representable", field=field) from error
    if not math.isfinite(result):
        raise ValidationError("value must be finite", field=field)
    return result


@dataclass(frozen=True, slots=True)
class Component:
    """Two-terminal linear element with current flowing positive to negative.

    Name prefix R/C/L/V/I selects resistance (Ohm), capacitance (F), inductance
    (H), voltage (V) or current (A). Values are finite real SI scalars; passive
    values must be positive. Identifiers are case insensitive. Invalid names,
    equal terminals or unsupported values raise ValidationError.
    """

    name: str
    positive: str
    negative: str
    value: float

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not re.fullmatch(
            r"[RCLVIrclvi][A-Za-z0-9_]*", self.name
        ):
            raise ValidationError("expected an R/C/L/V/I name", field="component.name")
        object.__setattr__(self, "name", self.name.upper())
        for terminal in ("positive", "negative"):
            node = getattr(self, terminal)
            if not isinstance(node, str) or not re.fullmatch(r"[A-Za-z0-9_]+", node):
                raise ValidationError("invalid node identifier", field=terminal)
            object.__setattr__(self, terminal, node.lower())
        if self.positive == self.negative:
            raise ValidationError("terminals must differ", field=self.name)
        value = _finite_scalar(self.value, field=self.name)
        if self.kind in "RCL" and value <= 0:
            raise ValidationError("R/C/L must be positive", field=self.name)
        object.__setattr__(self, "value", value)

    @property
    def kind(self) -> str:
        """Return the R/C/L/V/I discriminator inferred from the component name."""
        return self.name[0]


@dataclass(frozen=True, slots=True)
class Circuit:
    """Sorted owned components with reference node '0'; no implicit ground aliases.

    Empty circuits, duplicate names, missing ground and invalid components raise
    ValidationError. Analysis-specific topology is checked by the DC solver.
    """

    components: tuple[Component, ...]

    def __post_init__(self) -> None:
        components = tuple(self.components)
        if not components or any(not isinstance(c, Component) for c in components):
            raise ValidationError(
                "need nonempty Component sequence", field="components"
            )
        if len({c.name for c in components}) != len(components):
            raise ValidationError("duplicate component names", field="components")
        object.__setattr__(
            self, "components", tuple(sorted(components, key=lambda c: c.name))
        )
        if "0" not in self.nodes:
            raise ValidationError("missing reference node '0'", field="ground")

    @property
    def nodes(self) -> tuple[str, ...]:
        """Return sorted node identifiers, including reference node '0'."""
        return tuple(
            sorted({n for c in self.components for n in (c.positive, c.negative)})
        )

    @classmethod
    def from_netlist(cls, text: str) -> Circuit:
        """Parse the restricted R/C/L/V/I dialect; errors identify input lines."""
        from fem_toolbox.circuits.netlist import parse_netlist

        return parse_netlist(text)

    def to_netlist(self) -> str:
        """Serialize sorted elements as roundtrippable base-SI values, then .end."""
        return (
            "\n".join(
                f"{c.name} {c.positive} {c.negative} {c.value:.17g}"
                for c in self.components
            )
            + "\n.end\n"
        )
