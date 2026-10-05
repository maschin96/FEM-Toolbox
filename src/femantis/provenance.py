"""Deterministic, JSON-compatible provenance for simulation runs."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from hashlib import sha256
from math import isfinite
from types import MappingProxyType

from femantis import __version__
from femantis.errors import ProvenanceError

PROVENANCE_SCHEMA_VERSION = 1


def _freeze_json(value: object, *, path: str) -> object:
    """Validate and copy a JSON-compatible value into immutable containers."""
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not isfinite(value):
            raise ProvenanceError(
                field="solver_parameters", reason=f"{path} must be finite"
            )
        return value
    if isinstance(value, Mapping):
        items: list[tuple[str, object]] = []
        for key, nested in value.items():
            if not isinstance(key, str):
                raise ProvenanceError(
                    field="solver_parameters",
                    reason=f"{path} contains non-string key {key!r}",
                )
            items.append((key, _freeze_json(nested, path=f"{path}.{key}")))
        return MappingProxyType(dict(sorted(items)))
    if isinstance(value, (list, tuple)):
        return tuple(
            _freeze_json(nested, path=f"{path}[{index}]")
            for index, nested in enumerate(value)
        )
    raise ProvenanceError(
        field="solver_parameters",
        reason=f"{path} has unsupported type {type(value).__name__}",
    )


def _thaw_json(value: object) -> object:
    """Create plain JSON-ready dictionaries and lists from frozen values."""
    if isinstance(value, Mapping):
        return {str(key): _thaw_json(nested) for key, nested in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(nested) for nested in value]
    return value


@dataclass(frozen=True, slots=True)
class RunProvenance:
    """Minimal reproducibility record for one solver invocation.

    Args:
        input_sha256: Lowercase SHA-256 digest of the canonical model input.
        solver_name: Stable fully qualified solver identifier.
        solver_parameters: JSON-compatible parameters that affect the result.
        toolbox_version: Installed FEM Toolbox distribution version.

    The record deliberately excludes timestamps and host-specific paths so that
    identical inputs and solver settings yield identical serialized metadata.
    """

    input_sha256: str
    solver_name: str
    solver_parameters: Mapping[str, object]
    toolbox_version: str = field(default=__version__)
    schema_version: int = field(default=PROVENANCE_SCHEMA_VERSION, init=False)

    def __post_init__(self) -> None:
        if len(self.input_sha256) != 64 or any(
            character not in "0123456789abcdef" for character in self.input_sha256
        ):
            raise ProvenanceError(
                field="input_sha256",
                reason="expected a 64-character lowercase SHA-256 digest",
            )
        if not self.solver_name or self.solver_name.strip() != self.solver_name:
            raise ProvenanceError(
                field="solver_name",
                reason="expected a non-empty identifier without surrounding spaces",
            )
        if (
            not self.toolbox_version
            or self.toolbox_version.strip() != self.toolbox_version
        ):
            raise ProvenanceError(
                field="toolbox_version",
                reason="expected a non-empty version without surrounding spaces",
            )
        frozen_parameters = _freeze_json(
            self.solver_parameters, path="solver_parameters"
        )
        if not isinstance(frozen_parameters, Mapping):
            raise ProvenanceError(
                field="solver_parameters", reason="expected a mapping"
            )
        object.__setattr__(self, "solver_parameters", frozen_parameters)

    @classmethod
    def from_input_bytes(
        cls,
        model_input: bytes,
        *,
        solver_name: str,
        solver_parameters: Mapping[str, object],
        toolbox_version: str = __version__,
    ) -> RunProvenance:
        """Build provenance from canonical model-input bytes."""
        return cls(
            input_sha256=sha256(model_input).hexdigest(),
            solver_name=solver_name,
            solver_parameters=solver_parameters,
            toolbox_version=toolbox_version,
        )

    def to_dict(self) -> dict[str, object]:
        """Return a detached structure accepted by strict JSON serializers."""
        return {
            "schema_version": self.schema_version,
            "toolbox_version": self.toolbox_version,
            "input_sha256": self.input_sha256,
            "solver": {
                "name": self.solver_name,
                "parameters": _thaw_json(self.solver_parameters),
            },
        }
