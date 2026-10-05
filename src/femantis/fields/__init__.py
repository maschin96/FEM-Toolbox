"""Field-major DOFs, node-major and component-major within each field."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.assembly import indices
from femantis.errors import ValidationError


@dataclass(frozen=True)
class Field:
    """Active nodal field with ordered components/nodes and SI unit metadata.

    Node IDs are internal zero-based mesh indices; the caller selects only the
    participating nodes. Empty support is allowed. Names and unit must be nonempty.
    """

    name: str
    components: tuple[str, ...]
    nodes: tuple[int, ...]
    unit: str

    def __post_init__(self) -> None:
        components = tuple(self.components)
        nodes = tuple(int(n) for n in indices(self.nodes))
        if (
            not isinstance(self.name, str)
            or not self.name.strip()
            or not isinstance(self.unit, str)
            or not self.unit.strip()
            or not components
            or any(not isinstance(c, str) or not c.strip() for c in components)
            or len(set(components)) != len(components)
            or len(set(nodes)) != len(nodes)
        ):
            raise ValidationError(
                "invalid field name, unit, components or nodes", field="field"
            )
        object.__setattr__(self, "components", components)
        object.__setattr__(self, "nodes", nodes)


class DofMap:
    """Allocate exactly the supplied fields in input order, without padding."""

    def __init__(self, fields: Sequence[Field]) -> None:
        self._fields: dict[str, Field] = {}
        self._blocks: dict[str, slice] = {}
        self._lookup: dict[str, dict[int, int]] = {}
        offset = 0
        for f in fields:
            if f.name in self._fields:
                raise ValidationError("duplicate field name", field=f.name)
            self._fields[f.name] = f
            count = len(f.nodes) * len(f.components)
            self._blocks[f.name] = slice(offset, offset + count)
            self._lookup[f.name] = {
                node: offset + i * len(f.components) for i, node in enumerate(f.nodes)
            }
            offset += count
        self._size = offset

    @property
    def size(self) -> int:
        """Total number of active DOFs."""
        return self._size

    def block(self, field: str) -> slice:
        """Contiguous field slice, suitable for subsequent block solvers."""
        if field not in self._blocks:
            raise ValidationError("unknown field", field=field)
        return self._blocks[field]

    def indices(
        self, field: str, nodes: ArrayLike, components: Sequence[str] | None = None
    ) -> NDArray[np.int64]:
        """Flatten requested nodes/components in supplied order.

        Unknown fields, inactive nodes and unknown components raise ValidationError.
        """
        self.block(field)
        f = self._fields[field]
        selected = f.components if components is None else tuple(components)
        try:
            offsets = [f.components.index(c) for c in selected]
            return np.array(
                [
                    self._lookup[field][int(n)] + c
                    for n in indices(nodes)
                    for c in offsets
                ],
                dtype=np.int64,
            )
        except (KeyError, ValueError) as error:
            raise ValidationError(
                "inactive node or unknown component", field=field
            ) from error
