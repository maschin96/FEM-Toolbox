"""Exception hierarchy for explicit FEM Toolbox validation failures."""


class FEMToolboxError(Exception):
    """Base class for errors intentionally raised by the toolbox."""


class ValidationError(FEMToolboxError, ValueError):
    """Base class for invalid user or model data.

    Args:
        message: Human-readable failure description.
        field: Stable field identifier for programmatic error handling.
    """

    def __init__(self, message: str, *, field: str) -> None:
        self.field = field
        super().__init__(message)


class ShapeError(ValidationError):
    """An array or tensor has an incompatible shape."""

    def __init__(self, *, field: str, expected: str, actual: tuple[int, ...]) -> None:
        self.expected = expected
        self.actual = actual
        super().__init__(
            f"{field} has shape {actual}; expected {expected}", field=field
        )


class UnitError(ValidationError):
    """A physical quantity is not expressed in the required SI unit."""

    def __init__(self, *, field: str, expected: str, actual: str) -> None:
        self.expected = expected
        self.actual = actual
        super().__init__(
            f"{field} uses unit {actual!r}; expected SI unit {expected!r}",
            field=field,
        )


class ProvenanceError(ValidationError):
    """Simulation provenance is incomplete or cannot be serialized safely."""

    def __init__(self, *, field: str, reason: str) -> None:
        super().__init__(f"invalid provenance field {field!r}: {reason}", field=field)


class MeshValidationError(ValidationError):
    """Mesh topology or metadata violates the in-memory data contract."""

    def __init__(self, *, field: str, reason: str) -> None:
        super().__init__(f"invalid mesh field {field!r}: {reason}", field=field)


class ConnectivityError(MeshValidationError):
    """A cell or facet connectivity row is invalid.

    ``entity_id`` is the external entity ID when supplied and otherwise the
    zero-based internal row index. ``duplicate_of`` follows the same rule.
    """

    def __init__(
        self,
        *,
        entity_kind: str,
        entity_id: int,
        reason: str,
        duplicate_of: int | None = None,
    ) -> None:
        self.entity_kind = entity_kind
        self.entity_id = entity_id
        self.duplicate_of = duplicate_of
        detail = f"entity {entity_id}: {reason}"
        if duplicate_of is not None:
            detail += f"; duplicates entity {duplicate_of}"
        super().__init__(field=entity_kind, reason=detail)


class SelectionError(MeshValidationError):
    """A named volume or facet selection contains invalid internal indices."""

    def __init__(self, *, field: str, selection_name: str, reason: str) -> None:
        self.selection_name = selection_name
        label = repr(selection_name) if selection_name else "<empty>"
        super().__init__(field=field, reason=f"selection {label}: {reason}")


class ElementGeometryError(ValidationError):
    """Element coordinates or geometry settings are numerically invalid."""

    def __init__(self, *, field: str, reason: str) -> None:
        super().__init__(
            f"invalid element geometry field {field!r}: {reason}", field=field
        )


class DegenerateElementError(ElementGeometryError):
    """A Tet4 has insufficient dimensionless volume for stable operators."""

    def __init__(
        self,
        *,
        entity_id: int | None,
        signed_volume: float,
        scaled_volume: float,
        tolerance: float,
    ) -> None:
        self.entity_id = entity_id
        self.signed_volume = signed_volume
        self.scaled_volume = scaled_volume
        self.tolerance = tolerance
        label = "<unspecified>" if entity_id is None else str(entity_id)
        super().__init__(
            field="coordinates",
            reason=(
                f"Tet4 entity {label} is degenerate: scaled volume "
                f"{scaled_volume:.6e} <= tolerance {tolerance:.6e} "
                f"(signed volume {signed_volume:.6e} m^3)"
            ),
        )
