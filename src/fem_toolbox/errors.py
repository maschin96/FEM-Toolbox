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
