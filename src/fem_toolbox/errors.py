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


class BoundaryTopologyError(MeshValidationError):
    """A surface simplex is not a unique exterior volume face, or is missing."""

    def __init__(
        self,
        *,
        kind: str,
        entity_id: int | None,
        cell_id: int | None,
        groups: tuple[str, ...],
        nodes: tuple[int, int, int],
    ) -> None:
        self.kind = kind
        self.entity_id = entity_id
        self.cell_id = cell_id
        self.groups = groups
        self.nodes = nodes
        super().__init__(
            field="boundary_topology",
            reason=(
                f"{kind} face {nodes}; triangle ID {entity_id}, "
                f"cell ID {cell_id}, physical groups {groups}"
            ),
        )


class MaterialAssignmentError(ValidationError):
    """A named volume region does not yield one valid material per active cell."""

    def __init__(
        self,
        *,
        kind: str,
        field: str,
        reason: str,
        cell_id: int | None = None,
        regions: tuple[str, ...] = (),
    ) -> None:
        self.kind = kind
        self.cell_id = cell_id
        self.regions = regions
        super().__init__(
            f"material assignment for {field}: {kind}: {reason}; "
            f"cell ID {cell_id}, regions {regions}",
            field="material_assignment",
        )


class ModelSchemaError(ValidationError):
    """A versioned declarative model contains invalid or unsupported data."""

    def __init__(self, *, path: str, reason: str) -> None:
        self.path = path
        self.reason = reason
        super().__init__(f"invalid model {path}: {reason}", field=path)


class MeshIOError(FEMToolboxError):
    """A mesh could not be converted at an external-format boundary."""

    def __init__(self, *, operation: str, reason: str) -> None:
        self.operation = operation
        self.reason = reason
        super().__init__(f"mesh {operation} failed: {reason}")


class UnsupportedCellTypeError(MeshIOError):
    """An input mesh contains cell types outside the supported simplex families."""

    def __init__(self, cell_types: tuple[str, ...]) -> None:
        self.cell_types = cell_types
        names = ", ".join(cell_types)
        super().__init__(
            operation="import",
            reason=(
                f"unsupported cell type(s): {names}; expected only "
                "tetra/triangle or tetra10/triangle6"
            ),
        )


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


class ThermalNullspaceError(ValidationError):
    """An unanchored connected thermal component has a constant null mode.

    compatible identifies a balanced Neumann problem with an undetermined
    temperature offset; otherwise no steady solution exists. net_power and
    tolerance are in W; component is its zero-based connected-component index.
    """

    def __init__(
        self, *, component: int, compatible: bool, net_power: float, tolerance: float
    ) -> None:
        self.component = component
        self.compatible = compatible
        self.net_power = net_power
        self.tolerance = tolerance
        reason = (
            "compatible Neumann problem: temperature reference missing"
            if compatible
            else "incompatible Neumann power balance"
        )
        super().__init__(
            f"component {component}: {reason}; net power={net_power:.6e} W",
            field="thermal_nullspace",
        )


class ElectricGaugeError(ValidationError):
    """Unreferenced electrical component: compatible free gauge or charge defect.

    net_charge = volume charge minus prescribed outward D-flux, in C.
    tolerance is in C; component is the zero-based connected-component index.
    """

    def __init__(
        self, *, component: int, compatible: bool, net_charge: float, tolerance: float
    ) -> None:
        self.component = component
        self.compatible = compatible
        self.net_charge = net_charge
        self.tolerance = tolerance
        reason = (
            "compatible charge problem: potential gauge missing"
            if compatible
            else "incompatible charge balance without potential reference"
        )
        super().__init__(
            f"component {component}: {reason}; net charge={net_charge:.6e} C",
            field="electric_gauge",
        )
