"""Type-consumer contract: both import paths denote the same public types."""

from fem_toolbox.errors import ValidationError as LegacyError
from fem_toolbox.mesh import Mesh as LegacyMesh
from femantis.errors import ValidationError
from femantis.mesh import Mesh


def canonical_accepts_legacy(mesh: LegacyMesh) -> Mesh:
    return mesh


def legacy_accepts_canonical(mesh: Mesh) -> LegacyMesh:
    return mesh


def legacy_exception(error: ValidationError) -> LegacyError:
    return error
