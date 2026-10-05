"""Validated in-memory mesh data structures."""

from femantis.mesh.model import Mesh
from femantis.mesh.topology import (
    BoundaryFace,
    BoundaryTopology,
    inspect_boundary_topology,
)

__all__ = ("BoundaryFace", "BoundaryTopology", "Mesh", "inspect_boundary_topology")
