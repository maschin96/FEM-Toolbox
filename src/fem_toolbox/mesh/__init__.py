"""Validated in-memory mesh data structures."""

from fem_toolbox.mesh.model import Mesh
from fem_toolbox.mesh.topology import (
    BoundaryFace,
    BoundaryTopology,
    inspect_boundary_topology,
)

__all__ = ("BoundaryFace", "BoundaryTopology", "Mesh", "inspect_boundary_topology")
