"""Interfaces for physics-specific volume and surface integration adapters."""

from typing import Protocol

from femantis.assembly import VectorBlock
from femantis.fields import DofMap
from femantis.mesh import Mesh


class VolumeLoad(Protocol):
    """Integrate one Tet4 volume contribution, including its weak-form sign.

    cell is a zero-based mesh.tetra4 row. Returned values are integrated SI
    equation loads; the adapter owns quadrature, field and material conventions.
    """

    def integrate(self, mesh: Mesh, cell: int, dofs: DofMap) -> VectorBlock:
        """Return local global-DOF indices and integrated values."""
        ...


class SurfaceLoad(Protocol):
    """Integrate one Triangle3 load; facet is a zero-based mesh.triangle3 row.

    The implementation owns outward orientation and weak-form sign, as well as
    quadrature. No universal thermal/electric sign is imposed by assembly.
    """

    def integrate(self, mesh: Mesh, facet: int, dofs: DofMap) -> VectorBlock:
        """Return local global-DOF indices and integrated SI equation loads."""
        ...
