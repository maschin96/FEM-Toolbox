"""Deterministic Tet4 incidence and Triangle3 boundary preflight."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import BoundaryTopologyError
from fem_toolbox.mesh.model import Mesh

FloatArray = NDArray[np.float64]
FaceKey = tuple[int, int, int]


def _face_key(nodes: tuple[int, int, int]) -> FaceKey:
    ordered = sorted(nodes)
    return ordered[0], ordered[1], ordered[2]


@dataclass(frozen=True, slots=True)
class BoundaryFace:
    """Exterior face, oriented away from its unique Tet4 owner.

    ``cell_index`` and ``triangle_index`` are zero-based rows; IDs retain the
    imported entity numbers where present. ``triangle_index`` is absent for an
    exterior face with no Triangle3 entry.
    """

    nodes: FaceKey
    cell_index: int
    cell_id: int
    opposite_node: int
    area: float
    unit_normal: FloatArray
    triangle_index: int | None = None
    triangle_id: int | None = None
    groups: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class BoundaryTopology:
    """Preflight result with stable canonical order and triangle-row lookup."""

    facet_incidence: Mapping[FaceKey, tuple[int, ...]]
    exterior_faces: tuple[BoundaryFace, ...]
    triangle_faces: tuple[BoundaryFace, ...]


def inspect_boundary_topology(
    mesh: Mesh, *, require_complete: bool = False
) -> BoundaryTopology:
    """Validate each Triangle3 against the unique exterior Tet4 incidence.

    Set ``require_complete`` when every exterior Tet4 face must have a
    Triangle3 row. Partial boundary meshes remain valid by default. Named
    groups may overlap; conflicting boundary conditions belong to Model
    validation, where their physical meaning is known.
    """
    incidence: dict[FaceKey, list[BoundaryFace]] = {}
    for cell_index, cell in enumerate(mesh.tetra4):
        cell_id = (
            cell_index if mesh.tetra4_ids is None else int(mesh.tetra4_ids[cell_index])
        )
        geometry = Tet4Geometry(mesh.points[cell], entity_id=cell_id)
        for face in geometry.faces:
            nodes: FaceKey = (
                int(cell[face.local_nodes[0]]),
                int(cell[face.local_nodes[1]]),
                int(cell[face.local_nodes[2]]),
            )
            key = _face_key(nodes)
            owner = BoundaryFace(
                nodes=nodes,
                cell_index=cell_index,
                cell_id=cell_id,
                opposite_node=int(cell[face.opposite_node]),
                area=face.area,
                unit_normal=face.unit_normal,
            )
            incidence.setdefault(key, []).append(owner)

    assigned: dict[FaceKey, BoundaryFace] = {}
    triangle_faces: list[BoundaryFace] = []
    for triangle_index, triangle in enumerate(mesh.triangle3):
        key = _face_key((int(triangle[0]), int(triangle[1]), int(triangle[2])))
        owners = incidence.get(key, [])
        groups = tuple(
            sorted(
                name for name, rows in mesh.facet_sets.items() if triangle_index in rows
            )
        )
        triangle_id = (
            triangle_index
            if mesh.triangle3_ids is None
            else int(mesh.triangle3_ids[triangle_index])
        )
        if len(owners) != 1:
            kind = (
                "foreign"
                if not owners
                else "interior"
                if len(owners) == 2
                else "nonmanifold"
            )
            raise BoundaryTopologyError(
                kind=kind,
                entity_id=triangle_id,
                cell_id=owners[0].cell_id if owners else None,
                groups=groups,
                nodes=key,
            )
        owner = owners[0]
        mapped = BoundaryFace(
            nodes=owner.nodes,
            cell_index=owner.cell_index,
            cell_id=owner.cell_id,
            opposite_node=owner.opposite_node,
            area=owner.area,
            unit_normal=owner.unit_normal,
            triangle_index=triangle_index,
            triangle_id=triangle_id,
            groups=groups,
        )
        assigned[key] = mapped
        triangle_faces.append(mapped)

    exterior: list[BoundaryFace] = []
    for key in sorted(incidence):
        owners = incidence[key]
        if len(owners) > 2:
            raise BoundaryTopologyError(
                kind="nonmanifold",
                entity_id=None,
                cell_id=owners[0].cell_id,
                groups=(),
                nodes=key,
            )
        if len(owners) == 1:
            owner = assigned.get(key, owners[0])
            if require_complete and key not in assigned:
                raise BoundaryTopologyError(
                    kind="missing",
                    entity_id=None,
                    cell_id=owner.cell_id,
                    groups=(),
                    nodes=key,
                )
            exterior.append(owner)
    facet_incidence = MappingProxyType(
        {
            key: tuple(owner.cell_id for owner in incidence[key])
            for key in sorted(incidence)
        }
    )
    return BoundaryTopology(facet_incidence, tuple(exterior), tuple(triangle_faces))
