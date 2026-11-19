"""Deterministic simplex incidence and matching surface preflight."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.tet4 import DEFAULT_DEGENERACY_TOLERANCE
from fem_toolbox.errors import BoundaryTopologyError
from fem_toolbox.mesh.model import Mesh
from fem_toolbox.mesh.quadratic import edge_nodes, require_straight

FloatArray = NDArray[np.float64]
FaceKey = tuple[int, int, int]


def _face_key(nodes: tuple[int, ...]) -> FaceKey:
    ordered = sorted(nodes)
    return ordered[0], ordered[1], ordered[2]


@dataclass(frozen=True, slots=True)
class BoundaryFace:
    """Exterior face, oriented away from its unique tetrahedral owner.

    ``cell_index`` and ``triangle_index`` are zero-based rows; IDs retain the
    imported entity numbers where present. ``triangle_index`` is absent for an
    exterior face with no surface entry.
    """

    nodes: tuple[int, ...]
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


# Same local face order and degeneracy criterion as Tet4Geometry.
_FACE_NODES = np.array(((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)))
_EDGE_A = np.array((0, 0, 0, 1, 1, 2))
_EDGE_B = np.array((1, 2, 3, 2, 3, 3))


def _cell_id(mesh: Mesh, row: int) -> int:
    return row if mesh.cell_ids is None else int(mesh.cell_ids[row])


def _validate_linear_cells(mesh: Mesh) -> None:
    # Bound temporary coordinate arrays independently of the mesh size.
    for start in range(0, mesh.n_cells, 65536):
        points = mesh.points[mesh.cells[start : start + 65536]]
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            jacobians = (points[:, 1:] - points[:, :1]).transpose(0, 2, 1)
            determinants = np.linalg.det(jacobians)
            lengths = np.linalg.norm(points[:, _EDGE_B] - points[:, _EDGE_A], axis=2)
            longest = lengths.max(axis=1)
            scaled = 6.0 * np.sqrt(2.0) * np.abs(determinants / 6.0) / longest**3
        invalid = np.flatnonzero(
            ~np.isfinite(determinants)
            | ~np.isfinite(scaled)
            | (scaled <= DEFAULT_DEGENERACY_TOLERANCE * (1.0 + 1e-12))
        )
        if len(invalid):
            # Resolve roundoff at the tolerance boundary with the scalar oracle,
            # preserving its exception type and diagnostic values. Check every
            # candidate in input order, since a borderline one may still pass.
            for local_row in invalid:
                row = start + int(local_row)
                Tet4Geometry(
                    mesh.points[mesh.cells[row]], entity_id=_cell_id(mesh, row)
                )


def _exterior_geometry(
    mesh: Mesh, faces: NDArray[np.int64], rows: NDArray[np.int64]
) -> dict[int, BoundaryFace]:
    nodes = faces[rows].copy()
    points = mesh.points[nodes]
    cross = np.cross(points[:, 1] - points[:, 0], points[:, 2] - points[:, 0])
    norms = np.linalg.norm(cross, axis=1)
    normals = cross / norms[:, None]
    opposite = mesh.cells[rows // 4, rows % 4]
    flip = np.sum(normals * (mesh.points[opposite] - points.mean(axis=1)), axis=1) > 0
    nodes[flip] = nodes[flip][:, [0, 2, 1]]
    normals[flip] *= -1
    # Bytes-backed views preserve the public immutable-normal contract.
    normals = np.frombuffer(normals.tobytes(), dtype=np.float64).reshape(-1, 3)
    result: dict[int, BoundaryFace] = {}
    for index, row in enumerate(rows):
        cell_index = int(row) // 4
        oriented = tuple(int(node) for node in nodes[index])
        if mesh.cells.shape[1] == 10:
            edges = edge_nodes(mesh.cells[cell_index])
            oriented = (
                *oriented,
                *(
                    edges[tuple(sorted((oriented[a], oriented[b])))]
                    for a, b in ((0, 1), (1, 2), (2, 0))
                ),
            )
        result[int(row)] = BoundaryFace(
            nodes=oriented,
            cell_index=cell_index,
            cell_id=_cell_id(mesh, cell_index),
            opposite_node=int(opposite[index]),
            area=float(norms[index] / 2.0),
            unit_normal=normals[index],
        )
    return result


def inspect_boundary_topology(
    mesh: Mesh, *, require_complete: bool = False
) -> BoundaryTopology:
    """Validate each surface simplex against the unique exterior volume incidence.

    Set ``require_complete`` when every exterior volume face must have a
    surface row. Partial boundary meshes remain valid by default. Named
    groups may overlap; conflicting boundary conditions belong to Model
    validation, where their physical meaning is known.
    """
    if mesh.cells.shape[1] == 10:
        # Keep the established per-cell validation order for quadratic geometry.
        for cell_index, cell in enumerate(mesh.cells):
            Tet4Geometry(mesh.points[cell[:4]], entity_id=_cell_id(mesh, cell_index))
            require_straight(mesh.points, cell)
    else:
        _validate_linear_cells(mesh)

    # Stable lexicographic sorting retains input cell order among shared faces.
    faces = mesh.cells[:, _FACE_NODES].reshape(-1, 3)
    keys = np.sort(faces, axis=1)
    order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
    keys = keys[order]
    starts = np.r_[0, 1 + np.flatnonzero(np.any(keys[1:] != keys[:-1], axis=1))]
    ends = np.r_[starts[1:], len(keys)]
    if len(keys) == 0:
        starts = ends = np.empty(0, dtype=np.int64)
    incidence: dict[FaceKey, tuple[int, ...]] = {
        (int(key[0]), int(key[1]), int(key[2])): tuple(
            int(row) for row in order[start:end]
        )
        for key, start, end in zip(keys[starts], starts, ends, strict=True)
    }
    exterior_rows = order[starts[ends - starts == 1]]
    exterior_owners = _exterior_geometry(mesh, faces, exterior_rows)

    assigned: dict[FaceKey, BoundaryFace] = {}
    triangle_faces: list[BoundaryFace] = []
    for triangle_index, triangle in enumerate(mesh.facets):
        key = _face_key((int(triangle[0]), int(triangle[1]), int(triangle[2])))
        owners = incidence.get(key, ())
        groups = tuple(
            sorted(
                name for name, rows in mesh.facet_sets.items() if triangle_index in rows
            )
        )
        triangle_id = (
            triangle_index
            if mesh.facet_ids is None
            else int(mesh.facet_ids[triangle_index])
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
                cell_id=_cell_id(mesh, owners[0] // 4) if owners else None,
                groups=groups,
                nodes=key,
            )
        owner = exterior_owners[owners[0]]
        if len(triangle) == 6 and edge_nodes(triangle) != edge_nodes(
            np.array(owner.nodes, dtype=np.int64)
        ):
            raise BoundaryTopologyError(
                kind="nonconforming",
                entity_id=triangle_id,
                cell_id=owner.cell_id,
                groups=groups,
                nodes=key,
            )
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
    for key in incidence:
        owners = incidence[key]
        if len(owners) > 2:
            raise BoundaryTopologyError(
                kind="nonmanifold",
                entity_id=None,
                cell_id=_cell_id(mesh, owners[0] // 4),
                groups=(),
                nodes=key,
            )
        if len(owners) == 1:
            owner = assigned.get(key, exterior_owners[owners[0]])
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
            key: tuple(_cell_id(mesh, row // 4) for row in incidence[key])
            for key in incidence
        }
    )
    return BoundaryTopology(facet_incidence, tuple(exterior), tuple(triangle_faces))
