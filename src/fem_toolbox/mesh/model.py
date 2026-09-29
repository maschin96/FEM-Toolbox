"""Validated topology and selection storage for three-dimensional meshes."""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.errors import (
    ConnectivityError,
    MeshValidationError,
    SelectionError,
    ShapeError,
)
from fem_toolbox.mesh.quadratic import validate_edge_incidence

FloatArray = NDArray[np.float64]
IndexArray = NDArray[np.int64]
SelectionMap = Mapping[str, IndexArray]


def _frozen_array(array: NDArray[np.generic]) -> NDArray[np.generic]:
    """Copy an array onto immutable bytes-backed storage."""
    return np.frombuffer(array.tobytes(order="C"), dtype=array.dtype).reshape(
        array.shape
    )


def _points_array(value: ArrayLike) -> FloatArray:
    array = np.asarray(value)
    if array.ndim != 2 or array.shape[1:] != (3,):
        raise ShapeError(field="points", expected="(n, 3)", actual=array.shape)
    if array.dtype.kind not in "iuf":
        raise MeshValidationError(field="points", reason="expected real numbers")
    points = np.asarray(array, dtype=np.float64)
    if not np.isfinite(points).all():
        raise MeshValidationError(
            field="points", reason="all coordinates must be finite"
        )
    return _frozen_array(points).view(np.ndarray)


def _index_matrix(value: ArrayLike | None, *, field: str, width: int) -> IndexArray:
    if value is None:
        array = np.empty((0, width), dtype=np.int64)
    else:
        array = np.asarray(value)
    if array.ndim != 2 or array.shape[1:] != (width,):
        raise ShapeError(field=field, expected=f"(n, {width})", actual=array.shape)
    if array.dtype.kind not in "iu":
        raise MeshValidationError(field=field, reason="expected integer node indices")
    return _frozen_array(np.asarray(array, dtype=np.int64)).view(np.ndarray)


def _external_ids(
    value: ArrayLike | None, *, size: int, field: str
) -> IndexArray | None:
    if value is None:
        return None
    array = np.asarray(value)
    if array.shape != (size,):
        raise ShapeError(field=field, expected=f"({size},)", actual=array.shape)
    if array.dtype.kind not in "iu":
        raise MeshValidationError(field=field, reason="expected integer external IDs")
    ids = np.asarray(array, dtype=np.int64)
    seen: set[int] = set()
    for entity_id in ids:
        identifier = int(entity_id)
        if identifier in seen:
            raise MeshValidationError(
                field=field, reason=f"duplicate external ID {identifier}"
            )
        seen.add(identifier)
    return _frozen_array(ids).view(np.ndarray)


def _entity_tags(
    value: ArrayLike | None, *, size: int, field: str
) -> IndexArray | None:
    if value is None:
        return None
    array = np.asarray(value)
    if array.shape != (size,):
        raise ShapeError(field=field, expected=f"({size},)", actual=array.shape)
    if array.dtype.kind not in "iu":
        raise MeshValidationError(field=field, reason="expected integer tags")
    return _frozen_array(np.asarray(array, dtype=np.int64)).view(np.ndarray)


def _reported_id(ids: IndexArray | None, internal_index: int) -> int:
    return internal_index if ids is None else int(ids[internal_index])


def _validate_connectivity(
    connectivity: IndexArray,
    *,
    field: str,
    point_count: int,
    external_ids: IndexArray | None,
    selections: SelectionMap | None = None,
) -> None:
    canonical_rows: dict[tuple[int, ...], int] = {}
    for internal_index, row in enumerate(connectivity):
        entity_id = _reported_id(external_ids, internal_index)
        outside = np.flatnonzero((row < 0) | (row >= point_count))
        if outside.size:
            local_index = int(outside[0])
            node_index = int(row[local_index])
            raise ConnectivityError(
                entity_kind=field,
                entity_id=entity_id,
                reason=(
                    f"node index {node_index} at local position {local_index} is "
                    f"outside [0, {point_count})"
                ),
            )
        canonical = tuple(sorted(int(node) for node in row))
        if len(set(canonical)) != len(canonical):
            raise ConnectivityError(
                entity_kind=field,
                entity_id=entity_id,
                reason="contains a repeated node index",
            )
        duplicate_index = canonical_rows.get(canonical)
        if duplicate_index is not None:
            groups = (
                tuple(
                    sorted(
                        name
                        for name, rows in selections.items()
                        if internal_index in rows
                    )
                )
                if selections is not None
                else ()
            )
            raise ConnectivityError(
                entity_kind=field,
                entity_id=entity_id,
                reason=f"connectivity is topologically duplicated; physical groups {groups}",
                duplicate_of=_reported_id(external_ids, duplicate_index),
            )
        canonical_rows[canonical] = internal_index


def _selection_map(
    value: Mapping[str, ArrayLike] | None,
    *,
    field: str,
    entity_count: int,
) -> SelectionMap:
    if value is None:
        return MappingProxyType({})
    if not isinstance(value, Mapping):
        raise MeshValidationError(field=field, reason="expected a mapping")
    selections: dict[str, IndexArray] = {}
    for name, indices_value in value.items():
        if not isinstance(name, str) or not name or name.strip() != name:
            display_name = name if isinstance(name, str) else repr(name)
            raise SelectionError(
                field=field,
                selection_name=display_name,
                reason="name must be a non-empty string without surrounding spaces",
            )
        indices = np.asarray(indices_value)
        if indices.ndim != 1:
            raise SelectionError(
                field=field,
                selection_name=name,
                reason="indices must be one-dimensional",
            )
        if indices.dtype.kind not in "iu":
            raise SelectionError(
                field=field, selection_name=name, reason="indices must be integers"
            )
        normalized = np.asarray(indices, dtype=np.int64)
        if len(np.unique(normalized)) != len(normalized):
            raise SelectionError(
                field=field, selection_name=name, reason="contains duplicate indices"
            )
        outside = np.flatnonzero((normalized < 0) | (normalized >= entity_count))
        if outside.size:
            bad_index = int(normalized[int(outside[0])])
            raise SelectionError(
                field=field,
                selection_name=name,
                reason=f"entity index {bad_index} is outside [0, {entity_count})",
            )
        selections[name] = _frozen_array(normalized).view(np.ndarray)
    return MappingProxyType(selections)


class Mesh:
    """Immutable topology for one Tet4/Triangle3 or Tet10/Triangle6 family.

    Connectivity and selections always contain zero-based internal indices.
    Optional external IDs remain separate so importer-specific numbering cannot
    leak into numerical kernels. Inputs are copied into read-only canonical
    ``float64`` and ``int64`` arrays.
    """

    __slots__ = (
        "_cell_tags",
        "_facet_sets",
        "_facet_tags",
        "_point_ids",
        "_points",
        "_tetra4",
        "_tetra4_ids",
        "_tetra10",
        "_tetra10_ids",
        "_triangle3",
        "_triangle3_ids",
        "_triangle6",
        "_triangle6_ids",
        "_volume_sets",
    )

    def __init__(
        self,
        *,
        points: ArrayLike,
        tetra4: ArrayLike | None = None,
        triangle3: ArrayLike | None = None,
        tetra10: ArrayLike | None = None,
        triangle6: ArrayLike | None = None,
        tetra10_ids: ArrayLike | None = None,
        triangle6_ids: ArrayLike | None = None,
        point_ids: ArrayLike | None = None,
        tetra4_ids: ArrayLike | None = None,
        triangle3_ids: ArrayLike | None = None,
        cell_tags: ArrayLike | None = None,
        facet_tags: ArrayLike | None = None,
        volume_sets: Mapping[str, ArrayLike] | None = None,
        facet_sets: Mapping[str, ArrayLike] | None = None,
    ) -> None:
        points_array = _points_array(points)
        tetra4_array = _index_matrix(tetra4, field="tetra4", width=4)
        triangle3_array = _index_matrix(triangle3, field="triangle3", width=3)
        tetra10_array = _index_matrix(tetra10, field="tetra10", width=10)
        triangle6_array = _index_matrix(triangle6, field="triangle6", width=6)
        if (len(tetra4_array) or len(triangle3_array)) and (
            len(tetra10_array) or len(triangle6_array)
        ):
            raise MeshValidationError(
                field="cells",
                reason="hybrid Tet4/Tet10 or Triangle3/Triangle6 meshes are unsupported",
            )
        tetra10_ids_array = _external_ids(
            tetra10_ids, size=len(tetra10_array), field="tetra10_ids"
        )
        triangle6_ids_array = _external_ids(
            triangle6_ids, size=len(triangle6_array), field="triangle6_ids"
        )
        cell_count = len(tetra4_array) + len(tetra10_array)
        facet_count = len(triangle3_array) + len(triangle6_array)
        point_ids_array = _external_ids(
            point_ids, size=len(points_array), field="point_ids"
        )
        tetra4_ids_array = _external_ids(
            tetra4_ids, size=len(tetra4_array), field="tetra4_ids"
        )
        triangle3_ids_array = _external_ids(
            triangle3_ids, size=len(triangle3_array), field="triangle3_ids"
        )
        cell_tags_array = _entity_tags(cell_tags, size=cell_count, field="cell_tags")
        facet_tags_array = _entity_tags(
            facet_tags, size=facet_count, field="facet_tags"
        )
        volume_set_map = _selection_map(
            volume_sets, field="volume_sets", entity_count=cell_count
        )
        facet_set_map = _selection_map(
            facet_sets, field="facet_sets", entity_count=facet_count
        )
        _validate_connectivity(
            tetra4_array,
            field="tetra4",
            point_count=len(points_array),
            external_ids=tetra4_ids_array,
        )
        _validate_connectivity(
            triangle3_array,
            field="triangle3",
            point_count=len(points_array),
            external_ids=triangle3_ids_array,
            selections=facet_set_map,
        )

        for field, array, ids in (
            ("tetra10", tetra10_array, tetra10_ids_array),
            ("triangle6", triangle6_array, triangle6_ids_array),
        ):
            _validate_connectivity(
                array,
                field=field,
                point_count=len(points_array),
                external_ids=ids,
                selections=facet_set_map if field == "triangle6" else volume_set_map,
            )
            # The corner simplex identifies an entity even if its edge nodes differ.
            _validate_connectivity(
                array[:, : 4 if field == "tetra10" else 3],
                field=field,
                point_count=len(points_array),
                external_ids=ids,
                selections=facet_set_map if field == "triangle6" else volume_set_map,
            )
        validate_edge_incidence(tetra10_array, tetra10_ids_array)
        self._tetra10 = tetra10_array
        self._triangle6 = triangle6_array
        self._tetra10_ids = tetra10_ids_array
        self._triangle6_ids = triangle6_ids_array

        self._points = points_array
        self._tetra4 = tetra4_array
        self._triangle3 = triangle3_array
        self._point_ids = point_ids_array
        self._tetra4_ids = tetra4_ids_array
        self._triangle3_ids = triangle3_ids_array
        self._cell_tags = cell_tags_array
        self._facet_tags = facet_tags_array
        self._volume_sets = volume_set_map
        self._facet_sets = facet_set_map

    @classmethod
    def empty(cls) -> Mesh:
        """Create an empty mesh with canonical array shapes."""
        return cls(points=np.empty((0, 3), dtype=np.float64))

    @property
    def points(self) -> FloatArray:
        """Point coordinates in metres with shape ``(n_points, 3)``."""
        return self._points

    @property
    def tetra4(self) -> IndexArray:
        """Zero-based Tet4 node connectivity with shape ``(n_tetra4, 4)``."""
        return self._tetra4

    @property
    def triangle3(self) -> IndexArray:
        """Zero-based Triangle3 node connectivity with shape ``(n_triangle3, 3)``."""
        return self._triangle3

    @property
    def point_ids(self) -> IndexArray | None:
        """Optional importer-provided point IDs, separate from internal indices."""
        return self._point_ids

    @property
    def tetra4_ids(self) -> IndexArray | None:
        """Optional importer-provided Tet4 IDs."""
        return self._tetra4_ids

    @property
    def triangle3_ids(self) -> IndexArray | None:
        """Optional importer-provided Triangle3 IDs."""
        return self._triangle3_ids

    @property
    def cell_tags(self) -> IndexArray | None:
        """Optional integer region tag for every active volume cell."""
        return self._cell_tags

    @property
    def facet_tags(self) -> IndexArray | None:
        """Optional integer boundary tag for every active surface facet."""
        return self._facet_tags

    @property
    def volume_sets(self) -> SelectionMap:
        """Named active-volume selections containing zero-based internal cell indices."""
        return self._volume_sets

    @property
    def facet_sets(self) -> SelectionMap:
        """Named active-surface selections containing zero-based internal facet indices."""
        return self._facet_sets

    @property
    def n_points(self) -> int:
        """Number of mesh points."""
        return len(self._points)

    @property
    def n_tetra4(self) -> int:
        """Number of Tet4 cells."""
        return len(self._tetra4)

    @property
    def n_triangle3(self) -> int:
        """Number of Triangle3 facets."""
        return len(self._triangle3)

    @property
    def tetra10(self) -> IndexArray:
        """Tet10 rows: corners 0..3, edge nodes 01,12,20,03,13,23 (VTK order)."""
        return self._tetra10

    @property
    def triangle6(self) -> IndexArray:
        """Triangle6 rows: corners 0..2, edge nodes 01,12,20."""
        return self._triangle6

    @property
    def tetra10_ids(self) -> IndexArray | None:
        """Optional external Tet10 IDs."""
        return self._tetra10_ids

    @property
    def triangle6_ids(self) -> IndexArray | None:
        """Optional external Triangle6 IDs."""
        return self._triangle6_ids

    @property
    def n_tetra10(self) -> int:
        return len(self._tetra10)

    @property
    def n_triangle6(self) -> int:
        return len(self._triangle6)

    @property
    def cell_type(self) -> str:
        """Active volume family; empty meshes default to tetra4."""
        return "tetra10" if self.n_tetra10 or self.n_triangle6 else "tetra4"

    @property
    def facet_type(self) -> str:
        return "triangle6" if self.cell_type == "tetra10" else "triangle3"

    @property
    def cells(self) -> IndexArray:
        """Active volume rows addressed by cell_tags and volume_sets."""
        return self.tetra10 if self.cell_type == "tetra10" else self.tetra4

    @property
    def facets(self) -> IndexArray:
        """Active surface rows addressed by facet_tags and facet_sets."""
        return self.triangle6 if self.facet_type == "triangle6" else self.triangle3

    @property
    def cell_ids(self) -> IndexArray | None:
        return self.tetra10_ids if self.cell_type == "tetra10" else self.tetra4_ids

    @property
    def facet_ids(self) -> IndexArray | None:
        return (
            self.triangle6_ids if self.facet_type == "triangle6" else self.triangle3_ids
        )

    @property
    def n_cells(self) -> int:
        return len(self.cells)

    @property
    def n_facets(self) -> int:
        return len(self.facets)
