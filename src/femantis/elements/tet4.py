"""Geometry operators and quality measures for the linear tetrahedron."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import isfinite, sqrt

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.errors import (
    DegenerateElementError,
    ElementGeometryError,
    ShapeError,
)

FloatArray = NDArray[np.float64]
IndexArray = NDArray[np.int64]

DEFAULT_DEGENERACY_TOLERANCE = 1.0e-12

_REFERENCE_GRADIENTS = np.array(
    [
        [-1.0, -1.0, -1.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ],
    dtype=np.float64,
)
_FACE_NODES = (
    (1, 2, 3),
    (0, 2, 3),
    (0, 1, 3),
    (0, 1, 2),
)


def _frozen_float_array(value: ArrayLike, *, shape: tuple[int, ...]) -> FloatArray:
    array = np.asarray(value, dtype=np.float64)
    return np.frombuffer(array.tobytes(order="C"), dtype=np.float64).reshape(shape)


def _frozen_index_array(value: ArrayLike, *, shape: tuple[int, ...]) -> IndexArray:
    array = np.asarray(value, dtype=np.int64)
    return np.frombuffer(array.tobytes(order="C"), dtype=np.int64).reshape(shape)


def _coordinates_array(value: ArrayLike) -> FloatArray:
    raw = np.asarray(value)
    if raw.shape != (4, 3):
        raise ShapeError(field="coordinates", expected="(4, 3)", actual=raw.shape)
    if raw.dtype.kind not in "iuf":
        raise ElementGeometryError(field="coordinates", reason="expected real numbers")
    coordinates = np.asarray(raw, dtype=np.float64)
    if not np.isfinite(coordinates).all():
        raise ElementGeometryError(
            field="coordinates", reason="all coordinates must be finite"
        )
    return _frozen_float_array(coordinates, shape=(4, 3))


def _point_array(value: ArrayLike) -> FloatArray:
    raw = np.asarray(value)
    if raw.shape != (3,):
        raise ShapeError(field="point", expected="(3,)", actual=raw.shape)
    if raw.dtype.kind not in "iuf":
        raise ElementGeometryError(field="point", reason="expected real numbers")
    point = np.asarray(raw, dtype=np.float64)
    if not np.isfinite(point).all():
        raise ElementGeometryError(field="point", reason="coordinates must be finite")
    return point


@dataclass(frozen=True, slots=True)
class Tet4Face:
    """Geometry of one face, identified by local nodes and opposite node."""

    local_nodes: tuple[int, int, int]
    opposite_node: int
    area: float
    unit_normal: FloatArray
    centroid: FloatArray


@dataclass(frozen=True, slots=True)
class Tet4Quality:
    """Permutation-, rotation-, translation-, and scale-invariant quality data."""

    min_edge_length: float
    max_edge_length: float
    aspect_ratio: float
    scaled_volume: float
    condition_number: float


class Tet4Geometry:
    """Validated constant geometry operators for one Tet4 element.

    ``shape_gradients`` stores one physical gradient per node in a ``(4, 3)``
    array. ``gradient_matrix`` exposes the equivalent FEM convention ``G`` with
    shape ``(3, 4)``. The raw signed volume retains the input orientation;
    ``oriented_node_order`` and ``oriented_coordinates`` provide a positive
    orientation without mutating the input.
    """

    __slots__ = (
        "_coordinates",
        "_entity_id",
        "_faces",
        "_jacobian",
        "_oriented_coordinates",
        "_oriented_node_order",
        "_quality",
        "_shape_gradients",
        "_signed_volume",
        "_volume",
    )

    def __init__(
        self,
        coordinates: ArrayLike,
        *,
        entity_id: int | None = None,
        degeneracy_tolerance: float = DEFAULT_DEGENERACY_TOLERANCE,
    ) -> None:
        if not isfinite(degeneracy_tolerance) or degeneracy_tolerance <= 0.0:
            raise ElementGeometryError(
                field="degeneracy_tolerance",
                reason="expected a positive finite dimensionless value",
            )
        points = _coordinates_array(coordinates)
        jacobian = np.column_stack(
            (points[1] - points[0], points[2] - points[0], points[3] - points[0])
        )
        determinant = float(np.linalg.det(jacobian))
        if not isfinite(determinant):
            raise ElementGeometryError(
                field="coordinates", reason="Jacobian determinant is not finite"
            )
        signed_volume = determinant / 6.0
        volume = abs(signed_volume)
        edge_lengths = np.array(
            [
                np.linalg.norm(points[second] - points[first])
                for first, second in combinations(range(4), 2)
            ]
        )
        max_edge_length = float(edge_lengths.max())
        min_edge_length = float(edge_lengths.min())
        scaled_volume = (
            6.0 * sqrt(2.0) * volume / max_edge_length**3
            if max_edge_length > 0.0
            else 0.0
        )
        if not isfinite(scaled_volume) or scaled_volume <= degeneracy_tolerance:
            raise DegenerateElementError(
                entity_id=entity_id,
                signed_volume=signed_volume,
                scaled_volume=scaled_volume,
                tolerance=degeneracy_tolerance,
            )

        shape_gradients = np.linalg.solve(jacobian.T, _REFERENCE_GRADIENTS.T).T
        centered = points - points.mean(axis=0)
        condition_number = float(np.linalg.cond(centered))
        quality = Tet4Quality(
            min_edge_length=min_edge_length,
            max_edge_length=max_edge_length,
            aspect_ratio=max_edge_length / min_edge_length,
            scaled_volume=scaled_volume,
            condition_number=condition_number,
        )
        oriented_node_order = (0, 1, 2, 3) if determinant > 0.0 else (0, 2, 1, 3)

        self._coordinates = points
        self._entity_id = entity_id
        self._jacobian = _frozen_float_array(jacobian, shape=(3, 3))
        self._signed_volume = signed_volume
        self._volume = volume
        self._shape_gradients = _frozen_float_array(shape_gradients, shape=(4, 3))
        self._quality = quality
        self._oriented_node_order = oriented_node_order
        self._oriented_coordinates = _frozen_float_array(
            points[list(oriented_node_order)], shape=(4, 3)
        )
        self._faces: tuple[Tet4Face, ...] | None = None

    def _build_faces(self) -> tuple[Tet4Face, ...]:
        faces: list[Tet4Face] = []
        for opposite_node, local_nodes in enumerate(_FACE_NODES):
            face_points = self._coordinates[list(local_nodes)]
            cross_product = np.cross(
                face_points[1] - face_points[0], face_points[2] - face_points[0]
            )
            cross_norm = float(np.linalg.norm(cross_product))
            area = cross_norm / 2.0
            normal = cross_product / cross_norm
            centroid = face_points.mean(axis=0)
            oriented_nodes = local_nodes
            if np.dot(normal, self._coordinates[opposite_node] - centroid) > 0.0:
                normal = -normal
                oriented_nodes = (local_nodes[0], local_nodes[2], local_nodes[1])
            faces.append(
                Tet4Face(
                    local_nodes=oriented_nodes,
                    opposite_node=opposite_node,
                    area=area,
                    unit_normal=_frozen_float_array(normal, shape=(3,)),
                    centroid=_frozen_float_array(centroid, shape=(3,)),
                )
            )
        return tuple(faces)

    @property
    def coordinates(self) -> FloatArray:
        """Immutable input coordinates with shape ``(4, 3)`` in metres."""
        return self._coordinates

    @property
    def entity_id(self) -> int | None:
        """Optional mesh entity ID used in diagnostics."""
        return self._entity_id

    @property
    def signed_volume(self) -> float:
        """Oriented input volume in cubic metres."""
        return self._signed_volume

    @property
    def volume(self) -> float:
        """Positive element volume in cubic metres."""
        return self._volume

    @property
    def shape_gradients(self) -> FloatArray:
        """Constant physical gradients as rows ``[grad(N_0), ..., grad(N_3)]``."""
        return self._shape_gradients

    @property
    def gradient_matrix(self) -> FloatArray:
        """Constant scalar-field operator ``G`` with shape ``(3, 4)``."""
        return self._shape_gradients.T

    @property
    def quality(self) -> Tet4Quality:
        """Dimensionless and edge-length quality indicators."""
        return self._quality

    @property
    def oriented_node_order(self) -> tuple[int, int, int, int]:
        """Local permutation that yields a positive orientation."""
        return self._oriented_node_order

    @property
    def oriented_coordinates(self) -> FloatArray:
        """Coordinates permuted to positive orientation."""
        return self._oriented_coordinates

    @property
    def faces(self) -> tuple[Tet4Face, ...]:
        """Four outward faces, computed on first access and retained per instance.

        Volume-only operations do not construct faces. The unchanged geometric
        calculation uses owned immutable coordinates; no external invalidation
        or shared cache is needed. Concurrent first access may duplicate work.
        """
        if self._faces is None:
            self._faces = self._build_faces()
        return self._faces

    def shape_values(self, point: ArrayLike) -> FloatArray:
        """Evaluate the four affine shape functions at a physical point."""
        physical_point = _point_array(point)
        reference_coordinates = np.linalg.solve(
            self._jacobian, physical_point - self._coordinates[0]
        )
        values = np.concatenate(
            ([1.0 - float(reference_coordinates.sum())], reference_coordinates)
        )
        return _frozen_float_array(values, shape=(4,))

    def oriented_connectivity(self, connectivity: ArrayLike) -> IndexArray:
        """Return four node indices permuted consistently with the geometry."""
        raw = np.asarray(connectivity)
        if raw.shape != (4,):
            raise ShapeError(field="connectivity", expected="(4,)", actual=raw.shape)
        if raw.dtype.kind not in "iu":
            raise ElementGeometryError(
                field="connectivity", reason="expected integer node indices"
            )
        normalized = np.asarray(raw, dtype=np.int64)
        if len(np.unique(normalized)) != 4:
            raise ElementGeometryError(
                field="connectivity", reason="node indices must be unique"
            )
        return _frozen_index_array(
            normalized[list(self._oriented_node_order)], shape=(4,)
        )
