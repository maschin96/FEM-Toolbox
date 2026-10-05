"""Quadratic Lagrange interpolation on affine tetrahedra and planar triangles."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.elements.quadrature import (
    QuadratureRule,
    _rule,
    tetrahedron_quadrature,
    triangle_quadrature,
)
from femantis.elements.tet4 import (
    DEFAULT_DEGENERACY_TOLERANCE,
    Tet4Geometry,
    Tet4Quality,
    _frozen_float_array,
    _frozen_index_array,
    _point_array,
)
from femantis.errors import ElementGeometryError, ShapeError

FloatArray = NDArray[np.float64]


def _coordinates(value: ArrayLike, count: int) -> FloatArray:
    raw = np.asarray(value)
    if raw.shape != (count, 3):
        raise ShapeError(
            field="coordinates", expected=f"({count}, 3)", actual=raw.shape
        )
    if raw.dtype.kind not in "iuf" or not np.isfinite(raw).all():
        raise ElementGeometryError(
            field="coordinates", reason="expected finite real coordinates"
        )
    return _frozen_float_array(raw, shape=(count, 3))


def _values(barycentric: FloatArray, edges: tuple[tuple[int, int], ...]) -> FloatArray:
    return _frozen_float_array(
        np.concatenate(
            (
                barycentric * (2 * barycentric - 1),
                [4 * barycentric[a] * barycentric[b] for a, b in edges],
            )
        ),
        shape=(len(barycentric) + len(edges),),
    )


def _gradients(
    barycentric: FloatArray, gradients: FloatArray, edges: tuple[tuple[int, int], ...]
) -> FloatArray:
    return _frozen_float_array(
        np.vstack(
            (
                (4 * barycentric - 1)[:, None] * gradients,
                [
                    4 * (barycentric[a] * gradients[b] + barycentric[b] * gradients[a])
                    for a, b in edges
                ],
            )
        ),
        shape=(len(barycentric) + len(edges), 3),
    )


@dataclass(frozen=True, slots=True)
class Tet10Face:
    """Six outward-oriented local nodes, opposite corner and affine face data."""

    local_nodes: tuple[int, ...]
    opposite_node: int
    area: float
    unit_normal: FloatArray
    centroid: FloatArray


class Tet10Geometry:
    """Straight Tet10 with VTK node order: corners then 01,12,20,03,13,23.

    Coordinates are (10,3) in metres. Shape values and point-dependent physical
    gradients retain input ordering, including negative orientation. Curved or
    non-midpoint edges are rejected at relative edge tolerance 1e-10. Corner
    degeneracy uses the Tet4 dimensionless criterion. Evaluation outside the
    element is polynomial extrapolation, not an inclusion test.
    """

    __slots__ = ("_coordinates", "_edges", "_faces", "_linear", "_order")

    def __init__(
        self,
        coordinates: ArrayLike,
        *,
        entity_id: int | None = None,
        degeneracy_tolerance: float = DEFAULT_DEGENERACY_TOLERANCE,
    ) -> None:
        from femantis.mesh.quadratic import (
            TET10_EDGES,
            edge_nodes,
            reordered_tet10,
            require_straight,
        )

        self._coordinates = _coordinates(coordinates, 10)
        self._linear = Tet4Geometry(
            self._coordinates[:4],
            entity_id=entity_id,
            degeneracy_tolerance=degeneracy_tolerance,
        )
        nodes = np.arange(10, dtype=np.int64)
        require_straight(self._coordinates, nodes)
        self._edges = TET10_EDGES
        order = reordered_tet10(
            nodes, np.asarray(self._linear.oriented_node_order, dtype=np.int64)
        )
        self._order = tuple(int(i) for i in order)
        lookup = edge_nodes(nodes)
        faces = []
        for face in self._linear.faces:
            corners = face.local_nodes
            mids = tuple(
                lookup[tuple(sorted((corners[a], corners[b])))]
                for a, b in ((0, 1), (1, 2), (2, 0))
            )
            faces.append(
                Tet10Face(
                    (*corners, *mids),
                    face.opposite_node,
                    face.area,
                    face.unit_normal,
                    face.centroid,
                )
            )
        self._faces = tuple(faces)

    @property
    def coordinates(self) -> FloatArray:
        """Immutable (10,3) input coordinates in metres."""
        return self._coordinates

    @property
    def volume(self) -> float:
        """Positive volume in m³."""
        return self._linear.volume

    @property
    def signed_volume(self) -> float:
        """Signed volume of the input corner ordering in m³."""
        return self._linear.signed_volume

    @property
    def quality(self) -> Tet4Quality:
        """Affine corner-simplex quality, independent of interpolation degree."""
        return self._linear.quality

    @property
    def oriented_node_order(self) -> tuple[int, ...]:
        """Ten-node permutation giving positive corner orientation."""
        return self._order

    @property
    def oriented_coordinates(self) -> FloatArray:
        """Immutable coordinates consistently permuted to positive orientation."""
        return _frozen_float_array(self.coordinates[list(self._order)], shape=(10, 3))

    @property
    def faces(self) -> tuple[Tet10Face, ...]:
        """Four outward-oriented six-node faces."""
        return self._faces

    def oriented_connectivity(self, connectivity: ArrayLike) -> NDArray[np.int64]:
        """Permute ten unique integer indices with the geometry; retain IDs."""
        raw = np.asarray(connectivity)
        if raw.shape != (10,):
            raise ShapeError(field="connectivity", expected="(10,)", actual=raw.shape)
        if raw.dtype.kind not in "iu" or len(np.unique(raw)) != 10:
            raise ElementGeometryError(
                field="connectivity", reason="expected ten unique integer node indices"
            )
        return _frozen_index_array(raw[list(self._order)], shape=(10,))

    def shape_values(self, point: ArrayLike) -> FloatArray:
        """Ten dimensionless shape values at a physical (3,) point in metres."""
        return _values(self._linear.shape_values(point), self._edges)

    def shape_gradients(self, point: ArrayLike) -> FloatArray:
        """Physical gradients (10,3), one row per input node, in 1/m."""
        return _gradients(
            self._linear.shape_values(point), self._linear.shape_gradients, self._edges
        )

    def gradient_matrix(self, point: ArrayLike) -> FloatArray:
        """Scalar field gradient operator (3,10) in 1/m."""
        return self.shape_gradients(point).T

    def quadrature(self) -> QuadratureRule:
        """Physical points and positive m³ weights, exact through degree five."""
        rule = tetrahedron_quadrature()
        corners = self.coordinates[:4]
        points = corners[0] + rule.points @ (corners[1:] - corners[0])
        return _rule(points, rule.weights * (6 * self.volume), rule.degree)


class Triangle6Geometry:
    """Planar straight Triangle6: corners 0,1,2 then edge nodes 01,12,20.

    Coordinates are (6,3) in metres; input corner order defines the normal.
    Gradients are tangential vectors in ambient 3D. Off-plane evaluation is
    rejected at 1e-10 of the maximum edge length. In-plane extrapolation is
    allowed. Curved edges and relative area <= 1e-12 are rejected.
    """

    __slots__ = (
        "_area",
        "_coordinates",
        "_edges",
        "_gradients",
        "_inverse",
        "_normal",
        "_scale",
    )

    def __init__(self, coordinates: ArrayLike) -> None:
        from femantis.mesh.quadratic import TRI6_EDGES, require_straight

        self._coordinates = _coordinates(coordinates, 6)
        corners = self._coordinates[:3]
        jacobian = (corners[1:] - corners[0]).T
        cross = np.cross(jacobian[:, 0], jacobian[:, 1])
        measure = float(np.linalg.norm(cross))
        scale = max(
            float(np.linalg.norm(corners[b] - corners[a])) for a, b in TRI6_EDGES
        )
        if (
            not np.isfinite(measure)
            or scale == 0
            or measure / scale**2 <= DEFAULT_DEGENERACY_TOLERANCE
        ):
            raise ElementGeometryError(
                field="coordinates", reason="degenerate Triangle6"
            )
        require_straight(self._coordinates, np.arange(6, dtype=np.int64))
        q, r = np.linalg.qr(jacobian, mode="reduced")
        self._inverse = np.linalg.solve(r, q.T)
        self._gradients = np.vstack((-self._inverse.sum(axis=0), self._inverse))
        self._area, self._scale = measure / 2, scale
        self._normal = _frozen_float_array(cross / measure, shape=(3,))
        self._edges = TRI6_EDGES

    @property
    def coordinates(self) -> FloatArray:
        """Immutable (6,3) input coordinates in metres."""
        return self._coordinates

    @property
    def area(self) -> float:
        """Positive surface measure in m²."""
        return self._area

    @property
    def unit_normal(self) -> FloatArray:
        """Immutable unit normal following the input corner ordering."""
        return self._normal

    def _barycentric(self, point: ArrayLike) -> FloatArray:
        delta = _point_array(point) - self.coordinates[0]
        if abs(float(delta @ self.unit_normal)) > 1e-10 * self._scale:
            raise ElementGeometryError(
                field="point", reason="point lies outside Triangle6 plane"
            )
        uv = self._inverse @ delta
        return np.concatenate(([1 - float(uv.sum())], uv))

    def shape_values(self, point: ArrayLike) -> FloatArray:
        """Six dimensionless values at a physical (3,) point in the plane."""
        return _values(self._barycentric(point), self._edges)

    def shape_gradients(self, point: ArrayLike) -> FloatArray:
        """Tangential physical gradients (6,3) in 1/m, in input node order."""
        return _gradients(self._barycentric(point), self._gradients, self._edges)

    def quadrature(self) -> QuadratureRule:
        """Physical points and m² weights, exact through degree four."""
        rule = triangle_quadrature()
        corners = self.coordinates[:3]
        points = corners[0] + rule.points @ (corners[1:] - corners[0])
        return _rule(points, rule.weights * (2 * self.area), rule.degree)
