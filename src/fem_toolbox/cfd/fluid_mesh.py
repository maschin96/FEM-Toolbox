"""Experimental conforming 2D Quad4 fluid geometry and wall-resolution gates."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, NoReturn

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.spatial import cKDTree  # type: ignore[import-untyped]

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.errors import ModelSchemaError


def _fail(path: str, reason: str) -> NoReturn:
    raise ModelSchemaError(path=path, reason=reason)


def _keys(value: Any, required: set[str], optional: set[str], path: str) -> None:
    if not isinstance(value, Mapping):
        _fail(path, "expected a table")
    if required - value.keys() or value.keys() - required - optional:
        _fail(path, f"required keys {sorted(required)}; optional {sorted(optional)}")


def _float_array(value: ArrayLike, path: str) -> NDArray[np.float64]:
    try:
        if any(
            isinstance(x, (bool, np.bool_))
            for x in np.asarray(value, dtype=object).flat
        ):
            _fail(path, "boolean values are not real numbers")
        raw = np.asarray(value)
        if raw.dtype.kind not in "iuf" or raw.dtype.kind == "b":
            _fail(path, "expected real numbers")
        a = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ModelSchemaError(path=path, reason="invalid numeric array") from error
    if not np.all(np.isfinite(a)):
        _fail(path, "expected finite float64 values")
    return a


def _indices(
    value: ArrayLike, width: int | None, size: int, path: str
) -> NDArray[np.int64]:
    try:
        if any(
            isinstance(x, (bool, np.bool_))
            for x in np.asarray(value, dtype=object).flat
        ):
            _fail(path, "boolean values are not integer indices")
        raw = np.asarray(value)
        if raw.dtype.kind not in "iu" or (
            raw.size and (raw.min() < 0 or raw.max() >= size)
        ):
            _fail(path, f"expected integer indices in [0,{size})")
        a = np.asarray(raw, dtype=np.int64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ModelSchemaError(path=path, reason="invalid index array") from error
    if (width is None and a.ndim != 1) or (
        width is not None and (a.ndim != 2 or a.shape[1] != width)
    ):
        _fail(
            path,
            f"expected {'1D indices' if width is None else f'(n,{width}) indices'}",
        )
    if not a.size:
        _fail(path, "empty selections/connectivities are unsupported")
    if width is not None and any(len(set(row)) != width for row in a):
        _fail(path, "repeated node in entity")
    return a


def _freeze_float(a: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.frombuffer(a.tobytes(), dtype=np.float64).reshape(a.shape)


def _freeze_int(a: NDArray[np.int64]) -> NDArray[np.int64]:
    return np.frombuffer(a.tobytes(), dtype=np.int64).reshape(a.shape)


def _groups(value: Any, size: int, path: str) -> Mapping[str, NDArray[np.int64]]:
    if not isinstance(value, Mapping):
        _fail(path, "expected named selections")
    groups = {}
    for name, indices in value.items():
        if not isinstance(name, str) or not name.strip():
            _fail(path, "nonempty group names required")
        a = _indices(indices, None, size, f"{path}.{name}")
        if len(np.unique(a)) != len(a):
            _fail(f"{path}.{name}", "duplicate indices")
        groups[name] = _freeze_int(a)
    return MappingProxyType(groups)


def _cross(a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def _cell_geometry(
    points: NDArray[np.float64], cells: NDArray[np.int64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    areas = []
    centres = []
    with np.errstate(over="raise", invalid="raise", divide="raise", under="ignore"):
        try:
            for i, nodes in enumerate(cells):
                p = points[nodes]
                edges = np.roll(p, -1, axis=0) - p
                turns = _cross(edges, np.roll(edges, -1, axis=0))
                if np.any(turns <= 0):
                    _fail(
                        f"quad4[{i}]", "positive counterclockwise convex cell required"
                    )
                # Bilinear Jacobian at all corners and 2x2 Gauss points.
                g = 1 / np.sqrt(3)
                signs = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
                local = p - p[0]
                for xi, eta in [*signs, (-g, -g), (g, -g), (g, g), (-g, g)]:
                    dx = signs[:, 0] * (1 + signs[:, 1] * eta) / 4
                    dy = signs[:, 1] * (1 + signs[:, 0] * xi) / 4
                    if _cross(dx @ local, dy @ local) <= 0:
                        _fail(f"quad4[{i}]", "nonpositive bilinear Jacobian")
                a1 = float(_cross(local[1], local[2])) / 2
                a2 = float(_cross(local[2], local[3])) / 2
                area = a1 + a2
                centre = p[0] + (
                    a1 * (local[1] + local[2]) + a2 * (local[2] + local[3])
                ) / (3 * area)
                if (
                    area <= 0
                    or not np.isfinite(area)
                    or not np.all(np.isfinite(centre))
                ):
                    _fail(f"quad4[{i}]", "unrepresentable area/centroid")
                areas.append(area)
                centres.append(centre)
        except FloatingPointError as error:
            raise ModelSchemaError(
                path="quad4", reason="geometry exceeds float64 range"
            ) from error
    return np.array(areas), np.array(centres)


def _check_overlap(points: NDArray[np.float64], cells: NDArray[np.int64]) -> None:
    """Sweep bounding boxes then separating axes; touching faces are allowed."""
    polygons = points[cells]
    lower, upper = polygons.min(axis=1), polygons.max(axis=1)
    active: list[int] = []
    for raw_index in np.argsort(lower[:, 0]):
        i = int(raw_index)
        active = [j for j in active if upper[j, 0] > lower[i, 0]]
        for j in active:
            if upper[j, 1] <= lower[i, 1] or upper[i, 1] <= lower[j, 1]:
                continue
            a, b = polygons[i] - polygons[i, 0], polygons[j] - polygons[i, 0]
            separated = False
            for polygon in (a, b):
                edges = np.roll(polygon, -1, axis=0) - polygon
                normals = np.column_stack((edges[:, 1], -edges[:, 0]))
                normals /= np.linalg.norm(normals, axis=1)[:, None]
                ap, bp = a @ normals.T, b @ normals.T
                tolerance = (
                    64 * np.finfo(float).eps * max(np.max(np.abs(a)), np.max(np.abs(b)))
                )
                if np.any(
                    np.minimum(ap.max(axis=0), bp.max(axis=0))
                    - np.maximum(ap.min(axis=0), bp.min(axis=0))
                    <= tolerance
                ):
                    separated = True
                    break
            if not separated:
                _fail(f"quad4[{i}]", f"cell overlaps quad4[{j}]")
        active.append(i)


@dataclass(frozen=True, slots=True)
class WallQuality:
    """Distances in metres; boundary rows and growth refer to selected wall faces."""

    boundary_indices: NDArray[np.int64]
    wall_distance_m: NDArray[np.float64]
    first_cell_distance_m: NDArray[np.float64]
    layer_growth: NDArray[np.float64]

    def y_plus(
        self, *, density_kg_m3: float, viscosity_Pa_s: float, wall_shear_Pa: ArrayLike
    ) -> NDArray[np.float64]:
        """Evaluate actual wall shear, reject y+>1; no estimated shear default."""
        rho = _scalar(density_kg_m3, "density_kg_m3")
        mu = _scalar(viscosity_Pa_s, "viscosity_Pa_s")
        shear = _float_array(wall_shear_Pa, "wall_shear_Pa")
        if rho <= 0 or mu <= 0 or shear.shape != self.first_cell_distance_m.shape:
            _fail(
                "wall_shear_Pa",
                "positive SI fluid and one shear value per wall face required",
            )
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            values = (
                self.first_cell_distance_m * np.sqrt(np.abs(shear) / rho) / (mu / rho)
            )
        bad = np.flatnonzero((values > 1) | ~np.isfinite(values))
        if len(bad):
            i = int(bad[0])
            _fail(
                f"wall[{self.boundary_indices[i]}].y_plus",
                f"unresolved y_plus={values[i]}; maximum 1",
            )
        return _freeze_float(values)


@dataclass(frozen=True, slots=True, init=False)
class FluidMesh2D:
    """Owned SI Quad4 mesh; complete boundary rows, immutable named selections.

    points=(n,2) metres, quad4=(nc,4) CCW, line2=(nb,2). Indices are
    zero-based. Boundary row order is preserved; normals use cell ownership.
    Validates positive geometry, conformity, overlap and ADR-049 FV gates.
    This is a reference adapter, with no fluid solver or mesh motion.
    """

    points: NDArray[np.float64]
    quad4: NDArray[np.int64]
    line2: NDArray[np.int64]
    boundary_sets: Mapping[str, NDArray[np.int64]]
    volume_sets: Mapping[str, NDArray[np.int64]]
    cell_areas_m2: NDArray[np.float64]
    cell_aspect_ratio: NDArray[np.float64]
    cell_centres_m: NDArray[np.float64]
    faces: NDArray[np.int64]
    face_owners: NDArray[np.int64]
    face_neighbours: NDArray[np.int64]
    face_centres_m: NDArray[np.float64]
    face_normals: NDArray[np.float64]
    face_lengths_m: NDArray[np.float64]
    boundary_face_indices: NDArray[np.int64]
    nonorthogonality_deg: NDArray[np.float64]
    skewness: NDArray[np.float64]

    def __init__(
        self,
        points: ArrayLike,
        quad4: ArrayLike,
        line2: ArrayLike,
        boundary_sets: Mapping[str, ArrayLike],
        volume_sets: Mapping[str, ArrayLike] | None = None,
    ) -> None:
        p = _float_array(points, "points")
        if p.ndim != 2 or p.shape[1] != 2 or len(p) < 4:
            _fail("points", "expected at least four 2D points")
        if len(np.unique(p, axis=0)) != len(p):
            _fail("points", "duplicate coordinates; use shared node indices")
        cells = _indices(quad4, 4, len(p), "quad4")
        boundary = _indices(line2, 2, len(p), "line2")
        if len(cells) > 50000:
            _fail("quad4", "reference adapter limit is 50000 cells")
        if len(np.unique(cells)) != len(p):
            _fail("points", "unused nodes are unsupported")
        if len({tuple(sorted(row)) for row in cells}) != len(cells):
            _fail("quad4", "duplicate cells")
        areas, centres = _cell_geometry(p, cells)
        cell_edges = np.roll(p[cells], -1, axis=1) - p[cells]
        edge_lengths = np.linalg.norm(cell_edges, axis=2)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            aspect = edge_lengths.max(axis=1) / edge_lengths.min(axis=1)
        if not np.all(np.isfinite(aspect)):
            _fail("quad4", "unrepresentable edge aspect ratio")
        _check_overlap(p, cells)
        incidence: dict[tuple[int, int], list[tuple[int, int, int]]] = {}
        for i, row in enumerate(cells):
            for a, b in zip(row, np.roll(row, -1), strict=True):
                key = (min(int(a), int(b)), max(int(a), int(b)))
                incidence.setdefault(key, []).append((i, int(a), int(b)))
        for edge, records in incidence.items():
            if len(records) > 2 or (
                len(records) == 2 and records[0][1:] != records[1][1:][::-1]
            ):
                _fail("quad4", f"nonmanifold or inconsistent face {edge}")
        tree = cKDTree(p)
        for edge in incidence:
            a, b = p[list(edge)]
            delta = b - a
            length = float(np.linalg.norm(delta))
            tolerance = 64 * np.finfo(float).eps * length
            for node in tree.query_ball_point(a + delta / 2, length / 2 + tolerance):
                if node in edge:
                    continue
                t = float((p[node] - a) @ delta / (delta @ delta))
                if 0 < t < 1 and np.linalg.norm(p[node] - a - t * delta) <= tolerance:
                    _fail(
                        "quad4",
                        f"nonconforming face {edge}: node {node} lies inside edge",
                    )
        external = {e for e, owners in incidence.items() if len(owners) == 1}
        raw_boundary = [tuple(sorted(row)) for row in boundary]
        if len(set(raw_boundary)) != len(raw_boundary) or set(raw_boundary) != external:
            _fail(
                "line2",
                "requires exactly all exterior faces, no duplicates or interior faces",
            )
        # A conforming connected domain has one gauge in the initial contract.
        adjacent: dict[int, set[int]] = {i: set() for i in range(len(cells))}
        for records in incidence.values():
            if len(records) == 2:
                i, j = records[0][0], records[1][0]
                adjacent[i].add(j)
                adjacent[j].add(i)
        visited, queue = {0}, [0]
        while queue:
            for j in adjacent[queue.pop()] - visited:
                visited.add(j)
                queue.append(j)
        if len(visited) != len(cells):
            _fail("quad4", "one face-connected domain required")
        edges = list(incidence)
        faces = np.array(
            [owners[0][1:] for owners in incidence.values()], dtype=np.int64
        )
        owners = np.array([v[0][0] for v in incidence.values()], dtype=np.int64)
        neighbours = np.array(
            [v[1][0] if len(v) == 2 else -1 for v in incidence.values()], dtype=np.int64
        )
        delta = p[faces[:, 1]] - p[faces[:, 0]]
        lengths = np.linalg.norm(delta, axis=1)
        normals = np.column_stack((delta[:, 1], -delta[:, 0])) / lengths[:, None]
        midpoints = p[faces[:, 0]] + delta / 2
        direction = midpoints - centres[owners]
        inner = neighbours >= 0
        direction[inner] = centres[neighbours[inner]] - centres[owners[inner]]
        distance = np.einsum("ij,ij->i", direction, normals)
        if np.any(distance <= 0) or not np.all(np.isfinite(direction)):
            _fail("faces", "nonpositive owner-normal distance")
        angle = np.degrees(
            np.arccos(np.clip(distance / np.linalg.norm(direction, axis=1), -1, 1))
        )
        # Intersection of the centre-to-centre line with each face line.
        factor = np.einsum("ij,ij->i", midpoints - centres[owners], normals) / distance
        intersections = centres[owners] + factor[:, None] * direction
        skew = np.linalg.norm(midpoints - intersections, axis=1) / np.linalg.norm(
            direction, axis=1
        )
        if not np.all(np.isfinite(angle)) or not np.all(np.isfinite(skew)):
            _fail("faces", "unrepresentable quality measures")
        for bad_face in np.flatnonzero((angle > 60) | (skew > 0.25)):
            _fail(
                f"faces[{bad_face}]",
                f"FV quality: nonorthogonality={angle[bad_face]} deg, skewness={skew[bad_face]}",
            )
        edge_indices = {edge: i for i, edge in enumerate(edges)}
        ids = np.array([edge_indices[e] for e in raw_boundary], dtype=np.int64)
        for name, value in {
            "points": p,
            "cell_areas_m2": areas,
            "cell_aspect_ratio": aspect,
            "cell_centres_m": centres,
            "face_centres_m": midpoints,
            "face_normals": normals,
            "face_lengths_m": lengths,
            "nonorthogonality_deg": angle,
            "skewness": skew,
        }.items():
            object.__setattr__(self, name, _freeze_float(value))
        for name, indices in {
            "quad4": cells,
            "line2": boundary,
            "faces": faces,
            "face_owners": owners,
            "face_neighbours": neighbours,
            "boundary_face_indices": ids,
        }.items():
            object.__setattr__(self, name, _freeze_int(indices))
        object.__setattr__(
            self,
            "boundary_sets",
            _groups(boundary_sets, len(boundary), "boundary_sets"),
        )
        object.__setattr__(
            self,
            "volume_sets",
            _groups(
                volume_sets if volume_sets is not None else {},
                len(cells),
                "volume_sets",
            ),
        )

    @property
    def boundary_normals(self) -> NDArray[np.float64]:
        """Outward normals in the input line2 row order."""
        return _freeze_float(self.face_normals[self.boundary_face_indices])

    @property
    def max_nonorthogonality_deg(self) -> float:
        return float(self.nonorthogonality_deg.max())

    @property
    def max_skewness(self) -> float:
        return float(self.skewness.max())

    @classmethod
    def from_document(cls, document: Mapping[str, Any]) -> FluidMesh2D:
        _keys(
            document,
            {
                "fluid_mesh_schema_version",
                "length_unit",
                "points",
                "quad4",
                "line2",
                "boundary_sets",
            },
            {"volume_sets"},
            "fluid_mesh",
        )
        if (
            type(document["fluid_mesh_schema_version"]) is not int
            or document["fluid_mesh_schema_version"] != 1
        ):
            _fail("fluid_mesh_schema_version", "expected integer 1")
        if document["length_unit"] != "m":
            _fail("length_unit", "only explicit SI metres supported")
        return cls(
            document["points"],
            document["quad4"],
            document["line2"],
            document["boundary_sets"],
            document.get("volume_sets"),
        )

    @classmethod
    def from_file(cls, path: str | Path) -> FluidMesh2D:
        """Read UTF-8 schema-1 JSON; no implicit unit conversion or repair."""
        try:
            document = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise ModelSchemaError(path=str(path), reason=str(error)) from error
        return cls.from_document(document)

    def to_document(self) -> dict[str, Any]:
        """Return detached JSON-ready input data, preserving group order/indices."""
        return {
            "fluid_mesh_schema_version": 1,
            "length_unit": "m",
            "points": self.points.tolist(),
            "quad4": self.quad4.tolist(),
            "line2": self.line2.tolist(),
            "boundary_sets": {k: v.tolist() for k, v in self.boundary_sets.items()},
            "volume_sets": {k: v.tolist() for k, v in self.volume_sets.items()},
        }

    def wall_quality(self, wall_groups: Sequence[str]) -> WallQuality:
        """Compute exact segment distances and opposite-face layer growth.

        At least one named wall selection is required. For the selected
        structured Quad4 rows, trace opposite faces to the domain boundary;
        outward normal layer growth greater than 1.2 is rejected with wall ID.
        """
        if isinstance(wall_groups, str) or not wall_groups:
            _fail("walls", "named wall groups required")
        for name in wall_groups:
            if name not in self.boundary_sets:
                _fail("walls", f"unknown wall group {name!r}")
        ids = np.unique(
            np.concatenate([self.boundary_sets[name] for name in wall_groups])
        )
        distance = np.full(len(self.quad4), np.inf)
        first = []
        growth = []
        face_lookup = {tuple(sorted(row)): i for i, row in enumerate(self.faces)}
        for boundary_id in ids:
            face = int(self.boundary_face_indices[boundary_id])
            a, b = self.points[self.faces[face]]
            edge = b - a
            relative = self.cell_centres_m - a
            t = np.clip(relative @ edge / (edge @ edge), 0, 1)
            distance = np.minimum(
                distance, np.linalg.norm(relative - t[:, None] * edge, axis=1)
            )
            cell = int(self.face_owners[face])
            normal = self.face_normals[face]
            first.append(
                float((self.face_centres_m[face] - self.cell_centres_m[cell]) @ normal)
            )
            ratios, previous, seen = [1.0], None, set()
            while cell not in seen:
                seen.add(cell)
                row = self.quad4[cell]
                edge_ids = [
                    face_lookup[tuple(sorted((int(row[k]), int(row[(k + 1) % 4]))))]
                    for k in range(4)
                ]
                opposite = edge_ids[(edge_ids.index(face) + 2) % 4]
                thickness = abs(
                    float(
                        (self.face_centres_m[opposite] - self.face_centres_m[face])
                        @ normal
                    )
                )
                if thickness <= 0 or not np.isfinite(thickness):
                    _fail(
                        f"wall[{boundary_id}]", "nonpositive projected layer thickness"
                    )
                if previous is not None:
                    ratios.append(thickness / previous)
                previous = thickness
                owner, neighbour = (
                    int(self.face_owners[opposite]),
                    int(self.face_neighbours[opposite]),
                )
                next_cell = neighbour if owner == cell else owner
                if next_cell < 0:
                    break
                cell, face = next_cell, opposite
            else:
                _fail(f"wall[{boundary_id}]", "cyclic wall-normal layer trace")
            maximum = max(ratios)
            if maximum > 1.2 + 64 * np.finfo(float).eps:
                _fail(
                    f"wall[{boundary_id}].growth", f"layer growth {maximum} exceeds 1.2"
                )
            growth.append(maximum)
        if np.any(distance <= 0) or not np.all(np.isfinite(distance)):
            _fail("wall_distance_m", "positive finite cell-to-wall distances required")
        return WallQuality(
            _freeze_int(ids),
            _freeze_float(distance),
            _freeze_float(np.array(first)),
            _freeze_float(np.array(growth)),
        )
