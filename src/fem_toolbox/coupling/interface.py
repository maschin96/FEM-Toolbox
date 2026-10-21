"""Conservative P1 overlap quadrature on extruded planar fluid segments.

Coordinates are SI metres, traction Pa, integrated quadrature/nodal forces N.
The same interpolation matrix and its transpose define virtual-work pairing.
No nearest-neighbour assignment, extrapolation or geometry filling is used.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]

from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.errors import ValidationError
from fem_toolbox.mesh import Mesh
from fem_toolbox.mesh.topology import inspect_boundary_topology

FloatArray = NDArray[np.float64]


def _fail(message: str) -> None:
    raise ValidationError(message, field="fsi_interface")


def _real(value: ArrayLike, shape: tuple[int, ...], name: str) -> FloatArray:
    if any(
        isinstance(v, (bool, np.bool_)) for v in np.asarray(value, dtype=object).flat
    ):
        _fail(f"{name}: boolean values are not real input")
    raw = np.asarray(value)
    if raw.shape != shape or raw.dtype.kind not in "iuf":
        _fail(f"{name}: expected finite real array {shape}")
    result = np.asarray(raw, dtype=float)
    if not np.all(np.isfinite(result)):
        _fail(f"{name}: expected finite real array {shape}")
    return result


def _positive(value: float, name: str) -> float:
    number = _real([value], (1,), name)[0]
    if number <= 0:
        _fail(f"{name}: expected positive finite scalar")
    return float(number)


def _freeze(value: FloatArray) -> FloatArray:
    return np.frombuffer(value.tobytes(), dtype=value.dtype).reshape(value.shape)


def _area(poly: FloatArray) -> float:
    if len(poly) < 3:
        return 0.0
    # Translate before the shoelace sum to avoid large-origin cancellation.
    p = poly - poly[0]
    return float(
        np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1)) / 2
    )


def _cross(a: FloatArray, b: FloatArray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def _clip(poly: FloatArray, boundary: FloatArray) -> FloatArray:
    """Sutherland--Hodgman intersection, convex CCW clipping polygon."""
    result = poly
    for a, b in zip(boundary, np.roll(boundary, -1, axis=0), strict=True):
        if not len(result):
            break
        edge = b - a
        output = []
        for p, q in zip(result, np.roll(result, -1, axis=0), strict=True):
            dp, dq = _cross(edge, p - a), _cross(edge, q - a)
            pin, qin = dp >= 0, dq >= 0
            if pin:
                output.append(p)
            if pin != qin:
                output.append(p + dp / (dp - dq) * (q - p))
        result = np.asarray(output, dtype=float).reshape(-1, 2)
    return result


@dataclass(frozen=True, slots=True)
class TransferAudit:
    """Raw conservation errors (N, N m, J) and normalized project gate."""

    force_error_N: FloatArray
    moment_error_N_m: FloatArray
    virtual_work_error_J: float
    maximum_relative_error: float


@dataclass(frozen=True, slots=True)
class InterfaceLoads:
    """Net structural nodal load, separated fluid/external contributions."""

    nodal_forces_N: FloatArray
    fluid_force_N: FloatArray
    external_pressure_force_N: FloatArray
    audit: TransferAudit


@dataclass(frozen=True, slots=True)
class InterfaceMap:
    """Owned quadrature and scalar P1 matrix W; vector H=kron(W,I3).

    Build using :func:`build_interface_map`. The matrix property returns a copy;
    public numerical arrays use immutable bytes-backed storage. Arbitrary 3D
    virtual displacements may be used for work audits. Physical mesh motion
    additionally requires a span-invariant planar trace via the endpoint method.
    """

    quadrature_points_m: FloatArray
    quadrature_areas_m2: FloatArray
    segment_indices: NDArray[np.int64]
    structure_points_m: FloatArray
    segments_m: FloatArray
    fluid_normals: FloatArray
    width_m: float
    span_origin_m: float
    length_scale_m: float
    maximum_relative_coverage_error: float
    _weights: sparse.csr_matrix
    _triangles: NDArray[np.int64]

    @property
    def interpolation_matrix(self) -> sparse.csr_matrix:
        """Return an independent vector interpolation matrix H (3nq,3ns)."""
        return sparse.kron(self._weights, sparse.eye(3), format="csr")

    def interpolate(self, displacement_m: ArrayLike) -> FloatArray:
        """P1 structural displacement [nnode,3] to quadrature [nq,3], metres."""
        u = _real(displacement_m, self.structure_points_m.shape, "displacement_m")
        return _freeze(np.asarray(self._weights @ u))

    def transfer_forces(self, integrated_forces_N: ArrayLike) -> FloatArray:
        """Apply H transpose to integrated quadrature forces [nq,3] in N."""
        force = _real(
            integrated_forces_N, self.quadrature_points_m.shape, "integrated_forces_N"
        )
        result = np.asarray(self._weights.T @ force)
        if not np.all(np.isfinite(result)):
            _fail("force transfer overflow")
        return _freeze(result)

    def audit(
        self,
        integrated_forces_N: ArrayLike,
        virtual_displacement_m: ArrayLike | None = None,
    ) -> TransferAudit:
        """Check force, origin moments and optional virtual work, unchanged gates.

        Relative gate is 1e-10, with absolute floors 1e-12 N, N m and J.
        Normalization includes absolute summed contributions to handle cancellation.
        """
        f = _real(
            integrated_forces_N, self.quadrature_points_m.shape, "integrated_forces_N"
        )
        nodal = self.transfer_forces(f)
        force_error = nodal.sum(axis=0) - f.sum(axis=0)
        mq = np.cross(self.quadrature_points_m, f)
        mn = np.cross(self.structure_points_m, nodal)
        moment_error = mn.sum(axis=0) - mq.sum(axis=0)
        u = (
            np.zeros_like(nodal)
            if virtual_displacement_m is None
            else _real(virtual_displacement_m, nodal.shape, "virtual_displacement_m")
        )
        wn, wq = u * nodal, self.interpolate(u) * f
        work_error = float(wn.sum() - wq.sum())
        scales = [
            max(float(np.abs(f).sum()), float(np.abs(nodal).sum())),
            max(float(np.abs(mq).sum()), float(np.abs(mn).sum())),
            max(float(np.abs(wn).sum()), float(np.abs(wq).sum())),
        ]
        errors = [
            float(np.max(abs(force_error))),
            float(np.max(abs(moment_error))),
            abs(work_error),
        ]
        # Absolute floors correspond to an effective relative denominator 0.01.
        relative = max(
            error / max(scale, 0.01)
            for error, scale in zip(errors, scales, strict=True)
        )
        if not np.isfinite(relative) or relative > 1e-10:
            _fail(
                f"conservation/work gate failed: {relative=}, force={force_error}, moment={moment_error}, work={work_error}"
            )
        return TransferAudit(
            _freeze(force_error), _freeze(moment_error), work_error, relative
        )

    def loads_from_fluid_traction(
        self, traction_Pa: ArrayLike, *, external_pressure_Pa: float = 0.0
    ) -> InterfaceLoads:
        """Convert sigma_f*n_f [nsegment,2] to net structural load.

        t_structure=-t_fluid-p_external*n_f; pressure is absolute relative to
        the same datum as p_external. External pressure and fluid forces remain
        separate. Piecewise constant FV traction is integrated on overlap areas.
        """
        t = _real(traction_Pa, (len(self.segments_m), 2), "traction_Pa")
        p = _real([external_pressure_Pa], (1,), "external_pressure_Pa")[0]
        rows = self.segment_indices
        fluid = np.zeros_like(self.quadrature_points_m)
        external = np.zeros_like(fluid)
        fluid[:, :2] = -t[rows] * self.quadrature_areas_m2[:, None]
        external[:, :2] = (
            -p * self.fluid_normals[rows] * self.quadrature_areas_m2[:, None]
        )
        total = fluid + external
        return InterfaceLoads(
            self.transfer_forces(total),
            _freeze(fluid.sum(axis=0)),
            _freeze(external.sum(axis=0)),
            self.audit(total),
        )

    def segment_endpoint_displacements(self, displacement_m: ArrayLike) -> FloatArray:
        """Return [nsegment,2 endpoints,2] supported planar displacement.

        Reject span dependence, out-of-plane motion and P1 kinks not resolved
        by the fluid segments. This keeps the reduced trace and adjoint loads
        consistent with exactly the same H, without nearest-point interpolation.
        """
        u = _real(displacement_m, self.structure_points_m.shape, "displacement_m")
        tolerance = 1e-8 * self.length_scale_m
        nodes = np.unique(self._triangles)
        if np.max(abs(u[nodes, 2])) > tolerance:
            _fail("unsupported out-of-plane interface displacement")
        result = np.empty_like(self.segments_m)
        for i, segment in enumerate(self.segments_m):
            tangent = (segment[1] - segment[0]) / np.linalg.norm(
                segment[1] - segment[0]
            )
            for end, point in enumerate(segment):
                samples = []
                for z in self.span_origin_m + self.width_m * np.array([0.0, 0.5, 1.0]):
                    target = np.r_[point, z]
                    values = []
                    for triangle in self._triangles:
                        xyz = self.structure_points_m[triangle]
                        basis = np.column_stack([xyz[1] - xyz[0], xyz[2] - xyz[0]])
                        uv, _, _, _ = np.linalg.lstsq(
                            basis, target - xyz[0], rcond=None
                        )
                        bary = np.r_[1 - uv.sum(), uv]
                        if (
                            np.min(bary) >= -1e-10
                            and np.linalg.norm(xyz.T @ bary - target)
                            <= 1e-10 * self.length_scale_m
                        ):
                            values.append(bary @ u[triangle, :2])
                    if not values:
                        _fail(f"missing endpoint coverage at segment {i}")
                    if (
                        np.max(np.linalg.norm(np.asarray(values) - values[0], axis=1))
                        > tolerance
                    ):
                        _fail("inconsistent structural trace at common edge")
                    samples.append(values[0])
                if (
                    np.max(np.linalg.norm(np.asarray(samples) - samples[1], axis=1))
                    > tolerance
                ):
                    _fail("unsupported span-dependent interface displacement")
                result[i, end] = samples[1]
            rows = np.flatnonzero(self.segment_indices == i)
            length = np.linalg.norm(segment[1] - segment[0])
            tcoord = (
                (self.quadrature_points_m[rows, :2] - segment[0]) @ tangent
            ) / length
            linear = (1 - tcoord[:, None]) * result[i, 0] + tcoord[:, None] * result[
                i, 1
            ]
            actual = self.interpolate(u)[rows, :2]
            if np.max(np.linalg.norm(actual - linear, axis=1)) > tolerance:
                _fail(
                    "unsupported span dependence or unresolved piecewise-linear trace"
                )
        return _freeze(result)


def build_interface_map(
    *,
    segments_m: ArrayLike,
    fluid_normals: ArrayLike,
    width_m: float,
    structure: Mesh,
    facet_set: str,
    length_scale_m: float,
    span_origin_m: float = 0.0,
    max_quadrature_points: int = 200000,
) -> InterfaceMap:
    """Build exact degree-two overlap quadrature and conservative P1 map.

    Fluid segments [n,2,2] extrude over [span_origin,span_origin+width].
    Structure must have named exterior Triangle3 facets owned by Tet4 cells;
    outward structural normals must oppose the supplied fluid unit normals.
    Geometric tolerance is 1e-10 Lref, relative area gate 1e-10. Gaps, double
    coverage and unsupported geometry fail explicitly. No mesh generation/I/O.
    """
    raw = np.asarray(segments_m)
    if raw.ndim != 3 or raw.shape[1:] != (2, 2) or len(raw) == 0:
        _fail("segments_m: expected nonempty (n,2,2)")
    segments = _real(segments_m, raw.shape, "segments_m")
    normals = _real(fluid_normals, (len(segments), 2), "fluid_normals")
    width = _positive(width_m, "width_m")
    scale = _positive(length_scale_m, "length_scale_m")
    origin = float(_real([span_origin_m], (1,), "span_origin_m")[0])
    if (
        isinstance(max_quadrature_points, bool)
        or not isinstance(max_quadrature_points, int)
        or max_quadrature_points <= 0
    ):
        _fail("positive integer quadrature budget required")
    if len(structure.tetra10) or len(structure.triangle6):
        _fail("interface requires Tet4/Triangle3")
    if (
        facet_set not in structure.facet_sets
        or len(structure.facet_sets[facet_set]) == 0
    ):
        _fail("nonempty named structural interface facet set required")
    selected = structure.facet_sets[facet_set]
    topology = inspect_boundary_topology(structure)
    triangles = structure.triangle3[selected]
    xyz = structure.points[triangles]
    lengths = np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1)
    if np.any(lengths <= 1e-10 * scale) or not np.allclose(
        np.linalg.norm(normals, axis=1), 1, atol=1e-12, rtol=0
    ):
        _fail("nonzero segments and unit fluid normals required")
    tangent = (segments[:, 1] - segments[:, 0]) / lengths[:, None]
    if np.max(abs(np.einsum("ij,ij->i", tangent, normals))) > 1e-12:
        _fail("fluid normal is not perpendicular to segment")
    tol = 1e-10 * scale
    if len(segments) * len(triangles) > 2000000 or len(segments) ** 2 > 2000000:
        _fail("interface candidate-pair resource budget exhausted")
    for i, segment in enumerate(segments):
        for other in segments[i + 1 :]:
            if np.max(abs((other - segment[0]) @ normals[i])) <= tol:
                interval = (other - segment[0]) @ tangent[i]
                overlap = min(lengths[i], float(interval.max())) - max(
                    0.0, float(interval.min())
                )
                if overlap > tol:
                    _fail(f"double fluid interface coverage at segment {i}")
    covered = np.zeros(len(triangles))
    triangle_area = np.array(
        [topology.triangle_faces[int(row)].area for row in selected]
    )
    qpoints, qareas, qsegments, weights, columns = [], [], [], [], []
    maximum_error = 0.0
    polygon_pairs = 0
    gauss = np.array(
        [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]]
    )
    for face, (segment, normal, direction, length) in enumerate(
        zip(segments, normals, tangent, lengths, strict=True)
    ):
        rectangle = np.array([[0, 0], [length, 0], [length, width], [0, width]])
        polygons: list[FloatArray] = []
        face_area = 0.0
        for j, triangle in enumerate(xyz):
            distance = (triangle[:, :2] - segment[0]) @ normal
            if np.max(abs(distance)) > tol:
                continue
            local = np.column_stack(
                [(triangle[:, :2] - segment[0]) @ direction, triangle[:, 2] - origin]
            )
            polygon = _clip(local if _area(local) > 0 else local[::-1], rectangle)
            area = _area(polygon)
            if area <= 1e-14 * min(triangle_area[j], length * width):
                continue
            ns = topology.triangle_faces[int(selected[j])].unit_normal
            if np.linalg.norm(ns + np.r_[normal, 0]) > 1e-10:
                _fail(
                    f"wrong interface orientation at fluid segment {face}, structural facet {selected[j]}"
                )
            polygon_pairs += len(polygons)
            if polygon_pairs > 2000000:
                _fail("interface overlap-polygon resource budget exhausted")
            for previous in polygons:
                if abs(_area(_clip(polygon, previous))) > 1e-12 * min(
                    area, abs(_area(previous))
                ):
                    _fail(
                        f"double structural interface coverage at fluid segment {face}"
                    )
            polygons.append(polygon)
            face_area += area
            covered[j] += area
            inverse = np.linalg.inv(np.vstack([local.T, np.ones(3)]))
            for k in range(1, len(polygon) - 1):
                piece = polygon[[0, k, k + 1]]
                piece_area = abs(_area(piece))
                if piece_area == 0:
                    continue
                for point in gauss @ piece:
                    bary = inverse @ np.r_[point, 1]
                    if np.min(bary) < -1e-10:
                        _fail(
                            "overlap quadrature extrapolates outside structural triangle"
                        )
                    qpoints.append(
                        np.r_[segment[0] + point[0] * direction, origin + point[1]]
                    )
                    qareas.append(piece_area / 3)
                    qsegments.append(face)
                    weights.extend(bary)
                    columns.extend(triangles[j])
                    if len(qpoints) > max_quadrature_points:
                        _fail(
                            f"interface quadrature budget {max_quadrature_points} exhausted"
                        )
        error = abs(face_area - length * width) / (length * width)
        maximum_error = max(maximum_error, error)
        if error > 1e-10:
            _fail(
                f"incomplete or double interface coverage at fluid segment {face}: relative area error {error}"
            )
    errors = abs(covered - triangle_area) / triangle_area
    if np.max(errors) > 1e-10:
        _fail(
            f"incomplete or double fluid coverage of structural facet {selected[int(np.argmax(errors))]}"
        )
    maximum_error = max(maximum_error, float(np.max(errors)))
    points = np.asarray(qpoints)
    matrix = sparse.csr_matrix(
        (weights, (np.repeat(np.arange(len(points)), 3), columns)),
        shape=(len(points), len(structure.points)),
    )
    if (
        np.max(abs(matrix @ structure.points - points)) > tol
        or np.max(abs(matrix @ np.ones(len(structure.points)) - 1)) > 1e-12
    ):
        _fail("interface does not reproduce coordinates/constants")
    index = np.asarray(qsegments, dtype=np.int64)
    frozen_index = np.frombuffer(index.tobytes(), dtype=np.int64)
    frozen_tri = np.frombuffer(triangles.tobytes(), dtype=np.int64).reshape(
        triangles.shape
    )
    return InterfaceMap(
        _freeze(points),
        _freeze(np.asarray(qareas)),
        frozen_index,
        _freeze(structure.points),
        _freeze(segments),
        _freeze(normals),
        width,
        origin,
        scale,
        maximum_error,
        matrix,
        frozen_tri,
    )


def build_domain_interface(
    domain: FluidDomain,
    *,
    wall_name: str,
    structure: Mesh,
    facet_set: str,
    length_scale_m: float,
    span_origin_m: float = 0.0,
    max_quadrature_points: int = 200000,
) -> InterfaceMap:
    """Build an overlap map from a named explicitly FSI-marked no-slip wall.

    Segment order equals domain.mesh.boundary_sets[wall_name]. Tractions must
    be selected from the solver's wall-row order into that same segment order.
    Width is taken from FluidDomain, never guessed from structural coordinates.
    """
    bc = domain.boundaries.get(wall_name)
    if bc is None or bc["type"] != "wall" or not bc.get("fsi", False):
        _fail("interface requires a named explicitly FSI-marked wall")
    rows = domain.mesh.boundary_sets[wall_name]
    return build_interface_map(
        segments_m=domain.mesh.points[domain.mesh.line2[rows]],
        fluid_normals=domain.mesh.boundary_normals[rows],
        width_m=domain.width_m,
        structure=structure,
        facet_set=facet_set,
        length_scale_m=length_scale_m,
        span_origin_m=span_origin_m,
        max_quadrature_points=max_quadrature_points,
    )
