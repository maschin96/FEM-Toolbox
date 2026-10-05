"""Conservative SI Quad-FV primitives; no flow solve or boundary orchestration."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import NoReturn

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import coo_matrix, csr_matrix  # type: ignore[import-untyped]

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import (
    FluidMesh2D,
    _float_array,
    _freeze_float,
    _freeze_int,
    _indices,
)
from fem_toolbox.errors import ValidationError


def _error(field: str, reason: str) -> NoReturn:
    raise ValidationError(f"{field}: {reason}", field=field)


def _validated_mesh(value: object) -> FluidMesh2D:
    if not isinstance(value, FluidMesh2D):
        _error("mesh", "validated FluidMesh2D required")
    return value


def _validated_domain(value: object) -> FluidDomain:
    if not isinstance(value, FluidDomain):
        _error("domain", "validated FluidDomain required")
    return value


def _validated_pairs(value: object) -> Mapping[int, int]:
    if not isinstance(value, Mapping):
        _error("periodic_pairs", "reciprocal boundary-row mapping required")
    result: dict[int, int] = {}
    for a, b in value.items():
        if type(a) is not int or type(b) is not int:
            _error("periodic_pairs", "integer boundary rows required")
        result[a] = b
    return result


def _finite(value: NDArray[np.float64], field: str) -> NDArray[np.float64]:
    if not np.all(np.isfinite(value)):
        _error(field, "operation exceeds float64 range")
    return _freeze_float(value)


@dataclass(frozen=True, slots=True)
class FaceGradient:
    """Full corrected gradient and its compact/nonorthogonal normal parts.

    Scalar: gradient[nface,2], normal parts[nface]. Vector: gradient[nface,k,2],
    normal parts[nface,k]. Gradient last axis is spatial x/y.
    """

    gradient: NDArray[np.float64]
    compact_normal: NDArray[np.float64]
    nonorthogonal_normal: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class StressFlux:
    """sigma=-pI+mu*(grad U+grad U^T); traction is outward from fluid.

    stress_Pa[nface,2,2], traction_Pa[nface,2], face_force_N[nface,2],
    cell_force_N[ncell,2]=integrated div(sigma). No body/convection terms.
    Force on a structure is the negative fluid traction times face area.
    """

    stress_Pa: NDArray[np.float64]
    traction_Pa: NDArray[np.float64]
    face_force_N: NDArray[np.float64]
    cell_force_N: NDArray[np.float64]


@dataclass(frozen=True, slots=True, init=False)
class QuadFVGeometry:
    """Immutable operators on a validated mesh, width in metres.

    All face arrays follow mesh.faces; boundary arrays follow mesh.line2.
    Periodic values must already represent the same periodic field (pressure
    means its periodic part); affine pressure drops belong in the body force.
    One periodic face is active, its partner stores the opposite flux.
    No PDE solver, pressure gauge elimination or physical BC defaults.
    """

    mesh: FluidMesh2D
    width_m: float
    cell_volumes_m3: NDArray[np.float64]
    face_areas_m2: NDArray[np.float64]
    neighbours: NDArray[np.int64]
    active_faces: NDArray[np.int64]
    periodic_partners: Mapping[int, int]
    _direction: NDArray[np.float64]
    _distance: NDArray[np.float64]
    _fraction: NDArray[np.float64]
    _boundary_row: NDArray[np.int64]
    _gx: csr_matrix
    _gy: csr_matrix

    def __init__(
        self,
        mesh: FluidMesh2D,
        *,
        width_m: float,
        periodic_pairs: Mapping[int, int] | None = None,
        max_cells: int = 50000,
    ) -> None:
        mesh = _validated_mesh(mesh)
        n = len(mesh.quad4)
        if type(max_cells) is not int or max_cells < 1 or n > max_cells:
            _error("max_cells", "positive integer resource limit exceeded")
        width = _scalar(width_m, "width_m")
        if width <= 0:
            _error("width_m", "positive SI width required")
        volumes = mesh.cell_areas_m2 * width
        areas = mesh.face_lengths_m * width
        if (
            not np.all(np.isfinite(volumes))
            or not np.all(np.isfinite(areas))
            or min(volumes.min(), areas.min()) <= 0
        ):
            _error("width_m", "unrepresentable SI volumes/areas")
        owner = mesh.face_owners
        neighbour = mesh.face_neighbours.copy()
        direction = mesh.face_centres_m - mesh.cell_centres_m[owner]
        inner = neighbour >= 0
        direction[inner] = (
            mesh.cell_centres_m[neighbour[inner]] - mesh.cell_centres_m[owner[inner]]
        )
        pairs: dict[int, int] = {}
        boundary_row = np.full(len(owner), -1, dtype=np.int64)
        boundary_row[mesh.boundary_face_indices] = np.arange(len(mesh.line2))
        active = np.ones(len(owner), dtype=bool)
        raw = _validated_pairs(periodic_pairs if periodic_pairs is not None else {})
        tolerance = 1e-10 * float(np.max(np.ptp(mesh.points, axis=0)))
        for a, b in raw.items():
            if (
                type(a) is not int
                or type(b) is not int
                or not 0 <= a < len(mesh.line2)
                or not 0 <= b < len(mesh.line2)
                or a == b
                or raw.get(b) != a
            ):
                _error(
                    "periodic_pairs", "distinct valid reciprocal boundary rows required"
                )
            face, partner = (
                int(mesh.boundary_face_indices[a]),
                int(mesh.boundary_face_indices[b]),
            )
            shift = mesh.face_centres_m[partner] - mesh.face_centres_m[face]
            x, y = mesh.points[mesh.line2[a]] + shift, mesh.points[mesh.line2[b]]
            if (
                not (
                    np.max(np.abs(x - y)) <= tolerance
                    or np.max(np.abs(x - y[::-1])) <= tolerance
                )
                or np.linalg.norm(mesh.face_normals[face] + mesh.face_normals[partner])
                > 1e-10
            ):
                _error(
                    f"periodic_pairs[{a}]",
                    "congruent translated segments with opposite normals required",
                )
            neighbour[face] = owner[partner]
            direction[face] = (
                mesh.cell_centres_m[owner[partner]]
                - shift
                - mesh.cell_centres_m[owner[face]]
            )
            pairs[face] = partner
            if face > partner:
                active[face] = False
        distance = np.einsum("ij,ij->i", direction, mesh.face_normals)
        if not np.all(np.isfinite(direction)) or np.any(distance <= 0):
            _error("faces", "positive resolvable owner-normal distances required")
        fraction = (
            np.einsum(
                "ij,ij->i",
                mesh.face_centres_m - mesh.cell_centres_m[owner],
                mesh.face_normals,
            )
            / distance
        )
        if np.any(fraction <= 0) or np.any(fraction[neighbour >= 0] >= 1):
            _error("faces", "face plane must separate neighbouring centres")
        for name, value in {
            "cell_volumes_m3": volumes,
            "face_areas_m2": areas,
            "_direction": direction,
            "_distance": distance,
            "_fraction": fraction,
        }.items():
            object.__setattr__(self, name, _finite(value, name))
        object.__setattr__(self, "mesh", mesh)
        object.__setattr__(self, "width_m", width)
        object.__setattr__(self, "neighbours", _freeze_int(neighbour))
        object.__setattr__(self, "active_faces", _freeze_int(np.flatnonzero(active)))
        object.__setattr__(self, "_boundary_row", _freeze_int(boundary_row))
        object.__setattr__(self, "periodic_partners", MappingProxyType(pairs))
        self._assemble_gradient()

    @classmethod
    def from_domain(
        cls, domain: FluidDomain, *, max_cells: int = 50000
    ) -> QuadFVGeometry:
        """Use validated domain width and boundary-row periodic partners only."""
        domain = _validated_domain(domain)
        return cls(
            domain.mesh,
            width_m=domain.width_m,
            periodic_pairs=domain.periodic_partners,
            max_cells=max_cells,
        )

    def _assemble_gradient(self) -> None:
        """Weighted least squares, w=1/|delta|²; two equations per cell."""
        cells, targets, deltas = [], [], []
        n = len(self.mesh.quad4)
        for face in self.active_faces:
            p, q = int(self.mesh.face_owners[face]), int(self.neighbours[face])
            cells.append(p)
            targets.append(q if q >= 0 else n + int(self._boundary_row[face]))
            deltas.append(self._direction[face])
            if q >= 0:
                cells.append(q)
                targets.append(p)
                deltas.append(-self._direction[face])
        cell = np.array(cells, dtype=np.int64)
        delta = np.asarray(deltas)
        squared = np.einsum("ij,ij->i", delta, delta)
        if np.any(squared <= 0) or not np.all(np.isfinite(squared)):
            _error("gradient", "unresolved squared centre distances")
        weight = 1 / squared
        normal = np.zeros((n, 2, 2))
        np.add.at(
            normal, cell, delta[:, :, None] * delta[:, None, :] * weight[:, None, None]
        )
        # Anisotropic cells are supported, but unrecoverable float64 fits are not.
        if not np.all(np.isfinite(normal)) or np.any(np.linalg.cond(normal) > 1e12):
            _error("gradient", "rank deficient/unresolved least-squares stencil")
        coefficients = (
            np.einsum("ijk,ik->ij", np.linalg.inv(normal)[cell], delta)
            * weight[:, None]
        )
        if not np.all(np.isfinite(coefficients)):
            _error("gradient", "unrepresentable least-squares coefficients")
        for axis, name in enumerate(("_gx", "_gy")):
            data = np.concatenate((coefficients[:, axis], -coefficients[:, axis]))
            matrix = coo_matrix(
                (data, (np.tile(cell, 2), np.concatenate((targets, cell)))),
                shape=(n, n + len(self.mesh.line2)),
            ).tocsr()
            object.__setattr__(self, name, matrix)

    def _field(self, value: ArrayLike, count: int, field: str) -> NDArray[np.float64]:
        a = _float_array(value, field)
        if (
            a.ndim not in (1, 2)
            or a.shape[0] != count
            or (a.ndim == 2 and not a.shape[1])
        ):
            _error(field, f"expected ({count},) or ({count},k)")
        return a

    def _fields(
        self, values: ArrayLike, boundary_values: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        a = self._field(values, len(self.mesh.quad4), "cell_values")
        b = self._field(boundary_values, len(self.mesh.line2), "boundary_values")
        if a.shape[1:] != b.shape[1:]:
            _error("boundary_values", "components must match cell values")
        return a, b

    def gradient(
        self, values: ArrayLike, boundary_values: ArrayLike
    ) -> NDArray[np.float64]:
        """Linear-exact cell gradient from centre and prescribed face samples.

        Scalar input (nc,)/(nb,) -> (nc,2); k components -> (nc,k,2).
        All boundary rows are supplied; periodic rows are validated but ignored.
        """
        a, b = self._fields(values, boundary_values)
        combined = np.concatenate((a, b))
        return _finite(
            np.stack((self._gx @ combined, self._gy @ combined), axis=-1), "gradient"
        )

    def _interpolate(
        self,
        a: NDArray[np.float64],
        b: NDArray[np.float64],
        gradient: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        owner, neighbour = self.mesh.face_owners, self.neighbours
        inner = neighbour >= 0
        offset = self.mesh.face_centres_m - self.mesh.cell_centres_m[owner]
        result: NDArray[np.float64] = a[owner] + np.einsum(
            "i...j,ij->i...", gradient[owner], offset
        )
        displaced = offset[inner] - self._direction[inner]
        other = a[neighbour[inner]] + np.einsum(
            "i...j,ij->i...", gradient[neighbour[inner]], displaced
        )
        weight = self._fraction[inner].reshape((-1,) + (1,) * (a.ndim - 1))
        result[inner] = (1 - weight) * result[inner] + weight * other
        result[~inner] = b[self._boundary_row[~inner]]
        for face, partner in self.periodic_partners.items():
            if face < partner:
                result[partner] = result[face]
        return result

    def face_values(
        self, values: ArrayLike, boundary_values: ArrayLike
    ) -> NDArray[np.float64]:
        """Interpolate centre-to-face extrapolations; exact for affine fields."""
        a, b = self._fields(values, boundary_values)
        return _finite(self._interpolate(a, b, self.gradient(a, b)), "face_values")

    def _average_gradient(self, gradient: NDArray[np.float64]) -> NDArray[np.float64]:
        result = gradient[self.mesh.face_owners].copy()
        inner = self.neighbours >= 0
        weight = self._fraction[inner].reshape((-1,) + (1,) * (gradient.ndim - 1))
        result[inner] = (1 - weight) * result[inner] + weight * gradient[
            self.neighbours[inner]
        ]
        return result

    def face_gradient(
        self, values: ArrayLike, boundary_values: ArrayLike
    ) -> FaceGradient:
        """Over-relaxed directional correction with the same geometry as pressure.

        g_face = g_bar + n*(delta_phi - g_bar dot delta_x)/(delta_x dot n).
        Normal part is delta_phi/d_n - g_bar dot(delta_x-d_n*n)/d_n.
        """
        a, b = self._fields(values, boundary_values)
        owner, neighbour = self.mesh.face_owners, self.neighbours
        inner = neighbour >= 0
        difference = np.empty_like(a[owner])
        difference[inner] = a[neighbour[inner]] - a[owner[inner]]
        difference[~inner] = b[self._boundary_row[~inner]] - a[owner[~inner]]
        average = self._average_gradient(self.gradient(a, b))
        dn = self._distance.reshape((-1,) + (1,) * (a.ndim - 1))
        normal = self.mesh.face_normals
        compact = difference / dn
        tangent = self._direction - self._distance[:, None] * normal
        correction = -np.einsum("i...j,ij->i...", average, tangent) / dn
        current_normal = np.einsum("i...j,ij->i...", average, normal)
        full = average + (compact + correction - current_normal)[
            ..., None
        ] * normal.reshape((len(normal),) + (1,) * (a.ndim - 1) + (2,))
        for face, partner in self.periodic_partners.items():
            if face < partner:
                full[partner] = full[face]
                compact[partner], correction[partner] = (
                    -compact[face],
                    -correction[face],
                )
        return FaceGradient(
            _finite(full, "face_gradient"),
            _finite(compact, "compact_normal"),
            _finite(correction, "nonorthogonal_normal"),
        )

    def cell_balance(self, face_flux: ArrayLike) -> NDArray[np.float64]:
        """Integrated outward balance, scalar or k-vector flux per mesh face.

        Periodic alias must carry the exact negative master flux. No division
        by volume: mass flux kg/s -> kg/s; face force N -> cell force N.
        """
        flux = self._field(face_flux, len(self.mesh.faces), "face_flux")
        for face, partner in self.periodic_partners.items():
            if face < partner and not np.array_equal(flux[partner], -flux[face]):
                _error(
                    f"face_flux[{partner}]",
                    "periodic partner must be opposite master flux",
                )
        result = np.zeros((len(self.mesh.quad4), *flux.shape[1:]))
        active = self.active_faces
        np.add.at(result, self.mesh.face_owners[active], flux[active])
        inner = active[self.neighbours[active] >= 0]
        np.add.at(result, self.neighbours[inner], -flux[inner])
        return _finite(result, "cell_balance")

    def _positive(
        self, value: ArrayLike, count: int, field: str
    ) -> NDArray[np.float64]:
        a = _float_array(value, field)
        if a.ndim == 0:
            a = np.full(count, float(a))
        if a.shape != (count,) or np.any(a <= 0):
            _error(field, f"positive scalar or ({count},) SI values required")
        return a

    def stress_flux(
        self,
        velocity_m_s: ArrayLike,
        boundary_velocity_m_s: ArrayLike,
        *,
        viscosity_Pa_s: ArrayLike,
        boundary_viscosity_Pa_s: ArrayLike,
        pressure_Pa: ArrayLike | None = None,
        boundary_pressure_Pa: ArrayLike | None = None,
    ) -> StressFlux:
        """Full symmetric Newton stress with linear-exact face mu and pressure.

        U:(nc,2), Ub:(nb,2). mu: positive scalar/(nc,), mub:(nb,).
        Both pressure arrays are optional together, default zero viscous-only.
        No scalar Laplacian replacement; cross/shear terms retained.
        """
        a, b = self._fields(velocity_m_s, boundary_velocity_m_s)
        if a.shape != (len(self.mesh.quad4), 2):
            _error("velocity_m_s", "two velocity components required")
        mu = self._positive(viscosity_Pa_s, len(a), "viscosity_Pa_s")
        mub = self._positive(boundary_viscosity_Pa_s, len(b), "boundary_viscosity_Pa_s")
        mf = self.face_values(mu, mub)
        if np.any(mf <= 0):
            _error("viscosity_Pa_s", "face reconstruction yields nonpositive viscosity")
        g = self.face_gradient(a, b).gradient
        sigma = mf[:, None, None] * (g + np.swapaxes(g, 1, 2))
        if pressure_Pa is not None and boundary_pressure_Pa is not None:
            p, pb = self._fields(pressure_Pa, boundary_pressure_Pa)
            if p.ndim != 1:
                _error("pressure_Pa", "scalar pressure per cell required")
            sigma -= self.face_values(p, pb)[:, None, None] * np.eye(2)
        elif pressure_Pa is not None or boundary_pressure_Pa is not None:
            _error(
                "pressure_Pa", "cell and boundary pressure must be provided together"
            )
        traction = np.einsum("ijk,ik->ij", sigma, self.mesh.face_normals)
        force = traction * self.face_areas_m2[:, None]
        for face, partner in self.periodic_partners.items():
            if face < partner:
                force[partner], traction[partner], sigma[partner] = (
                    -force[face],
                    -traction[face],
                    sigma[face],
                )
        return StressFlux(
            _finite(sigma, "stress_Pa"),
            _finite(traction, "traction_Pa"),
            _finite(force, "face_force_N"),
            self.cell_balance(force),
        )

    def _pressure_data(
        self, diagonal: ArrayLike, density: float, boundary_indices: ArrayLike
    ) -> tuple[NDArray[np.float64], NDArray[np.bool_], float]:
        a = self._positive(diagonal, len(self.mesh.quad4), "momentum_diagonal_kg_s")
        rho = _scalar(density, "density_kg_m3")
        if rho <= 0:
            _error("density_kg_m3", "positive density required")
        try:
            raw_indices = np.asarray(boundary_indices)
        except (ValueError, TypeError) as error:
            raise ValidationError(
                "invalid pressure boundary indices", field="pressure_boundary_indices"
            ) from error
        if raw_indices.ndim != 1:
            _error("pressure_boundary_indices", "one-dimensional index array required")
        ids = (
            _indices(
                boundary_indices,
                None,
                len(self.mesh.line2),
                "pressure_boundary_indices",
            )
            if raw_indices.size
            else np.empty(0, dtype=np.int64)
        )
        if len(np.unique(ids)) != len(ids):
            _error("pressure_boundary_indices", "duplicate boundary indices")
        opened = np.zeros(len(self.mesh.faces), dtype=bool)
        opened[self.mesh.boundary_face_indices[ids]] = True
        if np.any(opened & (self.neighbours >= 0)):
            _error(
                "pressure_boundary_indices",
                "periodic face cannot also prescribe pressure",
            )
        d = self.cell_volumes_m3 / a
        df = d[self.mesh.face_owners].copy()
        inner = self.neighbours >= 0
        w = self._fraction[inner]
        df[inner] = (1 - w) * df[inner] + w * d[self.neighbours[inner]]
        if not np.all(np.isfinite(df)) or np.any(df <= 0):
            _error("momentum_diagonal_kg_s", "unrepresentable V/a pressure mobility")
        return df, opened, rho

    def pressure_correction_matrix(
        self,
        *,
        momentum_diagonal_kg_s: ArrayLike,
        density_kg_m3: float,
        pressure_boundary_indices: ArrayLike = (),
    ) -> csr_matrix:
        """Positive compact -div(rho*d*grad(pc)) block (kg/(Pa s)).

        d=V/a from the *unrelaxed integrated* momentum diagonal (kg/s).
        No gauge: closed/periodic block has exactly the constant null mode.
        Selected boundary line rows mean prescribed pc=0; other outer fluxes
        fixed. Nonorthogonal correction is explicit through FaceGradient.
        Returned sparse matrix is a detached mutable object.
        """
        df, opened, rho = self._pressure_data(
            momentum_diagonal_kg_s, density_kg_m3, pressure_boundary_indices
        )
        coefficient = rho * self.face_areas_m2 * df / self._distance
        if not np.all(np.isfinite(coefficient)) or np.any(coefficient <= 0):
            _error("pressure_correction", "unrepresentable compact coefficients")
        rows: list[int] = []
        cols: list[int] = []
        data: list[float] = []
        for face in self.active_faces:
            p, q, c = (
                int(self.mesh.face_owners[face]),
                int(self.neighbours[face]),
                float(coefficient[face]),
            )
            if q >= 0:
                rows.extend((p, p, q, q))
                cols.extend((p, q, p, q))
                data.extend((c, -c, -c, c))
            elif opened[face]:
                rows.append(p)
                cols.append(p)
                data.append(c)
        n = len(self.mesh.quad4)
        return coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()

    def rhie_chow_flux(
        self,
        velocity_m_s: ArrayLike,
        boundary_velocity_m_s: ArrayLike,
        pressure_Pa: ArrayLike,
        boundary_pressure_Pa: ArrayLike,
        *,
        momentum_diagonal_kg_s: ArrayLike,
        density_kg_m3: float,
        pressure_boundary_indices: ArrayLike = (),
    ) -> NDArray[np.float64]:
        """Common face mass flux rho*A*(Uf + df*(gbar_p-gface_p)) dot n.

        U (nc,2), p (nc,), explicit face boundary samples. Fixed outer fluxes
        use prescribed U; only selected pressure boundaries receive RC coupling.
        No relaxation applied here. Caller must supply consistent unrelaxed a.
        """
        u, ub = self._fields(velocity_m_s, boundary_velocity_m_s)
        p, pb = self._fields(pressure_Pa, boundary_pressure_Pa)
        if u.shape != (len(self.mesh.quad4), 2) or p.ndim != 1:
            _error("rhie_chow", "two velocity components and scalar pressure required")
        df, opened, rho = self._pressure_data(
            momentum_diagonal_kg_s, density_kg_m3, pressure_boundary_indices
        )
        average = self._average_gradient(self.gradient(p, pb))
        normal_p = np.einsum(
            "ij,ij->i", self.face_gradient(p, pb).gradient, self.mesh.face_normals
        )
        normal_average = np.einsum("ij,ij->i", average, self.mesh.face_normals)
        speed = np.einsum("ij,ij->i", self.face_values(u, ub), self.mesh.face_normals)
        coupled = (self.neighbours >= 0) | opened
        speed[coupled] += df[coupled] * (normal_average[coupled] - normal_p[coupled])
        flux = rho * self.face_areas_m2 * speed
        for face, partner in self.periodic_partners.items():
            if face < partner:
                flux[partner] = -flux[face]
        return _finite(flux, "mass_flux_kg_s")
