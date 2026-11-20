"""Sparse MAC operators with conservative convection and pressure outlet."""

import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import bmat, coo_matrix, diags  # type: ignore[import-untyped]
from scipy.sparse.linalg import (  # type: ignore[import-untyped]
    MatrixRankWarning,
    spsolve,
)

from fem_toolbox.cfd.flow_model import CartesianFluidGrid, ChannelBoundary, ClosedWalls
from fem_toolbox.errors import ValidationError


@dataclass
class MACSystem:
    """Private assembly state: free velocity maps, boundary values and volumes."""

    grid: CartesianFluidGrid
    boundary: ChannelBoundary | ClosedWalls
    mu: float
    uid: NDArray[np.int64]
    vid: NDArray[np.int64]
    u_boundary: NDArray[np.float64]
    weights: NDArray[np.float64]
    diffusion: Any
    divergence: Any
    gradient: Any
    boundary_rhs: NDArray[np.float64]
    mass_rhs: NDArray[np.float64]
    edges: list[tuple[int, int, float]]
    contacts: list[tuple[int, float, float]]

    @property
    def spacing(self) -> tuple[float, float]:
        nx, ny = self.grid.cells
        return self.grid.size_m[0] / nx, self.grid.size_m[1] / ny

    @property
    def nu(self) -> int:
        return int(np.count_nonzero(self.uid >= 0))

    def fields(
        self, q: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        u = self.u_boundary.copy()
        v = np.zeros(self.vid.shape)
        u[self.uid >= 0] = q[self.uid[self.uid >= 0]]
        v[self.vid >= 0] = q[self.vid[self.vid >= 0]]
        return u, v


def assemble(
    grid: CartesianFluidGrid, boundary: ChannelBoundary | ClosedWalls, mu: float
) -> MACSystem:
    nx, ny = grid.cells
    hx, hy = grid.size_m[0] / nx, grid.size_m[1] / ny
    opened = isinstance(boundary, ChannelBoundary)
    uid = np.full((ny, nx + 1), -1, dtype=np.int64)
    vid = np.full((ny + 1, nx), -1, dtype=np.int64)
    us = (ny, nx if opened else nx - 1)
    nu = us[0] * us[1]
    nv = (ny - 1) * nx
    uid[:, 1 : (nx + 1 if opened else nx)] = np.arange(nu).reshape(us)
    vid[1:-1, :] = np.arange(nu, nu + nv).reshape(ny - 1, nx)
    u0 = np.zeros(uid.shape)
    if opened:
        assert isinstance(boundary, ChannelBoundary)
        if len(boundary.inlet_velocity_m_s) != ny:
            raise ValidationError("inlet profile must match ny", field="inlet")
        u0[:, 0] = boundary.inlet_velocity_m_s
        walls = (0.0, 0.0, 0.0, 0.0)
    else:
        assert isinstance(boundary, ClosedWalls)
        walls = boundary.tangential_velocities_m_s
    weights = np.full(nu + nv, hx * hy)
    if opened:
        weights[uid[:, -1]] = hx * hy / 2
    rhs = np.zeros(nu + nv)
    rows = []
    cols = []
    data = []
    edges = []
    contacts = []
    for ids, normal_axis in [(uid, 1), (vid, 0)]:
        for j, i in zip(*np.nonzero(ids >= 0), strict=True):
            k = int(ids[j, i])
            diagonal = 0.0
            for axis, h in [(0, hy), (1, hx)]:
                coefficient = mu / h / h
                if ids is uid and opened and i == nx and axis == 1:
                    coefficient *= 2
                for direction in (-1, 1):
                    nj, ni = (
                        j + (direction if axis == 0 else 0),
                        i + (direction if axis == 1 else 0),
                    )
                    inside = 0 <= nj < ids.shape[0] and 0 <= ni < ids.shape[1]
                    neighbor = int(ids[nj, ni]) if inside else -1
                    if neighbor >= 0:
                        diagonal += coefficient
                        rows.append(k)
                        cols.append(neighbor)
                        data.append(-coefficient)
                        if neighbor > k:
                            edges.append((k, neighbor, coefficient * weights[k]))
                    elif axis == normal_axis:
                        if ids is uid and opened and i == nx and direction == 1:
                            continue
                        value = float(u0[nj, ni]) if ids is uid else 0.0
                        diagonal += coefficient
                        rhs[k] += coefficient * value
                        contacts.append((k, value, coefficient * weights[k]))
                    else:
                        if ids is vid and opened and axis == 1 and direction == 1:
                            continue
                        value = (
                            walls[0 if direction == -1 else 1]
                            if ids is uid
                            else walls[2 if direction == -1 else 3]
                        )
                        diagonal += 2 * coefficient
                        rhs[k] += 2 * coefficient * value
                        contacts.append((k, value, 2 * coefficient * weights[k]))
            rows.append(k)
            cols.append(k)
            data.append(diagonal)
    a = coo_matrix((data, (rows, cols)), shape=(nu + nv, nu + nv)).tocsc()
    dr = []
    dc = []
    dv = []
    mass = np.zeros(nx * ny)
    for j in range(ny):
        for i in range(nx):
            cell = j * nx + i
            for ids, nj, ni, coefficient in [
                (uid, j, i, -1 / hx),
                (uid, j, i + 1, 1 / hx),
                (vid, j, i, -1 / hy),
                (vid, j + 1, i, 1 / hy),
            ]:
                k = int(ids[nj, ni])
                if k >= 0:
                    dr.append(cell)
                    dc.append(k)
                    dv.append(coefficient)
                elif ids is uid:
                    mass[cell] -= coefficient * u0[nj, ni]
    d = coo_matrix((dv, (dr, dc)), shape=(nx * ny, nu + nv)).tocsc()
    g = -diags(hx * hy / weights) @ d.T
    if not np.all(np.isfinite(a.data)) or not np.all(np.isfinite(g.data)):
        raise ValidationError("MAC coefficients exceed float64 range", field="grid")
    return MACSystem(
        grid, boundary, mu, uid, vid, u0, weights, a, d, g, rhs, mass, edges, contacts
    )


def _limited_correction(
    ids: NDArray[np.int64],
    q: NDArray[np.float64],
    j: int,
    i: int,
    axis: int,
    flux: float,
    boundary_values: NDArray[np.float64],
) -> float:
    """Minmod linear donor reconstruction; no extrapolated new extrema."""
    center = q[ids[j, i]]
    samples = []
    for direction in (-1, 1):
        nj = j + (direction if axis == 0 else 0)
        ni = i + (direction if axis == 1 else 0)
        if 0 <= nj < ids.shape[0] and 0 <= ni < ids.shape[1]:
            index = int(ids[nj, ni])
            samples.append(q[index] if index >= 0 else boundary_values[nj, ni])
        else:
            samples.append(center)
    left, right = center - samples[0], samples[1] - center
    if left * right <= 0:
        return 0.0
    slope = np.copysign(min(abs(left), abs(right)), left)
    return float(np.copysign(0.5, flux) * slope)


def convection(
    system: MACSystem, q: NDArray[np.float64], rho: float, scheme: str = "donor-cell"
) -> tuple[Any, NDArray[np.float64]]:
    """Freeze shared dual-face mass fluxes and upwind donor values."""
    u, v = system.fields(q)
    hx, hy = system.spacing
    nx, _ny = system.grid.cells
    opened = isinstance(system.boundary, ChannelBoundary)
    rows = []
    cols = []
    data = []
    constant = np.zeros(len(q))
    v_boundary = np.zeros(system.vid.shape)
    for ids, normal_axis in [(system.uid, 1), (system.vid, 0)]:
        for j, i in zip(*np.nonzero(ids >= 0), strict=True):
            k = int(ids[j, i])
            dx = hx / 2 if ids is system.uid and opened and i == nx else hx
            if ids is system.uid:
                west = (u[j, i - 1] + u[j, i]) / 2
                east = u[j, i] if i == nx else (u[j, i] + u[j, i + 1]) / 2
                south = v[j, i - 1] if i == nx else (v[j, i - 1] + v[j, i]) / 2
                north = (
                    v[j + 1, i - 1] if i == nx else (v[j + 1, i - 1] + v[j + 1, i]) / 2
                )
            else:
                west = (u[j - 1, i] + u[j, i]) / 2
                east = (u[j - 1, i + 1] + u[j, i + 1]) / 2
                south = (v[j - 1, i] + v[j, i]) / 2
                north = (v[j, i] + v[j + 1, i]) / 2
            for axis, direction, flux, h in [
                (1, -1, west, dx),
                (1, 1, east, dx),
                (0, -1, south, hy),
                (0, 1, north, hy),
            ]:
                if flux == 0:
                    continue
                nj, ni = (
                    j + (direction if axis == 0 else 0),
                    i + (direction if axis == 1 else 0),
                )
                inside = 0 <= nj < ids.shape[0] and 0 <= ni < ids.shape[1]
                neighbor = int(ids[nj, ni]) if inside else -1
                if scheme == "central":
                    if rho * abs(flux) * h / system.mu >= 2:
                        raise ValidationError(
                            "central convection requires face Peclet <2; refine grid or choose donor-cell",
                            field="convection",
                        )
                    coefficient = rho * direction * flux / h
                    if neighbor >= 0:
                        for donor in (k, neighbor):
                            rows.append(k)
                            cols.append(donor)
                            data.append(coefficient / 2)
                    elif opened and axis == 1 and direction == 1:
                        rows.append(k)
                        cols.append(k)
                        data.append(coefficient)
                    elif inside and axis == normal_axis:
                        rows.append(k)
                        cols.append(k)
                        data.append(coefficient / 2)
                        value = (
                            float(system.u_boundary[nj, ni])
                            if ids is system.uid
                            else 0.0
                        )
                        constant[k] += coefficient * value / 2
                    # Tangential inflow v at the physical inlet is prescribed 0.
                    continue
                donor = k if flux * direction >= 0 else neighbor
                coefficient = rho * direction * flux / h
                # At the pressure outlet, extrapolate transported velocity.
                if donor < 0 and opened and axis == 1 and direction == 1:
                    donor = k
                if donor >= 0:
                    rows.append(k)
                    cols.append(donor)
                    data.append(coefficient)
                    if scheme == "linear-upwind" and neighbor >= 0:
                        dj, di = (j, i) if donor == k else (nj, ni)
                        boundary_values = (
                            system.u_boundary if ids is system.uid else v_boundary
                        )
                        constant[k] += coefficient * _limited_correction(
                            ids, q, dj, di, axis, float(flux), boundary_values
                        )
                else:
                    value = (
                        float(system.u_boundary[nj, ni])
                        if ids is system.uid and inside and axis == normal_axis
                        else 0.0
                    )
                    constant[k] += coefficient * value
    matrix = coo_matrix((data, (rows, cols)), shape=(len(q), len(q))).tocsc()
    if not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(constant)):
        raise ValidationError(
            "convective flux exceeds float64 range", field="convection"
        )
    return matrix, constant


def linear_step(
    system: MACSystem, operator: Any, rhs: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    opened = isinstance(system.boundary, ChannelBoundary)
    d = system.divergence if opened else system.divergence[1:, :]
    g = system.gradient if opened else system.gradient[:, 1:]
    mass = system.mass_rhs if opened else system.mass_rhs[1:]
    matrix = bmat([[operator, g], [d, None]], format="csc")
    b = np.concatenate((rhs, mass))
    if not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(b)):
        raise ValidationError("nonfinite Oseen coefficients", field="solution")
    row = np.asarray(abs(matrix).max(axis=1).toarray()).ravel()
    scaled = diags(1 / row) @ matrix
    col = np.asarray(abs(scaled).max(axis=0).toarray()).ravel()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            z = (
                np.asarray(
                    spsolve((scaled @ diags(1 / col)).tocsc(), b / row),
                    dtype=np.float64,
                )
                / col
            )
    except (RuntimeError, ValueError, MatrixRankWarning) as error:
        raise ValidationError(
            "Oseen pressure system singular or unresolved", field="solution"
        ) from error
    if not np.all(np.isfinite(z)):
        raise ValidationError("Oseen solution exceeds float64 range", field="solution")
    q = z[: len(system.weights)]
    p = z[len(system.weights) :]
    return q, p if opened else np.concatenate(([0.0], p))
