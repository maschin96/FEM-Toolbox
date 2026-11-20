"""Monolithic pressure/velocity coupling on a uniform closed 2D MAC grid."""

import warnings
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import bmat, coo_matrix, diags  # type: ignore[import-untyped]
from scipy.sparse.linalg import (  # type: ignore[import-untyped]
    MatrixRankWarning,
    spsolve,
)

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.errors import ValidationError


def _owned(a: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.frombuffer(a.tobytes(), dtype=np.float64).reshape(a.shape)


@dataclass(frozen=True, slots=True)
class StokesResult:
    """Immutable SI fields; residuals include the omitted gauge-cell equation."""

    pressure_Pa: NDArray[np.float64]
    u_m_s: NDArray[np.float64]
    v_m_s: NDArray[np.float64]
    divergence_s_inv: NDArray[np.float64]
    cell_flow_balance_m2_s: NDArray[np.float64]
    momentum_x_residual_N_m3: NDArray[np.float64]
    momentum_y_residual_N_m3: NDArray[np.float64]
    dissipation_W_m: float
    power_balance_W_m: float
    relative_residual: float
    converged: bool


def _force(value: ArrayLike, shape: tuple[int, int], field: str) -> NDArray[np.float64]:
    try:
        a = np.asarray(value, dtype=np.float64)
    except (ValueError, TypeError, OverflowError) as error:
        raise ValidationError(f"invalid {field}", field=field) from error
    if a.ndim == 0:
        a = np.full(shape, float(a))
    if a.shape != shape or not np.all(np.isfinite(a)):
        raise ValidationError(
            f"{field} needs finite values of shape {shape}", field=field
        )
    return a.copy()


def _power_balance(
    u: NDArray[np.float64],
    v: NDArray[np.float64],
    forcing: NDArray[np.float64],
    velocity: NDArray[np.float64],
    mu: float,
    hx: float,
    hy: float,
    walls: tuple[float, ...],
) -> tuple[float, float, float, float]:
    """Reconstruct viscous dissipation, body work and wall work from fields."""
    volume = hx * hy
    with np.errstate(over="ignore", invalid="ignore"):
        # Each shared internal gradient contributes once; walls use half-spacing.
        diss = (
            mu
            * volume
            * (
                np.sum(np.diff(u, axis=1) ** 2) / hx**2
                + np.sum(np.diff(v, axis=0) ** 2) / hy**2
            )
        )
        diss += (
            mu
            * volume
            * (
                np.sum(np.diff(u[:, 1:-1], axis=0) ** 2) / hy**2
                + np.sum(np.diff(v[1:-1, :], axis=1) ** 2) / hx**2
            )
        )
        diss += (
            2
            * mu
            * volume
            * (
                np.sum((u[0, 1:-1] - walls[0]) ** 2) / hy**2
                + np.sum((u[-1, 1:-1] - walls[1]) ** 2) / hy**2
                + np.sum((v[1:-1, 0] - walls[2]) ** 2) / hx**2
                + np.sum((v[1:-1, -1] - walls[3]) ** 2) / hx**2
            )
        )
        wall_work = (
            2
            * mu
            * volume
            * (
                np.sum(walls[0] * (walls[0] - u[0, 1:-1])) / hy**2
                + np.sum(walls[1] * (walls[1] - u[-1, 1:-1])) / hy**2
                + np.sum(walls[2] * (walls[2] - v[1:-1, 0])) / hx**2
                + np.sum(walls[3] * (walls[3] - v[1:-1, -1])) / hx**2
            )
        )
        body_work = float(forcing @ velocity) * volume
        balance = float(body_work + wall_work - diss)
    return float(diss), float(wall_work), body_work, balance


def solve_stokes(
    *,
    cells: tuple[int, int],
    size_m: tuple[float, float],
    viscosity_Pa_s: float,
    force_x_N_m3: ArrayLike = 0.0,
    force_y_N_m3: ArrayLike = 0.0,
    tangential_wall_velocities_m_s: tuple[float, float, float, float] = (0, 0, 0, 0),
    pressure_reference_Pa: float = 0.0,
    relative_tolerance: float = 1e-10,
    max_unknowns: int = 50000,
) -> StokesResult:
    """Solve -mu Laplace(u)+grad(p)=f and div(u)=0 on a closed rectangle.

    cells=(nx,ny), each >=2; size_m=(Lx,Ly)>0. Normal wall velocity is zero.
    Tangential wall values are bottom/top/left/right in m/s. Forces are finite
    scalars or interior-face arrays (ny,nx-1)/(ny-1,nx) in N/m³.
    Pin p[0,0] to pressure_reference_Pa; this does not impose outlet pressure.
    Return MAC fields u[ny,nx+1], v[ny+1,nx], p[ny,nx], including zero normal
    boundaries. ValidationError diagnoses inputs, singular/range failures or
    failed residual/conservation gates. No convection or turbulence model.
    """
    if (
        not isinstance(cells, (tuple, list))
        or len(cells) != 2
        or any(type(n) is not int or n < 2 for n in cells)
    ):
        raise ValidationError("cells must be two integers >=2", field="cells")
    nx, ny = cells
    nu = ny * (nx - 1)
    nv = (ny - 1) * nx
    nvel = nu + nv
    if (
        type(max_unknowns) is not int
        or max_unknowns < 1
        or nvel + nx * ny - 1 > max_unknowns
    ):
        raise ValidationError("Stokes unknown limit exceeded", field="max_unknowns")
    if not isinstance(size_m, (tuple, list)) or len(size_m) != 2:
        raise ValidationError("need two domain lengths", field="size_m")
    lx, ly = (_scalar(a, "size_m") for a in size_m)
    mu = _scalar(viscosity_Pa_s, "viscosity_Pa_s")
    tol = _scalar(relative_tolerance, "relative_tolerance")
    reference = _scalar(pressure_reference_Pa, "pressure_reference_Pa")
    if min(lx, ly, mu, tol) <= 0 or tol >= 1:
        raise ValidationError(
            "positive lengths/viscosity and 0<tolerance<1 required", field="parameters"
        )
    if (
        not isinstance(tangential_wall_velocities_m_s, (tuple, list))
        or len(tangential_wall_velocities_m_s) != 4
    ):
        raise ValidationError("need bottom/top/left/right wall values", field="walls")
    walls = tuple(_scalar(a, "walls") for a in tangential_wall_velocities_m_s)
    hx, hy = lx / nx, ly / ny
    if (
        not np.all(np.isfinite([hx * hy, hx * hx, hy * hy]))
        or min(hx * hy, hx * hx, hy * hy) == 0
    ):
        raise ValidationError("unresolved cell area or squared spacing", field="size_m")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        cx, cy = float(np.float64(mu) / hx / hx), float(np.float64(mu) / hy / hy)
    if not np.all(np.isfinite([cx, cy, 1 / hx, 1 / hy])) or min(cx, cy) == 0:
        raise ValidationError("unresolved grid/viscosity scale", field="parameters")
    fx = _force(force_x_N_m3, (ny, nx - 1), "force_x_N_m3")
    fy = _force(force_y_N_m3, (ny - 1, nx), "force_y_N_m3")
    forcing = np.concatenate((fx.ravel(), fy.ravel()))
    rhs = forcing.copy()
    rows = []
    cols = []
    values = []
    # u is face-aligned in x and cell-centered in y; v reverses these roles.
    for offset, shape, normal_axis, boundary in [
        (0, (ny, nx - 1), 1, (walls[0], walls[1])),
        (nu, (ny - 1, nx), 0, (walls[2], walls[3])),
    ]:
        for j in range(shape[0]):
            for i in range(shape[1]):
                k = offset + j * shape[1] + i
                diagonal = 0.0
                for axis, coefficient in [(0, cy), (1, cx)]:
                    for direction in [-1, 1]:
                        neighbor = [j, i]
                        neighbor[axis] += direction
                        diagonal += coefficient
                        if 0 <= neighbor[axis] < shape[axis]:
                            rows.append(k)
                            cols.append(offset + neighbor[0] * shape[1] + neighbor[1])
                            values.append(-coefficient)
                        elif axis != normal_axis:
                            diagonal += coefficient
                            rhs[k] += (
                                2 * coefficient * boundary[0 if direction == -1 else 1]
                            )
                rows.append(k)
                cols.append(k)
                values.append(diagonal)
    a = coo_matrix((values, (rows, cols)), shape=(nvel, nvel)).tocsc()
    dr = []
    dc = []
    dv = []
    for j in range(ny):
        for i in range(nx):
            cell = j * nx + i
            if i > 0:
                dr.append(cell)
                dc.append(j * (nx - 1) + i - 1)
                dv.append(-1 / hx)
            if i < nx - 1:
                dr.append(cell)
                dc.append(j * (nx - 1) + i)
                dv.append(1 / hx)
            if j > 0:
                dr.append(cell)
                dc.append(nu + (j - 1) * nx + i)
                dv.append(-1 / hy)
            if j < ny - 1:
                dr.append(cell)
                dc.append(nu + j * nx + i)
                dv.append(1 / hy)
    d = coo_matrix((dv, (dr, dc)), shape=(nx * ny, nvel)).tocsc()
    gradient = -d.T
    reduced = d[1:, :]
    matrix = bmat([[a, -reduced.T], [reduced, None]], format="csc")
    # Solve pressure relative to the gauge: constant shifts have zero gradient.
    b = np.concatenate((rhs, np.zeros(nx * ny - 1)))
    if not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(b)):
        raise ValidationError(
            "Stokes coefficients exceed float64 range", field="parameters"
        )
    row_max = np.asarray(abs(matrix).max(axis=1).toarray()).ravel()
    scaled = diags(1 / row_max) @ matrix
    col_max = np.asarray(abs(scaled).max(axis=0).toarray()).ravel()
    scaled = scaled @ diags(1 / col_max)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", MatrixRankWarning)
            solution = (
                np.asarray(spsolve(scaled.tocsc(), b / row_max), dtype=np.float64)
                / col_max
            )
    except (MatrixRankWarning, RuntimeError, ValueError) as error:
        raise ValidationError(
            "Stokes system singular or unresolved", field="solution"
        ) from error
    velocity = solution[:nvel]
    pressure_relative = np.concatenate(([0.0], solution[nvel:]))
    momentum = np.asarray(a @ velocity + gradient @ pressure_relative - rhs)
    pressure = pressure_relative + reference
    u = np.zeros((ny, nx + 1))
    u[:, 1:-1] = velocity[:nu].reshape(ny, nx - 1)
    v = np.zeros((ny + 1, nx))
    v[1:-1, :] = velocity[nu:].reshape(ny - 1, nx)
    # Reconstruct mass directly from returned face fields, including gauge cell.
    divergence = np.diff(u, axis=1) / hx + np.diff(v, axis=0) / hy
    divergence = divergence.ravel()
    volume = hx * hy
    diss, wall_work, body_work, balance = _power_balance(
        u, v, forcing, velocity, mu, hx, hy, walls
    )
    if not all(
        np.all(np.isfinite(z)) for z in [solution, momentum, divergence, pressure]
    ) or not np.all(np.isfinite([diss, balance])):
        raise ValidationError("Stokes solution exceeds float64 range", field="solution")
    force_scale = max(
        float(np.max(np.abs(rhs))),
        float(np.max(np.abs(a @ velocity))),
        float(np.max(np.abs(gradient @ pressure_relative))),
    )
    velocity_scale = float(np.max(np.abs(velocity)))
    relative = float(np.max(np.abs(momentum))) / (force_scale or 1.0)
    if (
        np.max(np.abs(momentum)) > 1e-12 + tol * force_scale
        or np.max(np.abs(divergence)) > 1e-12 + tol * velocity_scale * (1 / hx + 1 / hy)
        or abs(balance)
        > 1e-12 + tol * (abs(body_work) + abs(float(wall_work)) + float(diss))
    ):
        raise ValidationError(
            "Stokes momentum/mass/power gates failed", field="solution"
        )
    return StokesResult(
        _owned(pressure.reshape(ny, nx)),
        _owned(u),
        _owned(v),
        _owned(divergence.reshape(ny, nx)),
        _owned(divergence.reshape(ny, nx) * volume),
        _owned(momentum[:nu].reshape(ny, nx - 1)),
        _owned(momentum[nu:].reshape(ny - 1, nx)),
        float(diss),
        balance,
        relative,
        True,
    )
