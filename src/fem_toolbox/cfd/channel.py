"""Conservative fully developed channel momentum operator in SI units."""

from dataclasses import dataclass
from numbers import Real

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import diags  # type: ignore[import-untyped]
from scipy.sparse.linalg import spsolve  # type: ignore[import-untyped]

from fem_toolbox.errors import ValidationError


def _scalar(value: object, field: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValidationError(f"{field} must be a real SI scalar", field=field)
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValidationError(f"{field} exceeds float64 range", field=field) from error
    if not np.isfinite(result):
        raise ValidationError(f"{field} must be finite", field=field)
    return result


def _owned(values: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.frombuffer(values.tobytes(), dtype=np.float64)


@dataclass(frozen=True, slots=True)
class ChannelResult:
    """Cell velocities and face shear; flow per unit span, power per wall area."""

    centers_m: NDArray[np.float64]
    velocity_m_s: NDArray[np.float64]
    shear_Pa: NDArray[np.float64]
    cell_momentum_residual_Pa: NDArray[np.float64]
    flow_rate_m2_s: float
    momentum_balance_Pa: float
    dissipation_W_m2: float
    power_balance_W_m2: float


def solve_channel(
    *,
    faces_m: ArrayLike,
    viscosity_Pa_s: float,
    pressure_gradient_Pa_m: float,
    wall_velocities_m_s: tuple[float, float] = (0.0, 0.0),
) -> ChannelResult:
    """Solve -(mu u')'=-dp/dx with two prescribed wall velocities.

    faces_m has shape (n+1,), strictly increasing finite positions in m.
    Constant viscosity is positive (Pa s); pressure gradient may have either
    sign (Pa/m). Wall velocities are ordered lower/upper (m/s).
    Reject invalid data or numerically unresolved conservation with
    ValidationError. This is a reduced verification operator, not a RANS solver.
    """
    mu = _scalar(viscosity_Pa_s, "viscosity_Pa_s")
    g = -_scalar(pressure_gradient_Pa_m, "pressure_gradient_Pa_m")
    if mu <= 0:
        raise ValidationError("viscosity must be positive", field="viscosity_Pa_s")
    if (
        not isinstance(wall_velocities_m_s, (tuple, list))
        or len(wall_velocities_m_s) != 2
    ):
        raise ValidationError("need two wall velocities", field="wall_velocities_m_s")
    walls = np.array([_scalar(v, "wall_velocities_m_s") for v in wall_velocities_m_s])
    try:
        faces = np.asarray(faces_m, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValidationError("invalid face coordinates", field="faces_m") from error
    if faces.ndim != 1 or len(faces) < 2 or not np.all(np.isfinite(faces)):
        raise ValidationError("need finite 1D face coordinates", field="faces_m")
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        widths = np.diff(faces)
        centers = faces[:-1] + widths / 2
        distances = np.diff(np.concatenate(([faces[0]], centers, [faces[-1]])))
        conductance = mu / distances
    if (
        np.any(widths <= 0)
        or np.any(distances <= 0)
        or not np.all(np.isfinite(conductance))
    ):
        raise ValidationError("invalid or unresolved cell spacing", field="faces_m")
    rhs = g * widths
    rhs[0] += conductance[0] * walls[0]
    rhs[-1] += conductance[-1] * walls[1]
    matrix = diags(
        [-conductance[1:-1], conductance[:-1] + conductance[1:], -conductance[1:-1]],
        [-1, 0, 1],
        format="csc",
    )
    with np.errstate(over="ignore", invalid="ignore"):
        velocity = np.asarray(spsolve(matrix, rhs), dtype=np.float64)
        jumps = np.diff(np.concatenate(([walls[0]], velocity, [walls[1]])))
        shear = conductance * jumps
        residual = np.diff(shear) + g * widths
        flow = float(velocity @ widths)
        dissipation = float(np.sum(conductance * jumps * jumps))
        work = g * flow + shear[-1] * walls[1] - shear[0] * walls[0]
        balance = float(work - dissipation)
    if not (
        all(np.all(np.isfinite(a)) for a in [velocity, shear, residual])
        and np.all(np.isfinite([flow, dissipation, balance]))
    ):
        raise ValidationError(
            "channel solution exceeds float64 range", field="solution"
        )
    force_scale = float(np.max(np.abs(shear))) + abs(g) * float(np.sum(widths))
    power_scale = abs(float(work)) + dissipation
    if (
        np.max(np.abs(residual)) > 1e-10 * force_scale
        or abs(balance) > 1e-10 * power_scale
    ):
        raise ValidationError("channel conservation unresolved", field="solution")
    return ChannelResult(
        _owned(centers),
        _owned(velocity),
        _owned(shear),
        _owned(residual),
        flow,
        float(np.sum(residual)),
        dissipation,
        balance,
    )
