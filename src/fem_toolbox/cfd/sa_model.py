"""SI point terms for SA-noft2 with the Allmaras S modification (TMR note 1c).

No trip, ft2, negative-SA extension, viscosity cap or implicit laminar switch.
The model r cap is exactly ten; a zero modified vorticity uses that same cap.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.fluid_mesh import _float_array, _freeze_float
from fem_toolbox.cfd.quad_fv import _error


def _sa_constants() -> dict[str, float]:
    """Exact fixed SA-noft2 constants for immutable run metadata; no calibration."""
    return dict(
        c_b1=0.1355,
        c_b2=0.622,
        sigma=2 / 3,
        kappa=0.41,
        c_v1=7.1,
        c_w2=0.3,
        c_w3=2,
        c_2=0.7,
        c_3=0.9,
        c_w1=0.1355 / 0.41**2 + 1.622 / (2 / 3),
    )


@dataclass(frozen=True, slots=True)
class SATerms:
    """Immutable cell arrays; sources in m²/s², diffusivities in m²/s.

    destruction_rate_s_inv times nu_tilde gives the nonnegative destruction
    source. It is suitable for the frozen implicit sink in the scalar solve.
    r and wall_function are dimensionless; modified vorticity has unit 1/s.
    """

    eddy_viscosity_m2_s: NDArray[np.float64]
    diffusivity_m2_s: NDArray[np.float64]
    production_m2_s2: NDArray[np.float64]
    destruction_m2_s2: NDArray[np.float64]
    gradient_source_m2_s2: NDArray[np.float64]
    destruction_rate_s_inv: NDArray[np.float64]
    modified_vorticity_s_inv: NDArray[np.float64]
    r: NDArray[np.float64]
    wall_function: NDArray[np.float64]


def sa_terms(
    nu_tilde_m2_s: ArrayLike,
    *,
    velocity_gradient_s_inv: ArrayLike,
    gradient_nu_tilde_m_s: ArrayLike,
    wall_distance_m: ArrayLike,
    kinematic_viscosity_m2_s: float,
) -> SATerms:
    """Evaluate the complete SA-noft2 source law at positive-distance cells.

    a is (ncell,), grad(U) is (ncell,2,2) with [component,derivative],
    grad(a) is (ncell,2), d is scalar or (ncell,). Molecular nu is a positive
    scalar. Negative/nonfinite a, d<=0, boolean data and unrepresentable terms
    raise ValidationError; no input or output is clipped to hide a failure.
    Wall Dirichlet values a=0 are imposed separately, never evaluated at d=0.
    """
    a = _float_array(nu_tilde_m2_s, "nu_tilde_m2_s")
    if a.ndim != 1 or not a.size or np.any(a < 0):
        _error("nu_tilde_m2_s", "nonempty nonnegative cell vector required")
    grad_u = _float_array(velocity_gradient_s_inv, "velocity_gradient_s_inv")
    grad_a = _float_array(gradient_nu_tilde_m_s, "gradient_nu_tilde_m_s")
    d = _float_array(wall_distance_m, "wall_distance_m")
    if d.ndim == 0:
        d = np.full(a.size, float(d))
    if grad_u.shape != (a.size, 2, 2) or grad_a.shape != (a.size, 2):
        _error("gradients", "expected (ncell,2,2) velocity and (ncell,2) SA gradients")
    if d.shape != a.shape or np.any(d <= 0):
        _error("wall_distance_m", "positive scalar or matching cell vector required")
    nu = _scalar(kinematic_viscosity_m2_s, "kinematic_viscosity_m2_s")
    if nu <= 0:
        _error("kinematic_viscosity_m2_s", "positive molecular viscosity required")
    try:
        with np.errstate(over="raise", divide="raise", invalid="raise"):
            chi = a / nu
            chi3 = chi**3
            fv1 = chi3 / (chi3 + 7.1**3)
            fv2 = 1 - chi / (1 + chi * fv1)
            omega = np.abs(grad_u[:, 0, 1] - grad_u[:, 1, 0])
            sbar = a * fv2 / (0.41**2 * d**2)
            modified = omega + sbar
            branch = sbar < -0.7 * omega
            o, s = omega[branch], sbar[branch]
            modified[branch] = o + o * (0.7**2 * o + 0.9 * s) / (-0.5 * o - s)
            r = np.full_like(a, 10)
            active = modified > 0
            r[active] = np.minimum(
                a[active] / (modified[active] * 0.41**2 * d[active] ** 2), 10
            )
            g = r + 0.3 * (r**6 - r)
            fw = g * (65 / (g**6 + 64)) ** (1 / 6)
            sink = (0.1355 / 0.41**2 + 1.622 / (2 / 3)) * fw * a / d**2
            fields = (
                a * fv1,
                (nu + a) / (2 / 3),
                0.1355 * modified * a,
                sink * a,
                (0.622 / (2 / 3)) * np.einsum("ij,ij->i", grad_a, grad_a),
                sink,
                modified,
                r,
                fw,
            )
            if any(not np.all(np.isfinite(field)) for field in fields):
                raise FloatingPointError("nonfinite SA point term")
    except FloatingPointError as error:
        _error("sa_terms", f"unrepresentable SI model state: {error}")
    return SATerms(*(_freeze_float(field) for field in fields))


def positive_sa_update(
    previous_m2_s: ArrayLike,
    candidate_m2_s: ArrayLike,
    *,
    relaxation: float = 0.5,
    max_attempts: int = 20,
) -> tuple[NDArray[np.float64], float, int]:
    """One common relaxed step, halved until finite and nonnegative everywhere.

    Returns immutable state, accepted alpha and number of halvings. No component
    is clipped. Failure identifies the last minimum and its zero-based cell.
    """
    old = _float_array(previous_m2_s, "previous_nu_tilde_m2_s")
    candidate = _float_array(candidate_m2_s, "candidate_nu_tilde_m2_s")
    if old.ndim != 1 or not old.size or old.shape != candidate.shape or np.any(old < 0):
        _error(
            "nu_tilde_m2_s",
            "matching cell vectors and nonnegative previous state required",
        )
    alpha = _scalar(relaxation, "sa_relaxation")
    if not 0 < alpha <= 1 or type(max_attempts) is not int or max_attempts < 1:
        _error(
            "sa_line_search",
            "relaxation in (0,1] and positive integer attempts required",
        )
    state = old.copy()
    for attempt in range(max_attempts):
        with np.errstate(over="ignore", invalid="ignore"):
            state = (1 - alpha) * old + alpha * candidate
        if np.all(np.isfinite(state)) and np.all(state >= 0):
            return _freeze_float(state), alpha, attempt
        alpha *= 0.5
    cell = int(np.argmin(state))
    _error(
        "sa_line_search",
        f"no admissible step after {max_attempts} attempts; cell {cell}: minimum={state[cell]:.9g} m²/s",
    )
