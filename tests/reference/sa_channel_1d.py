"""Independent continuous half-channel SA-noft2 BVP, without production CFD imports.

Balances integrate dU/dy=tau/(mu+mu_t) analytically in momentum. SciPy collocation
solves the scalar model with a wall and a centre-plane Neumann condition.
All formulas are independently written here; no Quad-FV/SA point routines are used.
"""

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.integrate import solve_bvp


@dataclass(frozen=True)
class ChannelReference:
    half_height_m: float
    density_kg_m3: float
    kinematic_viscosity_m2_s: float
    pressure_drive_N_m3: float
    friction_velocity_m_s: float
    solution: Any

    @property
    def re_tau(self):
        return (
            self.friction_velocity_m_s
            * self.half_height_m
            / self.kinematic_viscosity_m2_s
        )

    @property
    def bulk_velocity_m_s(self):
        return float(self.friction_velocity_m_s * self.solution.sol.integrate(0, 1)[0])

    @property
    def flow_per_span_m2_s(self):
        return 2 * self.half_height_m * self.bulk_velocity_m_s

    @property
    def wall_shear_Pa(self):
        return self.pressure_drive_N_m3 * self.half_height_m

    def evaluate(self, y_m):
        """Full-channel coordinates 0<=y<=2delta; return U and a in SI."""
        y = np.asarray(y_m, dtype=float)
        if (
            np.any(~np.isfinite(y))
            or np.any(y < 0)
            or np.any(y > 2 * self.half_height_m)
        ):
            raise ValueError("reference coordinate outside channel")
        wall_distance = np.minimum(y, 2 * self.half_height_m - y) / self.half_height_m
        values = self.solution.sol(wall_distance)
        velocity = self.friction_velocity_m_s * values[0]
        a = self.friction_velocity_m_s * self.half_height_m * values[1]
        # Enforce the exact Dirichlet trace, not a clamp of interior iterates.
        return np.where(wall_distance == 0, 0, velocity), np.where(
            wall_distance == 0, 0, a
        )


def solve_channel_reference(
    *,
    half_height_m=1.0,
    density_kg_m3=1.0,
    kinematic_viscosity_m2_s=1 / 587.19,
    pressure_drive_N_m3=1.0,
    tolerance=1e-8,
    initial_amplitude=1.0,
):
    """Nonzero turbulent branch; nondimensionalize by delta and u_tau.

    At the wall a~s*y. The limiting complete source cancellation gives
    d[(nu+a)*a']/dy=(a')² there, avoiding any evaluation of a/d at d=0.
    This wall limit is derived from the continuous equation, separately from FV.
    """
    parameters = [
        half_height_m,
        density_kg_m3,
        kinematic_viscosity_m2_s,
        pressure_drive_N_m3,
        tolerance,
        initial_amplitude,
    ]
    if not all(np.isfinite(p) and p > 0 for p in parameters):
        raise ValueError("positive finite reference parameters required")
    utau = np.sqrt(pressure_drive_N_m3 * half_height_m / density_kg_m3)
    nu = kinematic_viscosity_m2_s / (utau * half_height_m)

    def rhs(y, state):
        _, a, q = state
        chi = a / nu
        fv1 = chi**3 / (chi**3 + 7.1**3)
        up = (1 - y) / (nu + a * fv1)
        ap = q / (nu + a)
        omega = np.abs(up)
        interior = y > 0
        distance = np.where(interior, y, 1.0)
        fv2 = 1 - chi / (1 + chi * fv1)
        sbar = a * fv2 / (0.41**2 * distance**2)
        modified = omega + sbar
        branch = sbar < -0.7 * omega
        o, s = omega[branch], sbar[branch]
        modified[branch] = o + o * (0.49 * o + 0.9 * s) / (-0.5 * o - s)
        r = np.full_like(y, 10)
        active = modified > 0
        r[active] = np.minimum(
            a[active] / (modified[active] * 0.41**2 * distance[active] ** 2), 10
        )
        g = r + 0.3 * (r**6 - r)
        fw = g * (65 / (g**6 + 64)) ** (1 / 6)
        production = 0.1355 * modified * a
        destruction = (0.1355 / 0.41**2 + 1.622 / (2 / 3)) * fw * (a / distance) ** 2
        qp = -0.622 * ap**2 - (2 / 3) * (production - destruction)
        qp[~interior] = ap[~interior] ** 2
        return np.array([up, ap, qp])

    def boundary(left, right):
        return np.array([left[0], left[1], right[2]])

    y = np.r_[0, np.geomspace(1e-6, 1, 600)]
    a = initial_amplitude * 0.41 * y * (1 - y / 2)
    guess = np.array(
        [
            np.log1p(0.41 * y / nu) / 0.41,
            a,
            (nu + a) * initial_amplitude * 0.41 * (1 - y),
        ]
    )
    result = solve_bvp(rhs, boundary, y, guess, tol=tolerance, max_nodes=30000)
    if not result.success:
        raise RuntimeError(f"independent SA BVP failed: {result.message}")
    sampled = result.sol(np.linspace(0, 1, 2001))
    if np.min(sampled[1, 1:]) < 0 or np.max(sampled[1]) < nu:
        raise RuntimeError(
            "independent SA BVP did not find a nonnegative turbulent branch"
        )
    return ChannelReference(
        half_height_m,
        density_kg_m3,
        kinematic_viscosity_m2_s,
        pressure_drive_N_m3,
        float(utau),
        result,
    )
