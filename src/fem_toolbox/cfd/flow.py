"""Stationary incompressible Navier-Stokes reference with Oseen/Picard steps."""

import json
from dataclasses import dataclass

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.flow_model import (
    CartesianFluidGrid,
    ChannelBoundary,
    ClosedWalls,
    IncompressibleFluid,
)
from fem_toolbox.cfd.mac import MACSystem, assemble, convection, linear_step
from fem_toolbox.cfd.stokes import _force, _owned
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


@dataclass(frozen=True, slots=True)
class StationaryFlowResult:
    """Immutable MAC solution and independently reconstructed SI diagnostics.

    u[ny,nx+1], v[ny+1,nx], p/divergence[ny,nx]. Momentum residuals live
    on free velocity faces; mass/power quantities are per unit span.
    """

    pressure_Pa: NDArray[np.float64]
    u_m_s: NDArray[np.float64]
    v_m_s: NDArray[np.float64]
    divergence_s_inv: NDArray[np.float64]
    cell_mass_balance_kg_s: NDArray[np.float64]
    momentum_x_residual_N_m3: NDArray[np.float64]
    momentum_y_residual_N_m3: NDArray[np.float64]
    inlet_mass_flow_kg_s: float
    outlet_mass_flow_kg_s: float
    global_mass_balance_kg_s: float
    global_momentum_balance_N_m: tuple[float, float]
    dissipation_W_m: float
    convective_power_W_m: float
    pressure_power_W_m: float
    power_balance_W_m: float
    iterations: int
    residual_history: tuple[float, ...]
    converged: bool
    provenance: RunProvenance


def _physical_momentum_terms(
    system: MACSystem, q: NDArray[np.float64], p: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Reconstruct viscous flux and pressure gradient without sparse products."""
    viscous = np.zeros_like(q)
    for i, j, coefficient in system.edges:
        flux = coefficient * (q[i] - q[j])
        viscous[i] += flux
        viscous[j] -= flux
    for i, value, coefficient in system.contacts:
        viscous[i] += coefficient * (q[i] - value)
    viscous /= system.weights
    nx, ny = system.grid.cells
    hx, hy = system.spacing
    field = p.reshape(ny, nx)
    gx = np.zeros(system.uid.shape)
    gy = np.zeros(system.vid.shape)
    gx[:, 1:-1] = np.diff(field, axis=1) / hx
    if isinstance(system.boundary, ChannelBoundary):
        gx[:, -1] = -2 * field[:, -1] / hx
    gy[1:-1, :] = np.diff(field, axis=0) / hy
    gradient = np.concatenate((gx[system.uid >= 0], gy[system.vid >= 0]))
    return viscous, gradient


def _powers(
    system: MACSystem,
    q: NDArray[np.float64],
    gradient: NDArray[np.float64],
    forcing: NDArray[np.float64],
    conv: NDArray[np.float64],
) -> tuple[float, float, float, float]:
    """Difference-based viscous work, explicit boundary/pressure/convection work."""
    diss = sum(c * (q[i] - q[j]) ** 2 for i, j, c in system.edges)
    diss += sum(c * (q[i] - value) ** 2 for i, value, c in system.contacts)
    boundary_work = sum(c * value * (value - q[i]) for i, value, c in system.contacts)
    pressure_work = -float(q @ (system.weights * gradient))
    conv_work = float(q @ (system.weights * conv))
    body_work = float(q @ (system.weights * forcing))
    return (
        float(diss),
        conv_work,
        pressure_work,
        float(body_work + boundary_work + pressure_work - diss - conv_work),
    )


def _check_outlet(system: MACSystem, q: NDArray[np.float64], tolerance: float) -> None:
    if isinstance(system.boundary, ChannelBoundary):
        u, _ = system.fields(q)
        if np.min(u[:, -1]) < -tolerance:
            raise ValidationError(
                "pressure outlet backflow is unsupported", field="outlet"
            )


def solve_stationary_flow(
    grid: CartesianFluidGrid,
    fluid: IncompressibleFluid,
    boundary: ChannelBoundary | ClosedWalls,
    *,
    force_x_N_m3: ArrayLike = 0.0,
    force_y_N_m3: ArrayLike = 0.0,
    pressure_reference_Pa: float | None = None,
    relaxation: float = 0.7,
    convection_scheme: str = "donor-cell",
    max_iterations: int = 200,
    relative_tolerance: float = 1e-10,
    momentum_tolerance_N_m3: float = 1e-12,
    mass_tolerance_kg_s: float = 1e-12,
    power_tolerance_W_m: float = 1e-12,
    max_unknowns: int = 50000,
) -> StationaryFlowResult:
    """Solve conservative steady 2D incompressible flow in SI units.

    Select donor-cell, linear-upwind (minmod deferred correction), or central
    (face Peclet <2). Constant material, monolithic pressure coupling and
    underrelaxed Picard iteration. ClosedWalls use p[0,0] as gauge (default 0).
    ChannelBoundary uses prescribed outlet pressure and rejects any extra gauge.
    Forces: scalar or free-face arrays u:(ny,nx-1) closed/(ny,nx) channel,
    v:(ny-1,nx). Gates test final nonlinear momentum, local/global mass and work.
    Raise ValidationError on unsupported backflow, invalid data, unresolved
    scales, linear solve failure or nonconvergence (with last residual).
    """
    if (
        not isinstance(grid, CartesianFluidGrid)
        or not isinstance(fluid, IncompressibleFluid)
        or not isinstance(boundary, (ChannelBoundary, ClosedWalls))
    ):
        raise ValidationError(
            "validated grid/fluid/boundary objects required", field="model"
        )
    if convection_scheme not in ("donor-cell", "central", "linear-upwind"):
        raise ValidationError("unknown convection scheme", field="convection_scheme")
    opened = isinstance(boundary, ChannelBoundary)
    if opened and pressure_reference_Pa is not None:
        raise ValidationError(
            "outlet pressure excludes additional pressure gauge",
            field="pressure_reference_Pa",
        )
    reference = (
        boundary.outlet_pressure_Pa
        if isinstance(boundary, ChannelBoundary)
        else _scalar(
            0 if pressure_reference_Pa is None else pressure_reference_Pa,
            "pressure_reference_Pa",
        )
    )
    alpha = _scalar(relaxation, "relaxation")
    tol = _scalar(relative_tolerance, "relative_tolerance")
    atols = [
        _scalar(v, n)
        for v, n in [
            (momentum_tolerance_N_m3, "momentum_tolerance_N_m3"),
            (mass_tolerance_kg_s, "mass_tolerance_kg_s"),
            (power_tolerance_W_m, "power_tolerance_W_m"),
        ]
    ]
    if not 0 < alpha <= 1 or not 0 < tol < 1 or min(atols) < 0:
        raise ValidationError("invalid relaxation or tolerances", field="parameters")
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValidationError(
            "positive integer iteration limit required", field="iterations"
        )
    nx, ny = grid.cells
    count = (
        ny * (nx if opened else nx - 1) + (ny - 1) * nx + nx * ny - (0 if opened else 1)
    )
    if type(max_unknowns) is not int or max_unknowns < 1 or count > max_unknowns:
        raise ValidationError("flow unknown limit exceeded", field="max_unknowns")
    system = assemble(grid, boundary, fluid.viscosity_Pa_s)
    hx, hy = system.spacing
    rho = fluid.density_kg_m3
    fx = _force(force_x_N_m3, (ny, nx if opened else nx - 1), "force_x_N_m3")
    fy = _force(force_y_N_m3, (ny - 1, nx), "force_y_N_m3")
    force = np.concatenate((fx.ravel(), fy.ravel()))
    rhs = force + system.boundary_rhs
    q, p = linear_step(system, system.diffusion, rhs)
    _check_outlet(system, q, 1e-12)
    history = []
    for _iteration in range(1, max_iterations + 1):
        c, constant = convection(system, q, rho, convection_scheme)
        candidate, pressure = linear_step(system, system.diffusion + c, rhs - constant)
        q = (1 - alpha) * q + alpha * candidate
        p = (1 - alpha) * p + alpha * pressure
        _check_outlet(system, q, 1e-12)
        final_c, final_constant = convection(system, q, rho, convection_scheme)
        conv = np.asarray(final_c @ q + final_constant)
        visc, grad = _physical_momentum_terms(system, q, p)
        momentum = visc + conv + grad - force
        u, v = system.fields(q)
        divergence = np.diff(u, axis=1) / hx + np.diff(v, axis=0) / hy
        local_mass = rho * hx * hy * divergence
        incoming = float(rho * hy * np.sum(u[:, 0]))
        outgoing = float(rho * hy * np.sum(u[:, -1]))
        diss, conv_work, pressure_work, power_balance = _powers(
            system, q, grad, force, conv
        )
        scale = max(float(np.max(np.abs(z))) for z in (rhs, visc, conv, grad))
        residual = float(np.max(np.abs(momentum)))
        history.append(residual / (scale or 1.0))
        mass_scale = max(
            abs(incoming), abs(outgoing), rho * max(hx, hy) * float(np.max(np.abs(q)))
        )
        power_scale = (
            abs(diss)
            + abs(conv_work)
            + abs(pressure_work)
            + abs(float(q @ (system.weights * force)))
        )
        if not all(
            np.all(np.isfinite(z)) for z in [q, p, momentum, divergence]
        ) or not np.all(
            np.isfinite(
                [
                    diss,
                    conv_work,
                    pressure_work,
                    power_balance,
                    scale,
                    mass_scale,
                    power_scale,
                ]
            )
        ):
            raise ValidationError(
                "nonlinear flow exceeds float64 range", field="solution"
            )
        if (
            residual <= atols[0] + tol * scale
            and np.max(np.abs(local_mass)) <= atols[1] + tol * mass_scale
            and abs(outgoing - incoming) <= atols[1] + tol * mass_scale
            and abs(power_balance) <= atols[2] + tol * power_scale
        ):
            break
    else:
        raise ValidationError(
            f"flow did not converge after {max_iterations} iterations; last momentum={residual:.6g} N/m³, relative={history[-1]:.6g}",
            field="iterations",
        )
    canonical = dict(
        cells=grid.cells,
        size_m=grid.size_m,
        density_kg_m3=rho,
        viscosity_Pa_s=fluid.viscosity_Pa_s,
        boundary={
            "type": "channel",
            "inlet_velocity_m_s": boundary.inlet_velocity_m_s.tolist(),
            "outlet_pressure_Pa": reference,
        }
        if isinstance(boundary, ChannelBoundary)
        else {
            "type": "closed",
            "walls": boundary.tangential_velocities_m_s,
            "pressure_reference_Pa": reference,
        },
        force_x_N_m3=fx.tolist(),
        force_y_N_m3=fy.tolist(),
    )
    provenance = RunProvenance.from_input_bytes(
        json.dumps(canonical, sort_keys=True, allow_nan=False).encode(),
        solver_name="experimental-mac-picard",
        solver_parameters=dict(
            convection=convection_scheme,
            numpy_version=np.__version__,
            scipy_version=scipy.__version__,
            revision=1,
            relaxation=alpha,
            max_iterations=max_iterations,
            relative_tolerance=tol,
            absolute_tolerances=atols,
            max_unknowns=max_unknowns,
        ),
    )
    return StationaryFlowResult(
        _owned((p + reference).reshape(ny, nx)),
        _owned(u),
        _owned(v),
        _owned(divergence),
        _owned(local_mass),
        _owned(momentum[: system.nu].reshape(fx.shape)),
        _owned(momentum[system.nu :].reshape(fy.shape)),
        incoming,
        outgoing,
        outgoing - incoming,
        (
            float(momentum[: system.nu] @ system.weights[: system.nu]),
            float(momentum[system.nu :] @ system.weights[system.nu :]),
        ),
        diss,
        conv_work,
        pressure_work,
        power_balance,
        len(history),
        tuple(history),
        True,
        provenance,
    )
