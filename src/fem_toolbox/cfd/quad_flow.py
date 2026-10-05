"""Stationary laminar Quad-FV flow with Krylov-accelerated SIMPLE steps."""

from __future__ import annotations

import json
from dataclasses import dataclass
from hashlib import sha256
from typing import Any

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import bmat, diags  # type: ignore[import-untyped]
from scipy.sparse.linalg import (  # type: ignore[import-untyped]
    LinearOperator,
    gmres,
    splu,
)

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import _float_array, _freeze_float
from fem_toolbox.cfd.quad_fv import _error, _validated_domain
from fem_toolbox.cfd.quad_system import Map, QuadSystem
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


class QuadFlowConvergenceError(ValidationError):
    """Nonconvergence with immutable complete outer residual history."""

    def __init__(self, history: list[tuple[float, ...]]) -> None:
        self.residual_history = tuple(history)
        super().__init__(
            f"quad flow did not converge after {len(history)} iterations; last residuals={history[-1]}",
            field="iterations",
        )


@dataclass(frozen=True, slots=True)
class QuadFlowResult:
    """Immutable cell/face SI fields and freshly reconstructed diagnostics.

    U[ncell,2], p[ncell], mass_flux[nface] in mesh.faces order. Pressure is
    the periodic part plus gauge, if periodic driving is present. History
    rows contain momentum, continuity, global mass, global momentum, power.
    """

    velocity_m_s: NDArray[np.float64]
    pressure_Pa: NDArray[np.float64]
    mass_flux_kg_s: NDArray[np.float64]
    face_velocity_m_s: NDArray[np.float64]
    mass_residual_kg_m3_s: NDArray[np.float64]
    momentum_residual_N_m3: NDArray[np.float64]
    nonorthogonal_correction_residual_N_m3: NDArray[np.float64]
    global_mass_balance_kg_s: float
    global_momentum_balance_N: NDArray[np.float64]
    inlet_mass_flow_kg_s: float
    outlet_mass_flow_kg_s: float
    wall_traction_Pa: NDArray[np.float64]
    viscous_dissipation_W: float
    gradient_dissipation_W: float
    viscous_discretization_power_W: float
    convective_power_W: float
    numerical_transport_power_W: float
    pressure_boundary_power_W: float
    pressure_interpolation_power_W: float
    body_power_W: float
    viscous_boundary_power_W: float
    power_balance_W: float
    iterations: int
    simple_steps: int
    residual_history: tuple[tuple[float, ...], ...]
    converged: bool
    provenance: RunProvenance


def _field(
    value: ArrayLike, n: int, *, vector: bool, field: str
) -> NDArray[np.float64]:
    a = _float_array(value, field)
    shape = (n, 2) if vector else (n,)
    if (not vector and a.ndim == 0) or (vector and a.shape == (2,)):
        a = np.broadcast_to(a, shape).copy()
    if a.shape != shape:
        _error(field, f"expected {shape} SI values")
    return a


def _driving(domain: FluidDomain) -> NDArray[np.float64]:
    shifts, drops = [], []
    for name, bc in domain.boundaries.items():
        if bc["type"] == "periodic" and name < bc["partner"]:
            shifts.append(bc["translation_m"])
            drops.append(bc["pressure_drop_Pa"])
    if not shifts:
        return np.zeros(2)
    drive = np.linalg.lstsq(np.array(shifts), np.array(drops), rcond=None)[0]
    if not np.allclose(np.array(shifts) @ drive, drops, atol=1e-12, rtol=1e-12):
        _error("periodic", "incompatible pressure drops")
    return np.asarray(drive)


def _linear_simple(
    system: QuadSystem,
    momentum: Map,
    rhs: NDArray[np.float64],
    diagonal: NDArray[np.float64],
    pressure_flux: Map,
    *,
    max_steps: int,
    momentum_gate: float,
    continuity_gate: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], int]:
    """Sparse-LU momentum/compact-pressure solves form the SIMPLE preconditioner.

    GMRES accelerates these pressure/velocity corrections; it is not a hidden
    coupled LU solve. Full nonorthogonal terms remain in the Oseen operator.
    """
    n = system.n
    a = momentum.matrix
    g = system.pressure_momentum.matrix
    divergence = system.mass_velocity.matrix.copy().tolil()
    h = (system.balance @ pressure_flux.matrix).tolil()
    mass_rhs = -system.mass_velocity.constant - np.asarray(
        system.balance @ pressure_flux.constant
    )
    gauge = system.domain.pressure_reference_cell
    compact = system.fv.pressure_correction_matrix(
        momentum_diagonal_kg_s=diagonal,
        density_kg_m3=system.domain.fluid.density_kg_m3,
        pressure_boundary_indices=system.pressure_rows,
    ).tolil()
    if gauge is not None:
        divergence[gauge, :] = 0
        h[gauge, :] = 0
        h[gauge, gauge] = 1
        mass_rhs[gauge] = 0
        compact[gauge, :] = 0
        compact[gauge, gauge] = 1
    divergence, h = divergence.tocsr(), h.tocsr()
    operator = bmat([[a, g], [divergence, h]], format="csr")
    target = np.concatenate((rhs, mass_rhs))
    scale = np.concatenate(
        (
            np.repeat(1 / (system.fv.cell_volumes_m3 * momentum_gate), 2),
            1 / (system.fv.cell_volumes_m3 * continuity_gate),
        )
    )
    if gauge is not None:
        scale[2 * n + gauge] = 1
    scaled = diags(scale) @ operator
    try:
        lu_u, lu_p = splu(a.tocsc()), splu(compact.tocsc())
    except RuntimeError as error:
        _error("linear_solver", f"SIMPLE block factorization failed: {error}")
    steps = [0]

    def precondition(residual: NDArray[np.float64]) -> NDArray[np.float64]:
        raw = residual / scale
        prediction = np.asarray(lu_u.solve(raw[: 2 * n]))
        correction = np.asarray(lu_p.solve(raw[2 * n :] - divergence @ prediction))
        prediction -= np.asarray(g @ correction) / diagonal.ravel()
        return np.concatenate((prediction, correction))

    def count(_residual: float) -> None:
        steps[0] += 1

    pre = LinearOperator((3 * n, 3 * n), matvec=precondition, dtype=np.float64)
    state = np.zeros(3 * n)
    # Refinement includes the entire nonorthogonal operator, at most 20 passes.
    for _correction in range(20):
        residual = (target - operator @ state) * scale
        if np.max(np.abs(residual)) <= 0.5:
            break
        try:
            delta, info = gmres(
                scaled,
                residual,
                M=pre,
                rtol=0,
                atol=0.25,
                restart=80,
                maxiter=max_steps,
                callback=count,
                callback_type="legacy",
            )
        except (ValueError, RuntimeError) as error:
            _error("linear_solver", f"SIMPLE iteration failed: {error}")
        state += delta
        reconstructed = float(np.max(np.abs((target - operator @ state) * scale)))
        if info != 0 and reconstructed >= 0.99 * float(np.max(np.abs(residual))):
            _error(
                "simple_iterations",
                f"SIMPLE did not converge in {max_steps} corrections; last linear residual={np.max(np.abs((target - operator @ state) * scale)):.6g}",
            )
    else:
        _error(
            "nonorthogonal_iterations",
            "nonorthogonal correction failed after 20 refinements",
        )
    if not np.all(np.isfinite(state)):
        _error("solution", "SIMPLE exceeds float64 range")
    return state[: 2 * n], state[2 * n :], steps[0]


def _diagnostics(
    system: QuadSystem,
    u: NDArray[np.float64],
    p: NDArray[np.float64],
    flux: NDArray[np.float64],
    transported: NDArray[np.float64],
    force: NDArray[np.float64],
) -> dict[str, Any]:
    """Direct face contractions, not an LU residual or cached convergence flag."""
    fv = system.fv
    cell_u = u.reshape(system.n, 2)
    face_u = system.uf.apply(u).reshape(system.nf, 2)
    owner, neighbour = system.mesh.face_owners, fv.neighbours
    inner_mask = neighbour >= 0
    polynomial = system.poly_u.apply(u).reshape(system.n, 2, 9)
    face_gradient = np.zeros((system.nf, 2, 2))
    third_difference = np.zeros((system.nf, 2))
    owner_offset = system.mesh.face_centres_m - system.mesh.cell_centres_m[owner]
    fractions = np.where(inner_mask, fv._fraction, 0.0)
    for cells, weights, offsets in (
        (owner, 1 - fractions, owner_offset),
        (np.maximum(neighbour, 0), fractions, owner_offset - fv._direction),
    ):
        scale = system.poly_scale[cells]
        x, y = (offsets / scale[:, None]).T
        zero, one = np.zeros(system.nf), np.ones(system.nf)
        evaluation = np.stack(
            (
                np.stack(
                    (one, zero, x, y, zero, x * x / 2, x * y, y * y / 2, zero), axis=1
                ),
                np.stack(
                    (zero, one, zero, x, y, zero, x * x / 2, x * y, y * y / 2), axis=1
                ),
            ),
            axis=1,
        )
        face_gradient += (
            weights[:, None, None]
            * np.einsum("fkj,faj->fka", polynomial[cells], evaluation)
            / scale[:, None, None]
        )
        dx, dy = (fv._direction / scale[:, None]).T
        cubic = np.column_stack((dx**3, 3 * dx * dx * dy, 3 * dx * dy * dy, dy**3))
        third_difference += weights[:, None] * np.einsum(
            "fkj,fj->fk", polynomial[cells, :, 5:], cubic
        )
    difference = face_u - cell_u[owner]
    difference[inner_mask] = cell_u[neighbour[inner_mask]] - cell_u[owner[inner_mask]]
    third_difference[~inner_mask] = 0
    correction = (
        difference
        - np.einsum("ijk,ik->ij", face_gradient, fv._direction)
        - third_difference / 24
    ) / fv._distance[:, None]
    face_gradient += correction[:, :, None] * system.mesh.face_normals[:, None, :]
    face_gradient.ravel()[system.physical_gradient_mask] = (
        system.physical_gradient.apply(u)[system.physical_gradient_mask]
    )
    stress = system.mu_face[:, None, None] * (
        face_gradient + np.swapaxes(face_gradient, 1, 2)
    )
    visc = system.alias_flux(
        np.einsum("ijk,ik->ij", stress, system.mesh.face_normals)
        * fv.face_areas_m2[:, None]
    )
    pressure = system.alias_flux(
        system.pf.apply(p)[:, None]
        * system.mesh.face_normals
        * fv.face_areas_m2[:, None]
    )
    mapped_visc = system.alias_flux(system.viscous_force.apply(u).reshape(system.nf, 2))
    nonorthogonal_residual = (
        fv.cell_balance(visc - mapped_visc) / fv.cell_volumes_m3[:, None]
    )
    conv = system.alias_flux(flux[:, None] * transported.reshape(system.nf, 2))
    momentum = (
        fv.cell_balance(conv - visc + pressure) - force * fv.cell_volumes_m3[:, None]
    )
    mass = fv.cell_balance(flux)
    exterior = fv.neighbours < 0
    inner = fv.active_faces[fv.neighbours[fv.active_faces] >= 0]
    outer = fv.active_faces[exterior[fv.active_faces]]
    jump = cell_u[fv.neighbours[inner]] - cell_u[system.mesh.face_owners[inner]]
    dissipation = float(np.sum(jump * visc[inner]))
    boundary_jump = face_u[outer] - cell_u[system.mesh.face_owners[outer]]
    dissipation += float(np.sum(boundary_jump * visc[outer]))
    boundary_work = float(np.sum(face_u[outer] * visc[outer]))
    gradient = system.gu.apply(u).reshape(system.n, 2, 2)
    strain = (gradient + np.swapaxes(gradient, 1, 2)) / 2
    gradient_diss = float(
        np.sum(2 * system.mu_cell * fv.cell_volumes_m3 * np.sum(strain**2, axis=(1, 2)))
    )
    body_work = float(np.sum(cell_u * force * fv.cell_volumes_m3[:, None]))
    conv_work = float(np.sum(cell_u * fv.cell_balance(conv)))
    conv_boundary = float(np.sum(flux[outer] * np.sum(face_u[outer] ** 2, axis=1) / 2))
    pressure_work = -float(np.sum(cell_u * fv.cell_balance(pressure)))
    pf = system.pf.apply(p) + system.reference
    pressure_boundary = -float(
        pf[outer] @ flux[outer] / system.domain.fluid.density_kg_m3
    )
    pressure_interpolation = pressure_work - pressure_boundary
    power = (
        body_work
        + boundary_work
        + pressure_boundary
        + pressure_interpolation
        - dissipation
        - conv_work
    )
    wall_faces = system.mesh.boundary_face_indices[system.walls]
    traction = (visc[wall_faces] - pressure[wall_faces]) / fv.face_areas_m2[
        wall_faces, None
    ] - system.reference * system.mesh.face_normals[wall_faces]
    return dict(
        momentum=momentum / fv.cell_volumes_m3[:, None],
        nonorthogonal=nonorthogonal_residual,
        mass=mass / fv.cell_volumes_m3,
        global_momentum=momentum.sum(axis=0),
        global_mass=float(np.sum(flux[outer])),
        inlet=-float(np.sum(flux[system.kind == "inlet"])),
        outlet=float(np.sum(flux[system.kind == "pressure_outlet"])),
        dissipation=dissipation,
        gradient_dissipation=gradient_diss,
        viscous_discretization=dissipation - gradient_diss,
        body_work=body_work,
        boundary_work=boundary_work,
        convective_power=conv_boundary,
        numerical_transport=conv_work - conv_boundary,
        pressure_boundary=pressure_boundary,
        pressure_interpolation=pressure_interpolation,
        power=power,
        traction=traction,
        face_u=face_u,
    )


def _json(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict) or hasattr(value, "items"):
        return {k: _json(v) for k, v in value.items()}
    return value


def solve_quad_flow(
    domain: FluidDomain,
    *,
    velocity_scale_m_s: float,
    length_scale_m: float,
    force_N_m3: ArrayLike = (0.0, 0.0),
    viscosity_Pa_s: ArrayLike | None = None,
    boundary_viscosity_Pa_s: ArrayLike | None = None,
    pressure_reference_Pa: float | None = None,
    convection_scheme: str = "donor-cell",
    velocity_relaxation: float = 0.5,
    pressure_relaxation: float = 0.5,
    max_iterations: int = 500,
    max_simple_steps: int = 500,
    max_unknowns: int = 50000,
    model: str = "laminar",
) -> QuadFlowResult:
    """Solve explicit laminar incompressible Quad-FV flow (not SA/RANS/FSI).

    Scales U_ref,L_ref are mandatory positive SI. Force is (2,) or (nc,2).
    Optional positive mu/mub fields must be provided together, shapes (nc,)/(nb,)
    or scalar. FluidDomain supplies fixed rho, width and physical boundaries.
    Full variable symmetric stress, component pressure mobility, donor-cell,
    central (|Pe|<2) or bounded linear-upwind; no implicit scheme substitution.
    Sparse LU block solves with GMRES-accelerated SIMPLE pressure corrections.
    Final absolute gates are those of Doc/71, no relaxed tolerance fallback.
    ValidationError covers backflow, nonconvergence, scales and resource limits.
    """
    domain = _validated_domain(domain)
    n, nb = len(domain.mesh.quad4), len(domain.mesh.line2)
    if model != "laminar":
        _error("model", "only explicit laminar model implemented; SA is #87")
    if convection_scheme not in ("donor-cell", "central", "linear-upwind"):
        _error("convection_scheme", "unknown transport scheme")
    for field, value in (
        ("max_iterations", max_iterations),
        ("max_simple_steps", max_simple_steps),
        ("max_unknowns", max_unknowns),
    ):
        if type(value) is not int or value < 1:
            _error(field, "positive integer limit required")
    if 3 * n > max_unknowns:
        _error("max_unknowns", "flow unknown limit exceeded")
    us, ls = (
        _scalar(velocity_scale_m_s, "velocity_scale_m_s"),
        _scalar(length_scale_m, "length_scale_m"),
    )
    au, ap = (
        _scalar(velocity_relaxation, "velocity_relaxation"),
        _scalar(pressure_relaxation, "pressure_relaxation"),
    )
    if min(us, ls) <= 0 or not 0 < au <= 1 or not 0 < ap <= 1:
        _error("parameters", "positive SI scales and relaxation in (0,1] required")
    if domain.pressure_reference_cell is None and pressure_reference_Pa is not None:
        _error("pressure_reference_Pa", "outlet excludes an additional gauge")
    reference = _scalar(
        0 if pressure_reference_Pa is None else pressure_reference_Pa,
        "pressure_reference_Pa",
    )
    if (viscosity_Pa_s is None) != (boundary_viscosity_Pa_s is None):
        _error("viscosity_Pa_s", "cell and boundary viscosities required together")
    mu = _field(
        domain.fluid.viscosity_Pa_s if viscosity_Pa_s is None else viscosity_Pa_s,
        n,
        vector=False,
        field="viscosity_Pa_s",
    )
    mub = _field(
        domain.fluid.viscosity_Pa_s
        if boundary_viscosity_Pa_s is None
        else boundary_viscosity_Pa_s,
        nb,
        vector=False,
        field="boundary_viscosity_Pa_s",
    )
    if min(mu.min(), mub.min()) <= 0:
        _error("viscosity_Pa_s", "positive viscosity required")
    force = _field(force_N_m3, n, vector=True, field="force_N_m3") + _driving(domain)
    system = QuadSystem(domain, mu, mub)
    if domain.pressure_reference_cell is not None:
        system.reference = reference
    rho, area_ref = domain.fluid.density_kg_m3, domain.width_m * ls
    gates = np.array(
        [
            1e-8 * rho * us**2 / ls,
            1e-8 * rho * us / ls,
            1e-10 * rho * us * area_ref,
            1e-8 * rho * us**2 * area_ref,
            1e-6 * rho * us**3 * area_ref,
        ]
    )
    if np.any(gates <= 0) or not np.all(np.isfinite(gates)):
        _error("scales", "unrepresentable absolute residual scales")
    u, p = np.zeros(2 * n), np.zeros(n)
    flux = system.alias_flux(system.velocity_flux.apply(u))
    history: list[tuple[float, ...]] = []
    steps = 0
    for _iteration in range(max_iterations):
        transport, _ = system.transport(flux, u, convection_scheme)
        momentum = Map(
            (system.diffusion.matrix + transport.matrix).tocsr(),
            system.diffusion.constant + transport.constant,
        )
        diagonal = momentum.matrix.diagonal().reshape(n, 2)
        if np.any(diagonal <= 0) or not np.all(np.isfinite(diagonal)):
            _error("momentum_diagonal_kg_s", "nonpositive/unresolved momentum diagonal")
        pflux = system.pressure_flux(diagonal)
        rhs = (
            (force * system.fv.cell_volumes_m3[:, None]).ravel()
            - momentum.constant
            - system.pressure_momentum.constant
        )
        candidate, pressure, count = _linear_simple(
            system,
            momentum,
            rhs,
            diagonal,
            pflux,
            max_steps=max_simple_steps,
            momentum_gate=float(gates[0]),
            continuity_gate=float(gates[1]),
        )
        steps += count
        u = (1 - au) * u + au * candidate
        p = (1 - ap) * p + ap * pressure
        flux = system.alias_flux(system.velocity_flux.apply(u) + pflux.apply(p))
        # Rebuild convection, the momentum diagonal and RC flux on final iterands.
        final_transport, transported = system.transport(flux, u, convection_scheme)
        final_diagonal = (
            (system.diffusion.matrix + final_transport.matrix).diagonal().reshape(n, 2)
        )
        final_pflux = system.pressure_flux(final_diagonal)
        flux = system.alias_flux(system.velocity_flux.apply(u) + final_pflux.apply(p))
        _, transported = system.transport(flux, u, convection_scheme)
        diag = _diagnostics(system, u, p, flux, transported, force)
        row = (
            float(np.max(np.abs(diag["momentum"]))),
            float(np.max(np.abs(diag["mass"]))),
            abs(diag["global_mass"]),
            float(np.max(np.abs(diag["global_momentum"]))),
            abs(diag["power"]),
        )
        history.append(row)
        if (
            not np.all(np.isfinite(row))
            or not np.all(np.isfinite(u))
            or not np.all(np.isfinite(p))
        ):
            _error("solution", "nonlinear flow exceeds float64 range")
        outlet_faces = np.flatnonzero(system.kind == "pressure_outlet")
        bad = outlet_faces[flux[outlet_faces] < -1e-12 * rho * us * area_ref]
        if len(bad):
            face = int(bad[0])
            _error(
                f"outlet.face[{face}]",
                f"unsupported backflow mass_flux={flux[face]:.6g} kg/s",
            )
        if (
            np.all(np.array(row) <= gates)
            and np.max(np.abs(diag["nonorthogonal"])) <= gates[0]
        ):
            break
    else:
        raise QuadFlowConvergenceError(history)
    canonical = dict(
        mesh=domain.mesh.to_document(),
        width_m=domain.width_m,
        fluid=dict(density_kg_m3=rho, viscosity_Pa_s=domain.fluid.viscosity_Pa_s),
        boundaries=_json(domain.boundaries),
        pressure_reference_cell=domain.pressure_reference_cell,
        pressure_reference_Pa=system.reference,
        force_N_m3=force.tolist(),
        viscosity_Pa_s=mu.tolist(),
        boundary_viscosity_Pa_s=mub.tolist(),
    )
    provenance = RunProvenance.from_input_bytes(
        json.dumps(canonical, sort_keys=True, allow_nan=False).encode(),
        solver_name="quad-fv-simple-gmres",
        solver_parameters=dict(
            revision=1,
            geometry_sha256=sha256(
                json.dumps(canonical["mesh"], sort_keys=True, allow_nan=False).encode()
            ).hexdigest(),
            residual_history=history,
            model=model,
            convection=convection_scheme,
            velocity_scale_m_s=us,
            length_scale_m=ls,
            velocity_relaxation=au,
            pressure_relaxation=ap,
            max_iterations=max_iterations,
            max_simple_steps=max_simple_steps,
            max_unknowns=max_unknowns,
            nonorthogonal_refinement_limit=20,
            absolute_gates=gates.tolist(),
            numpy_version=np.__version__,
            scipy_version=scipy.__version__,
        ),
    )
    return QuadFlowResult(
        _freeze_float(u.reshape(n, 2)),
        _freeze_float(p + system.reference),
        _freeze_float(flux),
        _freeze_float(diag["face_u"]),
        _freeze_float(diag["mass"]),
        _freeze_float(diag["momentum"]),
        _freeze_float(diag["nonorthogonal"]),
        diag["global_mass"],
        _freeze_float(diag["global_momentum"]),
        diag["inlet"],
        diag["outlet"],
        _freeze_float(diag["traction"]),
        diag["dissipation"],
        diag["gradient_dissipation"],
        diag["viscous_discretization"],
        diag["convective_power"],
        diag["numerical_transport"],
        diag["pressure_boundary"],
        diag["pressure_interpolation"],
        diag["body_work"],
        diag["boundary_work"],
        diag["power"],
        len(history),
        steps,
        tuple(history),
        True,
        provenance,
    )
