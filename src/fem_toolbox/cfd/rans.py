"""Coupled stationary 2D SA-noft2/SIMPLE with fresh field and wall-resolution gates."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import ArrayLike, NDArray
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import _freeze_float
from fem_toolbox.cfd.quad_flow import (
    QuadFlowResult,
    _diagnostics,
    _driving,
    _field,
    _json,
    _linear_simple,
)
from fem_toolbox.cfd.quad_fv import _error, _validated_domain
from fem_toolbox.cfd.quad_system import Map, QuadSystem
from fem_toolbox.cfd.sa_model import _sa_constants, positive_sa_update
from fem_toolbox.cfd.sa_transport import _sa_operator, _SATransportSystem
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


class RansConvergenceError(ValidationError):
    """Failed joint iteration with immutable last fields and complete history."""

    def __init__(
        self,
        reason: str,
        history: list[tuple[float, ...]],
        u: NDArray[np.float64],
        p: NDArray[np.float64],
        a: NDArray[np.float64],
    ) -> None:
        self.residual_history = tuple(history)
        self.velocity_m_s = _freeze_float(u.reshape(-1, 2))
        self.pressure_Pa = _freeze_float(p)
        self.nu_tilde_m2_s = _freeze_float(a)
        super().__init__(
            f"RANS did not converge: {reason}; iterations={len(history)}, last={history[-1] if history else None}, minimum SA={a.min():.9g} m²/s at cell {a.argmin()}",
            field="rans_iterations",
        )


@dataclass(frozen=True, slots=True)
class RansResult:
    """Joint solution: flow, SA cell fields and actual wall-resolution diagnostics.

    History columns: max momentum N/m³, mass kg/(m³ s), SA m²/s²,
    global mass kg/s, max global momentum N, power W, accepted SA alpha.
    wall_y_plus follows domain.wall_quality.boundary_indices. No DNS/model or
    FSI qualification is implied by convergence of these discrete equations.
    """

    flow: QuadFlowResult
    nu_tilde_m2_s: NDArray[np.float64]
    eddy_viscosity_m2_s: NDArray[np.float64]
    effective_viscosity_Pa_s: NDArray[np.float64]
    boundary_effective_viscosity_Pa_s: NDArray[np.float64]
    sa_residual_m2_s2: NDArray[np.float64]
    sa_global_balance_m5_s2: float
    wall_y_plus: NDArray[np.float64]
    iterations: int
    residual_history: tuple[tuple[float, ...], ...]
    converged: bool
    provenance: RunProvenance


def _eddy(a: NDArray[np.float64], nu: float) -> NDArray[np.float64]:
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            chi3 = (a / nu) ** 3
            result = a * chi3 / (chi3 + 7.1**3)
            if not np.all(np.isfinite(result)):
                raise FloatingPointError("nonfinite eddy viscosity")
    except FloatingPointError as error:
        _error("nu_tilde_m2_s", f"unrepresentable SA viscosity state: {error}")
    return result


def _reconstruct(
    system: QuadSystem,
    u: NDArray[np.float64],
    p: NDArray[np.float64],
    flux: NDArray[np.float64],
    force: NDArray[np.float64],
    scheme: str,
) -> tuple[NDArray[np.float64], dict[str, Any]]:
    transport, _ = system.transport(flux, u, scheme)
    diagonal = (
        (system.diffusion.matrix + transport.matrix).diagonal().reshape(system.n, 2)
    )
    if np.any(diagonal <= 0) or not np.all(np.isfinite(diagonal)):
        _error("momentum_diagonal_kg_s", "nonpositive/unresolved RANS diagonal")
    pressure_flux = system.pressure_flux(diagonal)
    fresh_flux = system.alias_flux(
        system.velocity_flux.apply(u) + pressure_flux.apply(p)
    )
    _, transported = system.transport(fresh_flux, u, scheme)
    return fresh_flux, _diagnostics(system, u, p, fresh_flux, transported, force)


def _flow_result(
    system: QuadSystem,
    u: NDArray[np.float64],
    p: NDArray[np.float64],
    flux: NDArray[np.float64],
    diag: dict[str, Any],
    history: list[tuple[float, ...]],
    steps: int,
    provenance: RunProvenance,
) -> QuadFlowResult:
    return QuadFlowResult(
        velocity_m_s=_freeze_float(u.reshape(-1, 2)),
        pressure_Pa=_freeze_float(p + system.reference),
        mass_flux_kg_s=_freeze_float(flux),
        face_velocity_m_s=_freeze_float(diag["face_u"]),
        mass_residual_kg_m3_s=_freeze_float(diag["mass"]),
        momentum_residual_N_m3=_freeze_float(diag["momentum"]),
        nonorthogonal_correction_residual_N_m3=_freeze_float(diag["nonorthogonal"]),
        global_mass_balance_kg_s=diag["global_mass"],
        global_momentum_balance_N=_freeze_float(diag["global_momentum"]),
        inlet_mass_flow_kg_s=diag["inlet"],
        outlet_mass_flow_kg_s=diag["outlet"],
        wall_traction_Pa=_freeze_float(diag["traction"]),
        viscous_dissipation_W=diag["dissipation"],
        gradient_dissipation_W=diag["gradient_dissipation"],
        viscous_discretization_power_W=diag["viscous_discretization"],
        convective_power_W=diag["convective_power"],
        numerical_transport_power_W=diag["numerical_transport"],
        pressure_boundary_power_W=diag["pressure_boundary"],
        pressure_interpolation_power_W=diag["pressure_interpolation"],
        body_power_W=diag["body_work"],
        viscous_boundary_power_W=diag["boundary_work"],
        power_balance_W=diag["power"],
        iterations=len(history),
        simple_steps=steps,
        residual_history=tuple(
            (row[0], row[1], row[3], row[4], row[5]) for row in history
        ),
        converged=True,
        provenance=provenance,
    )


def solve_rans(
    domain: FluidDomain,
    *,
    initial_nu_tilde_m2_s: ArrayLike,
    velocity_scale_m_s: float,
    length_scale_m: float,
    force_N_m3: ArrayLike = (0.0, 0.0),
    pressure_reference_Pa: float | None = None,
    convection_scheme: str = "linear-upwind",
    velocity_relaxation: float = 0.5,
    pressure_relaxation: float = 0.5,
    sa_relaxation: float = 0.5,
    max_iterations: int = 500,
    max_simple_steps: int = 500,
    max_unknowns: int = 50000,
) -> RansResult:
    """Joint Picard iteration for physical wall-resolved SA-noft2 and momentum.

    Only explicit SA-noft2 is solved: at least one initial cell must have a>0,
    all must be nonnegative; no inferred laminar mode or viscosity cap. Width,
    physical BCs, positive SI scales and all inlet SA data come from the caller.
    The same conservative Rhie-Chow mass flux transports momentum and SA.
    Frozen nonnegative destruction is implicit; all final coefficients,
    stresses, sources and field/balance residuals are recomputed at the accepted
    joint state. Positivity backtracking reduces U/p/a updates together.
    Actual wall shear must satisfy y+<=1 and layer growth<=1.2. Failure retains
    the last complete joint history and fields; no wall function is substituted.
    """
    domain = _validated_domain(domain)
    n = len(domain.mesh.quad4)
    for name, value in (
        ("max_iterations", max_iterations),
        ("max_simple_steps", max_simple_steps),
        ("max_unknowns", max_unknowns),
    ):
        if type(value) is not int or value < 1:
            _error(name, "positive integer budget required")
    if 4 * n > max_unknowns:
        _error("max_unknowns", "joint U/p/SA resource limit exceeded")
    us, ls = (
        _scalar(velocity_scale_m_s, "velocity_scale_m_s"),
        _scalar(length_scale_m, "length_scale_m"),
    )
    au, ap, aa = (
        _scalar(v, k)
        for v, k in (
            (velocity_relaxation, "velocity_relaxation"),
            (pressure_relaxation, "pressure_relaxation"),
            (sa_relaxation, "sa_relaxation"),
        )
    )
    if min(us, ls) <= 0 or not all(0 < x <= 1 for x in (au, ap, aa)):
        _error("parameters", "positive SI scales and relaxations in (0,1] required")
    if convection_scheme not in ("donor-cell", "linear-upwind"):
        _error(
            "convection_scheme", "explicit donor-cell or limited linear-upwind required"
        )
    if domain.pressure_reference_cell is None and pressure_reference_Pa is not None:
        _error("pressure_reference_Pa", "outlet excludes additional gauge")
    reference = _scalar(
        0 if pressure_reference_Pa is None else pressure_reference_Pa,
        "pressure_reference_Pa",
    )
    a = _field(initial_nu_tilde_m2_s, n, vector=False, field="initial_nu_tilde_m2_s")
    if np.any(a < 0) or not np.any(a > 0):
        _error(
            "initial_nu_tilde_m2_s",
            "nonnegative field with an explicit positive turbulence start required",
        )
    initial = a.copy()
    rho, mu0 = domain.fluid.density_kg_m3, domain.fluid.viscosity_Pa_s
    nu, area_ref = mu0 / rho, domain.width_m * ls
    gates = np.array(
        [
            1e-8 * rho * us**2 / ls,
            1e-8 * rho * us / ls,
            1e-8 * us * nu / ls,
            1e-10 * rho * us * area_ref,
            1e-8 * rho * us**2 * area_ref,
            1e-6 * rho * us**3 * area_ref,
        ]
    )
    if not np.all(np.isfinite(gates)) or np.any(gates <= 0):
        _error("scales", "unrepresentable absolute joint gates")
    if np.any(domain.wall_quality.layer_growth > 1.2):
        _error("layer_growth", "wall-resolved SA requires layer growth <=1.2")
    force = _field(force_N_m3, n, vector=True, field="force_N_m3") + _driving(domain)
    scalar = _SATransportSystem(domain, None)
    distance = domain.wall_quality.wall_distance_m
    volumes = scalar.fv.cell_volumes_m3
    u, p = np.zeros(2 * n), np.zeros(n)
    history: list[tuple[float, ...]] = []
    inner_momentum_gates: list[float] = []
    steps = 0
    previous_linear: NDArray[np.float64] | None = None
    system = QuadSystem(
        domain,
        mu0 + rho * _eddy(a, nu),
        mu0 + rho * _eddy(scalar.boundary.apply(a), nu),
    )
    if domain.pressure_reference_cell is not None:
        system.reference = reference
    flux = system.alias_flux(system.velocity_flux.apply(u))
    zero_source = np.zeros(n)
    for _ in range(max_iterations):
        try:
            transport, _ = system.transport(flux, u, convection_scheme)
            momentum = Map(
                (system.diffusion.matrix + transport.matrix).tocsr(),
                system.diffusion.constant + transport.constant,
            )
            diagonal = momentum.matrix.diagonal().reshape(n, 2)
            if np.any(diagonal <= 0) or not np.all(np.isfinite(diagonal)):
                _error("momentum_diagonal_kg_s", "nonpositive/unresolved RANS diagonal")
            pressure_flux = system.pressure_flux(diagonal)
            rhs = (
                (force * volumes[:, None]).ravel()
                - momentum.constant
                - system.pressure_momentum.constant
            )
            # SA production depends on velocity gradients multiplied by a.
            # A momentum-only stopping scale can therefore stall the stricter
            # SA equation when a >> molecular nu. Allocate an inner error budget
            # with that amplitude ratio; all outer physical gates stay intact.
            inner_momentum_gate = (
                float(gates[0]) * 0.1 * float(nu / max(nu, float(a.max())))
            )
            if not np.isfinite(inner_momentum_gate) or inner_momentum_gate <= 0:
                _error("scales", "unrepresentable SA-aware inner momentum budget")
            inner_momentum_gates.append(inner_momentum_gate)
            prediction, pressure, count = _linear_simple(
                system,
                momentum,
                rhs,
                diagonal,
                pressure_flux,
                max_steps=max_simple_steps,
                momentum_gate=inner_momentum_gate,
                continuity_gate=float(gates[1]),
                initial_state=previous_linear,
            )
            previous_linear = np.concatenate((prediction, pressure))
            steps += count
            candidate_flux = system.alias_flux(
                system.velocity_flux.apply(prediction) + pressure_flux.apply(pressure)
            )
            grad_u = system.gu.apply(prediction).reshape(n, 2, 2)
            sa_operator, _, _, _ = _sa_operator(
                scalar,
                a,
                candidate_flux / rho,
                grad_u,
                distance,
                zero_source,
                convection_scheme,
            )
            candidate_a = np.asarray(
                splu(sa_operator.matrix.tocsc()).solve(-sa_operator.constant)
            )
            updated, alpha, _ = positive_sa_update(a, candidate_a, relaxation=aa)
            fraction = alpha / aa
            u = u + fraction * au * (prediction - u)
            p = p + fraction * ap * (pressure - p)
            a = updated
            mu = mu0 + rho * _eddy(a, nu)
            mub = mu0 + rho * _eddy(scalar.boundary.apply(a), nu)
            system = system._with_viscosity(mu, mub)
            if domain.pressure_reference_cell is not None:
                system.reference = reference
            flux, diag = _reconstruct(system, u, p, flux, force, convection_scheme)
            grad_u = system.gu.apply(u).reshape(n, 2, 2)
            final_sa, terms, _, _ = _sa_operator(
                scalar, a, flux / rho, grad_u, distance, zero_source, convection_scheme
            )
            sa_residual = final_sa.apply(a) / volumes
            row = (
                float(np.max(np.abs(diag["momentum"]))),
                float(np.max(np.abs(diag["mass"]))),
                float(np.max(np.abs(sa_residual))),
                abs(diag["global_mass"]),
                float(np.max(np.abs(diag["global_momentum"]))),
                abs(diag["power"]),
                alpha,
            )
            history.append(row)
            if not np.all(np.isfinite(row)):
                _error("solution", "nonfinite freshly reconstructed joint state")
            outlet = np.flatnonzero(system.kind == "pressure_outlet")
            bad = outlet[flux[outlet] < -1e-12 * rho * us * area_ref]
            if len(bad):
                _error(
                    f"outlet.face[{bad[0]}]",
                    f"unsupported backflow: {flux[bad[0]]:.9g} kg/s",
                )
            if (
                np.all(np.array(row[:6]) <= gates)
                and np.max(np.abs(diag["nonorthogonal"])) <= gates[0]
            ):
                ids = domain.wall_quality.boundary_indices
                lookup = {boundary: i for i, boundary in enumerate(system.walls)}
                traction = diag["traction"][[lookup[int(i)] for i in ids]]
                normals = domain.mesh.boundary_normals[ids]
                tangent = (
                    traction
                    - np.einsum("ij,ij->i", traction, normals)[:, None] * normals
                )
                y_plus = domain.wall_quality.y_plus(
                    density_kg_m3=rho,
                    viscosity_Pa_s=mu0,
                    wall_shear_Pa=np.linalg.norm(tangent, axis=1),
                )
                break
        except (ValidationError, RuntimeError, FloatingPointError) as error:
            raise RansConvergenceError(
                str(error), history, u, p + system.reference, a
            ) from error
    else:
        raise RansConvergenceError(
            f"joint iteration budget {max_iterations} exhausted",
            history,
            u,
            p + system.reference,
            a,
        )
    canonical = dict(
        mesh=domain.mesh.to_document(),
        width_m=domain.width_m,
        fluid=dict(density_kg_m3=rho, viscosity_Pa_s=mu0),
        boundaries=_json(domain.boundaries),
        pressure_reference_cell=domain.pressure_reference_cell,
        pressure_reference_Pa=system.reference,
        force_N_m3=force.tolist(),
        initial_nu_tilde_m2_s=initial.tolist(),
    )
    constants = _sa_constants()
    provenance = RunProvenance.from_input_bytes(
        json.dumps(canonical, sort_keys=True, allow_nan=False).encode(),
        solver_name="quad-fv-simple-sa-noft2",
        solver_parameters=dict(
            model="SA-noft2",
            geometry_reconstruction="shared_only_within_identical_domain",
            inner_momentum_gate_rule="0.1*outer_momentum_gate*nu/max(nu,max(a))",
            inner_momentum_gates_N_m3=inner_momentum_gates,
            s_modification="Allmaras-2012-TMR-note-1c",
            constants=constants,
            convection=convection_scheme,
            velocity_scale_m_s=us,
            length_scale_m=ls,
            velocity_relaxation=au,
            pressure_relaxation=ap,
            sa_relaxation=aa,
            positivity_attempts=20,
            max_iterations=max_iterations,
            max_simple_steps=max_simple_steps,
            max_unknowns=max_unknowns,
            absolute_gates=gates.tolist(),
            residual_history=history,
            wall_y_plus=y_plus.tolist(),
            numpy_version=np.__version__,
            scipy_version=scipy.__version__,
        ),
    )
    flow = _flow_result(system, u, p, flux, diag, history, steps, provenance)
    return RansResult(
        flow,
        _freeze_float(a),
        _freeze_float(terms.eddy_viscosity_m2_s),
        _freeze_float(mu),
        _freeze_float(mub),
        _freeze_float(sa_residual),
        float(sa_residual @ volumes),
        _freeze_float(y_plus),
        len(history),
        tuple(history),
        True,
        provenance,
    )
