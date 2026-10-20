"""Conservative stationary SA-noft2 scalar equation on the validated Quad-FV mesh."""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
import scipy  # type: ignore[import-untyped]
from numpy.typing import ArrayLike, NDArray
from scipy.sparse import diags  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import _float_array, _freeze_float
from fem_toolbox.cfd.quad_flow import _field, _json
from fem_toolbox.cfd.quad_fv import QuadFVGeometry, _error, _validated_domain
from fem_toolbox.cfd.quad_system import Map, affine_product, sparse
from fem_toolbox.cfd.sa_model import (
    SATerms,
    _sa_constants,
    positive_sa_update,
    sa_terms,
)
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


class SAConvergenceError(ValidationError):
    """Nonconvergence retains the last nonnegative state and raw residual history."""

    def __init__(
        self, history: list[float], state: NDArray[np.float64], reason: str
    ) -> None:
        self.residual_history = tuple(history)
        self.nu_tilde_m2_s = _freeze_float(state)
        super().__init__(
            f"SA did not converge: {reason}; last residual={history[-1] if history else 'not evaluated'} m²/s²; minimum={state.min():.9g} m²/s at cell {state.argmin()}",
            field="sa_iterations",
        )


@dataclass(frozen=True, slots=True)
class SATransportResult:
    """Immutable scalar solution for prescribed conservative flow; not a RANS solve.

    Arrays: a/residual/wall distance are cell fields, boundary a follows line2,
    and advective/diffusive fluxes follow mesh.faces (outward from owner).
    Integrated fluxes and global balance have unit m^5/s². History is the
    maximum absolute fresh cell residual in m²/s² at every accepted update.
    """

    nu_tilde_m2_s: NDArray[np.float64]
    boundary_nu_tilde_m2_s: NDArray[np.float64]
    wall_distance_m: NDArray[np.float64]
    residual_m2_s2: NDArray[np.float64]
    advective_flux_m5_s2: NDArray[np.float64]
    diffusive_flux_m5_s2: NDArray[np.float64]
    global_balance_m5_s2: float
    terms: SATerms
    iterations: int
    residual_history: tuple[float, ...]
    accepted_relaxations: tuple[float, ...]
    converged: bool
    provenance: RunProvenance


class _SATransportSystem:
    def __init__(
        self, domain: FluidDomain, prescribed: NDArray[np.float64] | None
    ) -> None:
        self.domain = domain
        self.fv = fv = QuadFVGeometry.from_domain(domain)
        self.mesh = mesh = domain.mesh
        self.n, self.nf, self.nb = len(mesh.quad4), len(mesh.faces), len(mesh.line2)
        owner, other = mesh.face_owners, fv.neighbours
        boundary_owner = owner[mesh.boundary_face_indices]
        kind = np.full(self.nb, "periodic", dtype="U20")
        constant = np.zeros(self.nb)
        for name, bc in domain.boundaries.items():
            ids = mesh.boundary_sets[name]
            kind[ids] = bc["type"]
            if bc["type"] == "inlet":
                constant[ids] = bc["nu_tilde_m2_s"]
        neumann = np.isin(kind, ["pressure_outlet", "symmetry", "periodic"])
        if prescribed is not None:
            constant, neumann = prescribed.copy(), np.zeros(self.nb, dtype=bool)
        rows = np.flatnonzero(neumann)
        self.boundary = Map(
            sparse(
                rows.tolist(),
                boundary_owner[rows].tolist(),
                [1.0] * len(rows),
                (self.nb, self.n),
            ),
            constant,
        )
        self.gx = Map(
            (fv._gx[:, : self.n] + fv._gx[:, self.n :] @ self.boundary.matrix).tocsr(),
            np.asarray(fv._gx[:, self.n :] @ constant),
        )
        self.gy = Map(
            (fv._gy[:, : self.n] + fv._gy[:, self.n :] @ self.boundary.matrix).tocsr(),
            np.asarray(fv._gy[:, self.n :] @ constant),
        )
        r, c, v = [], [], []
        for f in fv.active_faces:
            r.append(int(owner[f]))
            c.append(int(f))
            v.append(1.0)
            if other[f] >= 0:
                r.append(int(other[f]))
                c.append(int(f))
                v.append(-1.0)
        self.balance = sparse(r, c, v, (self.n, self.nf))
        face_owner = sparse(
            list(range(self.nf)), owner.tolist(), [1.0] * self.nf, (self.nf, self.n)
        )
        inner = np.flatnonzero(other >= 0)
        neighbour = sparse(
            inner.tolist(), other[inner].tolist(), [1.0] * len(inner), (self.nf, self.n)
        )
        outer = np.flatnonzero(other < 0)
        trace = sparse(
            outer.tolist(),
            fv._boundary_row[outer].tolist(),
            [1.0] * len(outer),
            (self.nf, self.nb),
        )
        difference = Map(
            (neighbour + trace @ self.boundary.matrix - face_owner).tocsr(),
            np.asarray(trace @ constant),
        )
        fractions = np.where(other >= 0, fv._fraction, 0)
        average = diags(1 - fractions) @ face_owner + diags(fractions) @ neighbour
        tangent = fv._direction / fv._distance[:, None] - mesh.face_normals
        matrix = (
            diags(1 / fv._distance) @ difference.matrix
            - diags(tangent[:, 0]) @ average @ self.gx.matrix
            - diags(tangent[:, 1]) @ average @ self.gy.matrix
        )
        constant_n = (
            difference.constant / fv._distance
            - tangent[:, 0] * (average @ self.gx.constant)
            - tangent[:, 1] * (average @ self.gy.constant)
        )
        mask = np.ones(self.nf)
        mask[mesh.boundary_face_indices[neumann]] = 0
        self.normal_gradient = Map((diags(mask) @ matrix).tocsr(), mask * constant_n)

    def gradient(self, a: NDArray[np.float64]) -> NDArray[np.float64]:
        return np.stack((self.gx.apply(a), self.gy.apply(a)), axis=1)

    def flux_maps(
        self,
        a: NDArray[np.float64],
        volume_flux: NDArray[np.float64],
        nu: float,
        scheme: str,
    ) -> tuple[Map, Map]:
        fv, mesh = self.fv, self.mesh
        owner, other = mesh.face_owners, fv.neighbours
        b = self.boundary.apply(a)
        fraction = np.where(other >= 0, fv._fraction, 0)
        face_a = (1 - fraction) * a[owner] + fraction * a[np.maximum(other, 0)]
        outer = np.flatnonzero(other < 0)
        face_a[outer] = b[fv._boundary_row[outer]]
        conductivity = (nu + face_a) / (2 / 3)
        diffusion = affine_product(
            diags(-conductivity * fv.face_areas_m2), self.normal_gradient
        )
        rows, cols, vals = [], [], []
        constants = np.zeros(self.nf)
        gradients = self.gradient(a)
        for f in range(self.nf):
            p, q = int(owner[f]), int(other[f])
            if q < 0:
                row = int(fv._boundary_row[f])
                # Both inward and outward physical traces obey their explicit BC.
                bc = self.boundary.matrix.getrow(row)
                rows.extend([f] * len(bc.indices))
                cols.extend(bc.indices.tolist())
                vals.extend(bc.data.tolist())
                constants[f] = self.boundary.constant[row]
                continue
            up = p if volume_flux[f] >= 0 else q
            rows.append(f)
            cols.append(up)
            vals.append(1.0)
            if scheme == "linear-upwind":
                offset = mesh.face_centres_m[f] - mesh.cell_centres_m[p]
                if up == q:
                    offset = offset - fv._direction[f]
                prediction = a[up] + float(gradients[up] @ offset)
                # A documented face-reconstruction limiter, never a state clamp.
                limited = min(max(prediction, min(a[p], a[q])), max(a[p], a[q]))
                constants[f] = limited - a[up]
        advective = affine_product(
            diags(volume_flux),
            Map(sparse(rows, cols, vals, (self.nf, self.n)), constants),
        )
        return advective, diffusion


def _sa_operator(
    system: _SATransportSystem,
    state: NDArray[np.float64],
    volume_flux: NDArray[np.float64],
    grad_u: NDArray[np.float64],
    distance: NDArray[np.float64],
    source: NDArray[np.float64],
    scheme: str,
) -> tuple[Map, SATerms, Map, Map]:
    """Rebuild every scalar coefficient/source on the supplied current state."""
    nu = system.domain.fluid.viscosity_Pa_s / system.domain.fluid.density_kg_m3
    terms = sa_terms(
        state,
        velocity_gradient_s_inv=grad_u,
        gradient_nu_tilde_m_s=system.gradient(state),
        wall_distance_m=distance,
        kinematic_viscosity_m2_s=nu,
    )
    adv, diff = system.flux_maps(state, volume_flux, nu, scheme)
    volumes = system.fv.cell_volumes_m3
    operator = Map(
        (
            system.balance @ (adv.matrix + diff.matrix)
            + diags(volumes * terms.destruction_rate_s_inv)
        ).tocsr(),
        np.asarray(system.balance @ (adv.constant + diff.constant))
        - volumes * (terms.production_m2_s2 + terms.gradient_source_m2_s2 + source),
    )
    return operator, terms, adv, diff


def solve_sa_transport(
    domain: FluidDomain,
    *,
    mass_flux_kg_s: ArrayLike,
    velocity_gradient_s_inv: ArrayLike,
    initial_nu_tilde_m2_s: ArrayLike,
    velocity_scale_m_s: float,
    length_scale_m: float,
    boundary_values_m2_s: ArrayLike | None = None,
    wall_distance_m: ArrayLike | None = None,
    source_m2_s2: ArrayLike = 0.0,
    convection_scheme: str = "linear-upwind",
    relaxation: float = 0.5,
    max_iterations: int = 500,
    max_cells: int = 50000,
) -> SATransportResult:
    """Solve SA for frozen conservative flow with full fresh nonlinear residuals.

    Boundary defaults: wall a=0, inlet explicit a, outlet/symmetry zero normal
    derivative and periodic paired trace. Optional boundary and distance fields
    are explicit manufactured-solution hooks, recorded as such in provenance;
    they do not qualify physical wall resolution or a coupled RANS solution.
    Destruction is a frozen nonnegative implicit sink. A common positivity line
    search halves the relaxed candidate, at most twenty attempts, without state
    clipping. Failure carries the latest history, state, minimum and cell ID.
    """
    domain = _validated_domain(domain)
    mesh = domain.mesh
    n, nb, nf = len(mesh.quad4), len(mesh.line2), len(mesh.faces)
    for name, value in (("max_iterations", max_iterations), ("max_cells", max_cells)):
        if type(value) is not int or value < 1:
            _error(name, "positive integer required")
    if n > max_cells:
        _error("max_cells", "SA cell resource budget exceeded")
    us, ls = (
        _scalar(velocity_scale_m_s, "velocity_scale_m_s"),
        _scalar(length_scale_m, "length_scale_m"),
    )
    relax = _scalar(relaxation, "sa_relaxation")
    if min(us, ls) <= 0 or not 0 < relax <= 1:
        _error("scales", "positive SI scales and relaxation in (0,1] required")
    if convection_scheme not in ("donor-cell", "linear-upwind"):
        _error("convection_scheme", "SA requires donor-cell or limited linear-upwind")
    a = _field(initial_nu_tilde_m2_s, n, vector=False, field="initial_nu_tilde_m2_s")
    if np.any(a < 0):
        _error("initial_nu_tilde_m2_s", "nonnegative cell values required")
    grad_u = _float_array(velocity_gradient_s_inv, "velocity_gradient_s_inv")
    if grad_u.shape != (n, 2, 2):
        _error("velocity_gradient_s_inv", "expected (ncell,2,2)")
    flux = _field(mass_flux_kg_s, nf, vector=False, field="mass_flux_kg_s")
    prescribed = (
        None
        if boundary_values_m2_s is None
        else _field(
            boundary_values_m2_s, nb, vector=False, field="boundary_values_m2_s"
        )
    )
    if prescribed is not None and np.any(prescribed < 0):
        _error("boundary_values_m2_s", "negative SA boundary value")
    distance = (
        domain.wall_quality.wall_distance_m
        if wall_distance_m is None
        else _field(wall_distance_m, n, vector=False, field="wall_distance_m")
    )
    if np.any(distance <= 0):
        _error("wall_distance_m", "positive distance required at every cell")
    source = _field(source_m2_s2, n, vector=False, field="source_m2_s2")
    nu = domain.fluid.viscosity_Pa_s / domain.fluid.density_kg_m3
    gate = 1e-8 * us * nu / ls
    if not np.isfinite(gate) or gate <= 0:
        _error("scales", "unrepresentable absolute SA gate")
    system = _SATransportSystem(domain, prescribed)
    fv = system.fv
    rho = domain.fluid.density_kg_m3
    mass_gate = 1e-8 * rho * us / ls
    if np.max(np.abs(fv.cell_balance(flux) / fv.cell_volumes_m3)) > mass_gate:
        _error("mass_flux_kg_s", "SA requires an already conservative flow field")
    for face, partner in fv.periodic_partners.items():
        if abs(flux[face] + flux[partner]) > 1e-12 * rho * us * domain.width_m * ls:
            _error("mass_flux_kg_s", "periodic face fluxes must be opposite")
    for name, bc in domain.boundaries.items():
        faces = mesh.boundary_face_indices[mesh.boundary_sets[name]]
        if (
            bc["type"] in ("wall", "symmetry")
            and np.max(np.abs(flux[faces])) > 1e-10 * rho * us * domain.width_m * ls
        ):
            _error("mass_flux_kg_s", f"nonzero normal flow through {name}")
    volume_flux = flux / rho
    volumes = fv.cell_volumes_m3
    history: list[float] = []
    alphas: list[float] = []

    def rebuild(state: NDArray[np.float64]) -> tuple[Map, SATerms, Map, Map]:
        return _sa_operator(
            system, state, volume_flux, grad_u, distance, source, convection_scheme
        )

    for _ in range(max_iterations):
        operator, _, _, _ = rebuild(a)
        try:
            candidate = np.asarray(
                splu(operator.matrix.tocsc()).solve(-operator.constant)
            )
        except (ValueError, RuntimeError) as error:
            raise SAConvergenceError(
                history, a, f"scalar factorization: {error}"
            ) from error
        try:
            a, alpha, _ = positive_sa_update(a, candidate, relaxation=relax)
        except ValidationError as error:
            raise SAConvergenceError(history, a, str(error)) from error
        final, terms, adv, diff = rebuild(a)
        residual = final.apply(a) / volumes
        maximum = float(np.max(np.abs(residual)))
        history.append(maximum)
        alphas.append(alpha)
        if not np.isfinite(maximum):
            raise SAConvergenceError(history, a, "nonfinite fresh residual")
        if maximum <= gate:
            break
    else:
        raise SAConvergenceError(
            history, a, f"iteration budget {max_iterations} exhausted"
        )
    canonical = dict(
        mesh=mesh.to_document(),
        width_m=domain.width_m,
        density_kg_m3=rho,
        viscosity_Pa_s=domain.fluid.viscosity_Pa_s,
        boundaries=_json(domain.boundaries),
        mass_flux_kg_s=flux.tolist(),
        velocity_gradient_s_inv=grad_u.tolist(),
        initial_nu_tilde_m2_s=np.asarray(initial_nu_tilde_m2_s).tolist(),
        boundary_values_m2_s=system.boundary.apply(a).tolist(),
        wall_distance_m=distance.tolist(),
        source_m2_s2=source.tolist(),
    )
    provenance = RunProvenance.from_input_bytes(
        json.dumps(canonical, sort_keys=True, allow_nan=False).encode(),
        solver_name="quad-fv-sa-noft2",
        solver_parameters=dict(
            model="SA-noft2",
            constants=_sa_constants(),
            s_modification="Allmaras-2012-TMR-note-1c",
            convection=convection_scheme,
            wall_distance_source="mesh_segments"
            if wall_distance_m is None
            else "explicit_verification_field",
            boundary_source="domain"
            if prescribed is None
            else "explicit_verification_field",
            velocity_scale_m_s=us,
            length_scale_m=ls,
            relaxation=relax,
            positivity_attempts=20,
            max_iterations=max_iterations,
            max_cells=max_cells,
            absolute_residual_gate_m2_s2=gate,
            residual_history=history,
            accepted_relaxations=alphas,
            numpy_version=np.__version__,
            scipy_version=scipy.__version__,
        ),
    )
    return SATransportResult(
        _freeze_float(a),
        _freeze_float(system.boundary.apply(a)),
        _freeze_float(distance),
        _freeze_float(residual),
        _freeze_float(adv.apply(a)),
        _freeze_float(diff.apply(a)),
        float(residual @ volumes),
        terms,
        len(history),
        tuple(history),
        tuple(alphas),
        True,
        provenance,
    )
