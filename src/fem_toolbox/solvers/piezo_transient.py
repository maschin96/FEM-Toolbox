"""Undamped average-acceleration Newmark reference for piezoelectric DAEs."""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.linear_backend import LinearBackend, ScipyBackend, solve_linear
from fem_toolbox.solvers.piezo_dynamics import _build_dynamic_system, _readonly
from fem_toolbox.time import integration_times


@dataclass(frozen=True)
class KinematicValues:
    """Prescribed arrays in active-node order: u [m], v [m/s], a [m/s²], phi [V].

    Mechanical arrays have shape (nactive,3), potential (nactive,). Entries
    outside the corresponding fixed DOFs must be zero. The caller supplies
    consistent derivatives; finite jumps in prescribed u or v are unsupported.
    """

    displacement: ArrayLike
    velocity: ArrayLike
    acceleration: ArrayLike
    potential: ArrayLike


@dataclass(frozen=True, eq=False)
class TransientPiezoelectricResult:
    """Owned time histories, SI units, active-node order, no damping.

    Mechanical fields and original residuals: (ntime,nactive,3); potential
    and electric residual: (ntime,nactive). Residual_u includes inertia and
    contains support reactions. Minus residual_phi gives conductor charge
    contributions; floating electrodes require summation over their nodes.
    Total energy [J] is .5*v'Mv + .5*u'Kuu*u + .5*phi'D*phi, conserved for
    undriven homogeneous boundaries. Residual metrics use condensed electrodes.
    """

    node_ids: NDArray[np.int64]
    times: NDArray[np.float64]
    displacement: NDArray[np.float64]
    velocity: NDArray[np.float64]
    acceleration: NDArray[np.float64]
    potential: NDArray[np.float64]
    residual_u: NDArray[np.float64]
    residual_phi: NDArray[np.float64]
    total_energy: NDArray[np.float64]
    relative_residual_u: NDArray[np.float64]
    relative_residual_phi: NDArray[np.float64]
    provenance: RunProvenance

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "node_ids", np.frombuffer(self.node_ids.tobytes(), dtype=np.int64)
        )
        for name in (
            "times",
            "displacement",
            "velocity",
            "acceleration",
            "potential",
            "residual_u",
            "residual_phi",
            "total_energy",
            "relative_residual_u",
            "relative_residual_phi",
        ):
            object.__setattr__(self, name, _readonly(getattr(self, name)))


def solve_transient_piezoelectric(
    mesh: Mesh,
    material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial],
    density: ArrayLike,
    conditions: Iterable[Dirichlet],
    times: ArrayLike,
    *,
    backend: LinearBackend | None = None,
    electrodes: Sequence[Electrode] | None = None,
    initial_displacement: ArrayLike | None = None,
    initial_velocity: ArrayLike | None = None,
    loads: Callable[[float], ArrayLike] | None = None,
    prescribed: Callable[[float], KinematicValues] | None = None,
    events: ArrayLike = (),
) -> TransientPiezoelectricResult:
    """Integrate with beta=1/4, gamma=1/2, no physical or numerical damping.

    ``backend`` is an experimental linear-backend hook; None selects SciPy.
    Reduction, physical residuals and result units remain this solver's contract.

    Topology uses homogeneous Dirichlet/electrode declarations as in modes.
    loads(t) returns integrated signed [force (N), electric RHS (C)] of shape
    (4*nactive,); for total floating conductor charge Q put -Q on one of its
    potential DOFs. Contributions on tied nodes sum. No implicit volume/area
    integration is performed. prescribed(t) supplies KinematicValues, including
    voltage values on fixed potential DOFs. Defaults are all zero.

    Initial u/v have shape (nactive,3) and must match supports to 1e-12 SI;
    omitted arrays equal the prescribed state. Electric potential and mechanical
    acceleration are solved consistently at t0, never independent initial data.
    Every load/voltage jump or measured-curve knot must be listed in events [s];
    interior events are inserted. Integration uses the left limit at an interval
    end; stored potential/acceleration are the right-continuous state. u and v
    remain continuous. Step/Ramp/TimeSeries may be used inside load callbacks;
    pass their knots as events. Smooth callbacks are sampled at the grid only.

    At most 2048 active DOFs and 2,000,000 stored scalar state entries; arrays
    are retained in memory. Invalid data and failed residuals raise ValidationError.
    """
    if mesh.n_tetra10:
        raise ValidationError(
            "Tet10 transient/harmonic piezo is not yet supported", field="mesh"
        )
    system = _build_dynamic_system(
        mesh, material, density, conditions, electrodes, backend
    )
    grid = integration_times(times)
    raw_events = np.asarray(events)
    ev = real_array(events, raw_events.shape)
    if ev.ndim != 1:
        raise ValidationError("events must be a vector", field="events")
    grid = np.unique(np.concatenate((grid, ev[(ev > grid[0]) & (ev < grid[-1])])))
    n = len(system.active)
    nu = 3 * n
    size = 4 * n
    if len(grid) * size > 2_000_000:
        raise ValidationError(
            "history exceeds 2000000 scalar state entries", field="times"
        )
    fu, fp = system.free_u, system.free_phi
    fixed_u = system.fixed[system.fixed < nu]
    fixed_original = np.flatnonzero(np.isin(system.index, system.fixed))
    fixed_mask = np.zeros(size, dtype=bool)
    fixed_mask[fixed_original] = True
    dielectric = -system.reduced[np.ix_(fp, fp)]
    if len(fp):
        electric_inverse = solve_linear(
            dielectric,
            np.eye(len(fp)),
            fields={"phi": np.arange(len(fp))},
            scale=1 / np.sqrt(np.diag(dielectric)),
            kind="positive_definite",
            backend=backend,
        ).values
    else:
        electric_inverse = np.empty((0, 0))
    selected_backend = ScipyBackend() if backend is None else backend
    mass_factor = selected_backend.factorize(
        system.free_mass, fields={"u": np.arange(len(fu))}, kind="positive_definite"
    )
    _, representatives = np.unique(system.index, return_index=True)
    samples = []

    def evaluate(
        t: float,
    ) -> tuple[
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
    ]:
        f = np.zeros(size) if loads is None else real_array(loads(t), (size,))
        g = np.zeros(size)
        v = np.zeros(nu)
        a = np.zeros(nu)
        if prescribed is not None:
            state = prescribed(t)
            if not isinstance(state, KinematicValues):
                raise ValidationError(
                    "prescribed must return KinematicValues", field="prescribed"
                )
            g[:nu] = real_array(state.displacement, (n, 3)).ravel()
            g[nu:] = real_array(state.potential, (n,))
            v = real_array(state.velocity, (n, 3)).ravel()
            a = real_array(state.acceleration, (n, 3)).ravel()
        if (
            np.any(g[~fixed_mask])
            or np.any(v[~fixed_mask[:nu]])
            or np.any(a[~fixed_mask[:nu]])
        ):
            raise ValidationError(
                "prescribed data outside fixed DOFs", field="prescribed"
            )
        reduced_g = g[representatives]
        rhs = (
            np.bincount(system.index, weights=f, minlength=len(system.reduced))
            - system.reduced @ reduced_g
        )
        effective = (
            rhs[fu]
            + system.electric_map.T @ rhs[fp]
            - system.mass[np.ix_(fu, fixed_u)] @ a[fixed_u]
        )
        samples.append([t, f.tolist(), g.tolist(), v.tolist(), a.tolist()])
        return f, g, v, a, rhs, effective

    f, g, gv, ga, rhs, effective = evaluate(float(grid[0]))
    u = (
        g[:nu].copy()
        if initial_displacement is None
        else real_array(initial_displacement, (n, 3)).ravel()
    )
    v = (
        gv.copy()
        if initial_velocity is None
        else real_array(initial_velocity, (n, 3)).ravel()
    )
    if np.any(np.abs(u[fixed_u] - g[fixed_u]) > 1e-12) or np.any(
        np.abs(v[fixed_u] - gv[fixed_u]) > 1e-12
    ):
        raise ValidationError(
            "inconsistent initial support displacement/velocity", field="initial"
        )
    u[fixed_u] = g[fixed_u]
    v[fixed_u] = gv[fixed_u]
    q = u[fu].copy()
    w = v[fu].copy()
    acc = mass_factor.solve(effective - system.stiffness @ q).values
    displacements = []
    velocities = []
    accelerations = []
    potentials = []
    residuals = []
    energies = []
    relative_u = []
    relative_phi = []
    for step, t in enumerate(grid):
        if step:
            dt = float(t - grid[step - 1])
            left = evaluate(float(np.nextafter(t, grid[step - 1])))
            right = evaluate(float(t))
            for before, after in ((left[1][:nu], right[1][:nu]), (left[2], right[2])):
                if np.any(
                    np.abs(before - after)
                    > 1e-12 + 1e-10 * np.maximum(np.abs(before), np.abs(after))
                ):
                    raise ValidationError(
                        "discontinuous prescribed displacement/velocity",
                        field="prescribed",
                    )
            predicted = q + dt * w + 0.25 * dt * dt * acc
            new_acc = solve_linear(
                system.free_mass + 0.25 * dt * dt * system.stiffness,
                left[5] - system.stiffness @ predicted,
                fields={"u": np.arange(len(fu))},
                kind="positive_definite",
                backend=backend,
            ).values
            q = predicted + 0.25 * dt * dt * new_acc
            w = w + 0.5 * dt * (acc + new_acc)
            f, g, gv, ga, rhs, effective = right
            acc = mass_factor.solve(effective - system.stiffness @ q).values
        reduced_state = np.zeros(len(system.reduced))
        reduced_state[system.index[fixed_original]] = g[fixed_original]
        reduced_state[fu] = q
        reduced_state[fp] = system.electric_map @ q - electric_inverse @ rhs[fp]
        state = reduced_state[system.index]
        velocity = gv.copy()
        velocity[fu] = w
        acceleration = ga.copy()
        acceleration[fu] = acc
        inertia = system.mass @ acceleration
        residual = system.matrix @ state - f
        residual[:nu] += inertia
        tied = np.bincount(
            system.index, weights=residual, minlength=len(system.reduced)
        )
        denominator = np.abs(system.matrix) @ np.abs(state) + np.abs(f)
        denominator[:nu] += np.abs(system.mass) @ np.abs(acceleration)
        den = np.bincount(
            system.index, weights=denominator, minlength=len(system.reduced)
        )
        relative = []
        for free in (fu, fp):
            numerator = float(np.linalg.norm(tied[free]))
            norm = float(np.linalg.norm(den[free]))
            relative.append(
                numerator / norm if norm else (0.0 if numerator == 0 else np.inf)
            )
        if (
            not np.isfinite(state).all()
            or not np.isfinite(acceleration).all()
            or not np.isfinite(relative).all()
            or max(relative) > 1e-9
        ):
            raise ValidationError(
                "nonfinite state or dynamic field residual exceeds 1e-9", field="solver"
            )
        energy = 0.5 * (
            velocity @ system.mass @ velocity
            + state[:nu] @ system.matrix[:nu, :nu] @ state[:nu]
            - state[nu:] @ system.matrix[nu:, nu:] @ state[nu:]
        )
        if not np.isfinite(energy):
            raise ValidationError("unrepresentable total energy", field="solver")
        displacements.append(state[:nu].reshape(n, 3))
        potentials.append(state[nu:])
        velocities.append(velocity.reshape(n, 3))
        accelerations.append(acceleration.reshape(n, 3))
        residuals.append(residual)
        energies.append(energy)
        relative_u.append(relative[0])
        relative_phi.append(relative[1])
    provenance = RunProvenance.from_input_bytes(
        system.canonical
        + json.dumps(
            {
                "times": grid.tolist(),
                "samples": samples,
                "initial_u": u.tolist(),
                "initial_v": v.tolist(),
            },
            sort_keys=True,
            allow_nan=False,
        ).encode(),
        solver_name="fem_toolbox.solvers.piezo_transient.solve_transient_piezoelectric",
        solver_parameters={
            "linear_backend": "scipy" if backend is None else backend.name,
            "beta": 0.25,
            "gamma": 0.5,
            "physical_damping": 0,
            "numerical_damping": 0,
            "electrical_mass": 0,
            "event_endpoint": "left integrate; right store",
        },
    )
    residual_array = np.asarray(residuals)
    return TransientPiezoelectricResult(
        system.active,
        grid,
        np.asarray(displacements),
        np.asarray(velocities),
        np.asarray(accelerations),
        np.asarray(potentials),
        residual_array[:, :nu].reshape(-1, n, 3),
        residual_array[:, nu:],
        np.asarray(energies),
        np.asarray(relative_u),
        np.asarray(relative_phi),
        provenance,
    )
