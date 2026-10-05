"""Complex, undamped piezoelectric harmonic reference with exp(+i omega t)."""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.assembly import real_array
from femantis.boundary import Dirichlet, Electrode
from femantis.errors import ValidationError
from femantis.materials import PiezoelectricMaterial
from femantis.mesh import Mesh
from femantis.provenance import RunProvenance
from femantis.solvers.linear_backend import LinearBackend, solve_linear
from femantis.solvers.piezo_dynamics import _build_dynamic_system


def _complex_array(value: ArrayLike, shape: tuple[int, ...]) -> NDArray[np.complex128]:
    a = np.asarray(value)
    if a.shape != shape or a.dtype.kind not in "iufc" or not np.isfinite(a).all():
        raise ValidationError(
            f"expected finite complex array of shape {shape}", field="values"
        )
    return np.array(a, dtype=np.complex128)


@dataclass(frozen=True, eq=False)
class HarmonicPiezoelectricResult:
    """Complex peak phasors (not RMS) with x(t)=Re(x_hat exp(+i omega t)).

    Nodal u/v/a are (nactive,3) in m, m/s, m/s²; potential (nactive,) in V.
    Original residuals contain mechanical support reactions [N] and signed
    electric residuals [C]; conductor charge is their negative nodal sum.
    Hz and rad/s are explicit scalar fields. Arrays are immutable snapshots.
    """

    node_ids: NDArray[np.int64]
    frequency_hz: float
    angular_frequency: float
    displacement: NDArray[np.complex128]
    velocity: NDArray[np.complex128]
    acceleration: NDArray[np.complex128]
    potential: NDArray[np.complex128]
    residual_u: NDArray[np.complex128]
    residual_phi: NDArray[np.complex128]
    relative_residual_u: float
    relative_residual_phi: float
    provenance: RunProvenance

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "node_ids", np.frombuffer(self.node_ids.tobytes(), dtype=np.int64)
        )
        for name in (
            "displacement",
            "velocity",
            "acceleration",
            "potential",
            "residual_u",
            "residual_phi",
        ):
            a = getattr(self, name)
            object.__setattr__(
                self,
                name,
                np.frombuffer(a.tobytes(), dtype=np.complex128).reshape(a.shape),
            )


def solve_harmonic_piezoelectric(
    mesh: Mesh,
    material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial],
    density: ArrayLike,
    conditions: Iterable[Dirichlet],
    frequency_hz: float,
    *,
    backend: LinearBackend | None = None,
    electrodes: Sequence[Electrode] | None = None,
    loads: ArrayLike | None = None,
    prescribed: ArrayLike | None = None,
) -> HarmonicPiezoelectricResult:
    """Solve one undamped harmonic state, including DC (frequency_hz=0).

    ``backend`` is an experimental linear-backend hook; None selects SciPy.
    Reduction, physical residuals and result units remain this solver's contract.

    Straight Tet4/Tet10 cells carry consistent mechanical mass, no electric mass.
    Topology uses homogeneous conditions/electrodes; loads and prescribed
    arrays are finite real/complex (4*nactive,) in [u,phi] DOF order. Loads
    are integrated force [N] and signed electric RHS [C], prescribed values
    are displacement [m] and potential [V] and zero outside fixed DOFs.
    Floating electrode charge Q is a -Q load on one tied potential DOF.
    Support acceleration is -omega²*u_c and its inertia is included.

    Additive complex infrastructure preserves LinearSystem's real contract.
    The dynamic Schur matrix is K_eff-omega² M; no damping is added to hide
    resonance. A diagonal energy scaling gives a dimensionless symmetric
    matrix; |lambda_min|<=1e-10 is rejected as resonance/ill-conditioning,
    even with zero forcing. No infinite response is returned. The dense
    2048-total-DOF bound and modal gauge/support contract apply.
    """
    hz = float(real_array(frequency_hz, ()))
    if hz < 0:
        raise ValidationError("frequency_hz must be nonnegative", field="frequency_hz")
    omega = 2 * np.pi * hz
    if not np.isfinite(omega) or omega > np.sqrt(np.finfo(np.float64).max):
        raise ValidationError("frequency is not representable", field="frequency_hz")
    system = _build_dynamic_system(
        mesh, material, density, conditions, electrodes, backend
    )
    n = len(system.active)
    nu = 3 * n
    size = 4 * n
    f = (
        np.zeros(size, dtype=np.complex128)
        if loads is None
        else _complex_array(loads, (size,))
    )
    g = (
        np.zeros(size, dtype=np.complex128)
        if prescribed is None
        else _complex_array(prescribed, (size,))
    )
    fixed_mask = np.isin(system.index, system.fixed)
    if np.any(g[~fixed_mask]):
        raise ValidationError(
            "prescribed values outside fixed DOFs", field="prescribed"
        )
    _, representatives = np.unique(system.index, return_index=True)
    reduced_g = g[representatives]
    rhs = np.zeros(len(system.reduced), dtype=np.complex128)
    np.add.at(rhs, system.index, f)
    rhs -= system.reduced @ reduced_g
    fu, fp = system.free_u, system.free_phi
    effective = (
        rhs[fu]
        + system.electric_map.T @ rhs[fp]
        + omega**2 * (system.mass @ g[:nu])[fu]
    )
    dynamic = system.stiffness - omega**2 * system.free_mass
    scale = 1 / np.sqrt(
        np.diag(system.stiffness) + omega**2 * np.diag(system.free_mass)
    )
    scaled = dynamic * scale[:, None] * scale[None, :]
    if (
        not np.isfinite(scaled).all()
        or np.min(np.abs(np.linalg.eigvalsh(scaled))) <= 1e-10
    ):
        raise ValidationError(
            "resonance or ill-conditioned harmonic system", field="solver"
        )
    q = solve_linear(
        dynamic,
        effective,
        fields={"u": np.arange(len(fu))},
        scale=scale,
        backend=backend,
    ).values
    reduced_state = reduced_g.copy()
    reduced_state[fu] = q
    if len(fp):
        dielectric = system.dielectric
        electric = solve_linear(
            dielectric,
            rhs[fp],
            fields={"phi": np.arange(len(fp))},
            scale=1 / np.sqrt(np.diag(dielectric)),
            kind="positive_definite",
            backend=backend,
        ).values
        reduced_state[fp] = system.electric_map @ q - electric
    state = reduced_state[system.index]
    velocity = 1j * omega * state[:nu]
    acceleration = -(omega**2) * state[:nu]
    residual = system.matrix @ state - f
    residual[:nu] += system.mass @ acceleration
    tied = np.zeros(len(system.reduced), dtype=np.complex128)
    np.add.at(tied, system.index, residual)
    denominator = np.abs(system.matrix) @ np.abs(state) + np.abs(f)
    denominator[:nu] += np.abs(system.mass) @ np.abs(acceleration)
    den = np.bincount(system.index, weights=denominator, minlength=len(system.reduced))
    relative = []
    for free in (fu, fp):
        numerator = float(np.linalg.norm(tied[free]))
        norm = float(np.linalg.norm(den[free]))
        relative.append(
            numerator / norm if norm else (0.0 if numerator == 0 else np.inf)
        )
    if (
        not np.isfinite(state).all()
        or not np.isfinite(relative).all()
        or max(relative) > 1e-9
    ):
        raise ValidationError(
            "nonfinite state or harmonic field residual exceeds 1e-9", field="solver"
        )
    provenance = RunProvenance.from_input_bytes(
        system.canonical
        + json.dumps(
            {
                "frequency_hz": hz,
                "loads_real": f.real.tolist(),
                "loads_imag": f.imag.tolist(),
                "prescribed_real": g.real.tolist(),
                "prescribed_imag": g.imag.tolist(),
            },
            sort_keys=True,
            allow_nan=False,
        ).encode(),
        solver_name="fem_toolbox.solvers.piezo_harmonic.solve_harmonic_piezoelectric",
        solver_parameters={
            "linear_backend": "scipy" if backend is None else backend.name,
            "phasor": "Re(x_hat exp(+i omega t))",
            "frequency_hz": hz,
            "angular_frequency": omega,
            "damping": 0,
            "electrical_mass": 0,
            "scaled_resonance_tolerance": 1e-10,
        },
    )
    return HarmonicPiezoelectricResult(
        system.active,
        hz,
        omega,
        state[:nu].reshape(n, 3),
        velocity.reshape(n, 3),
        acceleration.reshape(n, 3),
        state[nu:],
        residual[:nu].reshape(n, 3),
        residual[nu:],
        relative[0],
        relative[1],
        provenance,
    )
