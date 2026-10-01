"""Index-1 RLC backward Euler with consistent initialization and energy audit."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from types import MappingProxyType

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.circuits.dc import (
    _assemble,
    _linear_solve,
    _numerical,
    _preflight,
    _tolerance,
    solve_dc,
)
from fem_toolbox.circuits.mna import provenance, source_names, source_rhs, storage
from fem_toolbox.circuits.model import Circuit, Component, _finite_scalar
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.time import (
    Ramp,
    Step,
    TimeSeries,
    TimeSignal,
    integration_times,
    signal_value,
)


def _snapshot(a: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.frombuffer(a.tobytes(), dtype=np.float64).reshape(a.shape)


@dataclass(frozen=True, slots=True)
class TransientResult:
    """Owned full-history SI arrays [ntime], keyed by normalized names.

    Q_C in C, flux_L in Wb. Energy/work/loss/damping/balance arrays use J;
    initial per-step entries are zero. KCL and instantaneous power balance
    are evaluated at all frames including the consistent initial frame.
    Events use left source values. There is no implicit adaptive time control.
    """

    times: NDArray[np.float64]
    node_potentials_V: Mapping[str, NDArray[np.float64]]
    branch_voltages_V: Mapping[str, NDArray[np.float64]]
    branch_currents_A: Mapping[str, NDArray[np.float64]]
    branch_powers_W: Mapping[str, NDArray[np.float64]]
    kcl_residual_A: Mapping[str, NDArray[np.float64]]
    capacitor_charge_C: Mapping[str, NDArray[np.float64]]
    inductor_flux_Wb: Mapping[str, NDArray[np.float64]]
    stored_energy_J: NDArray[np.float64]
    source_work_J: NDArray[np.float64]
    resistor_loss_J: NDArray[np.float64]
    numerical_dissipation_J: NDArray[np.float64]
    energy_balance_J: NDArray[np.float64]
    power_balance_W: NDArray[np.float64]
    relative_residual: NDArray[np.float64]
    provenance: RunProvenance

    def __post_init__(self) -> None:
        for name in (
            "times",
            "stored_energy_J",
            "source_work_J",
            "resistor_loss_J",
            "numerical_dissipation_J",
            "energy_balance_J",
            "power_balance_W",
            "relative_residual",
        ):
            object.__setattr__(self, name, _snapshot(getattr(self, name)))
        for name in (
            "node_potentials_V",
            "branch_voltages_V",
            "branch_currents_A",
            "branch_powers_W",
            "kcl_residual_A",
            "capacitor_charge_C",
            "inductor_flux_Wb",
        ):
            object.__setattr__(
                self,
                name,
                MappingProxyType(
                    {
                        key: _snapshot(value)
                        for key, value in getattr(self, name).items()
                    }
                ),
            )


def _initial_values(
    circuit: Circuit, kind: str, specified: Mapping[str, float] | None
) -> dict[str, float]:
    values = {c.name: 0.0 for c in circuit.components if c.kind == kind}
    seen = set()
    for name, value in (specified or {}).items():
        if (
            not isinstance(name, str)
            or name.upper() not in values
            or name.upper() in seen
        ):
            raise ValidationError(
                "unknown or duplicate initial state", field="initial_state"
            )
        seen.add(name.upper())
        values[name.upper()] = _finite_scalar(value, field="initial_state")
    return values


def _profiles(
    circuit: Circuit, sources: Mapping[str, TimeSignal] | None
) -> dict[str, TimeSignal]:
    result: dict[str, TimeSignal] = {}
    for name, value in source_names(circuit, sources).items():
        if isinstance(value, TimeSeries):
            if value.unit != ("V" if name[0] == "V" else "A"):
                raise ValidationError(
                    "profile must use source SI unit V/A", field="sources"
                )
        elif isinstance(value, Ramp):
            if not math.isfinite(value.end - value.start):
                raise ValidationError(
                    "ramp time span not representable", field="sources"
                )
        elif not isinstance(value, Step):
            value = _finite_scalar(value, field="sources")
        result[name] = value
    return result


def _signal_document(signal: TimeSignal) -> object:
    if isinstance(signal, TimeSeries):
        return {
            "kind": "series",
            "times": signal.times.tolist(),
            "values": signal.values.tolist(),
            "unit": signal.unit,
            "interpolation": signal.interpolation,
            "extrapolation": signal.extrapolation,
        }
    if isinstance(signal, (Step, Ramp)):
        return {"kind": type(signal).__name__, **asdict(signal)}
    return signal


def _initial_state(
    circuit: Circuit,
    sources: dict[str, float],
    cv: dict[str, float],
    li: dict[str, float],
    nodes: dict[str, int],
    branches: dict[str, int],
) -> tuple[NDArray[np.float64], dict[str, float], float]:
    proxy, names = [], {}
    for index, c in enumerate(circuit.components):
        kind = "V" if c.kind == "C" else "I" if c.kind == "L" else c.kind
        name = f"{kind}{index}"
        value = (
            cv[c.name]
            if c.kind == "C"
            else li[c.name]
            if c.kind == "L"
            else sources[c.name]
            if c.kind in "VI"
            else c.value
        )
        proxy.append(Component(name, c.positive, c.negative, value))
        names[c.name] = name
    try:
        initial = solve_dc(Circuit(proxy))
    except ValidationError as error:
        raise ValidationError(
            f"inconsistent or unsupported initial state: {error}", field="initial_state"
        ) from error
    x = np.zeros(len(nodes) + len(branches))
    for n, index in nodes.items():
        x[index] = initial.node_potentials_V[n]
    currents = {
        c.name: initial.branch_currents_A[names[c.name]] for c in circuit.components
    }
    for name, index in branches.items():
        x[index] = currents[name]
    return x, currents, initial.relative_residual


def _limit(value: int, name: str) -> int:
    if type(value) is not int or value < 1:
        raise ValidationError("expected positive integer limit", field=name)
    return value


def solve_transient(
    circuit: Circuit,
    times: ArrayLike,
    *,
    sources: Mapping[str, TimeSignal] | None = None,
    capacitor_voltages_V: Mapping[str, float] | None = None,
    inductor_currents_A: Mapping[str, float] | None = None,
    relative_tolerance: float = 1e-10,
    current_tolerance_A: float = 1e-12,
    power_tolerance_W: float = 1e-12,
    energy_tolerance_J: float = 1e-12,
    max_steps: int = 100_000,
    max_history_bytes: int = 256 * 1024 * 1024,
) -> TransientResult:
    """Integrate regular index-1 RLC networks by first-order backward Euler.

    Missing capacitor voltages and inductor currents mean zero initial state.
    Profile/initial names are case insensitive; unsupported ideal constraint
    loops/cutsets fail during consistent initialization. Knots enrich the grid,
    with left event values and no implicit prehistory. Retains every frame;
    max_history_bytes limits numeric result payload, not process overhead.
    KCL, branch power, charge/flux and BE energy identities are audited. Failed
    steps raise ValidationError with step index/time rather than partial success.
    """
    rtol = _tolerance(relative_tolerance, "relative_tolerance")
    atol = _tolerance(current_tolerance_A, "current_tolerance_A")
    ptol = _tolerance(power_tolerance_W, "power_tolerance_W")
    etol = _tolerance(energy_tolerance_J, "energy_tolerance_J")
    steps = _limit(max_steps, "max_steps")
    memory = _limit(max_history_bytes, "max_history_bytes")
    profiles = _profiles(circuit, sources)
    grid = integration_times(times, profiles.values())
    count = len(grid)
    if count - 1 > steps:
        raise ValidationError("event-enriched grid exceeds max_steps", field="times")
    payload = (
        count
        * 8
        * (
            2 * len(circuit.nodes)
            + 3 * len(circuit.components)
            + sum(c.kind in "CL" for c in circuit.components)
            + 8
        )
    )
    if payload > memory:
        raise ValidationError(
            "retained numeric history exceeds byte limit", field="history"
        )
    cv = _initial_values(circuit, "C", capacitor_voltages_V)
    li = _initial_values(circuit, "L", inductor_currents_A)
    _preflight(circuit, capacitors=True, ideal_inductors=False)
    a, _, nodes, branches = _assemble(circuit)
    e = storage(circuit, nodes, branches)
    values = {
        name: _finite_scalar(
            signal_value(p, float(grid[0]), side="left"), field="sources"
        )
        for name, p in profiles.items()
    }
    x, initial_currents, initial_relative = _initial_state(
        circuit, values, cv, li, nodes, branches
    )
    potentials = {n: np.zeros(count) for n in circuit.nodes}
    voltages = {c.name: np.zeros(count) for c in circuit.components}
    currents = {c.name: np.zeros(count) for c in circuit.components}
    powers = {c.name: np.zeros(count) for c in circuit.components}
    kcl = {n: np.zeros(count) for n in circuit.nodes}
    charge = {name: np.zeros(count) for name in cv}
    flux = {name: np.zeros(count) for name in li}
    energy, work, loss, damping, balance, power_balance, relative = (
        np.zeros(count) for _ in range(7)
    )
    relative[0] = initial_relative
    for index, time in enumerate(grid):
        dt = float(time - grid[index - 1]) if index else 0.0
        try:
            values = {
                name: _finite_scalar(
                    signal_value(p, float(time), side="left"), field="sources"
                )
                for name, p in profiles.items()
            }
            if index:
                with np.errstate(over="ignore", invalid="ignore"):
                    operator = (a + e / dt).tocsc()
                    rhs = (
                        source_rhs(
                            circuit,
                            nodes,
                            branches,
                            {n: complex(v) for n, v in values.items()},
                        ).real
                        + e @ x / dt
                    )
                if not np.isfinite(operator.data).all() or not np.isfinite(rhs).all():
                    raise _numerical("nonfinite transient operator or source")
                x, relative[index] = _linear_solve(operator, rhs)
                if relative[index] > rtol:
                    raise _numerical("transient backward error exceeds tolerance")
            for n, row in nodes.items():
                potentials[n][index] = x[row]
            outgoing: dict[str, list[float]] = {n: [] for n in circuit.nodes}
            for c in circuit.components:
                v = float(potentials[c.positive][index] - potentials[c.negative][index])
                if not index:
                    i = initial_currents[c.name]
                elif c.kind == "R":
                    i = v / c.value
                elif c.kind == "C":
                    i = c.value * (v - voltages[c.name][index - 1]) / dt
                elif c.kind == "I":
                    i = values[c.name]
                else:
                    i = float(x[branches[c.name]])
                (
                    voltages[c.name][index],
                    currents[c.name][index],
                    powers[c.name][index],
                ) = v, i, v * i
                outgoing[c.positive].append(i)
                outgoing[c.negative].append(-i)
                if c.kind == "C":
                    charge[c.name][index] = c.value * v
                    energy[index] += 0.5 * c.value * v * v
                    if index:
                        dv = v - voltages[c.name][index - 1]
                        damping[index] += 0.5 * c.value * dv * dv
                elif c.kind == "L":
                    flux[c.name][index] = c.value * i
                    energy[index] += 0.5 * c.value * i * i
                    if index:
                        di = i - currents[c.name][index - 1]
                        damping[index] += 0.5 * c.value * di * di
                        defect = v - c.value * di / dt
                        scale = abs(v) + abs(c.value * di / dt)
                        if abs(defect) > 1e-12 + rtol * scale:
                            raise _numerical(
                                f"inductor flux balance failed at {c.name}"
                            )
                if index:
                    if c.kind in "VI":
                        work[index] -= dt * v * i
                    elif c.kind == "R":
                        loss[index] += dt * v * i
            for n, terms in outgoing.items():
                kcl[n][index] = math.fsum(terms)
                if abs(kcl[n][index]) > atol + rtol * math.fsum(abs(i) for i in terms):
                    raise _numerical(f"Kirchhoff balance failed at node {n}")
            pvals = [float(p[index]) for p in powers.values()]
            power_balance[index] = math.fsum(pvals)
            if abs(power_balance[index]) > ptol + rtol * math.fsum(
                abs(p) for p in pvals
            ):
                raise _numerical("instantaneous power balance failed")
            if index:
                balance[index] = (
                    energy[index]
                    - energy[index - 1]
                    - work[index]
                    + loss[index]
                    + damping[index]
                )
                scale = (
                    abs(energy[index])
                    + abs(energy[index - 1])
                    + abs(work[index])
                    + abs(loss[index])
                    + abs(damping[index])
                )
                if abs(balance[index]) > etol + rtol * scale:
                    raise _numerical("BE energy balance failed")
            scalars = [
                float(arr[index])
                for mapping in (
                    potentials,
                    voltages,
                    currents,
                    powers,
                    kcl,
                    charge,
                    flux,
                )
                for arr in mapping.values()
            ] + [
                float(arr[index])
                for arr in (
                    energy,
                    work,
                    loss,
                    damping,
                    balance,
                    power_balance,
                    relative,
                )
            ]
            if not np.isfinite(scalars).all():
                raise _numerical("nonfinite transient field or balance")
        except (ValidationError, OverflowError, FloatingPointError) as error:
            raise ValidationError(
                f"step {index}, time {time:.17g} s: {error}", field="transient.step"
            ) from error
    prov = provenance(
        circuit,
        "transient",
        {
            "method": "backward_euler",
            "relative_tolerance": rtol,
            "current_tolerance_A": atol,
            "power_tolerance_W": ptol,
            "energy_tolerance_J": etol,
            "max_steps": steps,
            "max_history_bytes": memory,
        },
        {
            "times": grid.tolist(),
            "sources": {n: _signal_document(p) for n, p in profiles.items()},
            "capacitor_voltages_V": cv,
            "inductor_currents_A": li,
        },
    )
    return TransientResult(
        grid,
        potentials,
        voltages,
        currents,
        powers,
        kcl,
        charge,
        flux,
        energy,
        work,
        loss,
        damping,
        balance,
        power_balance,
        relative,
        prov,
    )
