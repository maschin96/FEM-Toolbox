"""Complex peak-phasor RLC reference with explicit DC and resonance semantics."""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from numbers import Complex
from types import MappingProxyType

import numpy as np

from fem_toolbox.circuits.dc import (
    _assemble,
    _linear_solve,
    _numerical,
    _preflight,
    _tolerance,
)
from fem_toolbox.circuits.mna import provenance, source_names, source_rhs, storage
from fem_toolbox.circuits.model import Circuit, _finite_scalar
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


@dataclass(frozen=True, slots=True)
class ACResult:
    """Peak phasors and absorbed complex power; exp(+i omega t).

    At f>0, S=.5*v*conj(i) and stored_energy_J is cycle-averaged. At f=0,
    real DC values use S=v*i and instantaneous stored energy. Maps are owned.
    Q_C is in C and flux_L in Wb; power_balance_VA is the sum of branch S.
    """

    frequency_hz: float
    node_potentials_V: Mapping[str, complex]
    branch_voltages_V: Mapping[str, complex]
    branch_currents_A: Mapping[str, complex]
    branch_powers_VA: Mapping[str, complex]
    kcl_residual_A: Mapping[str, complex]
    capacitor_charge_C: Mapping[str, complex]
    inductor_flux_Wb: Mapping[str, complex]
    power_balance_VA: complex
    stored_energy_J: float
    relative_residual: float
    provenance: RunProvenance

    def __post_init__(self) -> None:
        for name in (
            "node_potentials_V",
            "branch_voltages_V",
            "branch_currents_A",
            "branch_powers_VA",
            "kcl_residual_A",
            "capacitor_charge_C",
            "inductor_flux_Wb",
        ):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))


def _phasor(value: object) -> complex:
    if isinstance(value, bool) or not isinstance(value, Complex):
        raise ValidationError(
            "source must be a finite complex peak scalar", field="sources"
        )
    try:
        result = complex(value)
    except OverflowError as error:
        raise ValidationError(
            "source phasor is not representable", field="sources"
        ) from error
    if not np.isfinite(result):
        raise ValidationError("nonfinite source phasor", field="sources")
    return result


def solve_ac(
    circuit: Circuit,
    frequency_hz: float,
    *,
    sources: Mapping[str, complex] | None = None,
    relative_tolerance: float = 1e-10,
    current_tolerance_A: float = 1e-12,
    power_tolerance_VA: float = 1e-12,
    resonance_tolerance: float = 1e-12,
) -> ACResult:
    """Solve RLC peak phasors with explicit source amplitudes in V/A.

    Defaults use real component values as peak source amplitudes. f=0 requires
    real sources and DC-unique topology. Singular/ill-conditioned LU pivots
    relative to the largest equilibrated pivot diagnose resonance. No damping
    is inserted. Invalid frequencies, overrides or failed physical balances
    raise ValidationError; topology failures use CircuitTopologyError.
    """
    hz = _finite_scalar(frequency_hz, field="frequency_hz")
    omega = 2 * np.pi * hz
    if hz < 0 or not math.isfinite(omega):
        raise ValidationError(
            "frequency must be nonnegative and representable", field="frequency_hz"
        )
    rtol = _tolerance(relative_tolerance, "relative_tolerance")
    atol = _tolerance(current_tolerance_A, "current_tolerance_A")
    ptol = _tolerance(power_tolerance_VA, "power_tolerance_VA")
    resonance = _tolerance(resonance_tolerance, "resonance_tolerance")
    values = {name: _phasor(v) for name, v in source_names(circuit, sources).items()}
    if hz == 0 and any(v.imag for v in values.values()):
        raise ValidationError("DC sources must be real", field="sources")
    _preflight(circuit, capacitors=hz > 0, ideal_inductors=hz == 0)
    a, _, nodes, branches = _assemble(circuit)
    e = storage(circuit, nodes, branches)
    with np.errstate(over="ignore", invalid="ignore"):
        operator = (a + 1j * omega * e).tocsc()
        rhs = source_rhs(circuit, nodes, branches, values)
    if not np.isfinite(operator.data).all() or not np.isfinite(rhs).all():
        raise _numerical("nonfinite AC operator or source vector")
    x, relative = _linear_solve(
        operator, rhs, pivot_tolerance=resonance, singular_kind="resonance"
    )
    if relative > rtol:
        raise _numerical("AC algebraic backward error exceeds tolerance")
    potentials = {"0": 0j, **{n: complex(x[i]) for n, i in nodes.items()}}
    voltages, currents, powers, charge, flux = {}, {}, {}, {}, {}
    outgoing: dict[str, list[complex]] = {n: [] for n in circuit.nodes}
    energies = []
    factor = 0.5 if hz else 1.0
    for c in circuit.components:
        v = potentials[c.positive] - potentials[c.negative]
        if c.kind == "R":
            i = v / c.value
        elif c.kind == "C":
            i = 1j * omega * c.value * v
            charge[c.name] = c.value * v
            energies.append(factor * 0.5 * c.value * abs(v) * abs(v))
        elif c.kind == "I":
            i = values[c.name]
        else:
            i = complex(x[branches[c.name]])
            if c.kind == "L":
                flux[c.name] = c.value * i
                energies.append(factor * 0.5 * c.value * abs(i) * abs(i))
        p = factor * v * i.conjugate()
        if not np.isfinite([v, i, p]).all():
            raise _numerical(f"nonfinite AC branch result at {c.name}")
        voltages[c.name], currents[c.name], powers[c.name] = v, i, p
        outgoing[c.positive].append(i)
        outgoing[c.negative].append(-i)
    kcl = {n: sum(v, 0j) for n, v in outgoing.items()}
    total_power = sum(powers.values(), 0j)
    energy = sum(energies)
    if not np.isfinite(
        [*kcl.values(), total_power, energy, *charge.values(), *flux.values()]
    ).all():
        raise _numerical("nonfinite AC balance or storage")
    for n, terms in outgoing.items():
        if abs(kcl[n]) > atol + rtol * sum(abs(i) for i in terms):
            raise _numerical(f"AC Kirchhoff balance failed at node {n}")
    if abs(total_power) > ptol + rtol * sum(abs(p) for p in powers.values()):
        raise _numerical("AC complex power balance failed")
    prov = provenance(
        circuit,
        "ac",
        {
            "phasors": "peak_exp_positive_iwt",
            "relative_tolerance": rtol,
            "current_tolerance_A": atol,
            "power_tolerance_VA": ptol,
            "resonance_tolerance": resonance,
        },
        {
            "frequency_hz": hz,
            "sources": {n: [v.real, v.imag] for n, v in values.items()},
        },
    )
    return ACResult(
        hz,
        potentials,
        voltages,
        currents,
        powers,
        kcl,
        charge,
        flux,
        total_power,
        energy,
        relative,
        prov,
    )
