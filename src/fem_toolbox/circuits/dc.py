"""Sparse modified nodal analysis and independent branch balance reconstruction."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.metadata import version
from types import MappingProxyType

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import coo_matrix, csc_matrix, diags  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.circuits.model import (
    CIRCUIT_SCHEMA_VERSION,
    Circuit,
    Component,
    _finite_scalar,
)
from fem_toolbox.errors import ValidationError
from fem_toolbox.provenance import RunProvenance


class CircuitTopologyError(ValidationError):
    """DC topology lacks a unique solution; kind is a stable diagnosis identifier."""

    def __init__(self, *, kind: str, reason: str) -> None:
        self.kind = kind
        super().__init__(f"{kind}: {reason}", field="circuit.topology")


@dataclass(frozen=True, slots=True)
class DCResult:
    """Immutable SI solution with positive-to-negative currents and absorbed power.

    Maps use normalized node/component names. KCL residual is outgoing current
    at every node, including ground. power_balance_W is the sum of absorbed
    powers. relative_residual is the max-norm backward error on unscaled MNA.
    Provenance records circuit input, algorithm revision and balance tolerances.
    """

    node_potentials_V: Mapping[str, float]
    branch_voltages_V: Mapping[str, float]
    branch_currents_A: Mapping[str, float]
    branch_powers_W: Mapping[str, float]
    kcl_residual_A: Mapping[str, float]
    power_balance_W: float
    relative_residual: float
    provenance: RunProvenance

    def __post_init__(self) -> None:
        for name in (
            "node_potentials_V",
            "branch_voltages_V",
            "branch_currents_A",
            "branch_powers_W",
            "kcl_residual_A",
        ):
            object.__setattr__(self, name, MappingProxyType(dict(getattr(self, name))))


def _numerical(reason: str) -> ValidationError:
    return ValidationError(reason, field="circuit.numerics")


def _preflight(circuit: Circuit) -> None:
    conducting: dict[str, set[str]] = {n: set() for n in circuit.nodes}
    ideal: dict[str, list[tuple[str, float, str]]] = {n: [] for n in circuit.nodes}
    for c in circuit.components:
        if c.kind in "RVL":
            conducting[c.positive].add(c.negative)
            conducting[c.negative].add(c.positive)
        if c.kind in "VL":
            voltage = c.value if c.kind == "V" else 0.0
            ideal[c.positive].append((c.negative, voltage, c.name))
            ideal[c.negative].append((c.positive, -voltage, c.name))
    reached = {"0"}
    pending = ["0"]
    while pending:
        for node in conducting[pending.pop()]:
            if node not in reached:
                reached.add(node)
                pending.append(node)
    floating = sorted(set(circuit.nodes) - reached)
    if floating:
        raise CircuitTopologyError(
            kind="floating", reason=f"no DC R/V/L path to ground for nodes {floating}"
        )

    # Potentials relative to each ideal-source component, independent of MNA.
    offsets: dict[str, float] = {}
    path_scale: dict[str, float] = {}
    redundant = False
    for root in circuit.nodes:
        if root in offsets:
            continue
        offsets[root], path_scale[root] = 0.0, 0.0
        pending = [root]
        vertices = 0
        edges: set[str] = set()
        while pending:
            node = pending.pop()
            vertices += 1
            for neighbor, drop, name in ideal[node]:
                edges.add(name)
                if neighbor not in offsets:
                    offsets[neighbor] = offsets[node] - drop
                    path_scale[neighbor] = path_scale[node] + abs(drop)
                    if not math.isfinite(offsets[neighbor]) or not math.isfinite(
                        path_scale[neighbor]
                    ):
                        raise _numerical(
                            "ideal-source voltage path is not representable"
                        )
                    pending.append(neighbor)
        redundant |= len(edges) >= vertices
    for c in circuit.components:
        if c.kind not in "VL":
            continue
        expected = c.value if c.kind == "V" else 0.0
        # Roundoff allowance is based on voltages along both independent paths,
        # without an arbitrary absolute threshold that hides small sources.
        scale = max(path_scale[c.positive], path_scale[c.negative], abs(expected))
        defect = (
            math.fsum(
                (
                    offsets[c.positive] / scale,
                    -offsets[c.negative] / scale,
                    -expected / scale,
                )
            )
            if scale
            else 0.0
        )
        if abs(defect) > 64 * np.finfo(float).eps:
            raise CircuitTopologyError(
                kind="inconsistent_voltage_loop",
                reason=f"contradictory voltage at {c.name}",
            )
    if redundant:
        raise CircuitTopologyError(
            kind="redundant_voltage_loop",
            reason="ideal voltage/inductor cycle has undetermined currents",
        )


def _assemble(
    circuit: Circuit,
) -> tuple[csc_matrix, NDArray[np.float64], dict[str, int], dict[str, int]]:
    nodes = {n: i for i, n in enumerate(n for n in circuit.nodes if n != "0")}
    branches = {
        c.name: len(nodes) + i
        for i, c in enumerate(c for c in circuit.components if c.kind in "VL")
    }
    size = len(nodes) + len(branches)
    rhs = np.zeros(size)
    rows: list[int] = []
    columns: list[int] = []
    data: list[float] = []

    def entry(row: int, column: int, value: float) -> None:
        rows.append(row)
        columns.append(column)
        data.append(value)

    for c in circuit.components:
        incidence = [
            (nodes[n], sign)
            for n, sign in ((c.positive, 1), (c.negative, -1))
            if n != "0"
        ]
        if c.kind == "R":
            conductance = 1.0 / c.value
            if not math.isfinite(conductance):
                raise _numerical(f"conductance at {c.name} is not representable")
            for row, a in incidence:
                for column, b in incidence:
                    entry(row, column, a * b * conductance)
        elif c.kind == "I":
            for row, sign in incidence:
                rhs[row] -= sign * c.value
        elif c.kind in "VL":
            index = branches[c.name]
            for row, sign in incidence:
                entry(row, index, sign)
                entry(index, row, sign)
            rhs[index] = c.value if c.kind == "V" else 0.0
    matrix = coo_matrix((data, (rows, columns)), shape=(size, size)).tocsc()
    if not np.isfinite(matrix.data).all() or not np.isfinite(rhs).all():
        raise _numerical("assembled DC system is not finite")
    return matrix, rhs, nodes, branches


def _linear_solve(
    matrix: csc_matrix, rhs: NDArray[np.float64]
) -> tuple[NDArray[np.float64], float]:
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        row_scale = 1.0 / np.asarray(abs(matrix).max(axis=1).toarray()).ravel()
        row_scaled = diags(row_scale) @ matrix
        column_scale = 1.0 / np.asarray(abs(row_scaled).max(axis=0).toarray()).ravel()
        scaled = (row_scaled @ diags(column_scale)).tocsc()
        scaled_rhs = row_scale * rhs
    if (
        not np.isfinite(row_scale).all()
        or not np.isfinite(column_scale).all()
        or not np.isfinite(scaled.data).all()
        or not np.isfinite(scaled_rhs).all()
    ):
        raise _numerical("DC equilibration is not representable")
    try:
        with np.errstate(over="ignore", invalid="ignore"):
            solution = column_scale * splu(scaled).solve(scaled_rhs)
    except RuntimeError as error:
        raise CircuitTopologyError(
            kind="singular", reason="numerically singular DC system"
        ) from error
    with np.errstate(over="ignore", invalid="ignore"):
        residual = matrix @ solution - rhs
        denominator = float(np.max(abs(matrix) @ abs(solution) + abs(rhs)))
    if (
        not np.isfinite(solution).all()
        or not np.isfinite(residual).all()
        or not math.isfinite(denominator)
    ):
        raise _numerical("DC solution or residual is not finite")
    absolute = float(np.max(abs(residual)))
    relative = absolute / denominator if denominator else absolute
    return solution, relative


def _tolerance(value: float, name: str) -> float:
    result = _finite_scalar(value, field=name)
    if result < 0:
        raise ValidationError(
            "tolerance must be a finite nonnegative scalar", field=name
        )
    return result


def solve_dc(
    circuit: Circuit,
    *,
    relative_tolerance: float = 1e-10,
    current_tolerance_A: float = 1e-12,
    power_tolerance_W: float = 1e-12,
    voltage_tolerance_V: float = 1e-12,
) -> DCResult:
    """Solve a unique linear DC circuit with sparse MNA, without FEM coupling.

    C is open and L an ideal short. Preflight rejects floating DC nodes and
    incompatible/redundant ideal-source loops with CircuitTopologyError. Invalid
    numeric ranges or failed balances raise ValidationError. Physical gates use
    absolute SI tolerance plus relative_tolerance times local current throughput,
    total absolute power or ideal-branch voltage scale, respectively. Unscaled
    algebraic backward error must also not exceed relative_tolerance.
    """
    if not isinstance(circuit, Circuit):
        raise ValidationError("expected Circuit", field="circuit")
    rtol = _tolerance(relative_tolerance, "relative_tolerance")
    atol = _tolerance(current_tolerance_A, "current_tolerance_A")
    ptol = _tolerance(power_tolerance_W, "power_tolerance_W")
    vtol = _tolerance(voltage_tolerance_V, "voltage_tolerance_V")
    _preflight(circuit)
    matrix, rhs, nodes, branches = _assemble(circuit)
    solution, relative = _linear_solve(matrix, rhs)
    if relative > rtol:
        raise _numerical("DC algebraic backward error exceeds tolerance")
    potentials = {"0": 0.0, **{n: float(solution[i]) for n, i in nodes.items()}}
    voltages: dict[str, float] = {}
    currents: dict[str, float] = {}
    powers: dict[str, float] = {}
    outgoing: dict[str, list[float]] = {n: [] for n in circuit.nodes}
    for c in circuit.components:
        voltage = potentials[c.positive] - potentials[c.negative]
        current = _branch_current(c, voltage, solution, branches)
        power = voltage * current
        if not all(math.isfinite(v) for v in (voltage, current, power)):
            raise _numerical(f"nonfinite branch result at {c.name}")
        if c.kind in "VL":
            expected = c.value if c.kind == "V" else 0.0
            scale = max(
                abs(potentials[c.positive]), abs(potentials[c.negative]), abs(expected)
            )
            if abs(voltage - expected) > vtol + rtol * scale:
                raise _numerical(f"voltage constraint failed at {c.name}")
        voltages[c.name], currents[c.name], powers[c.name] = voltage, current, power
        outgoing[c.positive].append(current)
        outgoing[c.negative].append(-current)
    try:
        kcl = {n: math.fsum(values) for n, values in outgoing.items()}
        for n, values in outgoing.items():
            if abs(kcl[n]) > atol + rtol * math.fsum(abs(v) for v in values):
                raise _numerical(f"Kirchhoff current balance failed at node {n}")
        power_balance = math.fsum(powers.values())
        if abs(power_balance) > ptol + rtol * math.fsum(
            abs(p) for p in powers.values()
        ):
            raise _numerical("DC absorbed power balance failed")
    except OverflowError as error:
        raise _numerical("DC balance sums are not representable") from error
    canonical = json.dumps(
        {
            "circuit_schema_version": CIRCUIT_SCHEMA_VERSION,
            "components": [
                [c.name, c.positive, c.negative, c.value] for c in circuit.components
            ],
        },
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("ascii")
    provenance = RunProvenance.from_input_bytes(
        canonical,
        solver_name="fem_toolbox.circuits.solve_dc",
        solver_parameters={
            "analysis": "dc",
            "method": "equilibrated_sparse_mna_lu",
            "revision": 1,
            "numpy": np.__version__,
            "scipy": version("scipy"),
            "relative_tolerance": rtol,
            "current_tolerance_A": atol,
            "power_tolerance_W": ptol,
            "voltage_tolerance_V": vtol,
        },
    )
    return DCResult(
        potentials, voltages, currents, powers, kcl, power_balance, relative, provenance
    )


def _branch_current(
    c: Component,
    voltage: float,
    solution: NDArray[np.float64],
    branches: dict[str, int],
) -> float:
    if c.kind == "R":
        return voltage / c.value
    if c.kind == "I":
        return c.value
    if c.kind == "C":
        return 0.0
    return float(solution[branches[c.name]])
