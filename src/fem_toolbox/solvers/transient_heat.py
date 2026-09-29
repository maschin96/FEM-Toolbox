"""Backward-Euler Tet4/Tet10 heat reference with storage and reservoir reactions."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from itertools import combinations
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]

from fem_toolbox.assembly import (
    LinearSystem,
    MatrixBlock,
    assemble,
    indices,
    real_array,
)
from fem_toolbox.boundary import Dirichlet, eliminate_dirichlet
from fem_toolbox.elements import Tet4Geometry, Tet10Geometry
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.physics.heat import (
    Convection,
    ScalarField,
)
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.linear_backend import LinearBackend, solve_linear
from fem_toolbox.solvers.thermal_elements import (
    Lumping,
    ThermalGeometry,
    capacity_matrix,
    conductance_matrix,
    contains,
    convection_terms,
    flux_load,
    source_load,
    temperature_gradient,
)
from fem_toolbox.time import TimeSignal, integration_times, signal_value

SpaceTimeField: TypeAlias = TimeSignal | Callable[[NDArray[np.float64], float], float]


@dataclass(frozen=True)
class TemperatureBC:
    """Mesh node indices with one prescribed history in K (not global DOFs)."""

    nodes: ArrayLike
    value: TimeSignal


@dataclass(frozen=True)
class TransientConvection:
    """Uniform facet coefficient W/(m² K) and ambient K, optionally time-varying."""

    coefficient: TimeSignal
    ambient: TimeSignal


@dataclass(frozen=True)
class HeatFrame:
    """A stored time state: temperature K, gradient K/m, heat_flux W/m².

    temperature/reactions/residual use active node order; element fields use
    input cell order. Tet4 element fields are constant (field_points=None).
    Tet10 fields are raw centroid samples; field_points (ncells,3) locates them
    in metres. Initial reactions/residual are zero (not evaluated).
    """

    mesh_sha256: str
    index: int
    time: float
    temperature: NDArray[np.float64]
    gradient: NDArray[np.float64]
    heat_flux: NDArray[np.float64]
    reactions: NDArray[np.float64]
    residual: NDArray[np.float64]
    field_points: NDArray[np.float64] | None = None

    def __post_init__(self) -> None:
        for name in (
            "temperature",
            "gradient",
            "heat_flux",
            "reactions",
            "residual",
            "field_points",
        ):
            a = getattr(self, name)
            if a is None:
                continue
            object.__setattr__(
                self, name, np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)
            )


@dataclass(frozen=True)
class HeatBalance:
    """Per-step powers in W; stored_energy_change in J relative to initial T."""

    time: float
    dt: float
    source_power: float
    prescribed_outflow: float
    convection_outflow: float
    dirichlet_inflow: float
    storage_power: float
    balance: float
    relative_residual: float
    stored_energy_change: float


@dataclass(frozen=True)
class TransientHeatResult:
    """Complete time grid; retained output frames; probes (ntime,nprobe) in K.

    Probes and balances are recorded at every integration step, independently
    of output_every. node_ids maps frame nodal arrays to zero-based mesh nodes.
    """

    mesh_sha256: str
    times: NDArray[np.float64]
    node_ids: NDArray[np.int64]
    frames: tuple[HeatFrame, ...]
    balances: tuple[HeatBalance, ...]
    probe_names: tuple[str, ...]
    probe_temperature: NDArray[np.float64]
    provenance: RunProvenance

    def __post_init__(self) -> None:
        for name in ("times", "node_ids", "probe_temperature"):
            a = getattr(self, name)
            object.__setattr__(
                self, name, np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)
            )


def heat_mesh_fingerprint(mesh: Mesh) -> str:
    """Hash ordered geometry/topology to prevent exporting fields on another mesh."""
    return hashlib.sha256(
        json.dumps(
            {
                "points": mesh.points.tolist(),
                "cells": mesh.cells.tolist(),
                "facets": mesh.facets.tolist(),
            },
            separators=(",", ":"),
        ).encode()
    ).hexdigest()


def _at(field: SpaceTimeField, time: float) -> ScalarField:
    # Profile objects are callable too; only other callables receive (x,t).
    from fem_toolbox.time import Ramp, Step, TimeSeries

    if isinstance(field, (Step, Ramp, TimeSeries)) or not callable(field):
        return signal_value(field, time, side="left")
    left_time = float(np.nextafter(time, -np.inf))
    return lambda x: field(x, left_time)


def _probes(
    mesh: Mesh,
    active: NDArray[np.int64],
    geometries: Sequence[ThermalGeometry],
    probes: Mapping[str, int | ArrayLike],
) -> tuple[tuple[str, ...], NDArray[np.float64]]:
    rows = []
    for name, location in sorted(probes.items()):
        if not isinstance(name, str) or not name.strip():
            raise ValidationError("probe name must be nonempty", field="probes")
        row = np.zeros(len(active))
        if isinstance(location, (int, np.integer)) and not isinstance(
            location, (bool, np.bool_)
        ):
            if location not in active:
                raise ValidationError("probe node is not active", field="probes")
            row[np.searchsorted(active, location)] = 1.0
        else:
            point = real_array(location, (3,))
            for cell, geometry in zip(mesh.cells, geometries, strict=True):
                weights = geometry.shape_values(point)
                if contains(geometry, point):
                    row[np.searchsorted(active, cell)] = weights
                    break
            else:
                raise ValidationError(
                    f"probe {name!r} outside tetrahedral domain", field="probes"
                )
        rows.append(row)
    return tuple(sorted(probes)), np.asarray(rows).reshape(len(rows), len(active))


def solve_transient_heat(
    mesh: Mesh,
    material: ThermalMaterial | Sequence[ThermalMaterial],
    times: ArrayLike,
    initial_temperature: ArrayLike,
    *,
    backend: LinearBackend | None = None,
    conditions: Sequence[TemperatureBC] = (),
    source: SpaceTimeField = 0.0,
    outward_fluxes: Mapping[int, SpaceTimeField] | None = None,
    convection: Mapping[int, TransientConvection] | None = None,
    probes: Mapping[str, int | ArrayLike] | None = None,
    lumped: bool = False,
    lumping: Lumping | None = None,
    output_every: int = 1,
    max_history_bytes: int = 256 * 1024 * 1024,
    on_frame: Callable[[HeatFrame], None] | None = None,
) -> TransientHeatResult:
    """Solve C Tdot + K T = f by backward Euler on an event-enriched time grid.

    ``backend`` is an experimental linear-backend hook; None selects SciPy.
    Reduction, physical residuals and result units remain this solver's contract.

    Straight Tet10 uses consistent capacity by default; lumping="hrz" selects
    conservative diagonal scaling. Tet10 rejects lumped=True (negative row-sum
    masses). Tet4 retains its boolean lumped option and rejects explicit lumping.
    HRZ conserves total storage but does not preserve consistent source moments.

    SI: time s, initial_temperature K (scalar or active-node vector), source
    W/m³, outward flux W/m². Spatial sources/fluxes may be f(x,t); explicitly
    include their discontinuities in times. Profile knots are inserted for
    Step/Ramp/TimeSeries. Loads use the left limit at the interval endpoint;
    Dirichlet histories are right-continuous. No artificial thermal reference
    is imposed on insulated components. Positive rho/c required per cell.

    The initial state must match t0 Dirichlet values (1e-10 K absolute
    tolerance); it is never projected silently. Full residuals include C dT/dt.
    Output includes first/final frames and every output_every step; optional
    on_frame permits incremental export. max_history_bytes bounds estimated
    retained numerical history, not factorization/workspace memory.
    """
    if not mesh.n_cells:
        raise ValidationError(
            "transient heat requires Tet4 or Tet10 cells", field="mesh"
        )
    if (
        type(output_every) is not int
        or output_every < 1
        or type(max_history_bytes) is not int
        or max_history_bytes < 1
    ):
        raise ValidationError(
            "positive output interval/history budget required", field="output"
        )
    materials = (
        [material] * mesh.n_cells
        if isinstance(material, ThermalMaterial)
        else list(material)
    )
    if len(materials) != mesh.n_cells or any(
        not isinstance(m, ThermalMaterial)
        or m.density is None
        or m.specific_heat is None
        for m in materials
    ):
        raise ValidationError(
            "one material with storage rho/c required per cell", field="material"
        )
    fluxes, robin = dict(outward_fluxes or {}), dict(convection or {})
    if any(not isinstance(bc, TemperatureBC) for bc in conditions):
        raise ValidationError("expected TemperatureBC", field="conditions")
    if any(not isinstance(bc, TransientConvection) for bc in robin.values()):
        raise ValidationError("expected TransientConvection", field="convection")
    from fem_toolbox.time import Ramp, Step, TimeSeries

    signals = [
        v
        for v in [source, *fluxes.values()]
        if not callable(v) or isinstance(v, (Step, Ramp, TimeSeries))
    ]
    signals.extend(b.value for b in conditions)
    signals.extend(v for b in robin.values() for v in (b.coefficient, b.ambient))
    grid = integration_times(times, signals)
    if mesh.n_tetra10:
        inspect_boundary_topology(mesh)
    active = np.unique(mesh.cells)
    size = len(active)
    geometries = [
        (Tet10Geometry if mesh.n_tetra10 else Tet4Geometry)(mesh.points[cell])
        for cell in mesh.cells
    ]
    cell_dofs = [np.searchsorted(active, cell) for cell in mesh.cells]
    names, probe_matrix = _probes(mesh, active, geometries, probes or {})
    output_indices = set(range(0, len(grid), output_every)) | {len(grid) - 1}
    estimated = 8 * (
        len(output_indices) * (3 * size + (9 if mesh.n_tetra10 else 6) * mesh.n_cells)
        + len(grid) * (len(names) + 12)
    )
    if estimated > max_history_bytes:
        raise ValidationError(
            f"history budget exceeded: estimated {estimated} bytes", field="output"
        )
    t0 = np.asarray(initial_temperature)
    initial = (
        np.full(size, float(real_array(t0, ())))
        if t0.ndim == 0
        else real_array(t0, (size,))
    )
    if np.any(initial < 0):
        raise ValidationError(
            "initial temperature must be nonnegative K", field="initial"
        )
    bc_dofs = []
    for bc in conditions:
        nodes = indices(bc.nodes, mesh.n_points)
        if not len(nodes) or not np.isin(nodes, active).all():
            raise ValidationError(
                "temperature boundary requires active nodes", field="conditions"
            )
        bc_dofs.append(np.searchsorted(active, nodes))
    fixed = (
        np.unique(np.concatenate(bc_dofs)) if bc_dofs else np.empty(0, dtype=np.int64)
    )
    free = np.setdiff1d(np.arange(size), fixed)

    def constraints(time: float) -> list[Dirichlet]:
        values = [signal_value(b.value, time) for b in conditions]
        if any(v < 0 for v in values):
            raise ValidationError(
                "negative prescribed Kelvin temperature", field="conditions"
            )
        return [
            Dirichlet(ids, np.full(len(ids), v))
            for ids, v in zip(bc_dofs, values, strict=True)
        ]

    # Validate all event/end values and overlaps before assembly and output.
    identity = LinearSystem(sparse.eye(size, format="csr"), np.zeros(size))
    for time in grid:
        constrained = eliminate_dirichlet(identity, constraints(float(time)))
        if time == grid[0] and not np.allclose(
            initial[fixed], constrained.rhs[fixed], rtol=0, atol=1e-10
        ):
            raise ValidationError(
                "initial temperature conflicts with Dirichlet data", field="initial"
            )
        for condition in robin.values():
            if (
                signal_value(condition.coefficient, float(time)) < 0
                or signal_value(condition.ambient, float(time)) < 0
            ):
                raise ValidationError(
                    "negative convection coefficient or ambient K", field="convection"
                )
    stiffness_blocks, capacity_blocks = [], []
    for ids, geometry, m in zip(cell_dofs, geometries, materials, strict=True):
        assert m.density is not None and m.specific_heat is not None
        stiffness_blocks.append(MatrixBlock(ids, ids, conductance_matrix(geometry, m)))
        capacity_blocks.append(
            MatrixBlock(
                ids,
                ids,
                capacity_matrix(
                    geometry, m.density, m.specific_heat, lumped=lumped, lumping=lumping
                ),
            )
        )
    stiffness = assemble(size, stiffness_blocks).matrix
    capacity = assemble(size, capacity_blocks).matrix
    face_counts = Counter(
        tuple(sorted(f)) for cell in mesh.cells for f in combinations(cell[:4], 3)
    )
    face_data = {}
    for face in set(fluxes) | set(robin):
        i = int(indices([face], mesh.n_facets)[0])
        nodes = mesh.facets[i]
        if face_counts[tuple(sorted(nodes[:3]))] != 1:
            raise ValidationError(
                "boundary must be exterior tetrahedral face", field="surface"
            )
        face_data[face] = (np.searchsorted(active, nodes), mesh.points[nodes])
    digest = hashlib.sha256()

    def record(value: object) -> None:
        digest.update(
            json.dumps(
                value, sort_keys=True, allow_nan=False, separators=(",", ":")
            ).encode()
        )
        digest.update(b"\n")

    record(
        {
            "points": mesh.points.tolist(),
            "cells": mesh.cells.tolist(),
            "facets": mesh.facets.tolist(),
            "initial": initial.tolist(),
            "materials": [
                {
                    "k": m.conductivity.tolist(),
                    "rho": m.density,
                    "c": m.specific_heat,
                    "source": asdict(m.source) if m.source else None,
                }
                for m in materials
            ],
            "times": grid.tolist(),
            "lumped": lumped,
            **({"lumping": lumping or "consistent"} if mesh.n_tetra10 else {}),
            "probe_names": names,
            "probe_weights": probe_matrix.tolist(),
            "output_every": output_every,
        }
    )
    frames: list[HeatFrame] = []
    balances: list[HeatBalance] = []
    measurements = np.empty((len(grid), len(names)))

    def save(
        index: int,
        temperature: NDArray[np.float64],
        reactions: NDArray[np.float64],
        residual: NDArray[np.float64],
    ) -> None:
        measurements[index] = probe_matrix @ temperature
        if index not in output_indices:
            return
        gradients = np.array(
            [
                temperature_gradient(g, temperature[ids])
                for ids, g in zip(cell_dofs, geometries, strict=True)
            ]
        )
        frame = HeatFrame(
            mesh_hash,
            index,
            float(grid[index]),
            temperature,
            gradients,
            np.array(
                [
                    -m.conductivity @ grad
                    for m, grad in zip(materials, gradients, strict=True)
                ]
            ),
            reactions,
            residual,
            mesh.points[mesh.cells[:, :4]].mean(axis=1) if mesh.n_tetra10 else None,
        )
        if on_frame is not None:
            on_frame(frame)
        frames.append(frame)

    mesh_hash = heat_mesh_fingerprint(mesh)
    temperature = initial.copy()
    save(0, temperature, np.zeros(size), np.zeros(size))
    for index in range(1, len(grid)):
        time, dt = float(grid[index]), float(grid[index] - grid[index - 1])
        source_vector, flux_vector, ambient_vector = (np.zeros(size) for _ in range(3))
        for ids, geometry in zip(cell_dofs, geometries, strict=True):
            np.add.at(source_vector, ids, source_load(geometry, _at(source, time)))
        for face, flux in fluxes.items():
            ids, points = face_data[face]
            np.add.at(flux_vector, ids, flux_load(points, _at(flux, time)))
        robin_blocks, robin_terms = [], []
        for face, robin_bc in robin.items():
            ids, points = face_data[face]
            h, load = convection_terms(
                points,
                Convection(
                    signal_value(robin_bc.coefficient, time, side="left"),
                    signal_value(robin_bc.ambient, time, side="left"),
                ),
            )
            robin_blocks.append(MatrixBlock(ids, ids, h))
            robin_terms.append((ids, h, load))
            np.add.at(ambient_vector, ids, load)
        k = stiffness + assemble(size, robin_blocks).matrix
        load = source_vector + flux_vector + ambient_vector
        original = LinearSystem(capacity / dt + k, load + capacity @ temperature / dt)
        bcs = constraints(time)
        constrained = eliminate_dirichlet(original, bcs)
        record(
            {
                "time": time,
                "source": source_vector.tolist(),
                "flux": flux_vector.tolist(),
                "robin": [
                    (ids.tolist(), h.tolist(), f.tolist()) for ids, h, f in robin_terms
                ],
                "fixed": fixed.tolist(),
                "values": constrained.rhs[fixed].tolist(),
            }
        )
        next_temperature = constrained.rhs.copy()
        if len(free):
            a = constrained.matrix.tocsr()[free][:, free].tocsc()
            diagonal = a.diagonal()
            if np.any(diagonal <= 0):
                raise ValidationError("nonpositive transient diagonal", field="solver")
            scale = 1 / np.sqrt(diagonal)
            next_temperature[free] = solve_linear(
                a,
                constrained.rhs[free],
                fields={"T": np.arange(len(free))},
                scale=scale,
                backend=backend,
            ).values
        if not np.isfinite(next_temperature).all():
            raise ValidationError("nonfinite transient temperature", field="solver")
        storage = capacity @ ((next_temperature - temperature) / dt)
        residual = np.asarray(storage + k @ next_temperature - load)
        if not np.isfinite(storage).all() or not np.isfinite(residual).all():
            raise ValidationError("nonfinite transient residual", field="solver")
        reaction = np.zeros(size)
        reaction[fixed] = residual[fixed]
        roundoff = abs(original.matrix) @ np.abs(next_temperature) + np.abs(
            original.rhs
        )
        denominator = float(np.linalg.norm(roundoff))
        relative = (
            float(np.linalg.norm(residual[free])) / denominator if denominator else 0.0
        )
        convection_power = sum(
            float((h @ next_temperature[ids] - f).sum()) for ids, h, f in robin_terms
        )
        source_power, prescribed, reservoir, stored = (
            float(source_vector.sum()),
            -float(flux_vector.sum()),
            float(reaction.sum()),
            float(storage.sum()),
        )
        balance = source_power + reservoir - prescribed - convection_power - stored
        tolerance = 1e-9 * (
            abs(source_power)
            + abs(reservoir)
            + abs(prescribed)
            + abs(convection_power)
            + abs(stored)
        ) + 100 * np.finfo(float).eps * float(roundoff.sum())
        if relative > 1e-9 or abs(balance) > tolerance:
            raise ValidationError(
                "transient residual or energy balance failed", field="solver"
            )
        balances.append(
            HeatBalance(
                time,
                dt,
                source_power,
                prescribed,
                convection_power,
                reservoir,
                stored,
                balance,
                relative,
                float((capacity @ (next_temperature - initial)).sum()),
            )
        )
        temperature = next_temperature
        save(index, temperature, reaction, residual)
    provenance = RunProvenance.from_input_bytes(
        digest.digest(),
        solver_name="fem_toolbox.solvers.transient_heat.solve_transient_heat",
        solver_parameters={
            "linear_backend": "scipy" if backend is None else backend.name,
            "input_schema": "transient-heat-discrete-v1",
            "method": "backward_euler",
            "capacity": (lumping or "consistent")
            if mesh.n_tetra10
            else ("lumped" if lumped else "consistent"),
            "field_location": "cell_centroid" if mesh.n_tetra10 else "element_constant",
            "output_every": output_every,
            "load_endpoint": "left_limit",
            "dirichlet_endpoint": "right",
            "max_history_bytes": max_history_bytes,
            "residual_tolerance": 1e-9,
        },
    )
    return TransientHeatResult(
        mesh_hash,
        grid,
        active,
        tuple(frames),
        tuple(balances),
        names,
        measurements,
        provenance,
    )
