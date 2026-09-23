"""Strict additive schema 3 for modal, transient and harmonic piezo models."""

from __future__ import annotations

import json
from collections.abc import Mapping
from typing import Any, cast

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.assembly import indices, real_array
from fem_toolbox.dynamic_result import PiezoDynamicResult
from fem_toolbox.errors import ModelSchemaError, ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model, _mapping, _number, _required, _unit
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes
from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
from fem_toolbox.solvers.piezo_transient import (
    KinematicValues,
    solve_transient_piezoelectric,
)
from fem_toolbox.solvers.piezoelectric import PiezoelectricResult
from fem_toolbox.time import TimeSignal, integration_times, signal_value
from fem_toolbox.transient_model import _profile


def _vector_quantity(raw: object, unit: str, path: str) -> NDArray[np.float64]:
    table = _mapping(raw, path, {"values", "unit"})
    _unit(table, unit, path)
    values = np.asarray(_required(table, "values", path))
    if values.ndim != 1 or not values.size:
        raise ModelSchemaError(path=path, reason="expected nonempty vector")
    return real_array(values, values.shape)


def _phasor(raw: object, path: str) -> complex:
    table = _mapping(raw, path, {"real", "imag"})
    return complex(
        _number(table.get("real", 0), path + ".real"),
        _number(table.get("imag", 0), path + ".imag"),
    )


class PiezoDynamicModel:
    """Internal schema adapter, reached through Model and Model.from_file."""

    def __init__(self, mesh: Mesh, document: Mapping[str, Any]) -> None:
        try:
            self._initialize(mesh, document)
        except ValidationError as error:
            if isinstance(error, ModelSchemaError):
                raise
            raise ModelSchemaError(path="analysis", reason=str(error)) from error

    def _initialize(self, mesh: Mesh, document: Mapping[str, Any]) -> None:
        table = _mapping(
            document,
            "document",
            {
                "schema_version",
                "units",
                "physics",
                "mesh",
                "base",
                "density",
                "analysis",
            },
        )
        if (
            type(table.get("schema_version")) is not int
            or table["schema_version"] != 3
            or table.get("units") != "SI"
        ):
            raise ModelSchemaError(
                path="schema_version",
                reason="dynamic piezo requires schema 3 and SI units",
            )
        self.document = document
        self.physics = table["physics"]
        raw_base = dict(_mapping(_required(table, "base", "document"), "base"))
        if (
            raw_base.get("schema_version") != 1
            or raw_base.get("physics") != "piezoelectric"
        ):
            raise ModelSchemaError(
                path="base", reason="expected stationary schema-1 piezoelectric model"
            )
        if "mesh" in raw_base:
            raise ModelSchemaError(
                path="base.mesh", reason="declare mesh only at schema-3 root"
            )
        if "mesh" in table:
            raw_base["mesh"] = table["mesh"]
        self.base = Model(mesh=mesh, document=raw_base)
        self.materials = cast(
            tuple[PiezoelectricMaterial, ...], self.base.assignment.cell_materials
        )
        density = _mapping(
            _required(table, "density", "document"), "density", {"value", "unit"}
        )
        _unit(density, "kg/m^3", "density")
        raw = np.asarray(_required(density, "value", "density"))
        self.density = real_array(raw, raw.shape)
        if self.density.shape not in ((), (mesh.n_tetra4,)) or np.any(
            self.density <= 0
        ):
            raise ModelSchemaError(
                path="density",
                reason="positive scalar or one density per Tet4 required",
            )
        keys = {
            "piezo_modes": {"count"},
            "transient_piezoelectric": {
                "times",
                "load_factor",
                "prescribed",
                "initial_displacement",
                "initial_velocity",
            },
            "harmonic_piezoelectric": {"frequencies", "load_factor", "prescribed"},
        }
        self.analysis = _mapping(
            _required(table, "analysis", "document"), "analysis", keys[self.physics]
        )
        self.n = len(np.unique(mesh.tetra4))
        self.nu = 3 * self.n
        for bc in self.base._conditions:
            if np.any(bc.values):
                raise ModelSchemaError(
                    path="base.boundary",
                    reason="base conditions must be homogeneous; use analysis.prescribed",
                )
        for electrode in self.base._electrodes:
            if electrode.potential not in (None, 0) or electrode.charge not in (
                None,
                0,
            ):
                raise ModelSchemaError(
                    path="base.boundary", reason="base electrodes must be homogeneous"
                )
        self.signals: list[TimeSignal] = []
        self.entries: list[tuple[NDArray[np.int64], str, Any]] = []
        if self.physics == "piezo_modes":
            count = self.analysis.get("count", 6)
            if type(count) is not int or count < 1:
                raise ModelSchemaError(
                    path="analysis.count", reason="positive integer required"
                )
        elif self.physics == "transient_piezoelectric":
            self.grid = _vector_quantity(
                _required(self.analysis, "times", "analysis"), "s", "analysis.times"
            )
            self.factor = _profile(
                self.analysis.get("load_factor", 1.0), "1", "analysis.load_factor"
            )
            self.signals.append(self.factor)
            for key in ("initial_displacement", "initial_velocity"):
                if key in self.analysis:
                    real_array(self.analysis[key], (self.n, 3))
        else:
            self.grid = _vector_quantity(
                _required(self.analysis, "frequencies", "analysis"),
                "Hz",
                "analysis.frequencies",
            )
            if np.any(self.grid < 0) or np.any(np.diff(self.grid) <= 0):
                raise ModelSchemaError(
                    path="analysis.frequencies",
                    reason="strictly increasing nonnegative Hz required",
                )
            self.harmonic_factor = _phasor(
                self.analysis.get("load_factor", {"real": 1}), "analysis.load_factor"
            )
        raw_entries = self.analysis.get("prescribed", [])
        if not isinstance(raw_entries, list):
            raise ModelSchemaError(path="analysis.prescribed", reason="expected list")
        occupied: set[int] = set()
        for i, raw_entry in enumerate(raw_entries):
            path = f"analysis.prescribed[{i}]"
            entry = _mapping(
                raw_entry,
                path,
                {"dofs", "unit", "value", "coefficients", "real", "imag"},
            )
            ids = indices(_required(entry, "dofs", path), 4 * self.n)
            if (
                not len(ids)
                or len(np.unique(ids)) != len(ids)
                or occupied.intersection(ids.tolist())
            ):
                raise ModelSchemaError(
                    path=path, reason="nonempty disjoint unique DOFs required"
                )
            occupied.update(ids.tolist())
            mechanical = bool(np.all(ids < self.nu))
            if not mechanical and not np.all(ids >= self.nu):
                raise ModelSchemaError(
                    path=path, reason="cannot mix displacement and voltage"
                )
            _unit(entry, "m" if mechanical else "V", path)
            if self.physics == "harmonic_piezoelectric":
                _mapping(entry, path, {"dofs", "unit", "real", "imag"})
                self.entries.append(
                    (
                        ids,
                        "phasor",
                        _phasor(
                            {k: entry[k] for k in ("real", "imag") if k in entry}, path
                        ),
                    )
                )
            elif mechanical:
                _mapping(entry, path, {"dofs", "unit", "coefficients"})
                coefficients = np.asarray(_required(entry, "coefficients", path))
                if coefficients.ndim != 1 or not len(coefficients):
                    raise ModelSchemaError(
                        path=path,
                        reason="polynomial coefficients in ascending powers of seconds required",
                    )
                self.entries.append(
                    (ids, "polynomial", real_array(coefficients, coefficients.shape))
                )
            else:
                _mapping(entry, path, {"dofs", "unit", "value"})
                signal = _profile(_required(entry, "value", path), "V", path + ".value")
                self.entries.append((ids, "voltage", signal))
                self.signals.append(signal)
        if self.physics == "transient_piezoelectric":
            self.grid = integration_times(self.grid, self.signals)

    def _prescribed(self, time: float) -> KinematicValues:
        u = np.zeros(self.nu)
        v = u.copy()
        a = u.copy()
        phi = np.zeros(self.n)
        for ids, kind, data in self.entries:
            if kind == "voltage":
                phi[ids - self.nu] = signal_value(data, time)
            else:
                u[ids] = np.polynomial.polynomial.polyval(time, data)
                v[ids] = np.polynomial.polynomial.polyval(
                    time, np.polynomial.polynomial.polyder(data)
                )
                a[ids] = np.polynomial.polynomial.polyval(
                    time, np.polynomial.polynomial.polyder(data, 2)
                )
        return KinematicValues(
            u.reshape(-1, 3), v.reshape(-1, 3), a.reshape(-1, 3), phi
        )

    def solve(self) -> PiezoDynamicResult:
        base_result = cast(PiezoelectricResult, self.base.solve().field_result)
        force = np.concatenate(
            (base_result.applied_force.ravel(), base_result.applied_electric_load)
        )
        args = (self.base.mesh, self.materials, self.density, self.base._conditions)
        if self.physics == "piezo_modes":
            if np.any(force):
                raise ModelSchemaError(
                    path="base.loads",
                    reason="modal perturbation model requires zero loads",
                )
            solution: Any = solve_piezoelectric_modes(
                *args,
                count=self.analysis.get("count", 6),
                electrodes=self.base._electrodes,
            )
            hashes = [solution.provenance.input_sha256]
        elif self.physics == "transient_piezoelectric":
            solution = solve_transient_piezoelectric(
                *args,
                self.grid,
                electrodes=self.base._electrodes,
                loads=lambda t: force * signal_value(self.factor, t),
                prescribed=self._prescribed,
                initial_displacement=self.analysis.get("initial_displacement"),
                initial_velocity=self.analysis.get("initial_velocity"),
            )
            hashes = [solution.provenance.input_sha256]
        else:
            g = np.zeros(4 * self.n, dtype=complex)
            for ids, _, value in self.entries:
                g[ids] = value
            solution = tuple(
                solve_harmonic_piezoelectric(
                    *args,
                    float(hz),
                    electrodes=self.base._electrodes,
                    loads=force * self.harmonic_factor,
                    prescribed=g,
                )
                for hz in self.grid
            )
            hashes = [state.provenance.input_sha256 for state in solution]
        provenance = RunProvenance.from_input_bytes(
            json.dumps(
                {"document": self.document, "solver_hashes": hashes},
                sort_keys=True,
                allow_nan=False,
            ).encode(),
            solver_name="fem_toolbox.model.Model.solve." + self.physics,
            solver_parameters={
                "model_schema_version": 3,
                "units": "SI",
                "field_solver_hashes": hashes,
            },
        )
        return PiezoDynamicResult(self.base.mesh, self.materials, solution, provenance)
