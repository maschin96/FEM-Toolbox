"""Additive model schema 2 for transient heat; schema 1 remains stationary."""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

import numpy as np

from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.materials import ThermalMaterial, assign_region_materials
from fem_toolbox.materials.catalog import thermal_material
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.model import _mapping, _number, _required, _unit
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.transient_heat import (
    TemperatureBC,
    TransientConvection,
    TransientHeatResult,
    solve_transient_heat,
)
from fem_toolbox.time import Ramp, Step, TimeSeries, TimeSignal, integration_times


def _profile(raw: object, unit: str, path: str) -> TimeSignal:
    if not isinstance(raw, Mapping):
        return _number(raw, path)
    kind = raw.get("kind")
    keys = {
        "step": {"kind", "time", "before", "after"},
        "ramp": {"kind", "start", "end", "start_value", "end_value"},
        "series": {"kind", "times", "values", "interpolation", "extrapolation", "unit"},
    }
    if not isinstance(kind, str) or kind not in keys:
        raise ModelSchemaError(
            path=path, reason="expected step, ramp or series (CSV requires from_file)"
        )
    table = _mapping(raw, path, keys[kind])

    def number(key: str) -> float:
        return _number(_required(table, key, path), f"{path}.{key}")

    if kind == "step":
        return Step(number("time"), number("before"), number("after"))
    if kind == "ramp":
        return Ramp(
            number("start"), number("end"), number("start_value"), number("end_value")
        )
    if table.get("unit", unit) != unit:
        raise ModelSchemaError(
            path=path, reason="series unit conflicts with field unit"
        )
    return TimeSeries(
        _required(table, "times", path),
        _required(table, "values", path),
        interpolation=table.get("interpolation", "linear"),
        extrapolation=table.get("extrapolation", "error"),
        unit=unit,
    )


def _quantity(raw: object, unit: str, path: str) -> TimeSignal:
    table = _mapping(raw, path, {"value", "unit"})
    _unit(table, unit, path)
    return _profile(_required(table, "value", path), unit, path + ".value")


def resolve_csv_profiles(value: Any, root: Path) -> Any:
    """Embed CSV samples in a detached document; only files under model root."""
    if isinstance(value, list):
        return [resolve_csv_profiles(v, root) for v in value]
    if not isinstance(value, dict):
        return value
    if value.get("kind") != "csv":
        return {k: resolve_csv_profiles(v, root) for k, v in value.items()}
    table = _mapping(
        value,
        "csv",
        {
            "kind",
            "path",
            "time_column",
            "value_column",
            "unit",
            "delimiter",
            "interpolation",
            "extrapolation",
        },
    )
    name = _required(table, "path", "csv")
    if not isinstance(name, str) or not name or Path(name).is_absolute():
        raise ModelSchemaError(path="csv.path", reason="expected relative CSV path")
    path = (root / name).resolve()
    if not path.is_relative_to(root):
        raise ModelSchemaError(path="csv.path", reason="path escapes model directory")
    curve = TimeSeries.from_csv(
        path,
        value_column=_required(table, "value_column", "csv"),
        unit=_required(table, "unit", "csv"),
        time_column=table.get("time_column", "time_s"),
        delimiter=table.get("delimiter", ","),
        interpolation=table.get("interpolation", "linear"),
        extrapolation=table.get("extrapolation", "error"),
    )
    return {
        "kind": "series",
        "times": curve.times.tolist(),
        "values": curve.values.tolist(),
        "unit": curve.unit,
        "interpolation": curve.interpolation,
        "extrapolation": curve.extrapolation,
    }


class TransientHeatModel:
    """Implementation behind Model(physics='transient_heat', schema_version=2)."""

    def __init__(self, mesh: Mesh, document: Mapping[str, Any]) -> None:
        table = _mapping(
            document,
            "document",
            {
                "schema_version",
                "units",
                "physics",
                "mesh",
                "materials",
                "regions",
                "time",
                "initial_temperature",
                "loads",
                "boundary",
                "probes",
                "solver",
                "output",
            },
        )
        if (
            type(table.get("schema_version")) is not int
            or table["schema_version"] != 2
            or table.get("units") != "SI"
        ):
            raise ModelSchemaError(
                path="schema_version", reason="transient heat requires schema 2 and SI"
            )
        if "mesh" in table:
            mesh_spec = _mapping(
                table["mesh"], "mesh", {"path", "length_scale_to_m", "length_unit"}
            )
            if "length_unit" in mesh_spec:
                scales = {"m": 1.0, "mm": 1e-3}
                if (
                    mesh_spec["length_unit"] not in scales
                    or mesh_spec.get("length_scale_to_m")
                    != scales[mesh_spec["length_unit"]]
                ):
                    raise ModelSchemaError(
                        path="mesh", reason="inconsistent length unit and scale"
                    )
        self.mesh = mesh
        self.document = dict(table)
        inspect_boundary_topology(mesh)
        raw_materials = _mapping(_required(table, "materials", "document"), "materials")
        catalog = {}
        for name, raw in raw_materials.items():
            path = f"materials.{name}"
            entry = _mapping(raw, path)
            if "preset" in entry:
                _mapping(entry, path, {"preset"})
                catalog[name] = thermal_material(entry["preset"])
            else:
                _mapping(
                    entry, path, {"conductivity", "density", "specific_heat", "units"}
                )
                units = _mapping(
                    _required(entry, "units", path),
                    path + ".units",
                    {"conductivity", "density", "specific_heat"},
                )
                for key, unit in [
                    ("conductivity", "W/(m K)"),
                    ("density", "kg/m^3"),
                    ("specific_heat", "J/(kg K)"),
                ]:
                    if units.get(key) != unit:
                        raise ModelSchemaError(
                            path=path + ".units." + key, reason=f"expected {unit}"
                        )
                k = _required(entry, "conductivity", path)
                rho = _number(_required(entry, "density", path), path + ".density")
                cp = _number(
                    _required(entry, "specific_heat", path), path + ".specific_heat"
                )
                catalog[name] = (
                    ThermalMaterial.isotropic(
                        _number(k, path + ".conductivity"),
                        density=rho,
                        specific_heat=cp,
                    )
                    if np.asarray(k).ndim == 0
                    else ThermalMaterial(k, density=rho, specific_heat=cp)
                )
        regions = _mapping(_required(table, "regions", "document"), "regions")
        self.assignment = assign_region_materials(mesh, "T", regions, catalog)
        time = _mapping(
            _required(table, "time", "document"), "time", {"values", "unit"}
        )
        _unit(time, "s", "time")
        self.times = integration_times(_required(time, "values", "time"))
        initial = _mapping(
            _required(table, "initial_temperature", "document"),
            "initial_temperature",
            {"value", "unit"},
        )
        _unit(initial, "K", "initial_temperature")
        self.initial = _required(initial, "value", "initial_temperature")
        loads = _mapping(table.get("loads", {}), "loads", {"source"})
        self.source = (
            _quantity(loads["source"], "W/m^3", "loads.source")
            if "source" in loads
            else 0.0
        )
        boundary = _mapping(
            table.get("boundary", {}), "boundary", {"temperature", "flux", "convection"}
        )
        self.conditions = []
        self.fluxes: dict[int, TimeSignal] = {}
        self.convection: dict[int, TransientConvection] = {}
        for kind, entries in boundary.items():
            if not isinstance(entries, list):
                raise ModelSchemaError(
                    path="boundary." + kind, reason="expected array of tables"
                )
            for i, raw in enumerate(entries):
                path = f"boundary.{kind}[{i}]"
                keys = (
                    {"set", "nodes", "value", "unit"}
                    if kind == "temperature"
                    else {"set", "value", "unit"}
                    if kind == "flux"
                    else {
                        "set",
                        "coefficient",
                        "ambient_temperature",
                        "unit",
                        "ambient_unit",
                    }
                )
                entry = _mapping(raw, path, keys)
                faces = self._facets(entry["set"], path) if "set" in entry else None
                if kind == "temperature":
                    if ("nodes" in entry) == ("set" in entry):
                        raise ModelSchemaError(
                            path=path, reason="choose nodes or facet set"
                        )
                    nodes = (
                        entry["nodes"]
                        if faces is None
                        else np.unique(mesh.triangle3[faces])
                    )
                    _unit(entry, "K", path)
                    self.conditions.append(
                        TemperatureBC(
                            nodes, _profile(_required(entry, "value", path), "K", path)
                        )
                    )
                else:
                    if faces is None:
                        raise ModelSchemaError(path=path, reason="facet set required")
                    _unit(entry, "W/m^2" if kind == "flux" else "W/(m^2 K)", path)
                    if kind == "convection" and entry.get("ambient_unit") != "K":
                        raise ModelSchemaError(
                            path=path, reason="ambient_unit must be K"
                        )
                    for face in faces:
                        if kind == "flux":
                            if int(face) in self.fluxes:
                                raise ModelSchemaError(
                                    path=path, reason="duplicate flux assignment"
                                )
                            self.fluxes[int(face)] = _profile(
                                _required(entry, "value", path), "W/m^2", path
                            )
                        else:
                            if int(face) in self.convection:
                                raise ModelSchemaError(
                                    path=path, reason="duplicate convection assignment"
                                )
                            self.convection[int(face)] = TransientConvection(
                                _profile(
                                    _required(entry, "coefficient", path),
                                    "W/(m^2 K)",
                                    path,
                                ),
                                _profile(
                                    _required(entry, "ambient_temperature", path),
                                    "K",
                                    path,
                                ),
                            )
        self.probes = {}
        for name, raw in _mapping(table.get("probes", {}), "probes").items():
            entry = _mapping(raw, "probes." + name, {"node", "point"})
            if len(entry) != 1:
                raise ModelSchemaError(
                    path="probes." + name, reason="choose node or point"
                )
            self.probes[name] = next(iter(entry.values()))
        solver = _mapping(table.get("solver", {}), "solver", {"method", "lumped"})
        if (
            solver.get("method", "backward_euler") != "backward_euler"
            or type(solver.get("lumped", False)) is not bool
        ):
            raise ModelSchemaError(
                path="solver", reason="backward_euler and boolean lumped required"
            )
        self.lumped = solver.get("lumped", False)
        output = _mapping(
            table.get("output", {}), "output", {"every", "max_history_bytes"}
        )
        self.output_every = output.get("every", 1)
        self.max_history_bytes = output.get("max_history_bytes", 256 * 1024 * 1024)
        for value in (self.output_every, self.max_history_bytes):
            if type(value) is not int or value < 1:
                raise ModelSchemaError(
                    path="output", reason="positive integer required"
                )

    def _facets(self, name: Any, path: str) -> np.ndarray:
        if (
            not isinstance(name, str)
            or name not in self.mesh.facet_sets
            or not len(self.mesh.facet_sets[name])
        ):
            raise ModelSchemaError(path=path, reason="unknown or empty facet set")
        return self.mesh.facet_sets[name]

    def solve(self) -> TransientHeatResult:
        """Solve with model/document and material-source provenance included."""
        from dataclasses import replace

        result = solve_transient_heat(
            self.mesh,
            cast(list[ThermalMaterial], self.assignment.cell_materials),
            self.times,
            self.initial,
            conditions=self.conditions,
            source=self.source,
            outward_fluxes=self.fluxes,
            convection=self.convection,
            probes=self.probes,
            lumped=self.lumped,
            output_every=self.output_every,
            max_history_bytes=self.max_history_bytes,
        )
        provenance = RunProvenance.from_input_bytes(
            json.dumps(
                {
                    "document": self.document,
                    "assignment": self.assignment.to_provenance_dict(),
                    "discrete": result.provenance.to_dict(),
                },
                sort_keys=True,
                allow_nan=False,
            ).encode(),
            solver_name="fem_toolbox.model.Model.solve",
            solver_parameters={
                **dict(result.provenance.solver_parameters),
                "model_schema": 2,
                "materials": self.assignment.to_provenance_dict(),
                "resolved_document": self.document,
            },
        )
        return replace(result, provenance=provenance)
