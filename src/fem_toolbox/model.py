"""Strict versioned TOML/Python adapter for stationary field workflows."""

from __future__ import annotations

import json
import tomllib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, cast

import numpy as np

from fem_toolbox.assembly import real_array
from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.errors import ModelSchemaError, ValidationError
from fem_toolbox.fields import DofMap
from fem_toolbox.io import read_gmsh
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    MaterialAssignment,
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
    assign_region_materials,
)
from fem_toolbox.materials.regions import Material
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.physics.heat import Convection
from fem_toolbox.provenance import RunProvenance
from fem_toolbox.solvers.elasticity import (
    ElasticResult,
    mechanical_dofs,
    solve_elasticity,
)
from fem_toolbox.solvers.electrostatics import (
    ElectrostaticResult,
    electric_dofs,
    solve_electrostatics,
)
from fem_toolbox.solvers.heat import HeatResult, solve_heat, thermal_dofs
from fem_toolbox.solvers.piezoelectric import (
    PiezoelectricResult,
    piezoelectric_dofs,
    solve_piezoelectric,
)
from fem_toolbox.solvers.thermopiezoelectric import (
    ThermopiezoelectricResult,
    solve_thermopiezoelectric,
)

Physics = Literal[
    "heat", "elasticity", "electrostatics", "piezoelectric", "thermopiezoelectric"
]
FieldResult = (
    HeatResult
    | ElasticResult
    | ElectrostaticResult
    | PiezoelectricResult
    | ThermopiezoelectricResult
)
_FIELDS: dict[str, str] = {
    "heat": "T",
    "elasticity": "u",
    "electrostatics": "phi",
    "piezoelectric": "piezo",
    "thermopiezoelectric": "piezo",
}
_BOUNDARY_KEYS = {
    "dirichlet",
    "flux",
    "convection",
    "traction",
    "surface_charge",
    "electrodes",
}


def _mapping(
    value: object, path: str, allowed: set[str] | None = None
) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
        raise ModelSchemaError(path=path, reason="expected a table with string keys")
    if allowed is not None:
        unknown = sorted(set(value) - allowed)
        if unknown:
            raise ModelSchemaError(path=path, reason=f"unknown key {unknown[0]!r}")
    return value


def _required(table: Mapping[str, Any], key: str, path: str) -> Any:
    if key not in table:
        raise ModelSchemaError(path=f"{path}.{key}", reason="required key missing")
    return table[key]


def _number(value: object, path: str) -> float:
    try:
        return float(real_array(cast(Any, value), ()))
    except (ValidationError, TypeError, ValueError) as error:
        raise ModelSchemaError(
            path=path, reason="expected a finite SI number"
        ) from error


def _array(value: object, shape: tuple[int, ...], path: str) -> np.ndarray:
    try:
        return real_array(cast(Any, value), shape)
    except (ValidationError, TypeError, ValueError) as error:
        raise ModelSchemaError(
            path=path, reason=f"expected finite array {shape}"
        ) from error


def _unit(table: Mapping[str, Any], expected: str, path: str) -> None:
    actual = _required(table, "unit", path)
    if actual != expected:
        raise ModelSchemaError(
            path=f"{path}.unit", reason=f"expected {expected!r}, got {actual!r}"
        )


def _material(field: str, value: object, path: str) -> Material:
    table = _mapping(value, path)
    try:
        if field == "T":
            _mapping(table, path, {"conductivity", "unit"})
            _unit(table, "W/(m K)", path)
            raw = _required(table, "conductivity", path)
            return (
                ThermalMaterial.isotropic(_number(raw, f"{path}.conductivity"))
                if np.asarray(raw).ndim == 0
                else ThermalMaterial(_array(raw, (3, 3), f"{path}.conductivity"))
            )
        if field == "u":
            _mapping(table, path, {"stiffness", "young", "poisson", "unit"})
            _unit(table, "Pa", path)
            if "stiffness" in table and ("young" in table or "poisson" in table):
                raise ModelSchemaError(
                    path=path, reason="choose stiffness or young/poisson"
                )
            if "stiffness" in table:
                return ElasticMaterial(
                    _array(table["stiffness"], (6, 6), f"{path}.stiffness")
                )
            return ElasticMaterial.isotropic(
                _number(_required(table, "young", path), f"{path}.young"),
                _number(_required(table, "poisson", path), f"{path}.poisson"),
            )
        if field == "phi":
            _mapping(table, path, {"permittivity", "unit"})
            _unit(table, "F/m", path)
            raw = _required(table, "permittivity", path)
            return (
                DielectricMaterial.isotropic(_number(raw, f"{path}.permittivity"))
                if np.asarray(raw).ndim == 0
                else DielectricMaterial(_array(raw, (3, 3), f"{path}.permittivity"))
            )
        if field == "expansion":
            _mapping(
                table,
                path,
                {
                    "alpha",
                    "alpha_unit",
                    "principal_alpha",
                    "axes",
                    "reference_temperature",
                    "reference_temperature_unit",
                },
            )
            if table.get("alpha_unit") != "1/K":
                raise ModelSchemaError(
                    path=f"{path}.alpha_unit", reason="expected '1/K'"
                )
            if table.get("reference_temperature_unit") != "K":
                raise ModelSchemaError(
                    path=f"{path}.reference_temperature_unit", reason="expected 'K'"
                )
            reference = _number(
                _required(table, "reference_temperature", path),
                f"{path}.reference_temperature",
            )
            if "alpha" in table:
                if "principal_alpha" in table or "axes" in table:
                    raise ModelSchemaError(
                        path=path, reason="choose alpha or principal_alpha with axes"
                    )
                return ThermalExpansion(
                    _array(table["alpha"], (6,), f"{path}.alpha"), reference
                )
            return ThermalExpansion.from_principal(
                _array(
                    _required(table, "principal_alpha", path),
                    (3,),
                    f"{path}.principal_alpha",
                ),
                reference,
                _array(_required(table, "axes", path), (3, 3), f"{path}.axes"),
            )
        _mapping(
            table,
            path,
            {"stiffness", "stress_coupling", "permittivity", "axes", "units"},
        )
        units = _mapping(
            _required(table, "units", path),
            f"{path}.units",
            {"stiffness", "stress_coupling", "permittivity"},
        )
        for key, expected in {
            "stiffness": "Pa",
            "stress_coupling": "C/m^2",
            "permittivity": "F/m",
        }.items():
            if _required(units, key, f"{path}.units") != expected:
                raise ModelSchemaError(
                    path=f"{path}.units.{key}", reason=f"expected {expected}"
                )
        axes = (
            None
            if "axes" not in table
            else _array(table["axes"], (3, 3), f"{path}.axes")
        )
        return PiezoelectricMaterial(
            _array(_required(table, "stiffness", path), (6, 6), f"{path}.stiffness"),
            _array(
                _required(table, "stress_coupling", path),
                (3, 6),
                f"{path}.stress_coupling",
            ),
            _array(
                _required(table, "permittivity", path), (3, 3), f"{path}.permittivity"
            ),
            axes=axes,
        )
    except ValidationError as error:
        if isinstance(error, ModelSchemaError):
            raise
        raise ModelSchemaError(path=path, reason=str(error)) from error


@dataclass(frozen=True, slots=True)
class ModelResult:
    """A field result with schema-aware, deterministic model provenance."""

    physics: Physics
    field_result: FieldResult
    provenance: RunProvenance


class Model:
    """Validated stationary Model from a Gmsh-backed TOML file or Python mapping.

    ``document`` has schema version 1 and SI units. The Python constructor
    takes an already validated ``Mesh``; ``from_file`` imports the declared
    relative Gmsh path. No global settings or implicit material assignments.
    """

    def __init__(self, *, mesh: Mesh, document: Mapping[str, object]) -> None:
        if not isinstance(mesh, Mesh):
            raise ModelSchemaError(path="mesh", reason="expected Mesh")
        try:
            raw = json.loads(json.dumps(document, allow_nan=False, sort_keys=True))
        except (TypeError, ValueError) as error:
            raise ModelSchemaError(
                path="document", reason="expected JSON-safe finite TOML values"
            ) from error
        table = _mapping(
            raw,
            "document",
            {
                "schema_version",
                "units",
                "physics",
                "mesh",
                "materials",
                "regions",
                "boundary",
                "loads",
                "solver",
            },
        )
        version = _required(table, "schema_version", "document")
        if type(version) is not int or version != 1:
            raise ModelSchemaError(
                path="schema_version",
                reason="supported version is 1; migrate older documents explicitly",
            )
        if _required(table, "units", "document") != "SI":
            raise ModelSchemaError(path="units", reason="expected 'SI'")
        physics = _required(table, "physics", "document")
        if not isinstance(physics, str) or physics not in _FIELDS:
            raise ModelSchemaError(
                path="physics", reason=f"unsupported mode {physics!r}"
            )
        self.physics = cast(Physics, physics)
        self.mesh = mesh
        if "mesh" in table:
            _mapping(table["mesh"], "mesh", {"path", "length_scale_to_m"})
        fields = (
            ("T", "piezo", "expansion")
            if self.physics == "thermopiezoelectric"
            else (_FIELDS[self.physics],)
        )
        material_tables = _mapping(
            _required(table, "materials", "document"), "materials", set(fields)
        )
        region_tables = _mapping(
            _required(table, "regions", "document"), "regions", set(fields)
        )
        self.assignments: dict[str, MaterialAssignment] = {}
        for field in fields:
            catalog = _mapping(
                _required(material_tables, field, "materials"), f"materials.{field}"
            )
            materials = {
                name: _material(field, value, f"materials.{field}.{name}")
                for name, value in catalog.items()
            }
            region_map = _mapping(
                _required(region_tables, field, "regions"), f"regions.{field}"
            )
            if any(not isinstance(name, str) for name in region_map.values()):
                raise ModelSchemaError(
                    path=f"regions.{field}", reason="material names must be strings"
                )
            try:
                self.assignments[field] = assign_region_materials(
                    mesh, field, region_map, materials
                )
            except ValidationError as error:
                raise ModelSchemaError(
                    path=f"regions.{field}", reason=str(error)
                ) from error
        self.assignment = self.assignments[_FIELDS[self.physics]]
        self._document = raw
        self._boundary = _mapping(table.get("boundary", {}), "boundary", _BOUNDARY_KEYS)
        self._loads = _mapping(
            table.get("loads", {}), "loads", {"source", "body_force", "volume_charge"}
        )
        self._solver = _mapping(
            table.get("solver", {}),
            "solver",
            (
                {"block_scaling", "method", "electric_thermal_strain"}
                if self.physics == "thermopiezoelectric"
                else {"block_scaling", "method"}
                if self.physics == "piezoelectric"
                else set()
            ),
        )
        inspect_boundary_topology(mesh)
        self._conditions, self._thermal_conditions = self._parse_dirichlet()
        self._surface = self._parse_surfaces()
        self._electrodes = self._parse_electrodes()
        self._validate_loads()

    @classmethod
    def from_file(cls, path: str | Path) -> Model:
        """Load a UTF-8 TOML model and resolve its Gmsh path beside the file."""
        source = Path(path)
        if source.suffix.lower() != ".toml":
            raise ModelSchemaError(path="file", reason="expected .toml")
        try:
            with source.open("rb") as stream:
                document = tomllib.load(stream)
        except (OSError, tomllib.TOMLDecodeError) as error:
            raise ModelSchemaError(path="file", reason=str(error)) from error
        mesh_table = _mapping(
            _required(document, "mesh", "document"),
            "mesh",
            {"path", "length_scale_to_m"},
        )
        mesh_name = _required(mesh_table, "path", "mesh")
        if (
            not isinstance(mesh_name, str)
            or not mesh_name
            or Path(mesh_name).is_absolute()
        ):
            raise ModelSchemaError(
                path="mesh.path", reason="expected a relative Gmsh path"
            )
        root = source.resolve().parent
        mesh_path = (root / mesh_name).resolve()
        if not mesh_path.is_relative_to(root):
            raise ModelSchemaError(
                path="mesh.path", reason="path escapes model directory"
            )
        scale = _number(
            _required(mesh_table, "length_scale_to_m", "mesh"), "mesh.length_scale_to_m"
        )
        if scale <= 0:
            raise ModelSchemaError(
                path="mesh.length_scale_to_m", reason="must be positive"
            )
        return cls(
            mesh=read_gmsh(mesh_path, length_scale_to_m=scale), document=document
        )

    def to_document(self) -> dict[str, Any]:
        """Return a detached, JSON-safe copy for review or Python roundtrip."""
        return cast(
            dict[str, Any], json.loads(json.dumps(self._document, sort_keys=True))
        )

    def _dofs(self) -> DofMap:
        return {
            "heat": thermal_dofs,
            "elasticity": mechanical_dofs,
            "electrostatics": electric_dofs,
            "piezoelectric": piezoelectric_dofs,
        }[self.physics](self.mesh)

    def _nodes(self, entry: Mapping[str, Any], path: str) -> np.ndarray:
        if ("nodes" in entry) == ("set" in entry):
            raise ModelSchemaError(path=path, reason="choose nodes or facet set")
        if "set" in entry:
            name = entry["set"]
            if not isinstance(name, str) or name not in self.mesh.facet_sets:
                raise ModelSchemaError(
                    path=f"{path}.set", reason=f"unknown facet set {name!r}"
                )
            nodes = np.unique(self.mesh.triangle3[self.mesh.facet_sets[name]])
            if nodes.size == 0:
                raise ModelSchemaError(path=f"{path}.set", reason="empty facet set")
            return nodes
        value = entry["nodes"]
        array = np.asarray(value)
        if (
            array.ndim != 1
            or array.size == 0
            or array.dtype.kind not in "iu"
            or np.any(array < 0)
            or np.any(array >= self.mesh.n_points)
        ):
            raise ModelSchemaError(
                path=f"{path}.nodes",
                reason="expected active zero-based mesh node indices",
            )
        return array.astype(np.int64)

    def _parse_dirichlet(self) -> tuple[tuple[Dirichlet, ...], tuple[Dirichlet, ...]]:
        entries = self._boundary.get("dirichlet", [])
        if not isinstance(entries, list):
            raise ModelSchemaError(
                path="boundary.dirichlet", reason="expected array of tables"
            )
        conditions: list[Dirichlet] = []
        thermal_conditions: list[Dirichlet] = []
        dofs = None if self.physics == "thermopiezoelectric" else self._dofs()
        allowed_fields = (
            {"T", "u", "phi"}
            if self.physics == "thermopiezoelectric"
            else {"u", "phi"}
            if self.physics == "piezoelectric"
            else {_FIELDS[self.physics]}
        )
        for i, raw in enumerate(entries):
            path = f"boundary.dirichlet[{i}]"
            entry = _mapping(
                raw, path, {"field", "nodes", "set", "components", "value", "unit"}
            )
            field = _required(entry, "field", path)
            if not isinstance(field, str) or field not in allowed_fields:
                raise ModelSchemaError(
                    path=f"{path}.field", reason=f"inactive field {field!r}"
                )
            _unit(entry, {"u": "m", "T": "K", "phi": "V"}[field], path)
            nodes = self._nodes(entry, path)
            components = entry.get("components")
            if components is not None and (
                field != "u"
                or not isinstance(components, list)
                or not components
                or any(c not in ("x", "y", "z") for c in components)
            ):
                raise ModelSchemaError(
                    path=f"{path}.components",
                    reason="expected mechanical components x/y/z",
                )
            selected = None if components is None else tuple(components)
            try:
                field_dofs = (
                    (
                        thermal_dofs(self.mesh)
                        if field == "T"
                        else piezoelectric_dofs(self.mesh)
                    )
                    if self.physics == "thermopiezoelectric"
                    else dofs
                )
                assert field_dofs is not None
                ids = field_dofs.indices(field, nodes, selected)
            except ValidationError as error:
                raise ModelSchemaError(
                    path=f"{path}.nodes", reason=str(error)
                ) from error
            value = _required(entry, "value", path)
            per_node = (
                len(ids) // len(nodes)
                if len(nodes)
                else (len(selected) if selected else (3 if field == "u" else 1))
            )
            raw_value = np.asarray(value)
            if raw_value.ndim == 0:
                values = np.full(len(ids), _number(value, f"{path}.value"))
            elif raw_value.shape == (per_node,):
                values = np.tile(
                    _array(value, (per_node,), f"{path}.value"), len(nodes)
                )
            elif raw_value.shape == (len(ids),):
                values = _array(value, (len(ids),), f"{path}.value")
            else:
                raise ModelSchemaError(
                    path=f"{path}.value",
                    reason=f"expected scalar, ({per_node},) or ({len(ids)},)",
                )
            condition = Dirichlet(ids, values)
            if self.physics == "thermopiezoelectric" and field == "T":
                thermal_conditions.append(condition)
            else:
                conditions.append(condition)
        return tuple(conditions), tuple(thermal_conditions)

    def _facets(self, name: object, path: str) -> np.ndarray:
        if not isinstance(name, str) or name not in self.mesh.facet_sets:
            raise ModelSchemaError(path=path, reason=f"unknown facet set {name!r}")
        rows = self.mesh.facet_sets[name]
        if not len(rows):
            raise ModelSchemaError(path=path, reason="empty facet set")
        return rows

    def _parse_surfaces(self) -> dict[str, dict[int, Any]]:
        result: dict[str, dict[int, Any]] = {}
        contracts = {
            "flux": ({"heat", "thermopiezoelectric"}, "W/m^2"),
            "convection": ({"heat", "thermopiezoelectric"}, "W/(m^2 K)"),
            "traction": (
                {"elasticity", "piezoelectric", "thermopiezoelectric"},
                "N/m^2",
            ),
            "surface_charge": (
                {"electrostatics", "piezoelectric", "thermopiezoelectric"},
                "C/m^2",
            ),
        }
        for kind, (modes, unit) in contracts.items():
            entries = self._boundary.get(kind, [])
            if not isinstance(entries, list):
                raise ModelSchemaError(
                    path=f"boundary.{kind}", reason="expected array of tables"
                )
            if entries and self.physics not in modes:
                raise ModelSchemaError(
                    path=f"boundary.{kind}", reason="not valid for physics mode"
                )
            values: dict[int, Any] = {}
            for i, raw in enumerate(entries):
                path = f"boundary.{kind}[{i}]"
                allowed = (
                    {"set", "value", "unit"}
                    if kind != "convection"
                    else {
                        "set",
                        "coefficient",
                        "ambient_temperature",
                        "unit",
                        "ambient_unit",
                    }
                )
                entry = _mapping(raw, path, allowed)
                _unit(entry, unit, path)
                if kind == "convection":
                    if entry.get("ambient_unit") != "K":
                        raise ModelSchemaError(
                            path=f"{path}.ambient_unit", reason="expected 'K'"
                        )
                    value: Any = Convection(
                        _number(
                            _required(entry, "coefficient", path), f"{path}.coefficient"
                        ),
                        _number(
                            _required(entry, "ambient_temperature", path),
                            f"{path}.ambient_temperature",
                        ),
                    )
                elif kind == "traction":
                    value = _array(
                        _required(entry, "value", path), (3,), f"{path}.value"
                    )
                else:
                    value = _number(_required(entry, "value", path), f"{path}.value")
                for face in self._facets(_required(entry, "set", path), f"{path}.set"):
                    if int(face) in values:
                        raise ModelSchemaError(
                            path=path, reason=f"facet {face} assigned twice"
                        )
                    values[int(face)] = value
            result[kind] = values
        return result

    def _parse_electrodes(self) -> tuple[Electrode, ...]:
        entries = self._boundary.get("electrodes", [])
        if not isinstance(entries, list):
            raise ModelSchemaError(
                path="boundary.electrodes", reason="expected array of tables"
            )
        if entries and self.physics not in ("piezoelectric", "thermopiezoelectric"):
            raise ModelSchemaError(
                path="boundary.electrodes", reason="requires piezoelectric mode"
            )
        electrodes: list[Electrode] = []
        for i, raw in enumerate(entries):
            path = f"boundary.electrodes[{i}]"
            entry = _mapping(
                raw,
                path,
                {"name", "set", "potential", "potential_unit", "charge", "charge_unit"},
            )
            name = _required(entry, "name", path)
            if not isinstance(name, str) or not name.strip():
                raise ModelSchemaError(
                    path=f"{path}.name", reason="expected nonempty name"
                )
            face_set = _required(entry, "set", path)
            self._facets(face_set, f"{path}.set")
            if ("potential" in entry) == ("charge" in entry):
                raise ModelSchemaError(path=path, reason="choose potential or charge")
            if "potential" in entry:
                if entry.get("potential_unit") != "V" or "charge_unit" in entry:
                    raise ModelSchemaError(
                        path=f"{path}.potential_unit", reason="expected V"
                    )
                electrodes.append(
                    Electrode(
                        name,
                        face_set,
                        potential=_number(entry["potential"], f"{path}.potential"),
                    )
                )
            else:
                if entry.get("charge_unit") != "C" or "potential_unit" in entry:
                    raise ModelSchemaError(
                        path=f"{path}.charge_unit", reason="expected C"
                    )
                electrodes.append(
                    Electrode(
                        name,
                        face_set,
                        charge=_number(entry["charge"], f"{path}.charge"),
                    )
                )
        return tuple(electrodes)

    def _validate_loads(self) -> None:
        contracts = {
            "source": ({"heat", "thermopiezoelectric"}, "W/m^3", ()),
            "body_force": (
                {"elasticity", "piezoelectric", "thermopiezoelectric"},
                "N/m^3",
                (3,),
            ),
            "volume_charge": (
                {"electrostatics", "piezoelectric", "thermopiezoelectric"},
                "C/m^3",
                (),
            ),
        }
        for name, raw in self._loads.items():
            modes, unit, shape = contracts[name]
            if self.physics not in modes:
                raise ModelSchemaError(
                    path=f"loads.{name}", reason="not valid for physics mode"
                )
            table = _mapping(raw, f"loads.{name}", {"value", "unit"})
            _unit(table, unit, f"loads.{name}")
            _array(
                _required(table, "value", f"loads.{name}"), shape, f"loads.{name}.value"
            )
        if self.physics in ("piezoelectric", "thermopiezoelectric"):
            if "block_scaling" in self._solver and not isinstance(
                self._solver["block_scaling"], bool
            ):
                raise ModelSchemaError(
                    path="solver.block_scaling", reason="expected boolean"
                )
            if self._solver.get("method", "direct") != "direct":
                raise ModelSchemaError(
                    path="solver.method", reason="only direct is supported"
                )
        if (
            self.physics == "thermopiezoelectric"
            and "electric_thermal_strain" in self._solver
            and not isinstance(self._solver["electric_thermal_strain"], bool)
        ):
            raise ModelSchemaError(
                path="solver.electric_thermal_strain", reason="expected boolean"
            )

    def _load_value(self, name: str, default: Any) -> Any:
        return self._loads[name]["value"] if name in self._loads else default

    def solve(self) -> ModelResult:
        """Run the selected validated stationary reference solver."""
        mats = self.assignment.cell_materials
        if self.physics == "heat":
            result: FieldResult = solve_heat(
                self.mesh,
                cast(Sequence[ThermalMaterial], mats),
                self._conditions,
                source=self._load_value("source", 0.0),
                outward_fluxes=self._surface["flux"],
                convection=self._surface["convection"],
            )
        elif self.physics == "elasticity":
            result = solve_elasticity(
                self.mesh,
                cast(Sequence[ElasticMaterial], mats),
                self._conditions,
                body_force=self._load_value("body_force", (0.0, 0.0, 0.0)),
                tractions=self._surface["traction"],
            )
        elif self.physics == "electrostatics":
            result = solve_electrostatics(
                self.mesh,
                cast(Sequence[DielectricMaterial], mats),
                self._conditions,
                volume_charge=self._load_value("volume_charge", 0.0),
                surface_charges=self._surface["surface_charge"],
            )
        elif self.physics == "piezoelectric":
            result = solve_piezoelectric(
                self.mesh,
                cast(Sequence[PiezoelectricMaterial], mats),
                self._conditions,
                body_force=self._load_value("body_force", (0.0, 0.0, 0.0)),
                volume_charge=self._load_value("volume_charge", 0.0),
                tractions=self._surface["traction"],
                surface_charges=self._surface["surface_charge"],
                electrodes=self._electrodes,
                block_scaling=self._solver.get("block_scaling", True),
                method=self._solver.get("method", "direct"),
            )
        else:
            result = solve_thermopiezoelectric(
                self.mesh,
                cast(Sequence[ThermalMaterial], self.assignments["T"].cell_materials),
                cast(Sequence[PiezoelectricMaterial], mats),
                cast(
                    Sequence[ThermalExpansion],
                    self.assignments["expansion"].cell_materials,
                ),
                self._thermal_conditions,
                self._conditions,
                source=self._load_value("source", 0.0),
                outward_fluxes=self._surface["flux"],
                convection=self._surface["convection"],
                body_force=self._load_value("body_force", (0.0, 0.0, 0.0)),
                volume_charge=self._load_value("volume_charge", 0.0),
                tractions=self._surface["traction"],
                surface_charges=self._surface["surface_charge"],
                electrodes=self._electrodes,
                electric_thermal_strain=self._solver.get(
                    "electric_thermal_strain", False
                ),
                block_scaling=self._solver.get("block_scaling", True),
                method=self._solver.get("method", "direct"),
            )
        if isinstance(result, ThermopiezoelectricResult):
            solver_hashes: str | dict[str, str] = {
                "heat": result.heat.provenance.input_sha256,
                "piezoelectric": result.piezoelectric.provenance.input_sha256,
            }
            solver_names: str | list[str] = [
                result.heat.provenance.solver_name,
                result.piezoelectric.provenance.solver_name,
            ]
        else:
            solver_hashes = result.provenance.input_sha256
            solver_names = result.provenance.solver_name
        canonical_input = {
            "document": self._document,
            "mesh": {
                "points": self.mesh.points.tolist(),
                "tetra4": self.mesh.tetra4.tolist(),
                "triangle3": self.mesh.triangle3.tolist(),
            },
            "solver_input_sha256": solver_hashes,
        }
        if isinstance(result, ThermopiezoelectricResult):
            canonical_input["assignments"] = {
                field: assignment.to_provenance_dict()
                for field, assignment in sorted(self.assignments.items())
            }
        else:
            canonical_input["assignment"] = self.assignment.to_provenance_dict()
        canonical = json.dumps(
            canonical_input,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        provenance = RunProvenance.from_input_bytes(
            canonical,
            solver_name=f"fem_toolbox.model.Model.solve.{self.physics}",
            solver_parameters={
                "model_schema_version": 1,
                "units": "SI",
                "field_solver": solver_names,
                "field_solver_input_sha256": solver_hashes,
            },
        )
        return ModelResult(self.physics, result, provenance)
