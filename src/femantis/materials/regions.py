"""Resolve named volume sets to typed, traceable cell materials."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from typing import Any, TypeAlias

from femantis.errors import MaterialAssignmentError
from femantis.materials import ElasticMaterial, ThermalExpansion
from femantis.materials.electrostatics import DielectricMaterial
from femantis.materials.heat import ThermalMaterial
from femantis.materials.piezoelectric import PiezoelectricMaterial
from femantis.mesh import Mesh

Material: TypeAlias = (
    ThermalMaterial
    | ElasticMaterial
    | DielectricMaterial
    | PiezoelectricMaterial
    | ThermalExpansion
)

_FIELD_TYPES: dict[str, type[Material]] = {
    "T": ThermalMaterial,
    "u": ElasticMaterial,
    "phi": DielectricMaterial,
    "piezo": PiezoelectricMaterial,
    "expansion": ThermalExpansion,
}


def _material_data(material: Material) -> dict[str, Any]:
    if isinstance(material, ThermalMaterial):
        data: dict[str, Any] = {
            "conductivity": material.conductivity.tolist(),
            "unit": "W/(m K)",
        }
        if material.density is not None:
            data.update(
                density=material.density,
                specific_heat=material.specific_heat,
                density_unit="kg/m^3",
                specific_heat_unit="J/(kg K)",
            )
        if material.source is not None:
            data["source"] = asdict(material.source)
        return data
    if isinstance(material, DielectricMaterial):
        return {"permittivity": material.permittivity.tolist(), "unit": "F/m"}
    if isinstance(material, ElasticMaterial):
        return {"stiffness": material.stiffness.tolist(), "unit": "Pa"}
    if isinstance(material, PiezoelectricMaterial):
        return material.to_provenance_dict()
    return {
        "alpha": material.alpha.tolist(),
        "reference_temperature": material.reference_temperature,
        "axes_local_to_global": (
            None if material.source_axes is None else material.source_axes.tolist()
        ),
        "units": {"alpha": "1/K", "reference_temperature": "K"},
    }


@dataclass(frozen=True, slots=True)
class MaterialAssignment:
    """Per-cell materials and provenance in input volume-cell row order."""

    field: str
    cell_ids: tuple[int, ...]
    cell_materials: tuple[Material, ...]
    material_names: tuple[str, ...]
    region_names: tuple[tuple[str, ...], ...]

    def to_provenance_dict(self) -> dict[str, Any]:
        """Return JSON-safe effective cell mapping and SI material tensors."""
        catalog = {
            name: _material_data(material)
            for name, material in sorted(
                zip(self.material_names, self.cell_materials, strict=True)
            )
        }
        return {
            "field": self.field,
            "cells": [
                {"id": cell_id, "material": name, "regions": list(regions)}
                for cell_id, name, regions in zip(
                    self.cell_ids, self.material_names, self.region_names, strict=True
                )
            ],
            "catalog": catalog,
        }


def assign_region_materials(
    mesh: Mesh,
    field: str,
    regions: Mapping[str, str],
    materials: Mapping[str, Material],
) -> MaterialAssignment:
    """Map named volume sets onto all volume cells before field assembly.

    ``field`` is ``T``, ``u``, ``phi``, ``piezo`` or ``expansion``. Every
    region names one catalog entry of the matching immutable SI material type.
    Overlapping sets are allowed only when they name the same catalog entry.
    Unknown sets/materials, uncovered cells, conflicts and unused definitions
    raise :class:`MaterialAssignmentError` with the external volume-cell ID where
    applicable. Existing solvers accept ``cell_materials`` directly.
    """
    expected = _FIELD_TYPES.get(field)
    if expected is None:
        raise MaterialAssignmentError(
            kind="field", field=field, reason="unsupported active field"
        )
    if mesh.n_cells == 0:
        raise MaterialAssignmentError(
            kind="empty_mesh", field=field, reason="no active volume cells"
        )
    if not isinstance(regions, Mapping):
        raise MaterialAssignmentError(
            kind="regions_type", field=field, reason="expected a region mapping"
        )
    if not isinstance(materials, Mapping):
        raise MaterialAssignmentError(
            kind="materials_type", field=field, reason="expected a material mapping"
        )
    for name, material in materials.items():
        if not isinstance(name, str) or not name or name.strip() != name:
            raise MaterialAssignmentError(
                kind="material_name", field=field, reason=f"invalid name {name!r}"
            )
        if not isinstance(material, expected):
            raise MaterialAssignmentError(
                kind="material_type",
                field=field,
                reason=f"{name!r} must be {expected.__name__}",
            )
    for region, name in regions.items():
        if region not in mesh.volume_sets:
            raise MaterialAssignmentError(
                kind="unknown_region",
                field=field,
                reason=f"unknown set {region!r}",
                regions=(region,),
            )
        if name not in materials:
            raise MaterialAssignmentError(
                kind="unknown_material",
                field=field,
                reason=f"unknown material {name!r}",
                regions=(region,),
            )

    per_cell: list[list[tuple[str, str]]] = [[] for _ in range(mesh.n_cells)]
    for region, name in sorted(regions.items()):
        for cell_index in mesh.volume_sets[region]:
            per_cell[int(cell_index)].append((region, name))

    cell_ids: list[int] = []
    cell_materials: list[Material] = []
    material_names: list[str] = []
    region_names: list[tuple[str, ...]] = []
    for cell_index, entries in enumerate(per_cell):
        cell_id = (
            cell_index if mesh.cell_ids is None else int(mesh.cell_ids[cell_index])
        )
        if not entries:
            raise MaterialAssignmentError(
                kind="missing",
                field=field,
                reason="cell has no material",
                cell_id=cell_id,
            )
        names = {name for _, name in entries}
        if len(names) != 1:
            raise MaterialAssignmentError(
                kind="conflict",
                field=field,
                reason=f"different materials {sorted(names)}",
                cell_id=cell_id,
                regions=tuple(region for region, _ in entries),
            )
        name = entries[0][1]
        cell_ids.append(cell_id)
        cell_materials.append(materials[name])
        material_names.append(name)
        region_names.append(tuple(region for region, _ in entries))

    unused = sorted(set(materials) - set(material_names))
    if unused:
        raise MaterialAssignmentError(
            kind="unused_material", field=field, reason=f"unused definitions {unused}"
        )
    return MaterialAssignment(
        field,
        tuple(cell_ids),
        tuple(cell_materials),
        tuple(material_names),
        tuple(region_names),
    )
