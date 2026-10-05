"""Gmsh and VTU adapters implemented at the meshio infrastructure boundary."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from pathlib import Path
from typing import Any

import meshio  # type: ignore[import-untyped]
import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.elements import Tet4Geometry
from femantis.errors import MeshIOError, UnsupportedCellTypeError
from femantis.mesh import Mesh
from femantis.mesh.quadratic import reordered_tet10, require_straight

_INTERNAL_TO_MESHIO = {
    "tetra4": "tetra",
    "triangle3": "triangle",
    "tetra10": "tetra10",
    "triangle6": "triangle6",
}
_SUPPORTED_MESHIO_TYPES = frozenset(_INTERNAL_TO_MESHIO.values())
_RESERVED_DATA_PREFIX = "fem_toolbox:"
_MESHIO_METADATA_CELL_SETS = frozenset({"gmsh:bounding_entities"})


def _concatenated_cells(
    source: meshio.Mesh, cell_type: str, width: int
) -> NDArray[np.int64]:
    cells = source.cells_dict.get(cell_type)
    if cells is None:
        return np.empty((0, width), dtype=np.int64)
    return np.asarray(cells, dtype=np.int64)


def _physical_tags(source: meshio.Mesh, cell_type: str) -> NDArray[np.int64] | None:
    values = source.cell_data_dict.get("gmsh:physical", {}).get(cell_type)
    if values is None:
        return None
    return np.asarray(values, dtype=np.int64)


def _selection_sets(
    source: meshio.Mesh,
    *,
    cell_type: str,
    dimension: int,
    physical_tags: NDArray[np.int64] | None,
) -> dict[str, NDArray[np.int64]]:
    selections: dict[str, NDArray[np.int64]] = {}

    for name, cell_types in source.cell_sets_dict.items():
        if name in _MESHIO_METADATA_CELL_SETS:
            continue
        indices = cell_types.get(cell_type)
        if indices is not None:
            selections[name] = np.unique(np.asarray(indices, dtype=np.int64))

    if physical_tags is None:
        return selections

    for name, definition in source.field_data.items():
        data = np.asarray(definition)
        if data.shape != (2,) or int(data[1]) != dimension:
            continue
        indices = np.flatnonzero(physical_tags == int(data[0])).astype(np.int64)
        if not indices.size:
            continue
        existing = selections.get(name)
        selections[name] = (
            indices if existing is None else np.union1d(existing, indices)
        )

    return selections


def _orient_tetra4(
    points: NDArray[np.float64], connectivity: NDArray[np.int64]
) -> NDArray[np.int64]:
    oriented = connectivity.copy()
    for index, nodes in enumerate(connectivity):
        geometry = Tet4Geometry(points[nodes], entity_id=index)
        oriented[index] = geometry.oriented_connectivity(nodes)
    return oriented


def _metadata(
    source: meshio.Mesh, cell_type: str, name: str
) -> NDArray[np.int64] | None:
    data = source.cell_data_dict.get(f"fem_toolbox:{name}", {}).get(cell_type)
    present = source.cell_data_dict.get(f"fem_toolbox:has_{name}", {}).get(cell_type)
    if present is not None:
        if np.all(present == 0):
            return None
        if not np.all(present == 1) or data is None:
            raise MeshIOError(
                operation="import", reason=f"invalid presence metadata for {name}"
            )
    return None if data is None else np.asarray(data)


def _vtu_sets(source: meshio.Mesh, cell_type: str) -> dict[str, NDArray[np.int64]]:
    result = {}
    for key, by_type in source.cell_data_dict.items():
        if key.startswith("fem_toolbox:set:") and cell_type in by_type:
            name = key.removeprefix("fem_toolbox:set:")
            mask = np.asarray(by_type[cell_type])
            if mask.ndim != 1 or not np.isin(mask, [0, 1]).all():
                raise MeshIOError(
                    operation="import", reason=f"invalid selection mask {name}"
                )
            presence = source.cell_data_dict.get(
                f"fem_toolbox:set-present:{name}", {}
            ).get(cell_type)
            if presence is not None:
                if np.all(presence == 0):
                    if np.any(mask):
                        raise MeshIOError(
                            operation="import",
                            reason=f"absent selection {name} has members",
                        )
                    continue
                if not np.all(presence == 1):
                    raise MeshIOError(
                        operation="import", reason=f"invalid selection presence {name}"
                    )
            elif not np.any(mask):
                continue
            result[name] = np.flatnonzero(mask).astype(np.int64)
    return result


def _mesh_from_meshio(
    source: meshio.Mesh, *, length_scale_to_m: float, vtu: bool = False
) -> Mesh:
    unsupported = tuple(
        sorted(
            cell_type
            for cell_type, cells in source.cells_dict.items()
            if len(cells) and cell_type not in _SUPPORTED_MESHIO_TYPES
        )
    )
    if unsupported:
        raise UnsupportedCellTypeError(unsupported)
    points = np.asarray(source.points, dtype=np.float64)
    if points.ndim != 2 or points.shape[1:] != (3,):
        raise MeshIOError(
            operation="import",
            reason=f"point coordinates have shape {points.shape}; expected (n, 3)",
        )
    data: dict[str, Any] = {"points": points * length_scale_to_m}
    for internal, external in _INTERNAL_TO_MESHIO.items():
        width = {"tetra4": 4, "triangle3": 3, "tetra10": 10, "triangle6": 6}[internal]
        data[internal] = _concatenated_cells(source, external, width)
        if vtu:
            data[f"{internal}_ids"] = _metadata(source, external, "external_id")
    # Validate family, row indices and shared edge incidence before indexing points.
    preliminary = Mesh(**data)
    volume_type = _INTERNAL_TO_MESHIO[preliminary.cell_type]
    facet_type = _INTERNAL_TO_MESHIO[preliminary.facet_type]
    if vtu:
        data["point_ids"] = source.point_data.get("fem_toolbox:external_id")
        data["cell_tags"] = _metadata(source, volume_type, "physical_tag")
        data["facet_tags"] = _metadata(source, facet_type, "physical_tag")
        data["volume_sets"] = _vtu_sets(source, volume_type)
        data["facet_sets"] = _vtu_sets(source, facet_type)
    else:
        data["cell_tags"] = _physical_tags(source, volume_type)
        data["facet_tags"] = _physical_tags(source, facet_type)
        data["volume_sets"] = _selection_sets(
            source, cell_type=volume_type, dimension=3, physical_tags=data["cell_tags"]
        )
        data["facet_sets"] = _selection_sets(
            source, cell_type=facet_type, dimension=2, physical_tags=data["facet_tags"]
        )
        data["tetra4"] = _orient_tetra4(preliminary.points, preliminary.tetra4)
    if preliminary.n_tetra10:
        oriented = preliminary.tetra10.copy()
        for index, nodes in enumerate(preliminary.tetra10):
            require_straight(preliminary.points, nodes)
            geometry = Tet4Geometry(preliminary.points[nodes[:4]], entity_id=index)
            if not vtu:
                oriented[index] = reordered_tet10(
                    nodes, geometry.oriented_connectivity(nodes[:4])
                )
        data["tetra10"] = oriented
    for nodes in preliminary.triangle6:
        require_straight(preliminary.points, nodes)
    return Mesh(**data)


def read_gmsh(path: str | Path, *, length_scale_to_m: float = 1.0) -> Mesh:
    """Read one Tet4/Triangle3 or straight Tet10/Triangle6 Gmsh family.

    Coordinates are multiplied by ``length_scale_to_m`` and stored in metres.
    Physical names become immutable volume or facet selection sets. Gmsh 4.1
    entities belonging to several physical groups remain in every corresponding
    set, while the singular compatibility tags retain meshio's first physical tag.
    Negative volume orientations are normalized, including all Tet10 edge nodes.
    Curved quadratic geometry and degenerate volume cells are rejected.
    """
    if not isfinite(length_scale_to_m) or length_scale_to_m <= 0.0:
        raise MeshIOError(
            operation="import",
            reason="length_scale_to_m must be a positive finite scale factor",
        )
    source_path = Path(path)
    try:
        source = meshio.read(source_path, file_format="gmsh")
    except (OSError, meshio.ReadError) as error:
        raise MeshIOError(operation="import", reason=str(error)) from error
    return _mesh_from_meshio(source, length_scale_to_m=length_scale_to_m)


def read_vtu(path: str | Path) -> Mesh:
    """Read mesh topology and Toolbox metadata from VTU, without changing ordering.

    Restores optional IDs, tags and overlapping/empty sets written by write_vtu.
    Result arrays remain accessible through meshio; this function returns Mesh.
    """
    try:
        source = meshio.read(Path(path), file_format="vtu")
    except (OSError, meshio.ReadError) as error:
        raise MeshIOError(operation="import", reason=str(error)) from error
    return _mesh_from_meshio(source, length_scale_to_m=1.0, vtu=True)


def _result_array(
    value: ArrayLike,
    *,
    field: str,
    expected_rows: int,
) -> NDArray[np.generic]:
    array = np.asarray(value)
    if array.ndim == 0 or array.shape[0] != expected_rows:
        raise MeshIOError(
            operation="export",
            reason=(
                f"data field {field!r} has shape {array.shape}; expected first "
                f"dimension {expected_rows}"
            ),
        )
    if array.dtype.kind not in "biufc":
        raise MeshIOError(
            operation="export",
            reason=f"data field {field!r} must contain numeric values",
        )
    if not np.isfinite(array).all():
        raise MeshIOError(
            operation="export",
            reason=f"data field {field!r} contains non-finite values",
        )
    return array


def _active_cell_types(mesh: Mesh) -> list[tuple[str, str, NDArray[np.int64]]]:
    return [
        (internal, external, getattr(mesh, internal))
        for internal, external in _INTERNAL_TO_MESHIO.items()
        if len(getattr(mesh, internal))
    ]


def _export_cell_data(
    mesh: Mesh,
    active: list[tuple[str, str, NDArray[np.int64]]],
    supplied: Mapping[str, Mapping[str, ArrayLike]],
) -> dict[str, list[NDArray[np.generic]]]:
    result: dict[str, list[NDArray[np.generic]]] = {}
    expected_types = {internal for internal, _, _ in active}

    for name, values_by_type in supplied.items():
        if name.startswith(_RESERVED_DATA_PREFIX):
            raise MeshIOError(
                operation="export",
                reason=f"data field {name!r} uses the reserved fem_toolbox prefix",
            )
        actual_types = set(values_by_type)
        if actual_types != expected_types:
            missing = sorted(expected_types - actual_types)
            unknown = sorted(actual_types - expected_types)
            details = []
            if missing:
                details.append(f"missing {', '.join(missing)}")
            if unknown:
                details.append(f"unknown {', '.join(unknown)}")
            raise MeshIOError(
                operation="export",
                reason=f"cell data field {name!r}: {'; '.join(details)}",
            )
        arrays = [
            _result_array(
                values_by_type[internal], field=name, expected_rows=len(connectivity)
            )
            for internal, _, connectivity in active
        ]
        trailing_shapes = {array.shape[1:] for array in arrays}
        if len(trailing_shapes) != 1:
            raise MeshIOError(
                operation="export",
                reason=f"cell data field {name!r} has inconsistent component shapes",
            )
        result[name] = arrays

    selection_names = sorted(set(mesh.volume_sets) | set(mesh.facet_sets))
    for name in selection_names:
        arrays = []
        for internal, _, connectivity in active:
            mask = np.zeros(len(connectivity), dtype=np.uint8)
            selection = (
                mesh.volume_sets.get(name)
                if internal.startswith("tetra")
                else mesh.facet_sets.get(name)
            )
            if selection is not None:
                mask[selection] = 1
            arrays.append(mask)
        result[f"{_RESERVED_DATA_PREFIX}set:{name}"] = arrays
        result[f"{_RESERVED_DATA_PREFIX}set-present:{name}"] = [
            np.full(
                len(connectivity),
                int(
                    name
                    in (
                        mesh.volume_sets
                        if internal.startswith("tetra")
                        else mesh.facet_sets
                    )
                ),
                dtype=np.uint8,
            )
            for internal, _, connectivity in active
        ]

    for name, cell_value, facet_value in (
        ("physical_tag", mesh.cell_tags, mesh.facet_tags),
        ("external_id", mesh.cell_ids, mesh.facet_ids),
    ):
        if cell_value is None and facet_value is None:
            continue
        values: list[NDArray[np.generic]] = []
        presence: list[NDArray[np.generic]] = []
        for internal, _, connectivity in active:
            value = cell_value if internal.startswith("tetra") else facet_value
            values.append(
                np.zeros(len(connectivity), dtype=np.int64) if value is None else value
            )
            presence.append(
                np.full(len(connectivity), int(value is not None), dtype=np.uint8)
            )
        result[f"fem_toolbox:{name}"] = values
        result[f"fem_toolbox:has_{name}"] = presence

    return result


def write_vtu(
    path: str | Path,
    mesh: Mesh,
    *,
    point_data: Mapping[str, ArrayLike] | None = None,
    cell_data: Mapping[str, Mapping[str, ArrayLike]] | None = None,
    binary: bool = True,
) -> None:
    """Write a mesh and numeric result arrays to an unstructured VTK file.

    Every cell result field must provide arrays for all non-empty internal cell
    types, keyed by the active ``tetra4``/``triangle3`` or ``tetra10``/``triangle6`` family. Selection sets are
    exported as independent ``fem_toolbox:set:<name>`` masks, so overlapping
    memberships remain lossless in the VTU output.
    """
    active = _active_cell_types(mesh)
    if not active:
        raise MeshIOError(operation="export", reason="mesh contains no cells")

    if any(name.startswith(_RESERVED_DATA_PREFIX) for name in (point_data or {})):
        raise MeshIOError(
            operation="export", reason="point data uses the reserved fem_toolbox prefix"
        )
    exported_point_data = {
        name: _result_array(value, field=name, expected_rows=mesh.n_points)
        for name, value in (point_data or {}).items()
    }
    if mesh.point_ids is not None:
        exported_point_data["fem_toolbox:external_id"] = mesh.point_ids
    exported_cell_data = _export_cell_data(mesh, active, cell_data or {})
    output = meshio.Mesh(
        points=mesh.points,
        cells=[(external, connectivity) for _, external, connectivity in active],
        point_data=exported_point_data,
        cell_data=exported_cell_data,
    )
    try:
        meshio.write(Path(path), output, file_format="vtu", binary=binary)
    except (OSError, meshio.WriteError, ValueError) as error:
        raise MeshIOError(operation="export", reason=str(error)) from error
