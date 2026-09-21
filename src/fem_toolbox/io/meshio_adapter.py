"""Gmsh and VTU adapters implemented at the meshio infrastructure boundary."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from pathlib import Path

import meshio  # type: ignore[import-untyped]
import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import MeshIOError, UnsupportedCellTypeError
from fem_toolbox.mesh import Mesh

_INTERNAL_TO_MESHIO = {"tetra4": "tetra", "triangle3": "triangle"}
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


def _mesh_from_meshio(source: meshio.Mesh, *, length_scale_to_m: float) -> Mesh:
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
    tetra4 = _concatenated_cells(source, "tetra", 4)
    triangle3 = _concatenated_cells(source, "triangle", 3)
    cell_tags = _physical_tags(source, "tetra")
    facet_tags = _physical_tags(source, "triangle")

    scaled_points = points * length_scale_to_m
    oriented_tetra4 = _orient_tetra4(scaled_points, tetra4)
    return Mesh(
        points=scaled_points,
        tetra4=oriented_tetra4,
        triangle3=triangle3,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        volume_sets=_selection_sets(
            source,
            cell_type="tetra",
            dimension=3,
            physical_tags=cell_tags,
        ),
        facet_sets=_selection_sets(
            source,
            cell_type="triangle",
            dimension=2,
            physical_tags=facet_tags,
        ),
    )


def read_gmsh(path: str | Path, *, length_scale_to_m: float = 1.0) -> Mesh:
    """Read a Gmsh mesh containing only Tet4 and Triangle3 entities.

    Coordinates are multiplied by ``length_scale_to_m`` and stored in metres.
    Physical names become immutable volume or facet selection sets. Gmsh 4.1
    entities belonging to several physical groups remain in every corresponding
    set, while the singular compatibility tags retain meshio's first physical tag.
    Negative Tet4 orientations are normalized and degenerate elements are rejected.
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
    active: list[tuple[str, str, NDArray[np.int64]]] = []
    if mesh.n_tetra4:
        active.append(("tetra4", "tetra", mesh.tetra4))
    if mesh.n_triangle3:
        active.append(("triangle3", "triangle", mesh.triangle3))
    return active


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
                if internal == "tetra4"
                else mesh.facet_sets.get(name)
            )
            if selection is not None:
                mask[selection] = 1
            arrays.append(mask)
        result[f"{_RESERVED_DATA_PREFIX}set:{name}"] = arrays

    if mesh.cell_tags is not None or mesh.facet_tags is not None:
        tag_arrays: list[NDArray[np.generic]] = []
        for internal, _, connectivity in active:
            tags = mesh.cell_tags if internal == "tetra4" else mesh.facet_tags
            tag_arrays.append(
                np.zeros(len(connectivity), dtype=np.int64) if tags is None else tags
            )
        result[f"{_RESERVED_DATA_PREFIX}physical_tag"] = tag_arrays

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
    types, keyed by ``"tetra4"`` and/or ``"triangle3"``. Selection sets are
    exported as independent ``fem_toolbox:set:<name>`` masks, so overlapping
    memberships remain lossless in the VTU output.
    """
    active = _active_cell_types(mesh)
    if not active:
        raise MeshIOError(operation="export", reason="mesh contains no cells")

    exported_point_data = {
        name: _result_array(value, field=name, expected_rows=mesh.n_points)
        for name, value in (point_data or {}).items()
    }
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
