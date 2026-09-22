"""Boundary topology has independently known incidence and outward normals."""

from pathlib import Path

import meshio
import numpy as np
import pytest

from fem_toolbox.errors import BoundaryTopologyError, ConnectivityError
from fem_toolbox.io import read_gmsh
from fem_toolbox.mesh import Mesh, inspect_boundary_topology

POINTS = np.array(
    [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [2, 2, 2]],
    dtype=float,
)


def mesh_with_faces(
    faces: list[list[int]], cells: list[list[int]] | None = None
) -> Mesh:
    return Mesh(
        points=POINTS,
        tetra4=np.array(cells or [[0, 1, 2, 3]], dtype=int),
        triangle3=np.array(faces, dtype=int).reshape(-1, 3),
        tetra4_ids=np.arange(101, 101 + len(cells or [[0, 1, 2, 3]])),
        triangle3_ids=np.arange(201, 201 + len(faces)),
        facet_sets={
            "load": np.arange(len(faces)),
            "electrode": np.array([0] if faces else [], dtype=int),
        },
    )


def test_boundary_incidence_and_outward_geometry_ignore_node_order() -> None:
    for cell in ([0, 1, 2, 3], [0, 2, 1, 3]):
        topology = inspect_boundary_topology(mesh_with_faces([[0, 1, 2]], [cell]))
        assert len(topology.exterior_faces) == 4
        assert topology.facet_incidence[(0, 1, 2)] == (101,)
        face = topology.triangle_faces[0]
        assert face.cell_id == 101
        assert face.nodes == (0, 2, 1)
        np.testing.assert_allclose(face.unit_normal, [0, 0, -1])
        assert face.area == pytest.approx(0.5)


def test_overlapping_groups_are_preserved() -> None:
    topology = inspect_boundary_topology(mesh_with_faces([[0, 2, 1]]))
    assert topology.triangle_faces[0].groups == ("electrode", "load")


def test_shared_face_has_two_stable_owner_ids() -> None:
    mesh = mesh_with_faces([], [[0, 1, 2, 3], [0, 2, 1, 4]])
    topology = inspect_boundary_topology(mesh)
    assert topology.facet_incidence[(0, 1, 2)] == (101, 102)
    assert len(topology.exterior_faces) == 6


@pytest.mark.parametrize(
    ("faces", "cells", "kind"),
    [
        ([[0, 1, 2]], [[0, 1, 2, 3], [0, 2, 1, 4]], "interior"),
        ([[0, 1, 5]], None, "foreign"),
        (
            [[0, 1, 2]],
            [[0, 1, 2, 3], [0, 2, 1, 4], [0, 1, 2, 5]],
            "nonmanifold",
        ),
    ],
)
def test_invalid_triangle_reports_group_and_entity(
    faces: list[list[int]], cells: list[list[int]] | None, kind: str
) -> None:
    with pytest.raises(BoundaryTopologyError) as caught:
        inspect_boundary_topology(mesh_with_faces(faces, cells))
    assert caught.value.kind == kind
    assert caught.value.entity_id == 201
    assert "load" in str(caught.value)


def test_missing_exterior_triangle_is_reported_when_complete_cover_requested() -> None:
    with pytest.raises(BoundaryTopologyError) as caught:
        inspect_boundary_topology(mesh_with_faces([[0, 1, 2]]), require_complete=True)
    assert caught.value.kind == "missing"
    assert caught.value.cell_id == 101


def test_duplicate_surface_use_is_reported_even_with_distinct_semantic_groups() -> None:
    with pytest.raises(ConnectivityError) as caught:
        mesh_with_faces([[0, 1, 2], [0, 2, 1]])
    assert caught.value.entity_id == 202
    assert caught.value.duplicate_of == 201
    assert "load" in str(caught.value)


def test_existing_gmsh41_fixture_keeps_overlapping_groups() -> None:
    path = Path(__file__).parents[1] / "fixtures" / "overlapping_groups.msh"
    topology = inspect_boundary_topology(read_gmsh(path))
    assert topology.triangle_faces[0].groups == ("electrode", "ground")


def test_gmsh22_surface_preflight_reports_foreign_group(tmp_path: Path) -> None:
    path = tmp_path / "foreign.msh"
    meshio.write(
        path,
        meshio.Mesh(
            points=POINTS,
            cells=[
                ("tetra", np.array([[0, 1, 2, 3]])),
                ("triangle", np.array([[0, 1, 5]])),
            ],
            cell_data={"gmsh:physical": [np.array([1]), np.array([2])]},
            field_data={"load": np.array([2, 2])},
        ),
        file_format="gmsh22",
        binary=False,
    )
    with pytest.raises(BoundaryTopologyError) as caught:
        inspect_boundary_topology(read_gmsh(path))
    assert caught.value.kind == "foreign"
    assert "load" in str(caught.value)
