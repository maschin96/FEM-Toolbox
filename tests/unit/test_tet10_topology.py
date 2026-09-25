"""Quadratic topology, canonical edge order and explicit compatibility limits."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from fem_toolbox.errors import (
    BoundaryTopologyError,
    MeshValidationError,
    ValidationError,
)
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.materials.regions import assign_region_materials
from fem_toolbox.mesh import Mesh, inspect_boundary_topology

# VTK/meshio corners 0..3, then edges 01,12,20,03,13,23.
POINTS = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [0.5, 0, 0],
        [0.5, 0.5, 0],
        [0, 0.5, 0],
        [0, 0, 0.5],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
    ]
)
CELL = np.arange(10).reshape(1, 10)
FACES = np.array(
    [[0, 2, 1, 6, 5, 4], [0, 1, 3, 4, 8, 7], [0, 3, 2, 7, 9, 6], [1, 2, 3, 5, 9, 8]]
)


def quadratic(**kwargs):
    return Mesh(points=POINTS, tetra10=CELL, **kwargs)


def test_quadratic_storage_metadata_and_legacy_empty_shapes():
    nodes = CELL.copy()
    mesh = Mesh(
        points=POINTS,
        tetra10=nodes,
        tetra10_ids=[101],
        triangle6=FACES,
        triangle6_ids=[201, 202, 203, 204],
        cell_tags=[8],
        facet_tags=[1, 2, 3, 4],
        volume_sets={"solid": [0]},
        facet_sets={"port": [0, 1]},
    )
    nodes[:] = 0
    assert_array_equal(mesh.tetra10, CELL)
    assert mesh.tetra4.shape == (0, 4)
    assert mesh.triangle3.shape == (0, 3)
    assert mesh.n_tetra10 == mesh.n_cells == 1
    assert mesh.n_triangle6 == mesh.n_facets == 4
    assert mesh.cell_type == "tetra10"
    assert_array_equal(mesh.cells, CELL)
    assert_array_equal(mesh.cell_ids, [101])
    with pytest.raises(ValueError):
        mesh.tetra10[0, 0] = 7
    assert Mesh.empty().tetra10.shape == (0, 10)


@pytest.mark.parametrize(
    "extra", [dict(tetra4=[[0, 1, 2, 3]]), dict(triangle3=[[0, 1, 2]])]
)
def test_hybrid_families_are_explicitly_rejected(extra):
    with pytest.raises(MeshValidationError, match="hybrid"):
        quadratic(**extra)


def test_all_six_face_nodes_are_oriented_and_assigned():
    mesh = quadratic(
        triangle6=FACES,
        tetra10_ids=[101],
        triangle6_ids=[201, 202, 203, 204],
        facet_sets={"port": [0], "ground": [0]},
    )
    topology = inspect_boundary_topology(mesh, require_complete=True)
    face = topology.triangle_faces[0]
    assert face.nodes == (0, 2, 1, 6, 5, 4)
    assert face.cell_id == 101 and face.triangle_id == 201
    assert face.groups == ("ground", "port")
    assert_allclose(face.unit_normal, [0, 0, -1], rtol=0, atol=1e-14)
    assert_allclose(face.area, 0.5, rtol=0, atol=1e-14)
    assert len(topology.exterior_faces) == 4


def test_wrong_edge_assignment_on_surface_rejected():
    wrong = FACES[:1].copy()
    wrong[0, [3, 4]] = wrong[0, [4, 3]]
    with pytest.raises(BoundaryTopologyError, match="nonconforming"):
        inspect_boundary_topology(quadratic(triangle6=wrong))


def test_material_regions_address_quadratic_cell_rows():
    mesh = quadratic(tetra10_ids=[91], volume_sets={"bulk": [0], "overlap": [0]})
    material = ThermalMaterial.isotropic(2)
    result = assign_region_materials(
        mesh, "T", {"bulk": "a", "overlap": "a"}, {"a": material}
    )
    assert result.cell_ids == (91,)
    assert result.cell_materials == (material,)
    assert result.region_names == (("bulk", "overlap"),)


def test_tet4_solver_does_not_silently_drop_quadratic_cells():
    from fem_toolbox.solvers.heat import solve_heat

    with pytest.raises(ValidationError, match="requires Tet4"):
        solve_heat(quadratic(), ThermalMaterial.isotropic(2), [])


def test_model_rejects_quadratic_mesh_before_tet4_only_processing():
    from fem_toolbox.errors import ModelSchemaError
    from fem_toolbox.model import Model

    with pytest.raises(ModelSchemaError, match=r"Tet10.*not yet supported"):
        Model(mesh=quadratic(), document={})


def test_shared_quadratic_face_and_nonconforming_edge():
    from fem_toolbox.errors import ConnectivityError

    points = np.vstack(
        [
            POINTS,
            [[0, 0, -1], [0, 0, -0.5], [0.5, 0, -0.5], [0, 0.5, -0.5], [0.5, 0, 0]],
        ]
    )
    cells = np.vstack([CELL, [0, 2, 1, 10, 6, 5, 4, 11, 13, 12]])
    mesh = Mesh(points=points, tetra10=cells, tetra10_ids=[101, 102])
    topology = inspect_boundary_topology(mesh)
    assert len(topology.exterior_faces) == 6
    with pytest.raises(BoundaryTopologyError, match="interior"):
        inspect_boundary_topology(
            Mesh(points=points, tetra10=cells, triangle6=FACES[:1])
        )
    cells[1, 6] = 14
    with pytest.raises(ConnectivityError, match="nonconforming"):
        Mesh(points=points, tetra10=cells)


def test_curved_geometry_cannot_use_corner_geometry():
    points = POINTS.copy()
    points[4, 2] += 0.01
    with pytest.raises(MeshValidationError, match="curved"):
        inspect_boundary_topology(Mesh(points=points, tetra10=CELL))
