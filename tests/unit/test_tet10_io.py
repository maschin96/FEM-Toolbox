"""Independent Gmsh ordering and VTU mesh-metadata round trips."""

from pathlib import Path

import meshio
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from fem_toolbox.io import read_gmsh, write_vtu
from fem_toolbox.mesh import Mesh, inspect_boundary_topology

FIXTURE = Path(__file__).parents[1] / "fixtures" / "quadratic_groups.msh"


def test_gmsh41_node_order_overlapping_groups_and_vtu_roundtrip(tmp_path):
    from fem_toolbox.io import read_vtu

    mesh = read_gmsh(FIXTURE)
    assert_array_equal(mesh.tetra10, [np.arange(10)])
    assert_array_equal(mesh.triangle6, [[0, 2, 1, 6, 5, 4]])
    assert set(mesh.volume_sets) == {"ceramic", "active"}
    assert set(mesh.facet_sets) == {"ground", "electrode"}
    assert_allclose(
        inspect_boundary_topology(mesh).triangle_faces[0].unit_normal,
        [0, 0, -1],
        rtol=0,
        atol=1e-14,
    )
    output = tmp_path / "quadratic.vtu"
    write_vtu(output, mesh, cell_data={"value": {"tetra10": [3], "triangle6": [4]}})
    restored = read_vtu(output)
    assert_array_equal(restored.tetra10, mesh.tetra10)
    assert_array_equal(restored.triangle6, mesh.triangle6)
    assert_array_equal(restored.cell_tags, mesh.cell_tags)
    for name in mesh.volume_sets:
        assert_array_equal(restored.volume_sets[name], mesh.volume_sets[name])
    for name in mesh.facet_sets:
        assert_array_equal(restored.facet_sets[name], mesh.facet_sets[name])
    assert_array_equal(meshio.read(output).cell_data_dict["value"]["tetra10"], [3])


def test_gmsh22_negative_orientation_transports_midpoints(tmp_path):
    mesh = read_gmsh(FIXTURE)
    path = tmp_path / "negative.msh"
    negative = mesh.tetra10[:, [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]]
    meshio.write(
        path,
        meshio.Mesh(points=mesh.points * 1000, cells=[("tetra10", negative)]),
        file_format="gmsh22",
        binary=False,
    )
    result = read_gmsh(path, length_scale_to_m=1e-3)
    assert_array_equal(result.tetra10, mesh.tetra10)
    assert_allclose(result.points, mesh.points, rtol=0, atol=1e-14)


def test_vtu_preserves_optional_ids_and_absent_metadata(tmp_path):
    from fem_toolbox.io import read_vtu

    base = read_gmsh(FIXTURE)
    mesh = Mesh(
        points=base.points,
        tetra10=base.tetra10,
        triangle6=base.triangle6,
        tetra10_ids=[71],
        point_ids=np.arange(100, 110),
        facet_tags=[7],
    )
    path = tmp_path / "ids.vtu"
    write_vtu(path, mesh)
    result = read_vtu(path)
    assert_array_equal(result.point_ids, mesh.point_ids)
    assert_array_equal(result.tetra10_ids, [71])
    assert result.triangle6_ids is None
    assert result.cell_tags is None
    assert_array_equal(result.facet_tags, [7])


def test_vtu_legacy_family_empty_sets_and_reserved_metadata(tmp_path):
    import pytest

    from fem_toolbox.errors import MeshIOError
    from fem_toolbox.io import read_vtu

    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
        tetra4=[[0, 1, 2, 3]],
        triangle3=[[0, 2, 1]],
        volume_sets={"empty": np.array([], dtype=np.int64), "bulk": [0]},
        facet_sets={"wall": [0]},
        tetra4_ids=[42],
    )
    path = tmp_path / "linear.vtu"
    write_vtu(path, mesh)
    restored = read_vtu(path)
    assert set(restored.volume_sets) == {"empty", "bulk"}
    assert set(restored.facet_sets) == {"wall"}
    assert_array_equal(restored.tetra4_ids, [42])
    assert_array_equal(restored.tetra4, mesh.tetra4)
    with pytest.raises(MeshIOError, match="reserved"):
        write_vtu(path, mesh, point_data={"fem_toolbox:external_id": [1, 2, 3, 4]})


def test_import_rejects_curved_quadratic_geometry(tmp_path):
    import pytest

    from fem_toolbox.errors import MeshValidationError
    from fem_toolbox.io import read_vtu

    mesh = read_gmsh(FIXTURE)
    points = mesh.points.copy()
    points[4, 2] += 0.01
    path = tmp_path / "curved.vtu"
    meshio.write(path, meshio.Mesh(points, [("tetra10", mesh.tetra10)]))
    with pytest.raises(MeshValidationError, match="curved"):
        read_vtu(path)
