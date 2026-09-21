"""Contract tests for Gmsh import and VTU export adapters."""

from pathlib import Path

import meshio
import numpy as np
import pytest

from fem_toolbox.errors import MeshIOError, UnsupportedCellTypeError
from fem_toolbox.io import read_gmsh, write_vtu

FIXTURES = Path(__file__).parents[1] / "fixtures"


def test_read_gmsh_preserves_tet4_triangle3_and_overlapping_groups() -> None:
    mesh = read_gmsh(FIXTURES / "overlapping_groups.msh")

    np.testing.assert_allclose(
        mesh.points,
        np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 1.0, 1.0],
            ]
        ),
    )
    np.testing.assert_array_equal(mesh.tetra4, np.array([[0, 1, 2, 3], [1, 2, 3, 4]]))
    np.testing.assert_array_equal(mesh.triangle3, np.array([[0, 1, 2]]))
    assert mesh.cell_tags is not None
    assert mesh.facet_tags is not None
    assert set(mesh.volume_sets) == {"active", "ceramic"}
    assert set(mesh.facet_sets) == {"electrode", "ground"}
    np.testing.assert_array_equal(mesh.cell_tags, np.array([21, 21]))
    np.testing.assert_array_equal(mesh.facet_tags, np.array([11]))
    np.testing.assert_array_equal(mesh.volume_sets["ceramic"], np.array([0, 1]))
    np.testing.assert_array_equal(mesh.volume_sets["active"], np.array([0, 1]))
    np.testing.assert_array_equal(mesh.facet_sets["ground"], np.array([0]))
    np.testing.assert_array_equal(mesh.facet_sets["electrode"], np.array([0]))


def test_read_gmsh_rejects_unsupported_cell_type(tmp_path: Path) -> None:
    path = tmp_path / "unsupported.msh"
    meshio.write(
        path,
        meshio.Mesh(
            points=np.array(
                [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [1.0, 1.0, 0.0],
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                    [1.0, 0.0, 1.0],
                    [1.0, 1.0, 1.0],
                    [0.0, 1.0, 1.0],
                ]
            ),
            cells=[("hexahedron", np.arange(8).reshape(1, 8))],
        ),
        file_format="gmsh22",
        binary=False,
    )

    with pytest.raises(UnsupportedCellTypeError, match="hexahedron") as caught:
        read_gmsh(path)

    assert caught.value.cell_types == ("hexahedron",)


def test_read_gmsh_scales_to_metres_and_normalizes_tet4_orientation(
    tmp_path: Path,
) -> None:
    path = tmp_path / "negative_orientation.msh"
    points_mm = np.array(
        [
            [0.0, 0.0, 0.0],
            [1000.0, 0.0, 0.0],
            [0.0, 1000.0, 0.0],
            [0.0, 0.0, 1000.0],
        ]
    )
    meshio.write(
        path,
        meshio.Mesh(
            points=points_mm,
            cells=[("tetra", np.array([[0, 2, 1, 3]]))],
        ),
        file_format="gmsh22",
        binary=False,
    )

    mesh = read_gmsh(path, length_scale_to_m=1.0e-3)

    np.testing.assert_allclose(mesh.points, points_mm * 1.0e-3)
    np.testing.assert_array_equal(mesh.tetra4, np.array([[0, 1, 2, 3]]))


def test_vtu_roundtrip_preserves_topology_sets_and_result_data(
    tmp_path: Path,
) -> None:
    mesh = read_gmsh(FIXTURES / "overlapping_groups.msh")
    path = tmp_path / "result.vtu"
    point_data = {"temperature_K": np.linspace(293.15, 297.15, mesh.n_points)}
    cell_data = {
        "energy_density_J_m3": {
            "tetra4": np.array([1.25, 2.5]),
            "triangle3": np.array([0.0]),
        }
    }

    write_vtu(path, mesh, point_data=point_data, cell_data=cell_data)
    restored = meshio.read(path)

    np.testing.assert_allclose(restored.points, mesh.points)
    np.testing.assert_array_equal(restored.cells_dict["tetra"], mesh.tetra4)
    np.testing.assert_array_equal(restored.cells_dict["triangle"], mesh.triangle3)
    np.testing.assert_allclose(
        restored.point_data["temperature_K"], point_data["temperature_K"]
    )
    np.testing.assert_allclose(
        restored.cell_data_dict["energy_density_J_m3"]["tetra"],
        cell_data["energy_density_J_m3"]["tetra4"],
    )
    np.testing.assert_array_equal(
        restored.cell_data_dict["fem_toolbox:set:ceramic"]["tetra"],
        np.array([1, 1], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        restored.cell_data_dict["fem_toolbox:set:active"]["tetra"],
        np.array([1, 1], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        restored.cell_data_dict["fem_toolbox:set:ground"]["triangle"],
        np.array([1], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        restored.cell_data_dict["fem_toolbox:set:electrode"]["triangle"],
        np.array([1], dtype=np.uint8),
    )
    np.testing.assert_array_equal(
        restored.cell_data_dict["fem_toolbox:physical_tag"]["tetra"],
        mesh.cell_tags,
    )
    np.testing.assert_array_equal(
        restored.cell_data_dict["fem_toolbox:physical_tag"]["triangle"],
        mesh.facet_tags,
    )


def test_write_vtu_rejects_incomplete_cell_result_field(tmp_path: Path) -> None:
    mesh = read_gmsh(FIXTURES / "overlapping_groups.msh")

    with pytest.raises(MeshIOError, match="triangle3"):
        write_vtu(
            tmp_path / "result.vtu",
            mesh,
            cell_data={"stress_Pa": {"tetra4": np.ones(mesh.n_tetra4)}},
        )
