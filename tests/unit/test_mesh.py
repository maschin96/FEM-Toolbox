"""Contract tests for the validated in-memory mesh model."""

from collections.abc import Mapping

import numpy as np
import pytest

from fem_toolbox.errors import (
    ConnectivityError,
    MeshValidationError,
    SelectionError,
    ShapeError,
)
from fem_toolbox.mesh import Mesh


def test_empty_mesh_has_canonical_shapes() -> None:
    mesh = Mesh.empty()

    assert mesh.points.shape == (0, 3)
    assert mesh.tetra4.shape == (0, 4)
    assert mesh.triangle3.shape == (0, 3)
    assert mesh.n_points == 0
    assert mesh.n_tetra4 == 0
    assert mesh.n_triangle3 == 0
    assert mesh.point_ids is None
    assert mesh.tetra4_ids is None
    assert mesh.triangle3_ids is None
    assert mesh.cell_tags is None
    assert mesh.facet_tags is None
    assert mesh.volume_sets == {}
    assert mesh.facet_sets == {}


def test_mesh_preserves_external_ids_and_overlapping_selection_sets() -> None:
    points = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
    )
    tetra4 = np.array([[0, 1, 2, 3], [1, 2, 3, 4]], dtype=np.int64)
    triangle3 = np.array([[0, 2, 1], [1, 3, 2]], dtype=np.int64)
    volume_sets = {"ceramic": np.array([0, 1]), "electroded": np.array([1])}

    mesh = Mesh(
        points=points,
        tetra4=tetra4,
        triangle3=triangle3,
        point_ids=np.array([10, 20, 30, 40, 50]),
        tetra4_ids=np.array([101, 102]),
        triangle3_ids=np.array([201, 202]),
        cell_tags=np.array([7, 8]),
        facet_tags=np.array([12, 12]),
        volume_sets=volume_sets,
        facet_sets={"ground": np.array([0]), "electrode": np.array([0, 1])},
    )
    points[0, 0] = 99.0
    tetra4[0, 0] = 4
    volume_sets["ceramic"][0] = 1

    assert mesh.points[0, 0] == 0.0
    assert mesh.tetra4[0].tolist() == [0, 1, 2, 3]
    assert mesh.point_ids is not None
    assert mesh.point_ids.tolist() == [10, 20, 30, 40, 50]
    assert mesh.tetra4_ids is not None
    assert mesh.tetra4_ids.tolist() == [101, 102]
    assert mesh.triangle3_ids is not None
    assert mesh.triangle3_ids.tolist() == [201, 202]
    assert mesh.cell_tags is not None
    assert mesh.cell_tags.tolist() == [7, 8]
    assert mesh.facet_tags is not None
    assert mesh.facet_tags.tolist() == [12, 12]
    assert mesh.volume_sets["ceramic"].tolist() == [0, 1]
    assert mesh.volume_sets["electroded"].tolist() == [1]
    assert mesh.facet_sets["ground"].tolist() == [0]
    assert mesh.facet_sets["electrode"].tolist() == [0, 1]

    with pytest.raises(ValueError):
        mesh.points[0, 0] = 3.0
    with pytest.raises(ValueError):
        mesh.volume_sets["ceramic"][0] = 1
    with pytest.raises(TypeError):
        mesh.volume_sets["new"] = np.array([0])  # type: ignore[index]


def test_connectivity_is_zero_based_independently_of_external_point_ids() -> None:
    mesh = Mesh(
        points=np.eye(4, 3),
        tetra4=np.array([[0, 1, 2, 3]]),
        point_ids=np.array([11, 12, 13, 14]),
    )

    assert mesh.tetra4[0].tolist() == [0, 1, 2, 3]


@pytest.mark.parametrize(
    ("field", "value", "expected"),
    [
        ("points", np.zeros((4, 2)), "(n, 3)"),
        ("tetra4", np.zeros((1, 3), dtype=np.int64), "(n, 4)"),
        ("triangle3", np.zeros((1, 4), dtype=np.int64), "(n, 3)"),
    ],
)
def test_invalid_array_shape_raises_shape_error(
    field: str, value: np.ndarray, expected: str
) -> None:
    kwargs: dict[str, object] = {"points": np.zeros((4, 3))}
    kwargs[field] = value

    with pytest.raises(ShapeError) as caught:
        Mesh(**kwargs)  # type: ignore[arg-type]

    assert caught.value.field == field
    assert expected in str(caught.value)


@pytest.mark.parametrize(
    ("connectivity", "reason"),
    [
        (np.array([[0, 1, 2, 4]]), "outside"),
        (np.array([[0, 1, 1, 3]]), "repeated"),
    ],
)
def test_invalid_connectivity_reports_external_entity_id(
    connectivity: np.ndarray, reason: str
) -> None:
    with pytest.raises(ConnectivityError) as caught:
        Mesh(
            points=np.eye(4, 3),
            tetra4=connectivity,
            tetra4_ids=np.array([4711]),
        )

    assert caught.value.field == "tetra4"
    assert caught.value.entity_kind == "tetra4"
    assert caught.value.entity_id == 4711
    assert reason in str(caught.value)


def test_duplicate_connectivity_reports_both_external_entity_ids() -> None:
    with pytest.raises(ConnectivityError) as caught:
        Mesh(
            points=np.eye(4, 3),
            tetra4=np.array([[0, 1, 2, 3], [3, 2, 1, 0]]),
            tetra4_ids=np.array([41, 42]),
        )

    assert caught.value.entity_id == 42
    assert caught.value.duplicate_of == 41


@pytest.mark.parametrize(
    ("kwargs", "field"),
    [
        ({"points": np.array([[0.0, 0.0, np.nan]])}, "points"),
        (
            {
                "points": np.eye(4, 3),
                "tetra4": np.array([[0.0, 1.0, 2.0, 3.0]]),
            },
            "tetra4",
        ),
        (
            {
                "points": np.eye(4, 3),
                "point_ids": np.array([1, 2, 2, 4]),
            },
            "point_ids",
        ),
        (
            {
                "points": np.eye(4, 3),
                "tetra4": np.array([[0, 1, 2, 3]]),
                "cell_tags": np.array([1.0]),
            },
            "cell_tags",
        ),
    ],
)
def test_invalid_mesh_base_data_raises_specific_error(
    kwargs: dict[str, object], field: str
) -> None:
    with pytest.raises(MeshValidationError) as caught:
        Mesh(**kwargs)  # type: ignore[arg-type]

    assert caught.value.field == field


@pytest.mark.parametrize(
    ("volume_sets", "reason"),
    [
        ({"duplicated": np.array([0, 0])}, "duplicate"),
        ({"outside": np.array([1])}, "outside"),
        ({"": np.array([0])}, "non-empty"),
    ],
)
def test_invalid_volume_selection_is_rejected(
    volume_sets: Mapping[str, np.ndarray], reason: str
) -> None:
    with pytest.raises(SelectionError) as caught:
        Mesh(
            points=np.eye(4, 3),
            tetra4=np.array([[0, 1, 2, 3]]),
            volume_sets=volume_sets,
        )

    assert caught.value.field == "volume_sets"
    assert reason in str(caught.value)
