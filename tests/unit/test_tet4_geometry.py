"""Analytical and invariant tests for linear tetrahedron geometry."""

from itertools import permutations
from math import sqrt

import numpy as np
import pytest

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import DegenerateElementError, ElementGeometryError, ShapeError

RTOL = 2.0e-13
ATOL = 2.0e-14


@pytest.fixture
def reference_tetrahedron() -> np.ndarray:
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def test_reference_tetrahedron_has_analytical_volume_and_gradients(
    reference_tetrahedron: np.ndarray,
) -> None:
    geometry = Tet4Geometry(reference_tetrahedron)
    expected_gradients = np.array(
        [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )

    assert geometry.signed_volume == pytest.approx(1.0 / 6.0)
    assert geometry.volume == pytest.approx(1.0 / 6.0)
    np.testing.assert_allclose(geometry.shape_gradients, expected_gradients)
    np.testing.assert_allclose(geometry.gradient_matrix, expected_gradients.T)
    np.testing.assert_allclose(geometry.shape_gradients.sum(axis=0), 0.0, atol=ATOL)

    np.testing.assert_allclose(
        np.vstack([geometry.shape_values(point) for point in reference_tetrahedron]),
        np.eye(4),
        atol=ATOL,
    )
    np.testing.assert_allclose(
        geometry.shape_values(np.array([0.25, 0.25, 0.25])),
        np.full(4, 0.25),
        atol=ATOL,
    )
    point = np.array([0.2, 0.3, 0.1])
    values = geometry.shape_values(point)
    assert values.sum() == pytest.approx(1.0, abs=ATOL)
    np.testing.assert_allclose(values @ reference_tetrahedron, point, atol=ATOL)
    with pytest.raises(ValueError):
        geometry.coordinates[0, 0] = 2.0
    with pytest.raises(ValueError):
        geometry.faces[0].unit_normal[0] = 0.0


def test_gradient_and_shape_evaluation_never_use_explicit_inverse(
    reference_tetrahedron: np.ndarray, monkeypatch: pytest.MonkeyPatch
) -> None:
    def reject_inverse(*args: object, **kwargs: object) -> None:
        raise AssertionError("explicit inverse is forbidden")

    monkeypatch.setattr(np.linalg, "inv", reject_inverse)

    geometry = Tet4Geometry(reference_tetrahedron)
    geometry.shape_values(np.array([0.2, 0.3, 0.1]))


def test_negative_orientation_is_reported_and_corrected(
    reference_tetrahedron: np.ndarray,
) -> None:
    permuted = reference_tetrahedron[[0, 2, 1, 3]]
    geometry = Tet4Geometry(permuted)

    assert geometry.signed_volume == pytest.approx(-1.0 / 6.0)
    assert geometry.volume == pytest.approx(1.0 / 6.0)
    assert geometry.oriented_node_order == (0, 2, 1, 3)
    assert Tet4Geometry(geometry.oriented_coordinates).signed_volume > 0.0
    np.testing.assert_array_equal(
        geometry.oriented_connectivity(np.array([10, 12, 11, 13])),
        np.array([10, 11, 12, 13]),
    )


def test_faces_have_analytical_areas_and_outward_unit_normals(
    reference_tetrahedron: np.ndarray,
) -> None:
    geometry = Tet4Geometry(reference_tetrahedron)
    expected_areas = {0: sqrt(3.0) / 2.0, 1: 0.5, 2: 0.5, 3: 0.5}
    expected_normals = {
        0: np.array([1.0, 1.0, 1.0]) / sqrt(3.0),
        1: np.array([-1.0, 0.0, 0.0]),
        2: np.array([0.0, -1.0, 0.0]),
        3: np.array([0.0, 0.0, -1.0]),
    }

    assert len(geometry.faces) == 4
    for face in geometry.faces:
        face_points = geometry.coordinates[list(face.local_nodes)]
        face_centroid = face_points.mean(axis=0)
        opposite_vector = geometry.coordinates[face.opposite_node] - face_centroid
        oriented_cross_product = np.cross(
            face_points[1] - face_points[0], face_points[2] - face_points[0]
        )

        assert face.area == pytest.approx(expected_areas[face.opposite_node])
        np.testing.assert_allclose(
            face.unit_normal, expected_normals[face.opposite_node], atol=ATOL
        )
        np.testing.assert_allclose(
            oriented_cross_product / np.linalg.norm(oriented_cross_product),
            face.unit_normal,
            atol=ATOL,
        )
        assert np.linalg.norm(face.unit_normal) == pytest.approx(1.0)
        assert np.dot(face.unit_normal, opposite_vector) < 0.0


def test_regular_tetrahedron_has_unit_normalized_quality() -> None:
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, sqrt(3.0) / 2.0, 0.0],
            [0.5, sqrt(3.0) / 6.0, sqrt(2.0 / 3.0)],
        ]
    )

    quality = Tet4Geometry(coordinates).quality

    assert quality.min_edge_length == pytest.approx(1.0)
    assert quality.max_edge_length == pytest.approx(1.0)
    assert quality.aspect_ratio == pytest.approx(1.0)
    assert quality.scaled_volume == pytest.approx(1.0)
    assert quality.condition_number == pytest.approx(1.0)


def test_rigid_transformations_preserve_geometry_and_rotate_vectors() -> None:
    coordinates = np.array(
        [
            [0.2, -0.1, 0.3],
            [1.4, 0.2, 0.1],
            [-0.2, 1.1, 0.4],
            [0.1, 0.3, 1.6],
        ]
    )
    base = Tet4Geometry(coordinates)
    sample_point = np.array([0.3, 0.2, 0.5])
    random = np.random.default_rng(8128)

    for _ in range(16):
        orthogonal, _ = np.linalg.qr(random.normal(size=(3, 3)))
        if np.linalg.det(orthogonal) < 0.0:
            orthogonal[:, 0] *= -1.0
        translation = random.normal(size=3)
        transformed_coordinates = coordinates @ orthogonal.T + translation
        transformed = Tet4Geometry(transformed_coordinates)

        assert transformed.signed_volume == pytest.approx(base.signed_volume, rel=RTOL)
        assert transformed.quality.scaled_volume == pytest.approx(
            base.quality.scaled_volume, rel=RTOL
        )
        assert transformed.quality.condition_number == pytest.approx(
            base.quality.condition_number, rel=RTOL
        )
        np.testing.assert_allclose(
            transformed.shape_gradients,
            base.shape_gradients @ orthogonal.T,
            rtol=RTOL,
            atol=ATOL,
        )
        np.testing.assert_allclose(
            transformed.shape_values(sample_point @ orthogonal.T + translation),
            base.shape_values(sample_point),
            rtol=RTOL,
            atol=ATOL,
        )
        for base_face, transformed_face in zip(
            base.faces, transformed.faces, strict=True
        ):
            assert transformed_face.area == pytest.approx(base_face.area, rel=RTOL)
            np.testing.assert_allclose(
                transformed_face.unit_normal,
                base_face.unit_normal @ orthogonal.T,
                rtol=RTOL,
                atol=ATOL,
            )


def test_all_node_permutations_preserve_unoriented_geometry() -> None:
    coordinates = np.array(
        [
            [0.2, -0.1, 0.3],
            [1.4, 0.2, 0.1],
            [-0.2, 1.1, 0.4],
            [0.1, 0.3, 1.6],
        ]
    )
    base = Tet4Geometry(coordinates)
    point = coordinates.mean(axis=0)
    base_faces = {face.opposite_node: face for face in base.faces}

    for permutation in permutations(range(4)):
        permuted = Tet4Geometry(coordinates[list(permutation)])
        inversion_count = sum(
            permutation[i] > permutation[j] for i in range(4) for j in range(i + 1, 4)
        )
        sign = -1.0 if inversion_count % 2 else 1.0

        assert permuted.signed_volume == pytest.approx(
            sign * base.signed_volume, rel=RTOL
        )
        assert permuted.volume == pytest.approx(base.volume, rel=RTOL)
        assert permuted.quality.scaled_volume == pytest.approx(
            base.quality.scaled_volume, rel=RTOL
        )
        assert permuted.quality.condition_number == pytest.approx(
            base.quality.condition_number, rel=RTOL
        )
        np.testing.assert_allclose(
            permuted.shape_gradients,
            base.shape_gradients[list(permutation)],
            rtol=RTOL,
            atol=ATOL,
        )
        np.testing.assert_allclose(
            permuted.shape_values(point),
            base.shape_values(point)[list(permutation)],
            rtol=RTOL,
            atol=ATOL,
        )
        for face in permuted.faces:
            original_opposite_node = permutation[face.opposite_node]
            base_face = base_faces[original_opposite_node]
            assert face.area == pytest.approx(base_face.area, rel=RTOL)
            np.testing.assert_allclose(
                face.unit_normal, base_face.unit_normal, rtol=RTOL, atol=ATOL
            )
        assert Tet4Geometry(permuted.oriented_coordinates).signed_volume > 0.0


@pytest.mark.parametrize("scale", [1.0e-9, 1.0, 1.0e9])
def test_relative_degeneracy_check_is_scale_invariant(scale: float) -> None:
    nearly_coplanar = scale * np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.2, 0.2, 1.0e-14],
        ]
    )

    with pytest.raises(DegenerateElementError) as caught:
        Tet4Geometry(nearly_coplanar, entity_id=77)

    assert caught.value.entity_id == 77
    assert caught.value.scaled_volume <= caught.value.tolerance


def test_exactly_coplanar_tetrahedron_is_rejected() -> None:
    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ]
    )

    with pytest.raises(DegenerateElementError):
        Tet4Geometry(coordinates)


@pytest.mark.parametrize(
    ("coordinates", "error"),
    [
        (np.zeros((3, 3)), ShapeError),
        (np.full((4, 3), np.nan), ElementGeometryError),
    ],
)
def test_invalid_coordinates_are_rejected(
    coordinates: np.ndarray, error: type[Exception]
) -> None:
    with pytest.raises(error):
        Tet4Geometry(coordinates)


def test_invalid_degeneracy_tolerance_is_rejected(
    reference_tetrahedron: np.ndarray,
) -> None:
    with pytest.raises(ElementGeometryError):
        Tet4Geometry(reference_tetrahedron, degeneracy_tolerance=0.0)
