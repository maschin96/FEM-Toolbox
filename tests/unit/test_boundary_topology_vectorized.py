"""Parity with the frozen scalar oracle, including validation diagnostics."""

from dataclasses import fields
from unittest.mock import patch

import numpy as np
import pytest

from fem_toolbox.errors import BoundaryTopologyError
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from tests.reference.legacy_boundary_topology import (
    inspect_boundary_topology as legacy_inspect,
)
from tests.reference.quadratic_mesh import elevate


def assert_parity(mesh, **kwargs):
    expected = legacy_inspect(mesh, **kwargs)
    actual = inspect_boundary_topology(mesh, **kwargs)
    assert list(actual.facet_incidence.items()) == list(
        expected.facet_incidence.items()
    )
    for name in ("exterior_faces", "triangle_faces"):
        for old, new in zip(
            getattr(expected, name), getattr(actual, name), strict=True
        ):
            for field in fields(old):
                if field.name == "unit_normal":
                    np.testing.assert_allclose(
                        new.unit_normal, old.unit_normal, rtol=2e-14, atol=2e-14
                    )
                    with pytest.raises(ValueError):
                        new.unit_normal.setflags(write=True)
                elif field.name == "area":
                    assert new.area == pytest.approx(old.area, rel=2e-14, abs=0)
                else:
                    assert getattr(new, field.name) == getattr(old, field.name)


@pytest.mark.parametrize("seed", range(4))
@pytest.mark.parametrize("quadratic", [False, True])
def test_random_numbering_orientation_ids_and_overlapping_groups(
    cube_mesh, seed, quadratic
):
    mesh = cube_mesh(2, distort=True)
    if quadratic:
        mesh = elevate(mesh)
    rng = np.random.default_rng(seed)
    permutation = rng.permutation(len(mesh.points))
    points = np.empty_like(mesh.points)
    points[permutation] = mesh.points
    cells = permutation[mesh.cells[rng.permutation(mesh.n_cells)]]
    if not quadratic:
        cells = np.array([row[rng.permutation(4)] for row in cells])
    facets = permutation[mesh.facets[rng.permutation(mesh.n_facets)]]
    family = (
        dict(
            tetra10=cells,
            triangle6=facets,
            tetra10_ids=rng.permutation(mesh.n_cells) + 100,
            triangle6_ids=rng.permutation(mesh.n_facets) + 200,
        )
        if quadratic
        else dict(
            tetra4=cells,
            triangle3=facets,
            tetra4_ids=rng.permutation(mesh.n_cells) + 100,
            triangle3_ids=rng.permutation(mesh.n_facets) + 200,
        )
    )
    mesh = Mesh(
        points=points,
        **family,
        facet_sets={"all": np.arange(mesh.n_facets), "overlap": [0, 2]},
    )
    assert_parity(mesh, require_complete=True)


def test_empty_and_partial_boundary(cube_mesh):
    assert_parity(Mesh.empty())
    source = cube_mesh(2)
    assert_parity(Mesh(points=source.points, tetra4=source.cells))


@pytest.mark.parametrize(
    "height",
    [
        0.0,
        1e-14,
        1e-12,
        np.nextafter(2e-12, 0),
        2e-12,
        np.nextafter(2e-12, np.inf),
        1e-10,
    ],
)
def test_near_degenerate_cells_keep_original_diagnostics(height):
    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, height]],
        tetra4=[[0, 1, 2, 3]],
        tetra4_ids=[728],
    )
    try:
        legacy_inspect(mesh)
    except Exception as expected:
        with pytest.raises(type(expected)) as actual:
            inspect_boundary_topology(mesh)
        assert str(actual.value) == str(expected)
    else:
        assert_parity(mesh)


def test_nonmanifold_without_surface_keeps_diagnostic():
    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [0, 0, 2]],
        tetra4=[[0, 1, 2, 3], [0, 2, 1, 4], [0, 1, 2, 5]],
        tetra4_ids=[91, 73, 17],
    )
    with pytest.raises(BoundaryTopologyError) as expected:
        legacy_inspect(mesh)
    with pytest.raises(type(expected.value)) as actual:
        inspect_boundary_topology(mesh)
    assert str(actual.value) == str(expected.value)


def test_linear_fast_path_does_not_construct_per_cell_geometry(cube_mesh):
    mesh = cube_mesh(2)
    with patch("fem_toolbox.mesh.topology.Tet4Geometry", side_effect=AssertionError):
        topology = inspect_boundary_topology(mesh)
    assert len(topology.exterior_faces) == mesh.n_facets


@pytest.mark.parametrize("kind", ["foreign", "interior", "nonmanifold", "missing"])
def test_boundary_failures_keep_all_diagnostics(kind):
    points = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [0, 0, 2]]
    cells = [[0, 1, 2, 3]]
    if kind in ("interior", "nonmanifold"):
        cells.append([0, 2, 1, 4])
    if kind == "nonmanifold":
        cells.append([0, 1, 2, 5])
    triangle = [0, 1, 5] if kind == "foreign" else [0, 1, 2]
    mesh = Mesh(
        points=points,
        tetra4=cells,
        triangle3=[triangle],
        tetra4_ids=np.arange(len(cells)) + 981,
        triangle3_ids=[72],
        facet_sets={"load": [0], "overlap": [0]},
    )
    with pytest.raises(BoundaryTopologyError) as expected:
        legacy_inspect(mesh, require_complete=True)
    with pytest.raises(type(expected.value)) as actual:
        inspect_boundary_topology(mesh, require_complete=True)
    assert vars(actual.value) == vars(expected.value)
    assert str(actual.value) == str(expected.value)
