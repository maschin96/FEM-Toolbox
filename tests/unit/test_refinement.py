"""Conformity, inheritance and deterministic local refinement."""

import numpy as np
import pytest

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError
from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from fem_toolbox.mesh.refinement import refine_tet4


def test_refinement_preserves_regions_and_surfaces(cube_mesh):
    base = cube_mesh(2)
    tags = (base.points[base.tetra4].mean(axis=1)[:, 0] >= 0.5).astype(int)
    mesh = Mesh(
        points=base.points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        cell_tags=tags,
        facet_tags=np.arange(base.n_triangle3),
        volume_sets={"left": np.flatnonzero(tags == 0)},
        facet_sets={"wall": np.arange(base.n_triangle3)},
    )
    refined = refine_tet4(mesh, [0, 1])
    new = refined.mesh
    inspect_boundary_topology(new, require_complete=True)
    np.testing.assert_array_equal(new.cell_tags, tags[refined.cell_parents])
    np.testing.assert_array_equal(
        new.facet_tags, mesh.facet_tags[refined.facet_parents]
    )
    np.testing.assert_array_equal(
        new.volume_sets["left"], np.flatnonzero(new.cell_tags == 0)
    )
    assert len(new.facet_sets["wall"]) == new.n_triangle3
    for tag in (0, 1):
        before = sum(
            Tet4Geometry(mesh.points[c]).volume for c in mesh.tetra4[tags == tag]
        )
        after = sum(
            Tet4Geometry(new.points[c]).volume for c in new.tetra4[new.cell_tags == tag]
        )
        assert after == pytest.approx(before)
    repeat = refine_tet4(mesh, [1, 0])
    np.testing.assert_array_equal(new.tetra4, repeat.mesh.tetra4)
    assert new.n_points < refine_tet4(mesh, np.arange(mesh.n_tetra4)).mesh.n_points


def test_refinement_budget_and_invalid_marks(cube_mesh):
    mesh = cube_mesh(1)
    with pytest.raises(ValidationError, match="budget"):
        refine_tet4(mesh, [0], max_points=mesh.n_points)
    with pytest.raises(ValidationError):
        refine_tet4(mesh, [mesh.n_tetra4])
    assert refine_tet4(mesh, []).mesh is mesh
