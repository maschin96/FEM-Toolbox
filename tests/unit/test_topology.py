"""Independent invariants for filtered Tet4 compliance design."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ElasticMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.topology import LoadCase, TopologyProblem


def problem(mesh, **kwargs):
    dofs = mechanical_dofs(mesh)
    left = np.flatnonzero(mesh.points[:, 0] == 0)
    right = np.flatnonzero(mesh.points[:, 0] == 1)
    bc = (Dirichlet(dofs.indices("u", left), np.zeros(len(left) * 3)),)
    cases = [
        LoadCase("y", bc, point_loads={int(n): (0, -1 / len(right), 0) for n in right}),
        LoadCase("z", bc, point_loads={int(n): (0, 0, -1 / len(right)) for n in right}),
    ]
    return TopologyProblem(
        mesh,
        ElasticMaterial.isotropic(1000, 0.3),
        cases,
        filter_radius=0.6,
        weights=(0.4, 0.6),
        **kwargs,
    )


def test_full_material_matches_existing_mechanics(cube_mesh):
    mesh = cube_mesh(1)
    design = problem(mesh)
    state = design.evaluate(np.ones(mesh.n_cells))
    direct = solve_elasticity(
        mesh,
        design.materials,
        design.cases[0].conditions,
        point_loads=design.cases[0].point_loads,
    )
    assert_allclose(
        state.states[0].displacement, direct.displacement, rtol=2e-10, atol=1e-12
    )
    expected = sum(
        w * 2 * s.energy for w, s in zip(design.weights, state.states, strict=True)
    )
    assert_allclose(state.objective, expected, rtol=2e-10)
    assert state.volume_fraction == pytest.approx(1)


def test_constant_filter_and_fixed_regions(cube_mesh):
    design = problem(cube_mesh(1), fixed_solid=(0,), fixed_void=(5,))
    x = np.full(6, 0.5)
    x[0] = 1
    x[5] = 0
    state = design.evaluate(x)
    assert state.physical_density[0] == 1
    assert state.physical_density[5] == 0
    assert_allclose(state.objective_gradient[[0, 5]], 0, atol=0)
    assert_allclose(state.volume_gradient[[0, 5]], 0, atol=0)
    free = problem(cube_mesh(1)).evaluate(np.full(6, 0.4))
    assert_allclose(free.physical_density, 0.4, atol=1e-15)


@pytest.mark.parametrize("value", [float("nan"), -0.1, 1.1])
def test_invalid_density_fails(cube_mesh, value):
    with pytest.raises(ValidationError):
        problem(cube_mesh(1)).evaluate(np.full(6, value))


def test_nonzero_support_displacement_is_not_a_compliance_problem(cube_mesh):
    mesh = cube_mesh(1)
    with pytest.raises(ValidationError, match="homogeneous"):
        TopologyProblem(
            mesh,
            ElasticMaterial.isotropic(1000, 0.3),
            [LoadCase("bad", (Dirichlet([0], [0.01]),))],
            filter_radius=0.5,
        )
