"""A feasible improved design and a stopped computation are different states."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ValidationError
from fem_toolbox.topology_optimization import OptimizationOptions, optimize_topology
from tests.unit.test_topology import problem


def test_multiload_design_improves_and_repeats(cube_mesh):
    design = problem(cube_mesh(2), fixed_solid=(0,))
    options = OptimizationOptions(max_iterations=100)
    a = optimize_topology(design, 0.6, options=options)
    b = optimize_topology(design, 0.6, options=options)
    assert a.status == "converged", a.message
    assert a.final.objective < 0.99 * a.initial_objective
    assert a.final.volume_fraction <= 0.6 + 1e-8
    assert a.final.density[0] == 1
    assert_allclose(a.final.density, b.final.density, atol=1e-10, rtol=1e-10)
    assert a.history[0]["iteration"] == 0
    assert all("maximum_relative_residual" in row for row in a.history)


def test_evaluation_budget_is_not_convergence(cube_mesh):
    result = optimize_topology(
        problem(cube_mesh(1)), 0.5, options=OptimizationOptions(max_evaluations=1)
    )
    assert result.status == "evaluation_limit"
    assert result.evaluations == 1
    assert result.final.volume_fraction <= 0.5 + 1e-8


def test_infeasible_fixed_region_and_initial_design_fail(cube_mesh):
    with pytest.raises(ValidationError, match="infeasible"):
        optimize_topology(problem(cube_mesh(1), fixed_solid=(0, 1, 2, 3, 4, 5)), 0.5)
    with pytest.raises(ValidationError, match="initial"):
        optimize_topology(problem(cube_mesh(1)), 0.5, initial=np.ones(6))


def test_invalid_options_fail():
    for kwargs in (
        {"max_iterations": 0},
        {"max_evaluations": True},
        {"tolerance": float("nan")},
    ):
        with pytest.raises(ValidationError):
            OptimizationOptions(**kwargs)
