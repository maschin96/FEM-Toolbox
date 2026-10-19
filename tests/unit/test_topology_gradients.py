"""Finite differences re-solve states and never reuse the adjoint expression."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from tests.unit.test_topology import problem


@pytest.mark.parametrize("center", [0.02, 0.45, 0.98])
def test_multiload_filter_chain_and_step_study(cube_mesh, center):
    design = problem(cube_mesh(1), fixed_solid=(0,), fixed_void=(5,))
    x = np.full(6, center)
    x[0] = 1
    x[5] = 0
    exact = design.evaluate(x)
    errors = []
    for h in (1e-3, 1e-4, 1e-5):
        dc = np.zeros(6)
        dv = np.zeros(6)
        for j in design.free:
            plus = x.copy()
            minus = x.copy()
            plus[j] += h
            minus[j] -= h
            a, b = design.evaluate(plus), design.evaluate(minus)
            dc[j] = (a.objective - b.objective) / (2 * h)
            dv[j] = (a.volume_fraction - b.volume_fraction) / (2 * h)
        errors.append(
            np.linalg.norm(dc - exact.objective_gradient)
            / np.linalg.norm(exact.objective_gradient)
        )
        assert_allclose(dv, exact.volume_gradient, atol=2e-9, rtol=0)
    assert errors[-1] < 2e-4
    assert min(errors[1:]) < errors[0]


def test_distorted_unequal_volumes_directional_gradient(cube_mesh):
    design = problem(cube_mesh(2, distort=True), fixed_solid=(0,), fixed_void=(10,))
    x = np.linspace(0.25, 0.75, design.mesh.n_cells)
    x[0] = 1
    x[10] = 0
    direction = np.random.default_rng(78).normal(size=len(x))
    direction[[0, 10]] = 0
    direction /= np.linalg.norm(direction)
    exact = design.evaluate(x)
    h = 1e-5
    a, b = design.evaluate(x + h * direction), design.evaluate(x - h * direction)
    assert_allclose(
        (a.objective - b.objective) / (2 * h),
        exact.objective_gradient @ direction,
        rtol=2e-4,
        atol=1e-10,
    )
    assert_allclose(
        (a.volume_fraction - b.volume_fraction) / (2 * h),
        exact.volume_gradient @ direction,
        atol=2e-9,
        rtol=0,
    )
