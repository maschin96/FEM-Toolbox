"""Real PETSc solves: independent algebra, scalar fidelity and failure paths."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError
from fem_toolbox.solvers.linear_backend import solve_linear


def test_scaled_mixed_fields_and_reuse(petsc_backend):
    a = np.array([[4.0, 1.0], [1.0, -2.0]])
    factor = petsc_backend.factorize(
        a,
        fields={"u": [0], "phi": [1]},
        scale=[0.5, 2],
        pivot_tolerance=1e-12,
        relative_pivot=True,
    )
    a[:] = 99
    result = factor.solve([5, 8])
    np.testing.assert_allclose(result.values, [2, -3], rtol=0, atol=1e-13)
    np.testing.assert_allclose(
        factor.solve([10, 16]).values, [4, -6], rtol=0, atol=1e-13
    )
    assert result.report.backend.startswith("petsc")
    assert result.report.converged
    assert result.report.relative_residual < 1e-13
    assert set(result.report.field_relative_residuals) == {"u", "phi"}
    factor.close()
    factor.close()
    with pytest.raises(ValidationError, match="closed"):
        factor.solve([1, 2])


@pytest.mark.parametrize("complex_operator", [False, True])
def test_complex_multiple_rhs_without_scalar_loss(petsc_backend, complex_operator):
    a = np.array(
        [[2.0, 0.25], [0.25, 3.0]], dtype=complex if complex_operator else float
    )
    if complex_operator:
        a[0, 1] += 1j
        a[1, 0] -= 0.5j
    expected = np.array([[1 + 2j, -3j], [2 - 1j, 4 + 1j]])
    result = solve_linear(
        a, a @ expected, fields={"u": [0], "phi": [1]}, backend=petsc_backend
    )
    np.testing.assert_allclose(result.values, expected, rtol=0, atol=1e-13)
    assert np.iscomplexobj(result.values)
    assert result.report.relative_residual < 1e-13


def test_real_result_stays_real_on_complex_petsc_build(petsc_backend):
    result = solve_linear(
        [[2.0, 0], [0, 4.0]],
        [2.0, 8.0],
        fields={"T": [0, 1]},
        kind="positive_definite",
        backend=petsc_backend,
    )
    assert result.values.dtype == np.float64
    np.testing.assert_allclose(result.values, [1, 2], rtol=0, atol=1e-13)


def test_empty_system(petsc_backend):
    result = solve_linear(
        np.empty((0, 0)), [], fields={"u": [], "phi": []}, backend=petsc_backend
    )
    assert result.values.shape == (0,)
    assert result.report.relative_residual == 0


@pytest.mark.parametrize("a", [[[1.0, 1.0], [1.0, 1.0]], [[1.0, 0.0], [0.0, 1e-15]]])
def test_singularity_and_small_pivots_rejected(petsc_backend, a):
    with pytest.raises(ValidationError, match=r"PETSc|pivot|singular"):
        solve_linear(
            a,
            [0, 0],
            fields={"u": [0, 1]},
            pivot_tolerance=1e-12,
            backend=petsc_backend,
        )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"fields": {"u": [0]}},
        {"fields": {"u": [0, 1]}, "scale": [1, 0]},
        {"fields": {"u": [0, 1]}, "kind": "unknown"},
    ],
)
def test_invalid_contract(petsc_backend, kwargs):
    with pytest.raises(ValidationError):
        solve_linear(np.eye(2), [1, 2], backend=petsc_backend, **kwargs)


def test_invalid_spd_and_rhs(petsc_backend):
    with pytest.raises(ValidationError, match="positive"):
        solve_linear(
            np.diag([1.0, -1]),
            [1, 2],
            fields={"u": [0, 1]},
            kind="positive_definite",
            backend=petsc_backend,
        )
    factor = petsc_backend.factorize(np.eye(2), fields={"u": [0, 1]})
    with pytest.raises(ValidationError, match="RHS"):
        factor.solve([1, np.nan])
    factor.close()


def test_pure_imaginary_operator_and_no_scipy_fallback(petsc_backend, monkeypatch):
    from fem_toolbox.solvers import linear_backend

    def forbidden(*args, **kwargs):
        pytest.fail("PETSc must not delegate factorization to SciPy")

    monkeypatch.setattr(linear_backend, "splu", forbidden)
    a = np.array([[0, 2j], [3j, 0]])
    result = solve_linear(
        a, [-2 + 4j, -6 + 3j], fields={"u": [0, 1]}, backend=petsc_backend
    )
    np.testing.assert_allclose(result.values, [1 + 2j, 2 + 1j], rtol=0, atol=1e-13)


def test_no_rhs_columns(petsc_backend):
    factor = petsc_backend.factorize(np.eye(2), fields={"T": [0, 1]})
    assert factor.solve(np.empty((2, 0))).values.shape == (2, 0)
    factor.close()


def test_spd_hint_does_not_apply_general_lu_pivot_threshold(petsc_backend):
    result = solve_linear(
        [[2.0]],
        [6.0],
        fields={"phi": [0]},
        kind="positive_definite",
        pivot_tolerance=10.0,
        backend=petsc_backend,
    )
    np.testing.assert_allclose(result.values, [3.0], rtol=0, atol=1e-13)
