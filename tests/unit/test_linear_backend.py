"""Independent algebraic backend contract, including mixed units and phasors."""

import numpy as np
import pytest
from scipy import sparse

from fem_toolbox.errors import ValidationError
from fem_toolbox.solvers.linear_backend import ScipyBackend, solve_linear


def test_scaled_indefinite_system_and_field_residuals():
    # x=(2,-3) solves 4*x0+x1=5, x0-2*x1=8.
    a = sparse.csc_matrix([[4.0, 1.0], [1.0, -2.0]])
    result = solve_linear(
        a, [5.0, 8.0], fields={"u": [0], "phi": [1]}, scale=[0.5, 2.0]
    )
    np.testing.assert_allclose(result.values, [2, -3], rtol=0, atol=1e-14)
    assert result.report.converged
    assert result.report.reason == "direct_solve"
    assert set(result.report.field_relative_residuals) == {"u", "phi"}
    assert max(result.report.field_relative_residuals.values()) < 1e-14
    np.testing.assert_allclose(
        result.report.residual, a @ result.values - [5, 8], atol=0
    )
    with pytest.raises(ValueError):
        result.values[0] = 100


@pytest.mark.parametrize("kind", ["general", "positive_definite"])
def test_multiple_rhs_reuse_and_complex_rhs_for_real_factor(kind):
    a = np.array([[2.0, 0.0], [0.0, 4.0]])
    factor = ScipyBackend().factorize(a, fields={"phi": [0, 1]}, kind=kind)
    a[:] = 99  # factor owns operator snapshot
    first = factor.solve([2.0, 8.0])
    second = factor.solve([[2 + 4j, 6], [8 - 4j, 16]])
    np.testing.assert_allclose(first.values, [1, 2], rtol=0, atol=1e-14)
    np.testing.assert_allclose(
        second.values, [[1 + 2j, 3], [2 - 1j, 4]], rtol=0, atol=1e-14
    )
    assert second.report.relative_residual < 1e-14


def test_complex_operator_is_preserved():
    result = solve_linear(
        [[1 + 1j, 0], [0, 2 - 1j]], [2j, 5], fields={"u": [1], "phi": [0]}
    )
    np.testing.assert_allclose(result.values, [1 + 1j, 2 + 1j], rtol=0, atol=1e-14)


@pytest.mark.parametrize(
    "fields",
    [
        {"u": [0]},
        {"u": [0, 0], "phi": [1]},
        {"u": [0, 1], "phi": [1]},
        {"u": [0, 2]},
        {"u": [0.0, 1.0]},
        {"": [0, 1]},
    ],
)
def test_invalid_field_partition(fields):
    with pytest.raises(ValidationError, match="field"):
        solve_linear(np.eye(2), [1, 2], fields=fields)


@pytest.mark.parametrize("scale", [[1, 0], [1, -1], [1, np.nan], [1], [1, 1j]])
def test_invalid_scaling(scale):
    with pytest.raises(ValidationError, match="scal"):
        solve_linear(np.eye(2), [1, 2], fields={"u": [0, 1]}, scale=scale)


@pytest.mark.parametrize(
    "a,b",
    [
        ([[1, 0]], [1]),
        ([[1, np.nan], [0, 1]], [1, 2]),
        (np.eye(2), [1]),
        (np.eye(2), [1, np.inf]),
    ],
)
def test_invalid_operator_rhs(a, b):
    with pytest.raises(ValidationError):
        solve_linear(a, b, fields={"u": [0, 1]})


def test_singular_and_non_spd_fail_without_fallback():
    with pytest.raises(ValidationError, match="singular"):
        solve_linear([[1, 1], [1, 1]], [0, 0], fields={"u": [0, 1]})
    with pytest.raises(ValidationError, match="positive"):
        solve_linear(
            [[1, 0], [0, -1]], [1, 1], fields={"u": [0, 1]}, kind="positive_definite"
        )


def test_empty_system_and_empty_field():
    result = solve_linear(np.empty((0, 0)), [], fields={"u": [], "phi": []})
    assert result.values.shape == (0,)
    assert result.report.relative_residual == 0


def test_small_electric_field_error_is_not_hidden_by_large_mechanical_units(
    monkeypatch,
):
    from fem_toolbox.solvers import linear_backend

    class BadFactor:
        U = sparse.eye(2, format="csc")

        def solve(self, rhs):
            return np.array([1.0, 0.0])

    # Mock the numerical-library boundary: global residual is ~1e-12, but the
    # electric field is completely wrong. Per-field residuals must reject it.
    monkeypatch.setattr(linear_backend, "splu", lambda matrix: BadFactor())
    with pytest.raises(ValidationError, match="field residual"):
        solve_linear(np.diag([1.0, 1e-12]), [1.0, 1e-12], fields={"u": [0], "phi": [1]})


def test_backend_exception_is_not_replaced_by_scipy():
    class FailingBackend:
        name = "unavailable"

        def factorize(self, operator, **kwargs):
            raise ValidationError("selected backend unavailable", field="solver")

    with pytest.raises(ValidationError, match="selected backend unavailable"):
        solve_linear([[1.0]], [1.0], fields={"u": [0]}, backend=FailingBackend())
