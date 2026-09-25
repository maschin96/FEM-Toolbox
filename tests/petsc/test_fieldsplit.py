"""FieldSplit: real PETSc algebra, convergence failures and nullspace diagnostics."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ValidationError
from fem_toolbox.solvers.linear_backend import solve_linear


def iterative(petsc_backend, **kwargs):
    from fem_toolbox.solvers.petsc_backend import PetscFieldSplitBackend

    return PetscFieldSplitBackend(**kwargs)


def coupled():
    # Two different Schur eigenvalues, so A11 preconditioning is not exact.
    return np.array(
        [
            [4.0, 0.0, 1.0, 0.0],
            [0.0, 3.0, 0.0, 2.0],
            [1.0, 0.0, -2.0, 0.0],
            [0.0, 2.0, 0.0, -1.0],
        ]
    )


@pytest.mark.parametrize("complex_operator", [False, True])
def test_schur_reordered_indices_scaling_and_multiple_rhs(
    petsc_backend, complex_operator
):
    backend = iterative(petsc_backend)
    a = coupled().astype(complex if complex_operator else float)
    if complex_operator:
        a[0, 2] += 0.5j
        a[2, 0] -= 0.5j
    perm = [2, 0, 3, 1]
    a = a[np.ix_(perm, perm)]
    expected = np.array([[2 + 1j, 3j], [1 - 2j, 2], [-1j, 4], [3, 1 + 1j]])
    factor = backend.factorize(
        a, fields={"phi": [0, 2], "u": [1, 3]}, scale=[0.5, 2, 0.25, 3]
    )
    result = factor.solve(a @ expected)
    assert_allclose(result.values, expected, rtol=1e-11, atol=1e-12)
    assert result.report.method.startswith("fgmres_fieldsplit_schur")
    assert result.report.iterations > 0
    assert max(result.report.field_relative_residuals.values()) < 1e-11
    assert_allclose(
        factor.solve(a @ expected[:, 0]).values, expected[:, 0], rtol=1e-11, atol=1e-12
    )
    factor.close()


def test_iteration_limit_is_reported_without_fallback(petsc_backend, monkeypatch):
    from fem_toolbox.solvers import linear_backend

    def forbidden(*args, **kwargs):
        pytest.fail("no SciPy fallback")

    monkeypatch.setattr(linear_backend, "splu", forbidden)
    backend = iterative(petsc_backend, max_iterations=1)
    with pytest.raises(
        ValidationError,
        match=r"DIVERGED_MAX_IT.*reason=-3.*iterations=1.*max_iterations=1",
    ):
        solve_linear(
            coupled(),
            [1, 2, 3, 4],
            fields={"u": [0, 1], "phi": [2, 3]},
            backend=backend,
        )


def test_single_field_empty_field_and_empty_system(petsc_backend):
    backend = iterative(petsc_backend)
    result = solve_linear(
        np.diag([2.0, 3.0]), [4, 9], fields={"u": [0, 1], "phi": []}, backend=backend
    )
    assert_allclose(result.values, [2, 3], rtol=0, atol=1e-13)
    assert result.report.method == "fgmres_lu"
    empty = solve_linear(
        np.empty((0, 0)), [], fields={"u": [], "phi": []}, backend=backend
    )
    assert empty.report.iterations == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"rtol": 0},
        {"rtol": 1},
        {"rtol": np.nan},
        {"atol": -1},
        {"max_iterations": 0},
        {"max_iterations": True},
        {"restart": 0},
        {"restart": 2.5},
    ],
)
def test_invalid_iteration_options(petsc_backend, kwargs):
    with pytest.raises(ValidationError):
        iterative(petsc_backend, **kwargs)


def test_unconstrained_algebraic_block_diagnosed_even_with_zero_rhs(petsc_backend):
    backend = iterative(petsc_backend)
    with pytest.raises(ValidationError, match=r"singular|nullspace|pivot"):
        solve_linear(
            [[1.0, -1.0], [-1.0, 1.0]], [0, 0], fields={"phi": [0, 1]}, backend=backend
        )


def test_original_field_residual_can_reject_loose_ksp_convergence(petsc_backend):
    backend = iterative(petsc_backend, rtol=0.9)
    with pytest.raises(ValidationError, match="unscaled backend field residual"):
        solve_linear(
            coupled(),
            [1, 2, 3, 4],
            fields={"u": [0, 1], "phi": [2, 3]},
            backend=backend,
        )


@pytest.mark.parametrize("source", [0.0, 1.0])
def test_thermal_and_electric_nullspace_diagnostics(cube_mesh, petsc_backend, source):
    from fem_toolbox.errors import ElectricGaugeError, ThermalNullspaceError
    from fem_toolbox.materials import DielectricMaterial, ThermalMaterial
    from fem_toolbox.solvers.electrostatics import solve_electrostatics
    from fem_toolbox.solvers.heat import solve_heat

    backend = iterative(petsc_backend)
    mesh = cube_mesh(1)
    with pytest.raises(ThermalNullspaceError) as thermal:
        solve_heat(
            mesh, ThermalMaterial.isotropic(2), [], source=source, backend=backend
        )
    assert thermal.value.compatible == (source == 0)
    assert_allclose(thermal.value.net_power, source, rtol=0, atol=1e-14)
    with pytest.raises(ElectricGaugeError) as electric:
        solve_electrostatics(
            mesh,
            DielectricMaterial.isotropic(2),
            [],
            volume_charge=source,
            backend=backend,
        )
    assert electric.value.compatible == (source == 0)
    assert_allclose(electric.value.net_charge, source, rtol=0, atol=1e-14)


def test_unrestrained_rigid_modes_are_not_hidden_by_iterative_solve(
    cube_mesh, petsc_backend
):
    from fem_toolbox.materials import ElasticMaterial
    from fem_toolbox.solvers.elasticity import solve_elasticity

    with pytest.raises(ValidationError, match="unrestrained rigid modes"):
        solve_elasticity(
            cube_mesh(1),
            ElasticMaterial.isotropic(100, 0.25),
            [],
            backend=iterative(petsc_backend),
        )


def test_more_than_two_nonempty_fields_rejected(petsc_backend):
    with pytest.raises(ValidationError, match="at most two nonempty fields"):
        solve_linear(
            np.eye(3),
            [1, 2, 3],
            fields={"u": [0], "T": [1], "phi": [2]},
            backend=iterative(petsc_backend),
        )


def test_iteration_count_resets_and_failed_factor_is_reusable(petsc_backend):
    backend = iterative(petsc_backend, max_iterations=1)
    factor = backend.factorize(coupled(), fields={"u": [0, 1], "phi": [2, 3]})
    with pytest.raises(ValidationError, match="DIVERGED_MAX_IT"):
        factor.solve([1, 2, 3, 4])
    zero = factor.solve(np.zeros(4))
    assert_allclose(zero.values, 0, rtol=0, atol=0)
    assert zero.report.iterations == 0
    assert factor.solve(np.empty((4, 0))).report.iterations == 0
    factor.close()
    with pytest.raises(ValidationError, match="closed"):
        factor.solve(np.zeros(4))


def test_global_petsc_options_do_not_reconfigure_owned_subsolvers(petsc_backend):
    from petsc4py import PETSc

    options = PETSc.Options()
    key = "fieldsplit_u_pc_type"
    previous = options.getString(key) if options.hasName(key) else None
    options[key] = "not_a_valid_pc"
    try:
        result = solve_linear(
            coupled(),
            [1, 2, 3, 4],
            fields={"u": [0, 1], "phi": [2, 3]},
            backend=iterative(petsc_backend),
        )
        assert result.report.relative_residual < 1e-11
    finally:
        if previous is None:
            del options[key]
        else:
            options[key] = previous
