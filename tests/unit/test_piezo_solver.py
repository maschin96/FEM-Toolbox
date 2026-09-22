"""Independent small-system and error-path references for the piezo solver."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ElectricGaugeError, ValidationError
from fem_toolbox.materials import PiezoelectricMaterial, ThermalExpansion
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric


def single_tet():
    return Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], tetra4=[[0, 1, 2, 3]]
    )


def material(c=10.0, e=2.0, k=3.0):
    stiffness = np.diag([c, 20, 30, 4, 5, 6])
    coupling = np.zeros((3, 6))
    coupling[0, 0] = e
    return PiezoelectricMaterial(stiffness, coupling, np.diag([k, 4, 5]))


def two_free_dofs(mesh):
    dofs = piezoelectric_dofs(mesh)
    ux = int(dofs.indices("u", [1], ["x"])[0])
    phi = int(dofs.indices("phi", [1])[0])
    fixed = np.setdiff1d(np.arange(dofs.size), [ux, phi])
    return dofs, ux, phi, [Dirichlet(fixed, np.zeros(len(fixed)))]


def test_two_by_two_signed_saddle_reference():
    mesh = single_tet()
    dofs, ux, phi, conditions = two_free_dofs(mesh)
    assert dofs.block("u") == slice(0, 12)
    assert dofs.block("phi") == slice(12, 16)
    result = solve_piezoelectric(
        mesh,
        material(),
        conditions,
        point_loads={1: [1, 0, 0]},
        volume_charge=-2.4,
    )
    expected = np.linalg.solve(np.array([[10 / 6, 2 / 6], [2 / 6, -3 / 6]]), [1, 0.1])
    assert_allclose(result.displacement[1, 0], expected[0], rtol=2e-14)
    assert_allclose(result.potential[1], expected[1], rtol=2e-14)
    assert_allclose(result.residual_u.ravel()[ux], 0, atol=2e-15)
    assert_allclose(result.residual_phi[phi - 12], 0, atol=2e-15)
    assert result.relative_residual_u < 1e-13
    assert result.relative_residual_phi < 1e-13
    assert_allclose(result.force_balance, 0, atol=2e-15)
    assert_allclose(result.charge_balance, 0, atol=2e-15)
    assert result.provenance.solver_parameters["backend"] == "scipy.sparse.linalg.splu"


def test_field_block_scaling_preserves_solution_with_large_coefficient_ratio():
    mesh = single_tet()
    _, _, _, conditions = two_free_dofs(mesh)
    mat = material(c=1e6, e=10, k=1e-2)
    kwargs = dict(point_loads={1: [1, 0, 0]}, volume_charge=-0.24)
    scaled = solve_piezoelectric(mesh, mat, conditions, block_scaling=True, **kwargs)
    plain = solve_piezoelectric(mesh, mat, conditions, block_scaling=False, **kwargs)
    assert_allclose(scaled.displacement, plain.displacement, rtol=2e-12, atol=1e-20)
    assert_allclose(scaled.potential, plain.potential, rtol=2e-12, atol=1e-20)
    assert scaled.relative_residual_u < 1e-10
    assert scaled.relative_residual_phi < 1e-10
    assert scaled.mechanical_scale != scaled.electric_scale


def test_physical_si_block_scale_against_dense_reference():
    mesh = single_tet()
    _, _, _, conditions = two_free_dofs(mesh)
    result = solve_piezoelectric(
        mesh,
        material(c=1e11, e=1.0, k=1e-11),
        conditions,
        point_loads={1: [1, 0, 0]},
        volume_charge=-2.4e-11,
    )
    dense = np.array([[1e11 / 6, 1 / 6], [1 / 6, -1e-11 / 6]])
    expected = np.linalg.solve(dense, [1, 1e-12])
    assert_allclose(result.displacement[1, 0], expected[0], rtol=2e-14)
    assert_allclose(result.potential[1], expected[1], rtol=2e-14)
    assert_allclose(result.mechanical_scale, 1e11 / 6, rtol=1e-15)
    assert_allclose(result.electric_scale, 1e-11 / 6, rtol=1e-15)
    assert result.relative_residual_u < 1e-13
    assert result.relative_residual_phi < 1e-13


def test_missing_gauge_and_rigid_support_are_diagnosed():
    mesh = single_tet()
    dofs = piezoelectric_dofs(mesh)
    with pytest.raises(ElectricGaugeError):
        solve_piezoelectric(mesh, material(), [Dirichlet(np.arange(12), np.zeros(12))])
    with pytest.raises(ValidationError, match="rigid"):
        solve_piezoelectric(
            mesh,
            material(),
            [Dirichlet(dofs.indices("phi", np.arange(4)), np.zeros(4))],
        )


def test_cg_is_not_an_allowed_saddle_solver():
    mesh = single_tet()
    _, _, _, conditions = two_free_dofs(mesh)
    with pytest.raises(ValidationError, match="direct"):
        solve_piezoelectric(mesh, material(), conditions, method="cg")


def test_invalid_material_and_empty_mesh():
    mesh = single_tet()
    _, _, _, conditions = two_free_dofs(mesh)
    with pytest.raises(ValidationError, match="PiezoelectricMaterial"):
        solve_piezoelectric(mesh, [material(), material()], conditions)
    with pytest.raises(ValidationError, match="Tet4"):
        solve_piezoelectric(
            Mesh(points=mesh.points, tetra4=np.empty((0, 4), dtype=int)), material(), []
        )


def test_thermal_electric_contribution_requires_explicit_switch():
    mesh = single_tet()
    fixed = [Dirichlet(np.arange(16), np.zeros(16))]
    expansion = ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300)
    kwargs = dict(expansion=expansion, temperature=np.full(4, 310.0))
    default = solve_piezoelectric(mesh, material(), fixed, **kwargs)
    enabled = solve_piezoelectric(
        mesh, material(), fixed, electric_thermal_strain=True, **kwargs
    )
    assert_allclose(default.stress[0, 0], -1, atol=1e-14)
    assert_allclose(enabled.stress[0, 0], -1, atol=1e-14)
    assert_allclose(default.electric_displacement, 0, atol=1e-14)
    assert_allclose(enabled.electric_displacement[0], [-0.2, 0, 0], atol=1e-14)
    assert_allclose(default.reactions_phi, 0, atol=1e-14)
    assert_allclose(enabled.reactions_phi, [1 / 30, -1 / 30, 0, 0], atol=1e-14)
