"""Independent single-Tet4 references for coupled piezoelectric operators."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.piezoelectric import (
    piezoelectric_matrix,
    piezoelectric_response,
    piezoelectric_thermal_load,
)
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial

POINTS = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])


def material():
    stiffness = np.diag([10.0, 20, 30, 4, 5, 6])
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 2.0
    return PiezoelectricMaterial(stiffness, coupling, np.diag([3.0, 4, 5]))


def test_block_symmetry_reciprocity_and_signed_dielectric():
    geometry = Tet4Geometry(POINTS)
    matrix = piezoelectric_matrix(geometry, material())
    assert matrix.shape == (16, 16)
    assert_allclose(matrix, matrix.T, atol=2e-15)
    assert_allclose(matrix[12:, :12], matrix[:12, 12:].T, atol=2e-15)
    # Unit x-gradient from phi=[0,1,0,0] gives -V*kappa_xx on the x node.
    potential = np.array([0.0, 1, 0, 0])
    assert_allclose(potential @ matrix[12:, 12:] @ potential, -0.5, atol=1e-15)
    assert_allclose(matrix[:12, 12:] @ np.ones(4), 0, atol=2e-15)
    assert_allclose(matrix[12:, :12] @ np.tile([1.0, 2, 3], 4), 0, atol=2e-15)


def test_affine_coupled_patch_matches_constant_tractions_and_charge():
    geometry = Tet4Geometry(POINTS)
    a, b = 0.2, 0.3
    displacement = np.column_stack([a * POINTS[:, 0], np.zeros(4), np.zeros(4)])
    potential = b * POINTS[:, 0]
    nodal = np.concatenate([displacement.ravel(), potential])
    residual = piezoelectric_matrix(geometry, material()) @ nodal
    # sigma_xx=10a+2b=2.6 Pa; D_x=2a-3b=-0.5 C/m².
    expected = np.zeros(16)
    expected[[0, 3]] = [-2.6 / 6, 2.6 / 6]
    expected[[12, 13]] = [0.5 / 6, -0.5 / 6]
    assert_allclose(residual, expected, atol=1e-15)
    response = piezoelectric_response(geometry, material(), displacement, potential)
    assert_allclose(response.strain, [a, 0, 0, 0, 0, 0], atol=1e-15)
    assert_allclose(response.electric_field, [-b, 0, 0], atol=1e-15)
    assert_allclose(response.stress, [2.6, 0, 0, 0, 0, 0], atol=1e-15)
    assert_allclose(response.electric_displacement, [-0.5, 0, 0], atol=1e-15)
    assert_allclose(response.mechanical_energy, (5 * a * a) / 6, atol=1e-15)
    assert_allclose(response.coupling_energy, (2 * a * b) / 6, atol=1e-15)
    assert_allclose(response.dielectric_enthalpy, (-1.5 * b * b) / 6, atol=1e-15)


def test_energy_gradient_matches_matrix_with_central_differences():
    geometry = Tet4Geometry(POINTS)
    m = material()
    displacement = np.array(
        [[0, 0, 0], [0.02, 0.01, 0], [0, -0.01, 0.02], [0, 0, 0.03]]
    )
    potential = np.array([0.1, 0.4, -0.2, 0.05])
    nodal = np.concatenate([displacement.ravel(), potential])

    def energy(values):
        return piezoelectric_response(
            geometry, m, values[:12].reshape(4, 3), values[12:]
        ).total_enthalpy

    h = 1e-6
    for index in [0, 3, 7, 12, 13, 15]:
        plus, minus = nodal.copy(), nodal.copy()
        plus[index] += h
        minus[index] -= h
        derivative = (energy(plus) - energy(minus)) / (2 * h)
        assert_allclose(
            derivative, (piezoelectric_matrix(geometry, m) @ nodal)[index], atol=2e-10
        )


def test_explicit_thermal_strain_switch_changes_electric_rhs_and_d_field():
    geometry = Tet4Geometry(POINTS)
    m = material()
    thermal = np.array([0.1, 0, 0, 0, 0, 0])
    uncoupled = piezoelectric_thermal_load(geometry, m, thermal)
    coupled = piezoelectric_thermal_load(
        geometry, m, thermal, electric_thermal_strain=True
    )
    assert_allclose(uncoupled[[0, 3]], [-1 / 6, 1 / 6], atol=1e-15)
    assert_allclose(uncoupled[12:], 0, atol=0)
    assert_allclose(coupled[12:], [-0.2 / 6, 0.2 / 6, 0, 0], atol=1e-15)
    displacement = np.zeros((4, 3))
    potential = np.zeros(4)
    response = piezoelectric_response(
        geometry,
        m,
        displacement,
        potential,
        thermal_strain=thermal,
        electric_thermal_strain=True,
    )
    assert_allclose(response.stress, [-1, 0, 0, 0, 0, 0], atol=1e-15)
    assert_allclose(response.electric_displacement, [-0.2, 0, 0], atol=1e-15)
    assert_allclose(response.total_enthalpy, 0.5 * 10 * 0.1**2 / 6, atol=1e-15)


@pytest.mark.parametrize("electric_thermal_strain", [False, True])
def test_thermal_energy_derivative_matches_signed_thermal_rhs(electric_thermal_strain):
    geometry = Tet4Geometry(POINTS)
    m = material()
    thermal = np.array([0.1, 0, 0, 0, 0, 0])
    nodal = np.zeros(16)
    h = 1e-6
    thermal_rhs = piezoelectric_thermal_load(
        geometry, m, thermal, electric_thermal_strain=electric_thermal_strain
    )
    for index in [0, 3, 12, 13]:
        plus, minus = nodal.copy(), nodal.copy()
        plus[index] += h
        minus[index] -= h
        energy_plus = piezoelectric_response(
            geometry,
            m,
            plus[:12].reshape(4, 3),
            plus[12:],
            thermal_strain=thermal,
            electric_thermal_strain=electric_thermal_strain,
        ).total_enthalpy
        energy_minus = piezoelectric_response(
            geometry,
            m,
            minus[:12].reshape(4, 3),
            minus[12:],
            thermal_strain=thermal,
            electric_thermal_strain=electric_thermal_strain,
        ).total_enthalpy
        assert_allclose(
            (energy_plus - energy_minus) / (2 * h), -thermal_rhs[index], atol=2e-12
        )


def test_bad_shapes_and_missing_thermal_model_are_rejected():
    geometry = Tet4Geometry(POINTS)
    m = material()
    with pytest.raises(ValidationError):
        piezoelectric_response(geometry, m, np.zeros(12), np.zeros(4))
    with pytest.raises(ValidationError):
        piezoelectric_response(geometry, m, np.zeros((4, 3)), np.zeros(3))
    with pytest.raises(ValidationError):
        piezoelectric_thermal_load(geometry, m, np.zeros(3))
    with pytest.raises(ValidationError, match="thermal"):
        piezoelectric_response(
            geometry,
            m,
            np.zeros((4, 3)),
            np.zeros(4),
            electric_thermal_strain=True,
        )
