"""Independent SI and tensor-rotation checks for piezoelectric materials."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial


def test_strain_charge_conversion_and_roundtrip():
    # Independent scalar relations: C=1/s, e=d/s, kappa^S=kappa^T-d²/s.
    compliance = np.diag([0.5, 0.25, 0.2, 0.1, 0.125, 0.1])
    coupling = np.zeros((3, 6))
    coupling[2, 0] = 0.1
    permittivity_at_stress = np.diag([3.0, 4.0, 5.0])
    material = PiezoelectricMaterial.from_strain_charge(
        compliance, coupling, permittivity_at_stress
    )
    assert_allclose(material.stiffness, np.diag([2, 4, 5, 10, 8, 10]))
    assert_allclose(material.stress_coupling[2, 0], 0.2)
    assert_allclose(material.permittivity, np.diag([3, 4, 4.98]))
    s, d, k_t = material.as_strain_charge()
    assert_allclose(s, compliance, rtol=2e-14)
    assert_allclose(d, coupling, atol=2e-15)
    assert_allclose(k_t, permittivity_at_stress, rtol=2e-14)


def test_right_angle_rotation_and_inverse():
    axes = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
    stiffness = np.diag([10.0, 20, 30, 4, 5, 6])
    coupling = np.zeros((3, 6))
    coupling[2, 0] = 2.0
    permittivity = np.diag([3.0, 4, 5])
    material = PiezoelectricMaterial(stiffness, coupling, permittivity, axes=axes)
    assert_allclose(material.stiffness, np.diag([20, 10, 30, 5, 4, 6]))
    expected_coupling = np.zeros((3, 6))
    expected_coupling[2, 1] = 2.0
    assert_allclose(material.stress_coupling, expected_coupling, atol=1e-14)
    assert_allclose(material.permittivity, np.diag([4, 3, 5]))
    restored = PiezoelectricMaterial(
        material.stiffness,
        material.stress_coupling,
        material.permittivity,
        axes=axes.T,
    )
    assert_allclose(restored.stiffness, stiffness, atol=1e-13)
    assert_allclose(restored.stress_coupling, coupling, atol=1e-14)
    assert_allclose(restored.permittivity, permittivity, atol=1e-14)


def test_engineering_shear_rotation_preserves_energy_and_constitutive_response():
    angle = np.pi / 4
    axes = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    stiffness = np.diag([10.0, 20, 30, 4, 5, 6])
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.4
    material = PiezoelectricMaterial(stiffness, coupling, np.eye(3) * 2, axes=axes)
    # A local x extension becomes equal xx/yy and positive engineering xy shear.
    strain = np.array([0.5, 0.5, 0, 0, 0, 1.0])
    assert_allclose(0.5 * strain @ material.stiffness @ strain, 5.0, atol=1e-14)
    assert_allclose(material.stress_coupling @ strain, axes @ [0.4, 0, 0], atol=1e-14)


def test_general_rotation_matches_independent_tensor_constitutive_law():
    angle = np.pi / 6
    axes = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    stiffness = np.diag([11.0, 13, 17, 3, 4, 5])
    stiffness[0, 1] = stiffness[1, 0] = 1
    coupling = np.array(
        [[0.2, 0, 0, 0, 0.1, 0], [0, 0.3, 0, 0.2, 0, 0], [0, 0, 0.4, 0, 0, 0.1]]
    )
    permittivity = np.diag([2.0, 3, 4])
    material = PiezoelectricMaterial(stiffness, coupling, permittivity, axes=axes)

    strain_global = np.array(
        [[0.02, 0.03, -0.01], [0.03, -0.04, 0.02], [-0.01, 0.02, 0.01]]
    )
    field_global = np.array([0.5, -0.2, 0.7])
    strain_local = axes.T @ strain_global @ axes
    field_local = axes.T @ field_global
    local_voigt = np.array(
        [
            strain_local[0, 0],
            strain_local[1, 1],
            strain_local[2, 2],
            2 * strain_local[1, 2],
            2 * strain_local[0, 2],
            2 * strain_local[0, 1],
        ]
    )
    local_stress = stiffness @ local_voigt - coupling.T @ field_local
    local_stress_tensor = np.array(
        [
            [local_stress[0], local_stress[5], local_stress[4]],
            [local_stress[5], local_stress[1], local_stress[3]],
            [local_stress[4], local_stress[3], local_stress[2]],
        ]
    )
    expected_stress = axes @ local_stress_tensor @ axes.T
    global_voigt = np.array([0.02, -0.04, 0.01, 0.04, -0.02, 0.06])
    actual_stress = (
        material.stiffness @ global_voigt - material.stress_coupling.T @ field_global
    )
    actual_stress_tensor = np.array(
        [
            [actual_stress[0], actual_stress[5], actual_stress[4]],
            [actual_stress[5], actual_stress[1], actual_stress[3]],
            [actual_stress[4], actual_stress[3], actual_stress[2]],
        ]
    )
    assert_allclose(actual_stress_tensor, expected_stress, atol=2e-15)
    assert_allclose(
        material.stress_coupling @ global_voigt + material.permittivity @ field_global,
        axes @ (coupling @ local_voigt + permittivity @ field_local),
        atol=2e-15,
    )


def test_provenance_keeps_local_and_global_data_and_owns_inputs():
    stiffness = np.eye(6) * 10
    coupling = np.zeros((3, 6))
    permittivity = np.eye(3) * 2
    material = PiezoelectricMaterial(stiffness, coupling, permittivity)
    stiffness[:] = 0
    record = material.to_provenance_dict()
    assert record["source_form"] == "stress_charge"
    assert record["source"]["stiffness"][0][0] == 10
    assert record["global"]["stiffness"][0][0] == 10
    with pytest.raises(ValueError):
        material.stiffness[0, 0] = 0


@pytest.mark.parametrize(
    "stiffness,coupling,permittivity",
    [
        (np.eye(3), np.zeros((3, 6)), np.eye(3)),
        (np.eye(6), np.zeros((6, 3)), np.eye(3)),
        (np.eye(6), np.zeros((3, 6)), np.zeros((3, 3))),
        (np.eye(6) * -1, np.zeros((3, 6)), np.eye(3)),
        (np.eye(6) * np.nan, np.zeros((3, 6)), np.eye(3)),
        (np.eye(6) + np.triu(np.ones((6, 6)), 1), np.zeros((3, 6)), np.eye(3)),
    ],
)
def test_invalid_stress_charge_data(stiffness, coupling, permittivity):
    with pytest.raises(ValidationError):
        PiezoelectricMaterial(stiffness, coupling, permittivity)


def test_unstable_strain_charge_data_and_invalid_axes():
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 2
    with pytest.raises(ValidationError, match="permittivity"):
        PiezoelectricMaterial.from_strain_charge(np.eye(6), coupling, np.eye(3))
    with pytest.raises(ValidationError, match="axes"):
        PiezoelectricMaterial(np.eye(6), np.zeros((3, 6)), np.eye(3), axes=-np.eye(3))
