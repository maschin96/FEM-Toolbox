"""Independent point values, dimensions and admissibility for SA-noft2."""

import numpy as np
import pytest

from fem_toolbox.cfd.sa_model import positive_sa_update, sa_terms
from fem_toolbox.errors import ValidationError


def evaluate(a=1.0, nu=1.0, d=1.0, omega=2.0, gradient=(3.0, 4.0)):
    return sa_terms(
        [a],
        velocity_gradient_s_inv=[[[0, omega], [0, 0]]],
        gradient_nu_tilde_m_s=[gradient],
        wall_distance_m=d,
        kinematic_viscosity_m2_s=nu,
    )


def test_zero_state_and_zero_vorticity_are_finite_without_clipping():
    result = evaluate(a=0, omega=0, gradient=(0, 0))
    assert result.eddy_viscosity_m2_s[0] == 0
    assert result.production_m2_s2[0] == 0
    assert result.destruction_m2_s2[0] == 0
    assert result.gradient_source_m2_s2[0] == 0
    assert result.diffusivity_m2_s[0] == pytest.approx(1.5)
    assert result.r[0] == 10
    assert np.isfinite(result.wall_function[0])


def test_point_terms_against_independent_hand_algebra():
    result = evaluate()
    fv1 = 1 / (1 + 7.1**3)
    # chi=1 makes f_v2=f_v1/(1+f_v1); vorticity is exactly two.
    s = 2 + fv1 / ((1 + fv1) * 0.41**2)
    r = 1 / (s * 0.41**2)
    g = 0.7 * r + 0.3 * r**6
    fw = g * (65 / (g**6 + 64)) ** (1 / 6)
    cw1 = 0.1355 / 0.41**2 + 1.622 * 1.5
    assert result.eddy_viscosity_m2_s[0] == pytest.approx(fv1, rel=2e-14)
    assert result.modified_vorticity_s_inv[0] == pytest.approx(s, rel=2e-14)
    assert result.production_m2_s2[0] == pytest.approx(0.1355 * s, rel=2e-14)
    assert result.destruction_m2_s2[0] == pytest.approx(cw1 * fw, rel=2e-14)
    assert result.gradient_source_m2_s2[0] == pytest.approx(0.622 * 1.5 * 25)
    assert result.diffusivity_m2_s[0] == pytest.approx(3)


def test_allmaras_negative_sbar_branch_and_zero_omega_guard():
    result = evaluate(a=0.01, nu=0.001, omega=0.01, gradient=(0, 0))
    # Below the branch threshold, S remains positive without flooring/clipping.
    assert 0 < result.modified_vorticity_s_inv[0] < 0.003
    zero = evaluate(a=0.01, nu=0.001, omega=0, gradient=(0, 0))
    assert zero.modified_vorticity_s_inv[0] == 0
    assert zero.r[0] == 10
    assert zero.destruction_m2_s2[0] > 0


def test_rigid_rotation_uses_vorticity_and_not_strain_magnitude():
    result = sa_terms(
        [0.01],
        velocity_gradient_s_inv=[[[0, -2], [2, 0]]],
        gradient_nu_tilde_m_s=[[0, 0]],
        wall_distance_m=1,
        kinematic_viscosity_m2_s=0.001,
    )
    assert result.production_m2_s2[0] > 0
    assert result.modified_vorticity_s_inv[0] > 3.9


def test_dimension_scaling_and_immutable_outputs():
    base = evaluate()
    # Lengths x3 and times x2: nu x9/2, velocity gradients /2,
    # grad(nu_tilde) x3/2, all SA source terms x9/4.
    scaled = evaluate(a=4.5, nu=4.5, d=3, omega=1, gradient=(4.5, 6))
    for key in ("production_m2_s2", "destruction_m2_s2", "gradient_source_m2_s2"):
        np.testing.assert_allclose(getattr(scaled, key), 2.25 * getattr(base, key))
    np.testing.assert_allclose(
        scaled.eddy_viscosity_m2_s, 4.5 * base.eddy_viscosity_m2_s
    )
    with pytest.raises(ValueError):
        base.eddy_viscosity_m2_s[0] = 10
    with pytest.raises(ValueError):
        base.eddy_viscosity_m2_s.setflags(write=True)


@pytest.mark.parametrize(
    "change", [{"a": -1}, {"nu": 0}, {"d": 0}, {"a": np.nan}, {"a": 1e300}]
)
def test_unphysical_or_unrepresentable_states_are_diagnosed(change):
    with pytest.raises(ValidationError):
        evaluate(**change)


def test_shapes_and_boolean_data_are_rejected():
    for value in ([1, 2], [[1]], [True]):
        with pytest.raises(ValidationError):
            sa_terms(
                value,
                velocity_gradient_s_inv=np.zeros((1, 2, 2)),
                gradient_nu_tilde_m_s=[[0, 0]],
                wall_distance_m=1,
                kinematic_viscosity_m2_s=1,
            )


def test_positive_update_halves_one_common_step_without_clipping():
    old = np.array([1.0, 4.0])
    candidate = np.array([-4.0, 2.0])
    state, alpha, backtracks = positive_sa_update(old, candidate)
    assert alpha == 0.125
    assert backtracks == 2
    np.testing.assert_array_equal(state, old + alpha * (candidate - old))
    np.testing.assert_array_equal(old, [1, 4])
    with pytest.raises(ValidationError, match=r"cell.*0"):
        positive_sa_update([0, 1], [-1, 2])


def test_nonfinite_update_does_not_masquerade_as_convergence():
    with pytest.raises(ValidationError):
        positive_sa_update([1], [np.inf])
