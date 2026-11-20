"""Coupled SA/SIMPLE gates and explicit failure states on physical wall domains."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain
from fem_toolbox.cfd.quad_system import QuadSystem
from fem_toolbox.cfd.rans import RansConvergenceError, solve_rans
from fem_toolbox.errors import ValidationError


def channel(mesh, *, mu=0.1):
    return FluidDomain.from_document(
        {
            "fluid_domain_schema_version": 1,
            "mesh": "mesh.json",
            "width_m": 0.7,
            "fluid": {"density_kg_m3": 1, "viscosity_Pa_s": mu},
            "boundaries": {
                "bottom": {"type": "wall"},
                "top": {"type": "wall"},
                "left": {
                    "type": "periodic",
                    "partner": "right",
                    "translation_m": [1, 0],
                    "pressure_drop_Pa": 0.01,
                },
                "right": {
                    "type": "periodic",
                    "partner": "left",
                    "translation_m": [-1, 0],
                    "pressure_drop_Pa": -0.01,
                },
            },
            "pressure_reference_cell": 0,
        },
        mesh=mesh,
    )


def test_coupled_low_re_wall_channel_and_fresh_budgets(quad_mesh):
    domain = channel(quad_mesh(6, 8))
    result = solve_rans(
        domain, initial_nu_tilde_m2_s=0.3, velocity_scale_m_s=1, length_scale_m=1
    )
    assert result.converged
    assert result.iterations > 1
    assert np.all(result.nu_tilde_m2_s >= 0)
    assert np.all(result.eddy_viscosity_m2_s >= 0)
    assert np.max(np.abs(result.sa_residual_m2_s2)) <= 1e-9
    assert np.max(np.abs(result.flow.momentum_residual_N_m3)) <= 1e-8
    assert np.max(np.abs(result.flow.mass_residual_kg_m3_s)) <= 1e-8
    assert abs(result.flow.power_balance_W) <= 0.7e-6
    assert np.max(result.wall_y_plus) <= 1
    # In this low-Re limit SA diffuses to zero; no turbulent validation claim.
    y = domain.mesh.cell_centres_m[:, 1]
    np.testing.assert_allclose(
        result.flow.velocity_m_s[:, 0], 0.05 * y * (1 - y), atol=1e-7, rtol=0
    )
    # Outward fluid traction balances the independently known integrated drive.
    walls = domain.wall_quality.boundary_indices
    faces = domain.mesh.boundary_face_indices[walls]
    traction = result.flow.wall_traction_Pa
    np.testing.assert_allclose(
        (traction * (domain.mesh.face_lengths_m[faces] * 0.7)[:, None]).sum(axis=0),
        [-0.007, 0],
        atol=1e-8,
    )
    with pytest.raises(ValueError):
        result.nu_tilde_m2_s.setflags(write=True)


def test_outer_budget_failure_keeps_complete_history_and_fields(quad_mesh):
    with pytest.raises(RansConvergenceError) as caught:
        solve_rans(
            channel(quad_mesh(6, 8)),
            initial_nu_tilde_m2_s=0.3,
            velocity_scale_m_s=1,
            length_scale_m=1,
            max_iterations=1,
        )
    assert len(caught.value.residual_history) == 1
    assert caught.value.nu_tilde_m2_s.shape == (48,)
    assert np.all(caught.value.nu_tilde_m2_s >= 0)


@pytest.mark.parametrize("initial", [0, -1, np.nan, 1e300])
def test_rans_mode_needs_an_explicit_positive_start(quad_mesh, initial):
    with pytest.raises(ValidationError):
        solve_rans(
            channel(quad_mesh(6, 8)),
            initial_nu_tilde_m2_s=initial,
            velocity_scale_m_s=1,
            length_scale_m=1,
        )


def test_actual_wall_resolution_failure_is_not_reported_as_success(quad_mesh):
    # Force/width fix total wall shear; a coarse low-Re mesh has y+>1 even
    # though the discrete coupled equations themselves can converge.
    with pytest.raises(RansConvergenceError, match="y_plus") as caught:
        solve_rans(
            channel(quad_mesh(6, 8), mu=0.004),
            initial_nu_tilde_m2_s=0.012,
            velocity_scale_m_s=1,
            length_scale_m=1,
        )
    assert len(caught.value.residual_history) > 1


def test_pressure_gauge_only_shifts_absolute_pressure(quad_mesh):
    domain = channel(quad_mesh(6, 8))
    options = dict(initial_nu_tilde_m2_s=0.3, velocity_scale_m_s=1, length_scale_m=1)
    base = solve_rans(domain, **options)
    shifted = solve_rans(domain, pressure_reference_Pa=17, **options)
    np.testing.assert_allclose(
        shifted.flow.pressure_Pa - base.flow.pressure_Pa, 17, atol=1e-12
    )
    np.testing.assert_array_equal(shifted.flow.velocity_m_s, base.flow.velocity_m_s)
    np.testing.assert_array_equal(shifted.nu_tilde_m2_s, base.nu_tilde_m2_s)


def test_viscosity_update_reassembles_full_stress_without_mutating_previous_state(
    quad_mesh,
):
    mesh = quad_mesh(6, 8)
    domain = channel(mesh, mu=1)
    n, nb = len(mesh.quad4), len(mesh.line2)
    original = QuadSystem(domain, np.ones(n), np.ones(nb))
    y = mesh.cell_centres_m[:, 1]
    yb = mesh.face_centres_m[mesh.boundary_face_indices, 1]
    velocity = np.column_stack((y * (1 - y), np.zeros(n))).ravel()
    before = original.diffusion.apply(velocity).copy()
    updated = original._with_viscosity(1 + y, 1 + yb)
    result = (
        updated.diffusion.apply(velocity).reshape(n, 2)
        / updated.fv.cell_volumes_m3[:, None]
    )
    # Independently differentiated -d[(1+y)(1-2y)]/dy=1+4y.
    np.testing.assert_allclose(result[:, 0], 1 + 4 * y, atol=1e-10, rtol=0)
    np.testing.assert_allclose(result[:, 1], 0, atol=1e-10, rtol=0)
    np.testing.assert_array_equal(original.diffusion.apply(velocity), before)
    np.testing.assert_array_equal(original.mu_cell, np.ones(n))
    with pytest.raises(ValidationError):
        original._with_viscosity(-np.ones(n), np.ones(nb))
