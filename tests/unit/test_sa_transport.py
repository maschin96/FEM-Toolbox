"""Conservative scalar transport checked without a RANS momentum solution."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain
from fem_toolbox.cfd.sa_transport import SAConvergenceError, solve_sa_transport
from fem_toolbox.errors import ValidationError


def domain_for(mesh):
    return FluidDomain.from_document(
        {
            "fluid_domain_schema_version": 1,
            "mesh": "mesh.json",
            "width_m": 0.7,
            "fluid": {"density_kg_m3": 2, "viscosity_Pa_s": 0.02},
            "boundaries": {name: {"type": "wall"} for name in mesh.boundary_sets},
            "pressure_reference_cell": 0,
        },
        mesh=mesh,
    )


def test_affine_sa_diffusion_with_independent_exact_source(quad_mesh):
    mesh = quad_mesh(8, 8, shear=0.2)
    domain = domain_for(mesh)
    nu = 0.01
    exact = nu * (2 + mesh.cell_centres_m[:, 1])
    boundary = nu * (2 + mesh.face_centres_m[mesh.boundary_face_indices, 1])
    # chi>=2, Omega=0: S_tilde=0, r=10 and a/d=nu everywhere.
    # div(K grad a)=(nu**2)/sigma; full gradient-square source retained.
    g = 10 + 0.3 * (10**6 - 10)
    fw = g * (65 / (g**6 + 64)) ** (1 / 6)
    cw1 = 0.1355 / 0.41**2 + (1 + 0.622) / (2 / 3)
    source = nu**2 * (cw1 * fw - (1 + 0.622) / (2 / 3))
    result = solve_sa_transport(
        domain,
        mass_flux_kg_s=np.zeros(len(mesh.faces)),
        velocity_gradient_s_inv=np.zeros((len(mesh.quad4), 2, 2)),
        initial_nu_tilde_m2_s=nu,
        boundary_values_m2_s=boundary,
        wall_distance_m=2 + mesh.cell_centres_m[:, 1],
        source_m2_s2=source,
        velocity_scale_m_s=1,
        length_scale_m=1,
    )
    assert result.converged
    np.testing.assert_allclose(result.nu_tilde_m2_s, exact, atol=2e-8, rtol=0)
    assert np.max(np.abs(result.residual_m2_s2)) <= 1e-10
    assert abs(result.global_balance_m5_s2) <= 1e-10 * 0.7
    assert result.iterations > 1
    assert result.residual_history[-1] <= 1e-10
    assert min(result.nu_tilde_m2_s) >= 0
    assert (
        result.provenance.to_dict()["solver"]["parameters"]["wall_distance_source"]
        == "explicit_verification_field"
    )
    with pytest.raises(ValueError):
        result.nu_tilde_m2_s.setflags(write=True)


def test_zero_solution_is_explicitly_scalar_model_and_reports_fresh_residuals(
    quad_mesh,
):
    mesh = quad_mesh(4, 4)
    result = solve_sa_transport(
        domain_for(mesh),
        mass_flux_kg_s=np.zeros(len(mesh.faces)),
        velocity_gradient_s_inv=np.zeros((len(mesh.quad4), 2, 2)),
        initial_nu_tilde_m2_s=0,
        velocity_scale_m_s=1,
        length_scale_m=1,
    )
    np.testing.assert_array_equal(result.nu_tilde_m2_s, 0)
    np.testing.assert_array_equal(result.residual_m2_s2, 0)
    assert result.global_balance_m5_s2 == 0
    parameters = result.provenance.to_dict()["solver"]["parameters"]
    assert parameters["constants"]["kappa"] == 0.41
    assert parameters["constants"]["sigma"] == 2 / 3


def test_budget_failure_carries_latest_state_and_history(quad_mesh):
    mesh = quad_mesh(4, 4)
    with pytest.raises(SAConvergenceError) as caught:
        solve_sa_transport(
            domain_for(mesh),
            mass_flux_kg_s=np.zeros(len(mesh.faces)),
            velocity_gradient_s_inv=np.zeros((len(mesh.quad4), 2, 2)),
            initial_nu_tilde_m2_s=0.03,
            velocity_scale_m_s=1,
            length_scale_m=1,
            max_iterations=1,
        )
    assert len(caught.value.residual_history) == 1
    assert caught.value.nu_tilde_m2_s.shape == (len(mesh.quad4),)
    assert np.all(caught.value.nu_tilde_m2_s >= 0)


@pytest.mark.parametrize(
    "change",
    [
        {"velocity_scale_m_s": 0},
        {"max_iterations": True},
        {"convection_scheme": "magic"},
        {"wall_distance_m": 0},
    ],
)
def test_scalar_solver_invalid_settings_are_specific(quad_mesh, change):
    mesh = quad_mesh(4, 4)
    options = dict(
        mass_flux_kg_s=np.zeros(len(mesh.faces)),
        velocity_gradient_s_inv=np.zeros((len(mesh.quad4), 2, 2)),
        initial_nu_tilde_m2_s=0,
        velocity_scale_m_s=1,
        length_scale_m=1,
    )
    options.update(change)
    with pytest.raises(ValidationError):
        solve_sa_transport(domain_for(mesh), **options)
