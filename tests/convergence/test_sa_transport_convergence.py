"""Predefined SA MMS from Doc/71, with independently differentiated forcing."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain, solve_rans
from fem_toolbox.cfd.sa_transport import solve_sa_transport
from tests.reference.rans_channel_mesh import channel_domain


def exact_source(points, nu):
    """Differentiate nu*(1+sin(pi*x)*sin(pi*y)); never call production SA terms."""
    x, y = points.T
    sx, sy = np.sin(np.pi * x), np.sin(np.pi * y)
    a = nu * (1 + sx * sy)
    ax = nu * np.pi * np.cos(np.pi * x) * sy
    ay = nu * np.pi * sx * np.cos(np.pi * y)
    laplacian = -2 * np.pi**2 * nu * sx * sy
    chi = a / nu
    fv1 = chi**3 / (chi**3 + 7.1**3)
    fv2 = 1 - chi / (1 + chi * fv1)
    # Omega=0 reduces Allmaras note 1c to max(Sbar,0), including zero branch.
    modified = np.maximum(a * fv2 / (0.41**2 * (y + 1) ** 2), 0)
    r = np.full_like(a, 10)
    np.divide(a, modified * 0.41**2 * (y + 1) ** 2, out=r, where=modified > 0)
    r = np.minimum(r, 10)
    g = r + 0.3 * (r**6 - r)
    fw = g * (65 / (g**6 + 64)) ** (1 / 6)
    cw1 = 0.1355 / 0.41**2 + 1.622 / (2 / 3)
    # div(U*a), U=(1,0), minus RHS with independently analytic derivatives.
    source = ax - 0.1355 * modified * a + cw1 * fw * (a / (y + 1)) ** 2
    source -= ((nu + a) * laplacian + 1.622 * (ax * ax + ay * ay)) / (2 / 3)
    return a, source


@pytest.mark.cfd_convergence
def test_sa_positive_mms_four_meshes(quad_mesh, record_property):
    errors = []
    nu = 0.1
    # Uref=Lref=1 SI, nu=0.1 m²/s, width=0.7 m fixed before result evaluation.
    for n in (16, 32, 64, 128):
        print(f"SA MMS: starting {n}x{n}", flush=True)
        mesh = quad_mesh(n, n)
        domain = FluidDomain.from_document(
            {
                "fluid_domain_schema_version": 1,
                "mesh": "mesh.json",
                "width_m": 0.7,
                "fluid": {"density_kg_m3": 1, "viscosity_Pa_s": nu},
                "boundaries": {
                    "bottom": {"type": "wall"},
                    "top": {"type": "wall"},
                    "left": {
                        "type": "inlet",
                        "velocity_m_s": [1, 0],
                        "nu_tilde_m2_s": nu,
                    },
                    "right": {"type": "pressure_outlet", "pressure_Pa": 0},
                },
            },
            mesh=mesh,
        )
        expected, source = exact_source(mesh.cell_centres_m, nu)
        boundary, _ = exact_source(mesh.face_centres_m[mesh.boundary_face_indices], nu)
        flux = mesh.face_normals[:, 0] * mesh.face_lengths_m * 0.7
        result = solve_sa_transport(
            domain,
            mass_flux_kg_s=flux,
            velocity_gradient_s_inv=np.zeros((n * n, 2, 2)),
            initial_nu_tilde_m2_s=nu,
            boundary_values_m2_s=boundary,
            wall_distance_m=mesh.cell_centres_m[:, 1] + 1,
            source_m2_s2=source,
            velocity_scale_m_s=1,
            length_scale_m=1,
            convection_scheme="linear-upwind",
        )
        relative = float(
            np.linalg.norm(result.nu_tilde_m2_s - expected) / np.linalg.norm(expected)
        )
        errors.append(relative)
        assert np.min(result.nu_tilde_m2_s) >= 0
        assert np.max(np.abs(result.residual_m2_s2)) <= 1e-9
        print(
            f"SA MMS {n}: relative RMS={relative:.9g}, iterations={result.iterations}, residual={result.residual_history[-1]:.6g}",
            flush=True,
        )
        record_property(f"n{n}_relative_rms", relative)
        record_property(f"n{n}_iterations", result.iterations)
    orders = np.log2(np.array(errors[:-1]) / errors[1:])
    print(f"SA MMS orders: {orders}", flush=True)
    assert errors[-1] < 1e-2
    assert np.all(orders > 0.8)


@pytest.mark.cfd_convergence
def test_rans_predictor_accuracy_does_not_limit_fresh_sa_convergence():
    # The original momentum-only inner tolerance stalled the independently
    # reconstructed SA residual at ~1e-8, above its required ~3.2e-10 gate.
    nu = 0.95 / 587.19
    result = solve_rans(
        channel_domain(64, kinematic_viscosity_m2_s=nu),
        initial_nu_tilde_m2_s=3 * nu,
        velocity_scale_m_s=20,
        length_scale_m=1,
        sa_relaxation=0.1,
    )
    assert result.converged
    assert np.max(np.abs(result.sa_residual_m2_s2)) <= 1e-8 * 20 * nu
    assert np.max(np.abs(result.flow.momentum_residual_N_m3)) <= 4e-6
