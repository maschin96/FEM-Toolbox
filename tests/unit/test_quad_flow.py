"""Quad flow tests with independently prescribed SI channel references."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain, solve_quad_flow
from fem_toolbox.errors import ValidationError


def domain_for(mesh, *, periodic=False, rho=1, mu=1, outlet=0):
    boundaries = {"bottom": {"type": "wall"}, "top": {"type": "wall"}}
    if periodic:
        boundaries.update(
            left={
                "type": "periodic",
                "partner": "right",
                "translation_m": [1, 0],
                "pressure_drop_Pa": 2,
            },
            right={
                "type": "periodic",
                "partner": "left",
                "translation_m": [-1, 0],
                "pressure_drop_Pa": -2,
            },
        )
    else:
        b = mesh.face_centres_m[mesh.boundary_face_indices[mesh.boundary_sets["left"]]]
        boundaries.update(
            left={
                "type": "inlet",
                "velocity_m_s": np.column_stack(
                    (b[:, 1] * (1 - b[:, 1]), np.zeros(len(b)))
                ).tolist(),
                "nu_tilde_m2_s": 0,
            },
            right={"type": "pressure_outlet", "pressure_Pa": outlet},
        )
    doc = dict(
        fluid_domain_schema_version=1,
        mesh="unused.json",
        width_m=0.2,
        fluid={"density_kg_m3": rho, "viscosity_Pa_s": mu},
        boundaries=boundaries,
    )
    if periodic:
        doc["pressure_reference_cell"] = 0
    return FluidDomain.from_document(doc, mesh=mesh)


def test_quad_channel_fields_final_mass_and_momentum(quad_mesh):
    mesh = quad_mesh(8, 8)
    r = solve_quad_flow(domain_for(mesh), velocity_scale_m_s=1, length_scale_m=1)
    assert r.converged
    assert (
        np.max(
            np.abs(
                r.velocity_m_s[:, 0]
                - mesh.cell_centres_m[:, 1] * (1 - mesh.cell_centres_m[:, 1])
            )
        )
        < 0.01
    )
    assert np.max(np.abs(r.mass_residual_kg_m3_s)) < 1e-8
    assert np.max(np.abs(r.momentum_residual_N_m3)) < 1e-8
    assert abs(r.global_mass_balance_kg_s) < 2e-11
    assert np.max(np.abs(r.global_momentum_balance_N)) < 2e-9
    assert abs(r.power_balance_W) < 2e-7
    assert r.inlet_mass_flow_kg_s == pytest.approx(r.outlet_mass_flow_kg_s, rel=1e-9)
    parameters = r.provenance.solver_parameters
    assert len(parameters["geometry_sha256"]) == 64
    assert parameters["residual_history"] == r.residual_history
    with pytest.raises(ValueError):
        r.velocity_m_s.setflags(write=True)


def test_periodic_pressure_driven_channel_and_gauge(quad_mesh):
    mesh = quad_mesh(4, 16)
    a = solve_quad_flow(
        domain_for(mesh, periodic=True), velocity_scale_m_s=1, length_scale_m=1
    )
    b = solve_quad_flow(
        domain_for(mesh, periodic=True),
        velocity_scale_m_s=1,
        length_scale_m=1,
        pressure_reference_Pa=5,
    )
    y = mesh.cell_centres_m[:, 1]
    assert np.max(np.abs(a.velocity_m_s[:, 0] - y * (1 - y))) < 0.003
    np.testing.assert_allclose(b.velocity_m_s, a.velocity_m_s, atol=1e-11, rtol=0)
    np.testing.assert_allclose(b.pressure_Pa - a.pressure_Pa, 5, atol=1e-10, rtol=0)
    assert a.provenance.input_sha256 != b.provenance.input_sha256


def test_nonconvergence_and_limits(quad_mesh):
    d = domain_for(quad_mesh(4, 4))
    with pytest.raises(ValidationError, match="converg"):
        solve_quad_flow(d, velocity_scale_m_s=1, length_scale_m=1, max_iterations=1)
    with pytest.raises(ValidationError, match="limit"):
        solve_quad_flow(d, velocity_scale_m_s=1, length_scale_m=1, max_unknowns=1)


def test_component_momentum_diagonal_pressure_mobility(quad_mesh):
    from fem_toolbox.cfd import QuadFVGeometry

    fv = QuadFVGeometry(quad_mesh(2, 1), width_m=2)
    # x-face: d_x=(1/2,1/4), d_y=(1/8,1/16); only d_x enters n.D.n.
    matrix = fv.pressure_correction_matrix(
        momentum_diagonal_kg_s=[[2, 8], [4, 16]], density_kg_m3=3
    )
    np.testing.assert_allclose(
        matrix.toarray(), [[4.5, -4.5], [-4.5, 4.5]], atol=1e-13, rtol=0
    )


def test_sheared_periodic_flow_and_symmetry(quad_mesh):
    mesh = quad_mesh(4, 16, shear=np.tan(np.pi / 9), size_m=(1, 0.5))
    d = domain_for(mesh, periodic=True)
    doc = dict(
        fluid_domain_schema_version=1,
        mesh="unused.json",
        width_m=d.width_m,
        fluid={"density_kg_m3": 1, "viscosity_Pa_s": 1},
        pressure_reference_cell=0,
        boundaries={name: dict(bc) for name, bc in d.boundaries.items()},
    )
    doc["boundaries"]["top"] = {"type": "symmetry"}
    d = FluidDomain.from_document(doc, mesh=mesh)
    r = solve_quad_flow(d, velocity_scale_m_s=1, length_scale_m=1)
    y = mesh.cell_centres_m[:, 1]
    assert np.max(np.abs(r.velocity_m_s[:, 0] - y * (1 - y))) < 0.002
    assert np.max(np.abs(r.velocity_m_s[:, 1])) < 1e-9
    assert abs(r.power_balance_W) < 2e-7


@pytest.mark.parametrize(
    "options",
    [
        {"velocity_scale_m_s": 0},
        {"length_scale_m": False},
        {"velocity_relaxation": 1.1},
        {"pressure_relaxation": 0},
        {"max_iterations": True},
        {"max_simple_steps": 0},
        {"convection_scheme": "unknown"},
        {"model": "SA"},
        {"pressure_reference_Pa": 1},
        {"force_N_m3": np.ones(4)},
        {"force_N_m3": [float("nan"), 0]},
        {"viscosity_Pa_s": 1},
        {"viscosity_Pa_s": -1, "boundary_viscosity_Pa_s": 1},
    ],
)
def test_invalid_quad_flow_inputs(quad_mesh, options):
    d = domain_for(quad_mesh())
    args = dict(velocity_scale_m_s=1, length_scale_m=1)
    args.update(options)
    with pytest.raises(ValidationError):
        solve_quad_flow(d, **args)


def test_final_outlet_backflow_is_diagnosed(quad_mesh):
    mesh = quad_mesh(8, 8)
    force = np.column_stack((100 * (mesh.cell_centres_m[:, 1] - 0.5), np.zeros(64)))
    with pytest.raises(ValidationError, match="backflow"):
        solve_quad_flow(
            domain_for(mesh), velocity_scale_m_s=1, length_scale_m=1, force_N_m3=force
        )


def test_quad_stokes_channel_matches_independent_mac(quad_mesh):
    from fem_toolbox.cfd import (
        CartesianFluidGrid,
        ChannelBoundary,
        IncompressibleFluid,
        solve_stationary_flow,
    )

    n = 16
    mesh = quad_mesh(2 * n, n, size_m=(2, 1))
    quad = solve_quad_flow(
        domain_for(mesh, rho=1e-4), velocity_scale_m_s=1, length_scale_m=2
    )
    y = (np.arange(n) + 0.5) / n
    mac = solve_stationary_flow(
        CartesianFluidGrid((2 * n, n), (2, 1)),
        IncompressibleFluid(1e-4, 1),
        ChannelBoundary(y * (1 - y)),
    )
    cell_mac = np.column_stack(
        (
            (mac.u_m_s[:, :-1] + mac.u_m_s[:, 1:]).ravel() / 2,
            (mac.v_m_s[:-1] + mac.v_m_s[1:]).ravel() / 2,
        )
    )
    assert np.max(np.abs(quad.velocity_m_s - cell_mac)) < 0.003
    gp = np.diff(quad.pressure_Pa.reshape(n, 2 * n).mean(axis=0)) * n
    gm = np.diff(mac.pressure_Pa.mean(axis=0)) * n
    assert np.sqrt(np.mean((gp - gm) ** 2)) < 0.03


def test_closed_hydrostatic_patch_on_sheared_mesh(quad_mesh):
    mesh = quad_mesh(6, 6, shear=np.tan(np.pi / 9))
    d = domain_for(mesh, periodic=True)
    raw = dict(
        fluid_domain_schema_version=1,
        mesh="unused.json",
        width_m=0.2,
        fluid={"density_kg_m3": 1, "viscosity_Pa_s": 1},
        pressure_reference_cell=0,
        boundaries={name: {"type": "wall"} for name in d.boundaries},
    )
    r = solve_quad_flow(
        FluidDomain.from_document(raw, mesh=mesh),
        velocity_scale_m_s=1,
        length_scale_m=1,
        force_N_m3=[2, -3],
        velocity_relaxation=1,
        pressure_relaxation=1,
    )
    expected = mesh.cell_centres_m @ np.array([2.0, -3.0])
    expected -= expected[0]
    np.testing.assert_allclose(r.velocity_m_s, 0, atol=1e-10, rtol=0)
    np.testing.assert_allclose(r.pressure_Pa, expected, atol=1e-9, rtol=0)


def test_simple_step_matches_independent_coupled_lu(quad_mesh):
    from scipy.sparse import bmat
    from scipy.sparse.linalg import spsolve

    from fem_toolbox.cfd.quad_flow import _linear_simple
    from fem_toolbox.cfd.quad_system import QuadSystem

    d = domain_for(quad_mesh(4, 4))
    s = QuadSystem(d, np.ones(16), np.ones(16))
    diagonal = s.diffusion.matrix.diagonal().reshape(16, 2)
    pflux = s.pressure_flux(diagonal)
    rhs = -s.diffusion.constant
    u, p, count = _linear_simple(
        s,
        s.diffusion,
        rhs,
        diagonal,
        pflux,
        max_steps=500,
        momentum_gate=1e-8,
        continuity_gate=1e-8,
    )
    matrix = bmat(
        [
            [s.diffusion.matrix, s.pressure_momentum.matrix],
            [s.mass_velocity.matrix, s.balance @ pflux.matrix],
        ],
        format="csr",
    )
    exact = spsolve(
        matrix,
        np.concatenate(
            (rhs, -s.mass_velocity.constant - np.asarray(s.balance @ pflux.constant))
        ),
    )
    np.testing.assert_allclose(np.concatenate((u, p)), exact, atol=1e-10, rtol=0)
    assert count > 0


def test_complete_nonconvergence_history_is_attached(quad_mesh):
    with pytest.raises(ValidationError) as error:
        solve_quad_flow(
            domain_for(quad_mesh()),
            velocity_scale_m_s=1,
            length_scale_m=1,
            max_iterations=2,
        )
    assert len(error.value.residual_history) == 2
    assert len(error.value.residual_history[0]) == 5


@pytest.mark.parametrize("cubic", [False, True])
def test_polynomial_wall_normal_derivative_and_pressure_trace(quad_mesh, cubic):
    from fem_toolbox.cfd.quad_system import QuadSystem

    mesh = quad_mesh(4, 8)
    system = QuadSystem(domain_for(mesh, periodic=True), np.ones(32), np.ones(24))
    y = mesh.cell_centres_m[:, 1]
    u = np.column_stack((y * (1 - y) * (1 + y if cubic else 1), np.zeros(32)))
    gradient = system.guf.apply(u.ravel()).reshape(len(mesh.faces), 2, 2)
    walls = np.isin(system.kind, ["wall"])
    face_y = mesh.face_centres_m[walls, 1]
    expected = 1 - 3 * face_y**2 if cubic else 1 - 2 * face_y
    cell_gradient = system.gu.apply(u.ravel()).reshape(32, 2, 2)
    np.testing.assert_allclose(
        cell_gradient[:, 0, 1],
        1 - 3 * y**2 if cubic else 1 - 2 * y,
        atol=1e-11,
        rtol=0,
    )
    all_y = mesh.face_centres_m[:, 1]
    np.testing.assert_allclose(
        system.uf.apply(u.ravel()).reshape(-1, 2)[:, 0],
        all_y * (1 - all_y) * (1 + all_y if cubic else 1),
        atol=1e-11,
        rtol=0,
    )
    np.testing.assert_allclose(gradient[walls, 0, 1], expected, atol=1e-11, rtol=0)
    np.testing.assert_allclose(
        system.pf.apply(y ** (3 if cubic else 2))[walls],
        face_y ** (3 if cubic else 2),
        atol=1e-11,
        rtol=0,
    )


def test_cubic_diffusion_has_no_artificial_boundary_source(quad_mesh):
    from fem_toolbox.cfd.quad_system import QuadSystem

    mesh = quad_mesh(4, 8)
    system = QuadSystem(domain_for(mesh, periodic=True), np.ones(32), np.ones(24))
    y = mesh.cell_centres_m[:, 1]
    u = np.column_stack((y - y**3, np.zeros(32)))
    expected = np.column_stack((6 * y, np.zeros(32)))
    np.testing.assert_allclose(
        system.diffusion.apply(u.ravel()).reshape(32, 2)
        / system.fv.cell_volumes_m3[:, None],
        expected,
        atol=1e-10,
        rtol=0,
    )
