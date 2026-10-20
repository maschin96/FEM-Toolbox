"""Predefined SA wall-channel mesh, wall and physical-input sensitivity gates."""

import numpy as np
import pytest

from fem_toolbox.cfd import solve_rans
from tests.reference.rans_channel_mesh import channel_domain
from tests.reference.sa_channel_1d import solve_channel_reference

pytestmark = pytest.mark.cfd_convergence


def run_case(n, *, beta=4.5, viscosity_factor=1.0, relaxation=0.1):
    nu = viscosity_factor / 587.19
    domain = channel_domain(n, beta=beta, kinematic_viscosity_m2_s=nu)
    reference = solve_channel_reference(kinematic_viscosity_m2_s=nu)
    result = solve_rans(
        domain,
        initial_nu_tilde_m2_s=3 * nu,
        velocity_scale_m_s=20,
        length_scale_m=1,
        sa_relaxation=relaxation,
    )
    expected, _ = reference.evaluate(domain.mesh.cell_centres_m[:, 1])
    volume = domain.mesh.cell_areas_m2 * domain.width_m
    velocity = result.flow.velocity_m_s[:, 0]
    rms = float(
        np.sqrt(
            np.sum(volume * (velocity - expected) ** 2) / np.sum(volume * expected**2)
        )
    )
    bulk = float(velocity @ volume / volume.sum())
    walls = domain.wall_quality.boundary_indices
    faces = domain.mesh.boundary_face_indices[walls]
    force = (
        result.flow.wall_traction_Pa
        * (domain.mesh.face_lengths_m[faces] * domain.width_m)[:, None]
    ).sum(axis=0)
    relative_wall_force = float(abs(force[0] + 2) / 2)
    assert rms < 0.005
    assert relative_wall_force < 0.001
    assert result.wall_y_plus.max() <= 1
    assert np.all(result.nu_tilde_m2_s >= 0)
    assert np.max(abs(result.sa_residual_m2_s2)) <= 1e-8 * 20 * nu
    assert np.max(abs(result.flow.momentum_residual_N_m3)) <= 4e-6
    assert abs(result.flow.global_mass_balance_kg_s) <= 2e-9
    assert abs(result.flow.power_balance_W) <= 0.008
    return (
        domain,
        reference,
        result,
        {
            "relative_rms_u": rms,
            "bulk_velocity_m_s": bulk,
            "relative_wall_force": relative_wall_force,
            "maximum_y_plus": float(result.wall_y_plus.max()),
            "iterations": result.iterations,
        },
    )


def test_rans_channel_four_wall_meshes(record_property):
    bulks = []
    for n in (32, 64, 128, 256):
        print(f"RANS channel: starting {n} cells/half", flush=True)
        _, _, _, metrics = run_case(n)
        print(f"RANS channel {n}: {metrics}", flush=True)
        for name, value in metrics.items():
            record_property(f"n{n}_{name}", value)
        bulks.append(metrics["bulk_velocity_m_s"])
    assert abs(bulks[-1] / bulks[-2] - 1) < 0.01


@pytest.mark.parametrize("viscosity_factor", [0.95, 1.05])
def test_rans_channel_physical_viscosity_sensitivity(viscosity_factor, record_property):
    _, reference, _, metrics = run_case(64, viscosity_factor=viscosity_factor)
    assert abs(metrics["bulk_velocity_m_s"] / reference.bulk_velocity_m_s - 1) < 0.005
    for name, value in metrics.items():
        record_property(name, value)
    print(f"RANS physical nu factor {viscosity_factor}: {metrics}", flush=True)


def test_rans_channel_wall_distribution_sensitivity(record_property):
    _, _, _, baseline = run_case(128)
    _, _, _, changed = run_case(128, beta=5.0)
    relative = abs(changed["bulk_velocity_m_s"] / baseline["bulk_velocity_m_s"] - 1)
    record_property("relative_bulk_change", relative)
    assert relative < 0.01
    print(f"RANS wall distribution: beta4.5={baseline}; beta5={changed}", flush=True)


def test_rans_channel_relaxation_repeat(record_property):
    _, _, base, _ = run_case(64)
    _, _, repeat, _ = run_case(64, relaxation=0.08)
    delta = float(
        np.linalg.norm(repeat.flow.velocity_m_s - base.flow.velocity_m_s)
        / np.linalg.norm(base.flow.velocity_m_s)
    )
    assert delta < 1e-6
    record_property("relative_velocity_repeat_difference", delta)


def test_rans_small_si_channel_for_fsi_has_the_same_reference_branch(record_property):
    # Fixed before this scaling check: H=0.1 m, u_tau=0.01 m/s, same Re_tau.
    # Lref=H and Uref=20*u_tau; no change to the model or physical gates.
    nu = 0.0005 / 587.19
    domain = channel_domain(
        32, half_height_m=0.05, kinematic_viscosity_m2_s=nu, pressure_drive_N_m3=0.002
    )
    reference = solve_channel_reference(
        half_height_m=0.05, kinematic_viscosity_m2_s=nu, pressure_drive_N_m3=0.002
    )
    result = solve_rans(
        domain,
        initial_nu_tilde_m2_s=3 * nu,
        velocity_scale_m_s=0.2,
        length_scale_m=0.1,
        sa_relaxation=0.1,
    )
    expected, _ = reference.evaluate(domain.mesh.cell_centres_m[:, 1])
    weights = domain.mesh.cell_areas_m2
    relative = float(
        np.sqrt(
            np.sum(weights * (result.flow.velocity_m_s[:, 0] - expected) ** 2)
            / np.sum(weights * expected**2)
        )
    )
    assert relative < 0.005
    assert result.wall_y_plus.max() <= 1
    assert np.max(np.abs(result.sa_residual_m2_s2)) <= 1e-8 * 0.2 * nu / 0.1
    record_property("relative_rms_u", relative)
    record_property("maximum_y_plus", float(result.wall_y_plus.max()))
    print(
        f"RANS small SI channel: relative RMS={relative}, y+={result.wall_y_plus.max()}",
        flush=True,
    )
