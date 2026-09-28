"""Tet10 Newmark against exact scalar dynamics and rigid-translation inertia."""

import csv
import json

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.solvers.piezo_transient import (
    KinematicValues,
    solve_transient_piezoelectric,
)
from tests.unit.test_tet10_modes import modal_document, quadratic_oscillator


@pytest.mark.parametrize("short", [False, True])
def test_free_oscillator_energy_and_time_order(short):
    mesh, mat, bc = quadratic_oscillator(short)
    mass = 6 / 420
    stiffness = (10 if short else 10 + 4 / 3) / 10
    omega = np.sqrt(stiffness / mass)
    initial = np.zeros((10, 3))
    initial[1, 0] = 0.01
    velocity = np.zeros((10, 3))
    velocity[1, 0] = 0.02
    errors = []
    for n in (40, 80, 160):
        result = solve_transient_piezoelectric(
            mesh,
            mat,
            6,
            bc,
            np.linspace(0, 0.4, n + 1),
            initial_displacement=initial,
            initial_velocity=velocity,
        )
        errors.append(
            abs(
                result.displacement[-1, 1, 0]
                - (0.01 * np.cos(0.4 * omega) + 0.02 / omega * np.sin(0.4 * omega))
            )
        )
        assert_allclose(
            result.total_energy,
            0.5 * stiffness * 0.01**2 + 0.5 * mass * 0.02**2,
            rtol=3e-12,
            atol=0,
        )
        assert_allclose(
            result.potential[:, 1],
            0 if short else 2 / 3 * result.displacement[:, 1, 0],
            rtol=3e-12,
            atol=1e-14,
        )
        assert_allclose(
            result.acceleration[:, 1, 0],
            -(omega**2) * result.displacement[:, 1, 0],
            rtol=3e-12,
            atol=1e-13,
        )
        assert max(result.relative_residual_u) < 1e-12
        assert max(result.relative_residual_phi) < 1e-12
    assert 3.8 < errors[0] / errors[1] < 4.2
    assert 3.8 < errors[1] / errors[2] < 4.2


def test_moving_rigid_translation_with_independent_p2_body_load():
    mesh, mat, bc = quadratic_oscillator()
    fixed = [i for i in range(10) if i != 1]

    def prescribed(t):
        u = np.zeros((10, 3))
        v = u.copy()
        a = u.copy()
        u[fixed, 0] = 0.003 * t * t
        v[fixed, 0] = 0.006 * t
        a[fixed, 0] = 0.006
        return KinematicValues(u, v, a, np.zeros(10))

    force = np.zeros(40)
    # int N_corner=-V/20, int N_edge=V/5; V=1/6, rho=6.
    force[np.arange(10) * 3] = np.array([-0.05] * 4 + [0.2] * 6) * 0.006
    result = solve_transient_piezoelectric(
        mesh,
        mat,
        6,
        bc,
        np.linspace(0, 1, 41),
        loads=lambda t: force,
        prescribed=prescribed,
    )
    assert_allclose(
        result.displacement[:, 1, 0], 0.003 * result.times**2, rtol=0, atol=3e-15
    )
    assert_allclose(result.velocity[:, 1, 0], 0.006 * result.times, rtol=0, atol=3e-15)
    assert_allclose(result.acceleration[:, 1, 0], 0.006, rtol=0, atol=3e-14)
    assert_allclose(result.potential, 0, rtol=0, atol=3e-15)
    assert_allclose(result.residual_u, 0, rtol=0, atol=3e-14)


def transient_document(mesh, mat):
    document = modal_document(mesh, mat)
    document["physics"] = "transient_piezoelectric"
    document["base"]["loads"] = {"body_force": {"value": [6, 0, 0], "unit": "N/m^3"}}
    document["analysis"] = {
        "times": {"values": [0, 0.1, 0.2], "unit": "s"},
        "load_factor": {"kind": "step", "time": 0.05, "before": 0, "after": 1},
    }
    return document


def test_model_events_body_force_export_and_roundtrip(tmp_path):
    original, mat, _ = quadratic_oscillator()
    mesh = Mesh(
        points=original.points,
        tetra10=original.cells,
        triangle6=original.facets,
        volume_sets={"all": [0]},
    )
    model = Model(mesh=mesh, document=transient_document(mesh, mat))
    solved = model.solve()
    result = solved.field_result.solution
    assert_allclose(result.times, [0, 0.05, 0.1, 0.2], rtol=0, atol=0)
    assert_allclose(result.displacement[1], 0, rtol=0, atol=1e-15)
    assert_allclose(result.velocity[1], 0, rtol=0, atol=1e-15)
    assert_allclose(result.acceleration[1, 1, 0], -0.05 / (6 / 420), rtol=3e-14, atol=0)
    output = export_piezo_series(
        tmp_path / "series", solved, probes={"p": [0.2, 0.1, 0.1]}
    )
    report = json.loads((output / "manifest.json").read_text())
    assert report["schema_version"] == 2
    data = meshio.read(output / "frame_000003.vtu")
    assert_allclose(
        data.cell_data_dict["stress_at_centroid_real"]["tetra10"], 0, rtol=0, atol=2e-12
    )
    with (output / "probes.csv").open() as stream:
        row = list(csv.DictReader(stream))[-1]
    for name, values in [
        ("u", result.displacement),
        ("v", result.velocity),
        ("a", result.acceleration),
    ]:
        assert_allclose(
            float(row[name + "_x_real"]),
            -0.12 * values[-1, 1, 0],
            rtol=3e-14,
            atol=1e-14,
        )
    assert_allclose(
        Model(mesh=mesh, document=model.to_document())
        .solve()
        .field_result.solution.displacement,
        result.displacement,
        rtol=0,
        atol=0,
    )


@pytest.mark.parametrize("short", [False, True])
def test_electrode_mode_free_vibration(cube_mesh, short):
    from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes
    from tests.convergence.test_tet10_modal import bar

    mesh, mat, bc, electrodes = bar(cube_mesh, 2, short)
    mode = solve_piezoelectric_modes(
        mesh, mat, 7500, bc, electrodes=electrodes, count=1
    )
    initial = mode.displacement[0] * 1e-6
    dt = 0.04 / mode.angular_frequencies[0]
    result = solve_transient_piezoelectric(
        mesh,
        mat,
        7500,
        bc,
        dt * np.arange(21),
        electrodes=electrodes,
        initial_displacement=initial,
    )
    # Exact discrete Newmark phase on one eigenmode; tests mass/electric reconstruction.
    phase = np.arange(21) * 2 * np.arctan(mode.angular_frequencies[0] * dt / 2)
    assert_allclose(
        result.displacement,
        np.cos(phase)[:, None, None] * initial,
        rtol=2e-9,
        atol=1e-14,
    )
    assert_allclose(
        result.total_energy,
        0.5 * (mode.angular_frequencies[0] * 1e-6) ** 2,
        rtol=2e-9,
        atol=0,
    )
    top = np.flatnonzero(mesh.points[:, 2] == 0.02)
    assert_allclose(
        result.potential[:, top],
        np.broadcast_to(result.potential[:, top[:1]], (21, len(top))),
        rtol=0,
        atol=1e-12,
    )
    if short:
        assert_allclose(result.potential[:, top], 0, rtol=0, atol=1e-12)
    else:
        assert_allclose(result.residual_phi[:, top].sum(axis=1), 0, rtol=0, atol=1e-18)


def test_voltage_jump_keeps_mechanics_continuous():
    mesh, mat, bc = quadratic_oscillator(short=True)

    def prescribed(t):
        zero = np.zeros((10, 3))
        phi = np.zeros(10)
        phi[1] = float(t >= 0.05)
        return KinematicValues(zero, zero, zero, phi)

    result = solve_transient_piezoelectric(
        mesh, mat, 6, bc, [0, 0.1], prescribed=prescribed, events=[0.05]
    )
    assert_allclose(result.times, [0, 0.05, 0.1], rtol=0, atol=0)
    assert_allclose(result.displacement[1], 0, rtol=0, atol=1e-14)
    assert_allclose(result.velocity[1], 0, rtol=0, atol=1e-14)
    assert_allclose(result.potential[:, 1], [0, 1, 1], rtol=0, atol=1e-14)
    assert_allclose(
        result.acceleration[1, 1, 0], -(2 / 10) / (6 / 420), rtol=3e-14, atol=0
    )


def charged_bar(cube_mesh):
    from tests.convergence.test_tet10_modal import bar

    mesh, mat, bc, electrodes = bar(cube_mesh, 2, False)
    strain, gradient = 1e-4, 1.0
    stress = 1e9 * strain + 5 * gradient
    # Conductor charge is -integral D.n with n outward from the solid.
    charge = -(5 * strain - 1e-8 * gradient) * 0.01**2
    top = np.flatnonzero(mesh.points[:, 2] == 0.02)
    force = np.zeros(4 * mesh.n_points)
    for face in mesh.facets[np.all(mesh.points[mesh.facets, 2] == 0.02, axis=1)]:
        p = mesh.points[face[:3]]
        area = np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])) / 2
        # Triangle6 constant traction: corner integrals vanish; edge integrals A/3.
        force[3 * face[3:] + 2] += stress * area / 3
    force[3 * mesh.n_points + top[0]] = -charge
    initial = np.zeros((mesh.n_points, 3))
    initial[:, 2] = strain * mesh.points[:, 2]
    return mesh, mat, bc, electrodes, force, initial


def test_charged_electrode_affine_equilibrium(cube_mesh):
    mesh, mat, bc, electrodes, force, initial = charged_bar(cube_mesh)
    result = solve_transient_piezoelectric(
        mesh,
        mat,
        7500,
        bc,
        [0, 1e-6, 3e-6],
        electrodes=electrodes,
        loads=lambda t: force,
        initial_displacement=initial,
    )
    assert_allclose(
        result.displacement,
        np.broadcast_to(initial, result.displacement.shape),
        rtol=0,
        atol=2e-16,
    )
    assert_allclose(
        result.potential,
        np.broadcast_to(mesh.points[:, 2], result.potential.shape),
        rtol=0,
        atol=2e-11,
    )
    assert_allclose(result.velocity, 0, rtol=0, atol=2e-10)
    top = np.flatnonzero(mesh.points[:, 2] == 0.02)
    assert_allclose(result.residual_phi[:, top].sum(axis=1), 0, rtol=0, atol=1e-18)


def test_invalid_initial_support_and_history_limit():
    from fem_toolbox.errors import ValidationError

    mesh, mat, bc = quadratic_oscillator()
    with pytest.raises(ValidationError, match="initial"):
        solve_transient_piezoelectric(
            mesh, mat, 6, bc, [0, 1], initial_displacement=np.ones((10, 3))
        )
    with pytest.raises(ValidationError, match="2000000"):
        solve_transient_piezoelectric(mesh, mat, 6, bc, np.linspace(0, 1, 50001))

    def prescribed(t):
        zero = np.zeros((10, 3))
        u = zero.copy()
        u[0, 0] = float(t >= 0.5)
        return KinematicValues(u, zero, zero, np.zeros(10))

    with pytest.raises(ValidationError, match="discontinuous"):
        solve_transient_piezoelectric(
            mesh, mat, 6, bc, [0, 1], prescribed=prescribed, events=[0.5]
        )


def test_gmsh_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_newmark/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "history")], check=True
    )
    assert (tmp_path / "history/manifest.json").is_file()
