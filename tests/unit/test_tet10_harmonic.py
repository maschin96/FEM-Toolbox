"""Complex Tet10 response from independent scalar integrals and P2 export."""

import csv
import json

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.errors import ValidationError
from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
from tests.unit.test_tet10_modes import modal_document, quadratic_oscillator


@pytest.mark.parametrize("hz", [0.0, 0.2, 2.0])
def test_complex_corner_force_charge_and_dc(hz):
    mesh, mat, bc = quadratic_oscillator()
    force = np.zeros(40, dtype=complex)
    force[3] = 1 + 2j
    force[31] = 0.01 - 0.02j
    result = solve_harmonic_piezoelectric(mesh, mat, 6, bc, hz, loads=force)
    omega = 2 * np.pi * hz
    # int (dN1/dx)^2=1/10; int N1^2=1/420.
    q = (force[3] + 2 / 3 * force[31]) / ((10 + 4 / 3) / 10 - 6 / 420 * omega**2)
    phi = (0.2 * q - force[31]) / 0.3
    assert_allclose(result.displacement[1, 0], q, rtol=3e-14, atol=0)
    assert_allclose(result.potential[1], phi, rtol=3e-14, atol=0)
    assert_allclose(
        result.velocity, 1j * omega * result.displacement, rtol=3e-14, atol=1e-14
    )
    assert_allclose(
        result.acceleration, -(omega**2) * result.displacement, rtol=3e-14, atol=1e-14
    )
    assert result.relative_residual_u < 1e-12
    assert result.relative_residual_phi < 1e-12


def test_complex_voltage_moving_support_and_resonance():
    mesh, mat, bc = quadratic_oscillator(short=True)
    g = np.zeros(40, dtype=complex)
    g[31] = 2 - 3j
    g[6] = 1 + 0.5j
    omega = 1.7
    result = solve_harmonic_piezoelectric(
        mesh, mat, 6, bc, omega / (2 * np.pi), prescribed=g
    )
    # Corner cross-mass int N1*N2=1/2520; normal stiffness cross-term is zero.
    q = (-0.2 * g[31] + omega**2 * (6 / 2520) * g[6]) / (1 - omega**2 * (6 / 420))
    assert_allclose(result.displacement[1, 0], q, rtol=3e-14, atol=0)
    assert_allclose(result.potential[1], g[31], rtol=3e-14, atol=0)
    with pytest.raises(ValidationError, match="resonan"):
        solve_harmonic_piezoelectric(mesh, mat, 6, bc, np.sqrt(70) / (2 * np.pi))


def harmonic_document(mesh, mat):
    document = modal_document(mesh, mat)
    document["physics"] = "harmonic_piezoelectric"
    document["base"]["loads"] = {"body_force": {"value": [6, 0, 0], "unit": "N/m^3"}}
    document["analysis"] = {
        "frequencies": {"values": [0.0, 0.2, 2.0], "unit": "Hz"},
        "load_factor": {"real": 1, "imag": 2},
    }
    return document


def test_model_complex_export_and_quadratic_probes(tmp_path):
    original, mat, _ = quadratic_oscillator()
    mesh = Mesh(
        points=original.points,
        tetra10=original.cells,
        triangle6=original.facets,
        volume_sets={"all": [0]},
    )
    model = Model(mesh=mesh, document=harmonic_document(mesh, mat))
    solved = model.solve()
    states = solved.field_result.solution
    for state in states:
        q = (
            -0.05
            * (1 + 2j)
            / ((10 + 4 / 3) / 10 - 6 / 420 * state.angular_frequency**2)
        )
        assert_allclose(state.displacement[1, 0], q, rtol=3e-14, atol=0)
    output = export_piezo_series(
        tmp_path / "sweep", solved, probes={"p": [0.2, 0.1, 0.1]}
    )
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["cell_fields"]["location"] == "cell_centroid"
    data = meshio.read(output / "frame_000002.vtu")
    assert_allclose(
        data.cell_data_dict["stress_at_centroid_imag"]["tetra10"], 0, rtol=0, atol=2e-12
    )
    assert_allclose(
        data.point_data["u_imag"][1],
        states[-1].displacement[1].imag,
        rtol=3e-14,
        atol=1e-14,
    )
    with (output / "probes.csv").open() as stream:
        row = list(csv.DictReader(stream))[-1]
    value = -0.12 * states[-1].displacement[1, 0]
    for suffix, expected in [
        ("real", value.real),
        ("imag", value.imag),
        ("amplitude", abs(value)),
        ("phase_rad", np.angle(value)),
    ]:
        assert_allclose(float(row["u_x_" + suffix]), expected, rtol=3e-14, atol=1e-14)
    repeated = (
        Model(mesh=mesh, document=model.to_document()).solve().field_result.solution
    )
    assert_allclose(repeated[-1].displacement, states[-1].displacement, rtol=0, atol=0)


def test_complex_nonzero_centroid_fields(tmp_path):
    original, mat, _ = quadratic_oscillator()
    mesh = Mesh(
        points=original.points,
        tetra10=original.cells,
        triangle6=original.facets,
        volume_sets={"all": [0]},
    )
    document = harmonic_document(mesh, mat)
    document["base"]["loads"] = {}
    document["analysis"] = {
        "frequencies": {"values": [0.2], "unit": "Hz"},
        "prescribed": [
            {"dofs": [12], "unit": "m", "real": 1, "imag": 2},
            {"dofs": [34], "unit": "V", "real": 0.3, "imag": -0.4},
        ],
    }
    solved = Model(mesh=mesh, document=document).solve()
    output = export_piezo_series(
        tmp_path / "complex_fields", solved, probes={"p": [0.2, 0.1, 0.1]}
    )
    data = meshio.read(output / "frame_000000.vtu")
    # Edge 0--1 has grad N4=(0,-1,-1) at the centroid. Free N1 has zero gradient.
    g = 1 + 2j
    p = 0.3 - 0.4j
    expected = {
        "strain": np.array([0, 0, 0, 0, -g, -g]),
        "stress": np.array([0, 0, 0, 0, -5 * g, -6 * g]),
        "electric_field": np.array([0, p, p]),
        "electric_displacement": np.array([0, 4 * p, 5 * p]),
    }
    for name, value in expected.items():
        for suffix, part in [
            ("real", value.real),
            ("imag", value.imag),
            ("amplitude", np.abs(value)),
        ]:
            assert_allclose(
                data.cell_data_dict[name + "_at_centroid_" + suffix]["tetra10"][0],
                part,
                rtol=3e-14,
                atol=2e-12,
            )
    assert_allclose(data.point_data["u_phase_valid"][0], 0, rtol=0, atol=0)
    with (output / "probes.csv").open() as stream:
        row = next(csv.DictReader(stream))
    state = solved.field_result.solution[0]
    expected_phi = -0.12 * state.potential[1] + 0.48 * p
    assert_allclose(
        float(row["phi_real"]) + 1j * float(row["phi_imag"]),
        expected_phi,
        rtol=3e-14,
        atol=1e-14,
    )


def test_complex_charged_bar_dc(cube_mesh):
    from tests.unit.test_tet10_newmark import charged_bar

    mesh, mat, bc, electrodes, force, initial = charged_bar(cube_mesh)
    factor = 1 + 2j
    result = solve_harmonic_piezoelectric(
        mesh, mat, 7500, bc, 0, electrodes=electrodes, loads=force * factor
    )
    assert_allclose(result.displacement, initial * factor, rtol=0, atol=2e-16)
    assert_allclose(result.potential, mesh.points[:, 2] * factor, rtol=0, atol=5e-11)
    top = np.flatnonzero(mesh.points[:, 2] == 0.02)
    assert abs(result.residual_phi[top].sum()) < 1e-18


def test_time_domain_converges_to_complex_phasor():
    from fem_toolbox.solvers.piezo_transient import solve_transient_piezoelectric

    mesh, mat, bc = quadratic_oscillator()
    omega = 1.3
    force = np.zeros(40, dtype=complex)
    force[3] = 0.01 + 0.02j
    harmonic = solve_harmonic_piezoelectric(
        mesh, mat, 6, bc, omega / (2 * np.pi), loads=force
    )
    errors = []
    for n in (100, 200):
        result = solve_transient_piezoelectric(
            mesh,
            mat,
            6,
            bc,
            np.linspace(0, 1, n + 1),
            loads=lambda t: np.real(force * np.exp(1j * omega * t)),
            initial_displacement=harmonic.displacement.real,
            initial_velocity=harmonic.velocity.real,
        )
        exact = np.real(harmonic.displacement[1, 0] * np.exp(1j * omega * result.times))
        errors.append(np.max(np.abs(result.displacement[:, 1, 0] - exact)))
    assert errors[-1] < 1e-6
    assert 3.8 < errors[0] / errors[1] < 4.2


@pytest.mark.parametrize("frequency", [-1, np.nan, 1j, True, 1e200])
def test_invalid_frequency(frequency):
    mesh, mat, bc = quadratic_oscillator()
    with pytest.raises(ValidationError):
        solve_harmonic_piezoelectric(mesh, mat, 6, bc, frequency)


def test_gmsh_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_harmonic/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "sweep")], check=True
    )
    assert (tmp_path / "sweep/manifest.json").is_file()
