"""Independent Tet10 scalar oscillator, model and modal export checks."""

import csv
import json

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.model import Model
from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes
from tests.reference.quadratic_mesh import elevate
from tests.unit.test_piezo_modal import oscillator


def quadratic_oscillator(short=False):
    linear, mat, _ = oscillator(short)
    from fem_toolbox.mesh import Mesh

    mesh = elevate(
        Mesh(
            points=linear.points,
            tetra4=linear.cells,
            triangle3=[[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
        )
    )
    fixed = np.setdiff1d(np.arange(40), [3] if short else [3, 31])
    return mesh, mat, [Dirichlet(fixed, np.zeros(len(fixed)))]


@pytest.mark.parametrize("short", [False, True])
def test_exact_corner_mode_frequency_and_mass(short):
    mesh, mat, bc = quadratic_oscillator(short)
    result = solve_piezoelectric_modes(mesh, mat, 6, bc, count=1)
    # N1=x(2x-1): int N1²=1/420; int (4x-1)²=1/10 on unit tetrahedron.
    mass = 6 / 420
    omega = np.sqrt((10 if short else 10 + 4 / 3) / 10 / mass)
    assert_allclose(result.angular_frequencies, [omega], rtol=3e-14, atol=0)
    assert_allclose(result.displacement[0, 1, 0], 1 / np.sqrt(mass), rtol=3e-14, atol=0)
    assert_allclose(
        result.potential[0, 1],
        0 if short else 2 / 3 / np.sqrt(mass),
        rtol=3e-14,
        atol=1e-14,
    )
    assert result.displacement.shape == (1, 10, 3)
    assert max(result.relative_residual_u) < 1e-13
    assert max(result.relative_residual_phi) < 1e-13


def modal_document(mesh, mat):
    from tests.unit.test_tet10_piezoelectric_integration import document

    base = document()
    base["materials"]["piezo"]["p"].update(
        stiffness=mat.stiffness.tolist(),
        stress_coupling=mat.stress_coupling.tolist(),
        permittivity=mat.permittivity.tolist(),
    )
    base["regions"]["piezo"] = {"all": "p"}
    base["boundary"] = {
        "dirichlet": [
            {
                "field": "u",
                "nodes": [i for i in range(10) if i != 1],
                "value": 0,
                "unit": "m",
            },
            {
                "field": "u",
                "nodes": [1],
                "components": ["y", "z"],
                "value": 0,
                "unit": "m",
            },
            {
                "field": "phi",
                "nodes": [i for i in range(10) if i != 1],
                "value": 0,
                "unit": "V",
            },
        ]
    }
    return {
        "schema_version": 3,
        "units": "SI",
        "physics": "piezo_modes",
        "base": base,
        "density": {"value": [6], "unit": "kg/m^3"},
        "analysis": {"count": 1},
    }


def test_model_modal_export_centroid_and_quadratic_probe(tmp_path):
    from fem_toolbox.mesh import Mesh

    original, mat, _ = quadratic_oscillator()
    mesh = Mesh(
        points=original.points,
        tetra10=original.cells,
        triangle6=original.facets,
        volume_sets={"all": [0]},
    )
    model = Model(mesh=mesh, document=modal_document(mesh, mat))
    result = model.solve()
    modes = result.field_result.solution
    output = export_piezo_series(
        tmp_path / "modes", result, probes={"p": [0.2, 0.1, 0.1]}
    )
    report = json.loads((output / "manifest.json").read_text())
    assert report["schema_version"] == 2
    assert report["cell_fields"]["location"] == "cell_centroid"
    data = meshio.read(output / "frame_000000.vtu")
    assert_allclose(
        report["cell_fields"]["evaluation_points_m"],
        [[0.25, 0.25, 0.25]],
        rtol=0,
        atol=0,
    )
    assert_allclose(data.cell_data_dict["tetra_valid"]["tetra10"], 1)
    assert_allclose(data.cell_data_dict["tetra_valid"]["triangle6"], 0)
    assert_allclose(data.cell_data_dict["stress_at_centroid_real"]["triangle6"], 0)
    # At centroid x=1/4: derivative of N1 vanishes. Probe N1(.2)=-.12.
    assert_allclose(
        data.cell_data_dict["stress_at_centroid_real"]["tetra10"], 0, rtol=0, atol=2e-12
    )
    row = next(csv.DictReader((output / "probes.csv").open()))
    assert_allclose(
        float(row["u_x_real"]),
        -0.12 * modes.displacement[0, 1, 0],
        rtol=3e-14,
        atol=1e-14,
    )
    assert_allclose(
        float(row["phi_real"]), -0.12 * modes.potential[0, 1], rtol=3e-14, atol=1e-14
    )
    assert_allclose(
        Model(mesh=mesh, document=model.to_document())
        .solve()
        .field_result.solution.frequencies_hz,
        modes.frequencies_hz,
        rtol=0,
        atol=0,
    )


def test_two_modes_mass_orthogonality_and_orientation():
    from fem_toolbox.mesh import Mesh

    mesh, mat, _ = quadratic_oscillator()
    for reverse in (False, True):
        changed = Mesh(
            points=np.vstack([mesh.points, [2, 2, 2]]),
            tetra10=mesh.cells[:, [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]]
            if reverse
            else mesh.cells,
        )
        fixed = np.setdiff1d(np.arange(40), [3, 7])
        result = solve_piezoelectric_modes(
            changed, mat, [6], [Dirichlet(fixed, np.zeros(len(fixed)))], count=2
        )
        assert_allclose(
            result.angular_frequencies**2,
            42 * np.array([10, 20]) / 6,
            rtol=4e-14,
            atol=0,
        )
        values = result.displacement.reshape(2, 30)
        assert_allclose(6 / 420 * values @ values.T, np.eye(2), rtol=0, atol=4e-14)
        assert result.displacement.shape == (2, 10, 3)
        with pytest.raises(ValueError):
            result.displacement[0, 0, 0] = 1


def test_modal_validation_and_dynamic_scope(cube_mesh):
    from fem_toolbox.mesh import Mesh
    from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
    from fem_toolbox.solvers.piezo_transient import solve_transient_piezoelectric

    mesh, mat, bc = quadratic_oscillator()
    for density in (0, -1, np.nan, [1, 2]):
        with pytest.raises(ValidationError):
            solve_piezoelectric_modes(mesh, mat, density, bc, count=1)
    with pytest.raises(ValidationError, match="homogeneous"):
        solve_piezoelectric_modes(mesh, mat, 6, [*bc, Dirichlet([31], [1])], count=1)
    with pytest.raises(ValidationError, match="rigid"):
        solve_piezoelectric_modes(mesh, mat, 6, [], count=1)
    with pytest.raises(ValidationError, match="gauge"):
        solve_piezoelectric_modes(
            mesh,
            mat,
            6,
            [Dirichlet(np.setdiff1d(np.arange(30), [3]), np.zeros(29))],
            count=1,
        )
    with pytest.raises(ValidationError, match="2048"):
        solve_piezoelectric_modes(elevate(cube_mesh(4)), mat, 6, [], count=1)
    points = mesh.points.copy()
    points[4, 2] += 0.01
    with pytest.raises(ValidationError):
        solve_piezoelectric_modes(
            Mesh(points=points, tetra10=mesh.cells), mat, 6, bc, count=1
        )
    with pytest.raises(ValidationError, match="not yet supported"):
        solve_harmonic_piezoelectric(mesh, mat, 6, bc, 0)
    with pytest.raises(ValidationError, match="not yet supported"):
        solve_transient_piezoelectric(mesh, mat, 6, bc, [0, 1])


def test_gmsh_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_modes/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "modes")], check=True
    )
    assert (tmp_path / "modes/manifest.json").is_file()
