"""Production Tet10 heat solves, field locations, models and portable output."""

import json

import meshio
import numpy as np
import pytest

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.io.heat_series import export_heat_series
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.model import Model
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs
from fem_toolbox.solvers.transient_heat import (
    TemperatureBC,
    TransientConvection,
    solve_transient_heat,
)
from tests.reference.quadratic_mesh import elevate


def stationary_document():
    return {
        "schema_version": 1,
        "units": "SI",
        "physics": "heat",
        "materials": {"T": {"k": {"conductivity": 2.0, "unit": "W/(m K)"}}},
        "regions": {"T": {"body": "k"}},
        "loads": {"source": {"value": -4.0, "unit": "W/m^3"}},
        "boundary": {
            "dirichlet": [{"field": "T", "set": "left", "value": 300.0, "unit": "K"}],
            "flux": [{"set": "right", "value": -4.0, "unit": "W/m^2"}],
        },
    }


def test_stationary_quadratic_patch_and_centroid_export(cube_mesh, tmp_path):
    mesh = elevate(cube_mesh(2))
    xyz = mesh.points
    boundary = np.unique(mesh.facets)
    exact = 300 + xyz[:, 0] ** 2 + 0.5 * xyz[:, 1] ** 2
    result = solve_heat(
        mesh,
        ThermalMaterial.isotropic(2),
        [Dirichlet(thermal_dofs(mesh).indices("T", boundary), exact[boundary])],
        source=-6.0,
    )
    np.testing.assert_allclose(result.temperature, exact, rtol=0, atol=3e-12)
    centroids = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    gradient = np.column_stack(
        (2 * centroids[:, 0], centroids[:, 1], np.zeros(mesh.n_cells))
    )
    np.testing.assert_allclose(result.field_points, centroids, rtol=0, atol=1e-15)
    np.testing.assert_allclose(result.gradient, gradient, rtol=0, atol=3e-12)
    np.testing.assert_allclose(result.heat_flux, -2 * gradient, rtol=0, atol=6e-12)
    assert abs(result.power_balance) < 2e-11
    assert result.source_power == pytest.approx(-6, abs=2e-13)
    vtu, manifest = export_result_bundle(tmp_path / "quadratic", mesh, result)
    exported = meshio.read(vtu)
    np.testing.assert_allclose(
        exported.cell_data_dict["heat_flux_W_m2_at_centroid"]["tetra10"],
        result.heat_flux,
    )
    payload = json.loads(manifest.read_text())
    assert payload["bundle_schema_version"] == 2
    assert (
        payload["fields"]["tetra10"]["heat_flux_W_m2_at_centroid"]["location"]
        == "cell_centroid"
    )
    np.testing.assert_allclose(payload["field_evaluation_points_m"], centroids)


@pytest.mark.parametrize("robin", [False, True])
def test_triangle6_flux_robin_and_named_model_boundary(cube_mesh, robin):
    mesh = elevate(cube_mesh(1))
    raw = stationary_document()
    if robin:
        raw["boundary"].pop("flux")
        raw["boundary"]["convection"] = [
            {
                "set": "right",
                "coefficient": 2.0,
                "ambient_temperature": 303.0,
                "unit": "W/(m^2 K)",
                "ambient_unit": "K",
            }
        ]
    model = Model(mesh=mesh, document=raw)
    result = model.solve().field_result
    np.testing.assert_allclose(
        result.temperature, 300 + mesh.points[:, 0] ** 2, rtol=0, atol=4e-12
    )
    assert abs(result.power_balance) < 2e-11
    replay = Model(mesh=mesh, document=model.to_document()).solve().field_result
    np.testing.assert_array_equal(replay.temperature, result.temperature)


def test_transient_storage_probes_and_series(cube_mesh, tmp_path):
    mesh = elevate(cube_mesh(1))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    result = solve_transient_heat(
        mesh,
        material,
        [0, 0.1, 0.4, 1],
        300.0,
        source=12.0,
        probes={"inside": [0.2, 0.2, 0.2]},
    )
    np.testing.assert_allclose(result.frames[-1].temperature, 301, rtol=0, atol=4e-12)
    np.testing.assert_allclose(
        result.probe_temperature[:, 0], 300 + result.times, rtol=0, atol=4e-12
    )
    assert max(abs(b.balance) for b in result.balances) < 2e-10
    assert result.balances[-1].stored_energy_change == pytest.approx(12, abs=1e-10)
    directory = export_heat_series(tmp_path / "series", mesh, result)
    manifest = json.loads((directory / "manifest.json").read_text())
    assert manifest["schema_version"] == 2
    assert manifest["field_location"] == "cell_centroid"
    output = meshio.read(directory / manifest["frames"][-1]["file"])
    assert output.cells_dict["tetra10"].shape == mesh.cells.shape
    assert "heat_flux_W_m2_at_centroid" in output.cell_data_dict
    # Negative quadratic shape values inside a cell are legitimate, not outside.
    initial = 300 + mesh.points[:, 0] ** 2
    sampled = solve_transient_heat(
        mesh, material, [0, 0.01], initial, probes={"inside": [0.2, 0.2, 0.2]}
    )
    assert sampled.probe_temperature[0, 0] == pytest.approx(300.04, abs=1e-12)


def test_quadratic_lumping_requires_explicit_hrz(cube_mesh):
    mesh = elevate(cube_mesh(1))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    with pytest.raises(ValidationError, match="lump"):
        solve_transient_heat(mesh, material, [0, 1], 300.0, lumped=True)
    result = solve_transient_heat(
        mesh, material, [0, 1], 300.0, source=12.0, lumping="hrz"
    )
    # HRZ preserves total capacity, not the consistent nodal source moments.
    assert result.balances[-1].stored_energy_change == pytest.approx(12.0, abs=1e-10)
    unforced = solve_transient_heat(mesh, material, [0, 1], 300.0, lumping="hrz")
    np.testing.assert_allclose(
        unforced.frames[-1].temperature, 300.0, rtol=0, atol=4e-12
    )
    assert result.provenance.solver_parameters["capacity"] == "hrz"


def test_bundle_reader_supports_centroid_schema(cube_mesh, tmp_path):
    import runpy
    from pathlib import Path

    mesh = elevate(cube_mesh(1))
    result = Model(mesh=mesh, document=stationary_document()).solve()
    _, manifest = export_result_bundle(tmp_path / "result", mesh, result)
    reader = runpy.run_path(
        str(Path(__file__).parents[2] / "scripts/read_result_bundle.py")
    )
    assert reader["read_bundle"](manifest)["bundle_schema_version"] == 2


def test_transient_model_profiles_midpoints_and_preflight(cube_mesh):
    from tests.unit.test_transient_model import document

    mesh = elevate(cube_mesh(1))
    raw = document()
    result = Model(mesh=mesh, document=raw).solve().field_result
    np.testing.assert_allclose(result.frames[-1].temperature, 300.5, rtol=0, atol=3e-12)
    raw["solver"] = {"lumped": True}
    with pytest.raises(ValidationError, match="lump"):
        Model(mesh=mesh, document=raw)
    raw["solver"] = {"lumping": "hrz"}
    raw.pop("loads")
    raw["boundary"] = {
        "temperature": [
            {
                "set": "left",
                "value": {
                    "kind": "ramp",
                    "start": 0.0,
                    "end": 1.0,
                    "start_value": 300.0,
                    "end_value": 301.0,
                },
                "unit": "K",
            }
        ]
    }
    model = Model(mesh=mesh, document=raw)
    solved = model.solve().field_result
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    np.testing.assert_allclose(
        solved.frames[-1].temperature[left], 301.0, rtol=0, atol=3e-12
    )
    assert max(abs(b.balance) for b in solved.balances) < 1e-7
    no_reference = stationary_document()
    no_reference.pop("boundary")
    with pytest.raises(ValidationError, match="missing temperature reference"):
        Model(mesh=mesh, document=no_reference)


def test_quadratic_curvature_probe_and_export_rejections(cube_mesh, tmp_path):
    from dataclasses import replace

    from fem_toolbox.mesh import Mesh

    mesh = elevate(cube_mesh(1))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    bad_points = mesh.points.copy()
    bad_points[mesh.cells[0, 4], 2] += 0.01
    curved = Mesh(points=bad_points, tetra10=mesh.cells, triangle6=mesh.facets)
    with pytest.raises(ValidationError):
        solve_heat(curved, material, [])
    with pytest.raises(ValidationError, match="outside"):
        solve_transient_heat(
            mesh, material, [0, 1], 300.0, probes={"outside": [2, 2, 2]}
        )
    with pytest.raises(ValidationError, match="history budget"):
        solve_transient_heat(mesh, material, [0, 1], 300.0, max_history_bytes=1)
    result = Model(mesh=mesh, document=stationary_document()).solve().field_result
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(tmp_path / "bad", mesh, replace(result, field_points=None))
    for physics in (
        "thermopiezoelectric",
        "piezo_modes",
    ):
        with pytest.raises(ValidationError, match="not yet supported"):
            Model(mesh=mesh, document={"physics": physics})


@pytest.mark.parametrize("lumping", ["", "row_sum", False])
def test_invalid_quadratic_lumping_is_not_silently_defaulted(cube_mesh, lumping):
    mesh = elevate(cube_mesh(1))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    with pytest.raises(ValidationError, match="lumping"):
        solve_transient_heat(mesh, material, [0, 1], 300, lumping=lumping)


def test_quadratic_material_interface_preserves_normal_flux(cube_mesh):
    mesh = elevate(cube_mesh(2))
    centers = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    materials = [ThermalMaterial.isotropic(1 if x < 0.5 else 8) for x in centers[:, 0]]
    x = mesh.points[:, 0]
    # Unit-scale temperature-increment patch; 300 K offsets are covered above.
    exact = np.where(x <= 0.5, x, 0.5 + (x - 0.5) / 8)
    fixed = np.flatnonzero((x == 0) | (x == 1))
    result = solve_heat(mesh, materials, [Dirichlet(fixed, exact[fixed])])
    np.testing.assert_allclose(result.temperature, exact, rtol=0, atol=4e-12)
    np.testing.assert_allclose(
        result.heat_flux, np.tile([-1, 0, 0], (mesh.n_cells, 1)), rtol=0, atol=2e-11
    )
    assert abs(result.power_balance) < 2e-11


def test_tet10_toml_examples_use_real_gmsh_and_export(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_heat/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "example")], check=True
    )
    data = json.loads((tmp_path / "example" / "stationary.json").read_text())
    assert data["mesh"]["tetra10_cells"] == 6
    assert data["bundle_schema_version"] == 2
    series = json.loads(
        (tmp_path / "example" / "transient" / "manifest.json").read_text()
    )
    assert series["status"] == "complete"
    assert len(series["node_ids"]) == 27


@pytest.mark.parametrize("robin", [False, True])
def test_transient_quadratic_equilibrium_with_triangle6_boundaries(cube_mesh, robin):
    mesh = elevate(cube_mesh(1))
    exact = 300 + mesh.points[:, 0] ** 2
    right = mesh.facet_sets["right"]
    result = solve_transient_heat(
        mesh,
        ThermalMaterial.isotropic(2, density=3, specific_heat=4),
        [0, 0.1, 0.4],
        exact,
        source=-4,
        conditions=[
            TemperatureBC(np.unique(mesh.facets[mesh.facet_sets["left"]]), 300)
        ],
        outward_fluxes={} if robin else {int(face): -4 for face in right},
        convection={int(face): TransientConvection(2, 303) for face in right}
        if robin
        else {},
    )
    for frame in result.frames:
        np.testing.assert_allclose(frame.temperature, exact, rtol=0, atol=5e-12)
    for balance in result.balances:
        assert abs(balance.balance) < 2e-10
