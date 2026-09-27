"""Thermal eigenstrain integration, centroid stress and full energy on Tet10."""

import json
from dataclasses import replace

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.solvers.heat import thermal_dofs
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics
from tests.reference.quadratic_mesh import elevate


def test_blocked_quadratic_temperature_integrates_energy(cube_mesh):
    mesh = elevate(cube_mesh(1))
    material = ElasticMaterial.isotropic(100, 0.25)
    result = solve_elasticity(
        mesh,
        material,
        [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
        temperature=300 + mesh.points[:, 0] ** 2,
        expansion=ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
    )
    x = result.field_points[:, 0]
    assert_allclose(result.strain, 0, rtol=0, atol=0)
    assert_allclose(
        result.stress,
        -(x[:, None] ** 2) * np.array([1.2, 0.4, 0.4, 0, 0, 0]),
        rtol=0,
        atol=2e-13,
    )
    # 1/2 C11 alpha_x² integral x^4 = .5*120*.01²/5 J.
    assert result.energy == pytest.approx(0.0012, rel=0, abs=2e-15)
    assert_allclose(result.force_balance, 0, rtol=0, atol=2e-13)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-13)


def test_linear_temperature_free_expansion_is_quadratic(cube_mesh):
    mesh = elevate(cube_mesh(2))
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(100, 0.25),
        [Dirichlet(mechanical_dofs(mesh).indices("u", left), np.zeros(3 * len(left)))],
        temperature=300 + mesh.points[:, 0],
        expansion=ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
    )
    exact = np.zeros_like(mesh.points)
    exact[:, 0] = 0.005 * mesh.points[:, 0] ** 2
    assert_allclose(result.displacement, exact, rtol=0, atol=2e-14)
    assert_allclose(result.stress, 0, rtol=0, atol=2e-11)
    assert result.energy < 1e-24


def test_coupled_quadratic_fields_and_bundle(cube_mesh, tmp_path):
    mesh = elevate(cube_mesh(2))
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    result = solve_thermomechanics(
        mesh,
        ThermalMaterial.isotropic(2),
        ElasticMaterial.isotropic(100, 0),
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 1),
        [Dirichlet(thermal_dofs(mesh).indices("T", left), np.full(len(left), 1.0))],
        [Dirichlet(mechanical_dofs(mesh).indices("u", left), np.zeros(3 * len(left)))],
        source=-4,
        outward_fluxes={int(f): -4 for f in mesh.facet_sets["right"]},
        body_force=lambda p: [-2 + 2 * p[0], 0, 0],
        tractions={int(f): [1, 0, 0] for f in mesh.facet_sets["right"]},
    )
    assert_allclose(
        result.heat.temperature, 1 + mesh.points[:, 0] ** 2, rtol=0, atol=4e-12
    )
    exact = np.zeros_like(mesh.points)
    exact[:, 0] = 0.01 * mesh.points[:, 0] ** 2
    assert_allclose(result.elasticity.displacement, exact, rtol=0, atol=1e-13)
    x = result.elasticity.field_points[:, 0]
    stress = np.zeros((mesh.n_cells, 6))
    stress[:, 0] = 2 * x - x**2
    assert_allclose(result.elasticity.stress, stress, rtol=0, atol=2e-11)
    assert result.elasticity.energy == pytest.approx(0.008 / 3, rel=0, abs=2e-14)
    assert abs(result.heat.power_balance) < 2e-11
    assert_allclose(result.elasticity.force_balance, 0, rtol=0, atol=2e-12)
    assert_allclose(result.elasticity.moment_balance, 0, rtol=0, atol=2e-12)
    vtu, manifest = export_result_bundle(tmp_path / "coupled", mesh, result)
    data = json.loads(manifest.read_text())
    assert data["physics"] == "thermomechanics"
    assert data["bundle_schema_version"] == 2
    assert "heat_flux_W_m2_at_centroid" in data["fields"]["tetra10"]
    assert data["fields"]["tetra10"]["mechanical_energy_J"]["kind"] == "cell_integral"
    assert_allclose(
        meshio.read(vtu).cell_data_dict["stress_Pa_at_centroid"]["tetra10"],
        result.elasticity.stress,
    )
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(
            tmp_path / "bad",
            mesh,
            replace(
                result,
                elasticity=replace(
                    result.elasticity, field_points=result.elasticity.field_points + 1
                ),
            ),
        )

    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(
            tmp_path / "bad_heat",
            mesh,
            replace(result, heat=replace(result.heat, field_points=None)),
        )


def test_region_expansion_reference_temperatures_and_rotation(cube_mesh):
    mesh = elevate(cube_mesh(2))
    centers = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    material = ElasticMaterial.isotropic(100, 0.25)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    expansions = [
        ThermalExpansion.from_principal(
            [0.001, 0, 0], 300 if x < 0.5 else 310, rotation
        )
        for x in centers[:, 0]
    ]
    result = solve_elasticity(
        mesh,
        material,
        [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
        temperature=np.full(mesh.n_points, 325),
        expansion=expansions,
    )
    theta = np.where(centers[:, 0] < 0.5, 25, 15)
    expected = -theta[:, None] * np.array([0.04, 0.12, 0.04, 0, 0, 0])
    assert_allclose(result.stress, expected, rtol=0, atol=2e-13)
    assert result.energy == pytest.approx(
        0.5 * 120 * 1e-6 * (25**2 + 15**2) / 2, rel=0, abs=2e-15
    )
    with pytest.raises(ValidationError, match="one ThermalExpansion"):
        solve_elasticity(
            mesh,
            material,
            [],
            temperature=np.full(mesh.n_points, 325),
            expansion=expansions[:-1],
        )
    with pytest.raises(ValidationError, match="nonnegative"):
        solve_elasticity(
            mesh,
            material,
            [],
            temperature=np.full(mesh.n_points, -1),
            expansion=expansions,
        )


def test_thermal_reference_shift_preserves_quadratic_energy(cube_mesh):
    mesh = elevate(cube_mesh(1))
    results = []
    for reference in [1, 301]:
        results.append(
            solve_elasticity(
                mesh,
                ElasticMaterial.isotropic(100, 0.25),
                [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
                temperature=reference + mesh.points[:, 0] ** 2,
                expansion=ThermalExpansion([0.01, 0, 0, 0, 0, 0], reference),
            )
        )
    assert_allclose(
        results[0].element_energy, results[1].element_energy, rtol=0, atol=2e-16
    )
    assert_allclose(results[0].stress, results[1].stress, rtol=0, atol=2e-13)


def test_tet10_thermomechanics_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_thermomechanics/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "example")], check=True
    )
    data = json.loads((tmp_path / "example" / "coupled.json").read_text())
    assert data["physics"] == "thermomechanics"
    assert data["bundle_schema_version"] == 2


def test_sequential_backends_are_independent_and_recorded(cube_mesh):
    from fem_toolbox.solvers.linear_backend import ScipyBackend

    class RecordingBackend(ScipyBackend):
        def __init__(self, name):
            self.name = name
            self.fields = []

        def factorize(self, operator, **kwargs):
            self.fields.append(tuple(kwargs["fields"]))
            return super().factorize(operator, **kwargs)

    mesh = elevate(cube_mesh(1))
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    thermal, mechanical = RecordingBackend("thermal"), RecordingBackend("mechanical")
    result = solve_thermomechanics(
        mesh,
        ThermalMaterial.isotropic(2),
        ElasticMaterial.isotropic(100, 0),
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
        [Dirichlet(thermal_dofs(mesh).indices("T", left), np.full(len(left), 301))],
        [Dirichlet(mechanical_dofs(mesh).indices("u", left), np.zeros(3 * len(left)))],
        thermal_backend=thermal,
        mechanical_backend=mechanical,
    )
    assert thermal.fields == [("T",)]
    assert mechanical.fields == [("u",)]
    assert result.heat.provenance.solver_parameters["backend"] == "thermal"
    assert result.elasticity.provenance.solver_parameters["backend"] == "mechanical"
    expected = np.zeros_like(mesh.points)
    expected[:, 0] = 0.01 * mesh.points[:, 0]
    assert_allclose(result.elasticity.displacement, expected, rtol=0, atol=1e-13)
