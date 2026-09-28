"""Independent polynomial thermal-piezo patches and sequential transfer."""

import json
from dataclasses import replace

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.materials import ThermalExpansion, ThermalMaterial
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from fem_toolbox.solvers.thermopiezoelectric import solve_thermopiezoelectric
from tests.reference.quadratic_mesh import elevate
from tests.unit.test_tet10_piezoelectric_integration import material


@pytest.mark.parametrize("electric", [False, True])
def test_blocked_quadratic_temperature_and_integrated_enthalpy(cube_mesh, electric):
    mesh = elevate(cube_mesh(1))
    dofs = piezoelectric_dofs(mesh)
    nodes = np.arange(mesh.n_points)
    x = mesh.points[:, 0]
    result = solve_piezoelectric(
        mesh,
        material(),
        [
            Dirichlet(dofs.indices("u", nodes), np.zeros(3 * len(nodes))),
            Dirichlet(dofs.indices("phi", nodes), x * x),
        ],
        temperature=300 + x * x,
        expansion=ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
        electric_thermal_strain=electric,
    )
    centers = result.field_points[:, 0]
    stress = np.zeros_like(result.stress)
    stress[:, 0] = 4 * centers - centers**2
    displacement = np.zeros_like(result.electric_displacement)
    displacement[:, 0] = -6 * centers - electric * 0.02 * centers**2
    assert_allclose(result.stress, stress, rtol=0, atol=2e-12)
    assert_allclose(result.electric_displacement, displacement, rtol=0, atol=2e-12)
    assert_allclose(
        [
            result.mechanical_energy.sum(),
            result.coupling_energy.sum(),
            result.dielectric_enthalpy.sum(),
        ],
        [0.001, -0.01 * electric, -2],
        rtol=0,
        atol=2e-13,
    )


def coupled_case(mesh, electric=True, thermal_backend=None, piezoelectric_backend=None):
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    dofs = piezoelectric_dofs(mesh)
    return solve_thermopiezoelectric(
        mesh,
        ThermalMaterial.isotropic(2),
        material(),
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 1),
        [Dirichlet(left, np.ones(len(left)))],
        [
            Dirichlet(dofs.indices("u", left), np.zeros(3 * len(left))),
            Dirichlet(dofs.indices("phi", left), np.zeros(len(left))),
        ],
        source=-4,
        outward_fluxes={int(f): -4 for f in mesh.facet_sets["right"]},
        body_force=lambda p: [-6 + 2 * p[0], 0, 0],
        volume_charge=lambda p: -5.96 - electric * 0.04 * p[0],
        tractions={int(f): [5, 0, 0] for f in mesh.facet_sets["right"]},
        surface_charges={
            int(f): -5.96 - electric * 0.02 for f in mesh.facet_sets["right"]
        },
        electric_thermal_strain=electric,
        thermal_backend=thermal_backend,
        piezoelectric_backend=piezoelectric_backend,
    )


@pytest.mark.parametrize("electric", [False, True])
def test_sequential_quadratic_patch_and_export(cube_mesh, tmp_path, electric):
    mesh = elevate(cube_mesh(2))
    result = coupled_case(mesh, electric)
    x = mesh.points[:, 0]
    assert_allclose(result.heat.temperature, 1 + x * x, rtol=0, atol=4e-12)
    expected = np.zeros_like(mesh.points)
    expected[:, 0] = 0.01 * x * x
    piezo = result.piezoelectric
    assert_allclose(piezo.displacement, expected, rtol=0, atol=2e-13)
    assert_allclose(piezo.potential, x * x, rtol=0, atol=2e-12)
    center = piezo.field_points[:, 0]
    stress = np.zeros_like(piezo.stress)
    stress[:, 0] = 6 * center - center**2
    assert_allclose(piezo.stress, stress, rtol=0, atol=2e-11)
    # eps_m=.02x-.01x²; coupling strain subtracts .01x² only with switch.
    assert_allclose(
        [
            piezo.mechanical_energy.sum(),
            piezo.coupling_energy.sum(),
            piezo.dielectric_enthalpy.sum(),
        ],
        [0.008 / 3, 0.08 / 3 - 0.01 * electric, -2],
        rtol=0,
        atol=2e-12,
    )
    assert_allclose(piezo.force_balance, 0, rtol=0, atol=2e-11)
    assert_allclose(piezo.moment_balance, 0, rtol=0, atol=2e-11)
    assert abs(piezo.charge_balance) < 2e-11
    assert abs(result.heat.power_balance) < 2e-11
    _, report = export_result_bundle(tmp_path / "coupled", mesh, result)
    data = json.loads(report.read_text())
    assert data["physics"] == "thermopiezoelectric"
    assert data["bundle_schema_version"] == 2
    assert data["fields"]["tetra10"]["coupling_energy_J"]["kind"] == "cell_integral"
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(
            tmp_path / "bad",
            mesh,
            replace(result, heat=replace(result.heat, field_points=None)),
        )


def document(mesh, electric=True):
    from tests.unit.test_tet10_piezoelectric_integration import (
        document as piezo_document,
    )

    raw = piezo_document()
    raw["physics"] = "thermopiezoelectric"
    raw["materials"]["T"] = {"k": {"conductivity": 2, "unit": "W/(m K)"}}
    raw["materials"]["expansion"] = {
        "a": {
            "alpha": [0.01, 0, 0, 0, 0, 0],
            "alpha_unit": "1/K",
            "reference_temperature": 1,
            "reference_temperature_unit": "K",
        }
    }
    raw["regions"]["T"] = {"body": "k"}
    raw["regions"]["expansion"] = {"body": "a"}
    raw["boundary"]["dirichlet"] = [
        {"field": "u", "nodes": list(range(mesh.n_points)), "value": 0, "unit": "m"},
        {"field": "T", "set": "left", "value": 1, "unit": "K"},
        {"field": "T", "set": "right", "value": 2, "unit": "K"},
    ]
    raw["boundary"]["electrodes"][1]["charge"] = 0
    raw["solver"] = {"electric_thermal_strain": electric}
    return raw


@pytest.mark.parametrize("electric", [False, True])
def test_model_floating_electrode_thermal_voltage(cube_mesh, tmp_path, electric):
    from fem_toolbox.model import Model

    mesh = elevate(cube_mesh(2))
    model = Model(mesh=mesh, document=document(mesh, electric))
    solved = model.solve()
    result = solved.field_result
    assert_allclose(result.heat.temperature, 1 + mesh.points[:, 0], rtol=0, atol=4e-12)
    assert_allclose(
        result.piezoelectric.potential,
        -electric * 0.01 / 3 * mesh.points[:, 0] ** 2,
        rtol=0,
        atol=2e-13,
    )
    assert_allclose(result.piezoelectric.electric_displacement, 0, rtol=0, atol=2e-12)
    assert abs(result.piezoelectric.electrode_states[1].charge) < 2e-12
    replay = Model(mesh=mesh, document=model.to_document()).solve()
    assert_allclose(
        replay.field_result.piezoelectric.potential,
        result.piezoelectric.potential,
        rtol=0,
        atol=0,
    )
    export_result_bundle(tmp_path / "model", mesh, solved)


def test_gmsh_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_thermopiezoelectric/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "example")], check=True
    )
    assert (
        json.loads((tmp_path / "example/coupled.json").read_text())["physics"]
        == "thermopiezoelectric"
    )


def test_region_expansion_rotation_reference_shift_and_invalid_inputs(cube_mesh):
    mesh = elevate(cube_mesh(1))
    dofs = piezoelectric_dofs(mesh)
    nodes = np.arange(mesh.n_points)
    bc = [
        Dirichlet(dofs.indices("u", nodes), np.zeros(3 * len(nodes))),
        Dirichlet(dofs.indices("phi", nodes), np.zeros(len(nodes))),
    ]
    angle = 0.3
    rotation = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1],
        ]
    )
    alpha_tensor = rotation @ np.diag([0.01, 0.02, 0.03]) @ rotation.T
    alpha = np.array(
        [
            alpha_tensor[0, 0],
            alpha_tensor[1, 1],
            alpha_tensor[2, 2],
            2 * alpha_tensor[1, 2],
            2 * alpha_tensor[0, 2],
            2 * alpha_tensor[0, 1],
        ]
    )
    expansions = [ThermalExpansion(alpha, 300 + i) for i in range(mesh.n_cells)]
    result = solve_piezoelectric(
        mesh,
        material(),
        bc,
        temperature=np.full(mesh.n_points, 310),
        expansion=expansions,
        electric_thermal_strain=True,
    )
    eps = np.array([alpha * (10 - i) for i in range(mesh.n_cells)])
    assert_allclose(result.stress, -eps @ material().stiffness.T, rtol=0, atol=2e-12)
    assert_allclose(
        result.electric_displacement,
        -eps @ material().stress_coupling.T,
        rtol=0,
        atol=2e-12,
    )
    shifted = solve_piezoelectric(
        mesh,
        material(),
        bc,
        temperature=np.full(mesh.n_points, 410),
        expansion=[ThermalExpansion(alpha, 400 + i) for i in range(mesh.n_cells)],
        electric_thermal_strain=True,
    )
    assert_allclose(
        shifted.mechanical_energy, result.mechanical_energy, rtol=0, atol=2e-13
    )
    assert shifted.provenance.input_sha256 != result.provenance.input_sha256
    with pytest.raises(ValidationError, match="together"):
        solve_piezoelectric(mesh, material(), bc, temperature=np.ones(mesh.n_points))
    with pytest.raises(ValidationError, match="nonnegative"):
        solve_piezoelectric(
            mesh,
            material(),
            bc,
            temperature=-np.ones(mesh.n_points),
            expansion=expansions,
        )
    with pytest.raises(ValidationError, match="per cell"):
        solve_piezoelectric(
            mesh,
            material(),
            bc,
            temperature=np.ones(mesh.n_points),
            expansion=expansions[:1],
        )


def test_separate_backend_routing(cube_mesh):
    from fem_toolbox.solvers.linear_backend import ScipyBackend

    class RecordingBackend(ScipyBackend):
        def __init__(self, name):
            self.name = name
            self.fields = []

        def factorize(self, operator, **kwargs):
            self.fields.append(tuple(kwargs["fields"]))
            return super().factorize(operator, **kwargs)

    thermal, piezo = RecordingBackend("thermal"), RecordingBackend("piezo")
    result = coupled_case(elevate(cube_mesh(1)), True, thermal, piezo)
    assert thermal.fields == [("T",)]
    assert piezo.fields == [("u", "phi")]
    assert result.heat.provenance.solver_parameters["backend"] == "thermal"
    assert result.piezoelectric.provenance.solver_parameters["backend"] == "piezo"


def test_free_expansion_has_zero_mechanical_strain_and_electric_field(cube_mesh):
    mesh = elevate(cube_mesh(2))
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    dofs = piezoelectric_dofs(mesh)
    x = mesh.points[:, 0]
    result = solve_piezoelectric(
        mesh,
        material(),
        [
            Dirichlet(dofs.indices("u", left), np.zeros(3 * len(left))),
            Dirichlet(dofs.indices("phi", left), np.zeros(len(left))),
        ],
        temperature=300 + x,
        expansion=ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
        electric_thermal_strain=True,
    )
    expected = np.zeros_like(mesh.points)
    expected[:, 0] = 0.005 * x * x
    assert_allclose(result.displacement, expected, rtol=0, atol=2e-13)
    assert_allclose(result.potential, 0, rtol=0, atol=2e-13)
    assert_allclose(result.stress, 0, rtol=0, atol=2e-11)
    assert result.mechanical_energy.sum() < 1e-24
