"""Model and direct solver entry points represent the same stationary problems."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from fem_toolbox.solvers.thermopiezoelectric import solve_thermopiezoelectric

ROOT = Path(__file__).parents[2]


@pytest.mark.parametrize(
    ("script", "artifact"),
    [
        ("elastic_tension.py", "elastic_tension.vtu"),
        ("steady_heat.py", "steady_heat.vtu"),
        ("parallel_plate.py", "parallel_plate.vtu"),
    ],
)
def test_v03_direct_solver_examples_still_run(
    tmp_path: Path, script: str, artifact: str
) -> None:
    completed = subprocess.run(
        [sys.executable, str(ROOT / "examples" / script)],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT / "src")},
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert (tmp_path / artifact).is_file()


def _mesh() -> Mesh:
    return Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        volume_sets={"body": np.array([0])},
    )


def _document(physics: str) -> dict:
    heat = {"conductivity": 2.0, "unit": "W/(m K)"}
    elastic = {"young": 10.0, "poisson": 0.0, "unit": "Pa"}
    dielectric = {"permittivity": 4.0, "unit": "F/m"}
    piezo = {
        "stiffness": (np.eye(6) * 10).tolist(),
        "stress_coupling": np.array([[2, 0, 0, 0, 0, 0], [0] * 6, [0] * 6]).tolist(),
        "permittivity": (np.eye(3) * 4).tolist(),
        "units": {
            "stiffness": "Pa",
            "stress_coupling": "C/m^2",
            "permittivity": "F/m",
        },
    }
    expansion = {
        "alpha": [0.01, 0, 0, 0, 0, 0],
        "alpha_unit": "1/K",
        "reference_temperature": 300.0,
        "reference_temperature_unit": "K",
    }
    materials = {
        "heat": {"T": {"m": heat}},
        "elasticity": {"u": {"m": elastic}},
        "electrostatics": {"phi": {"m": dielectric}},
        "piezoelectric": {"piezo": {"m": piezo}},
        "thermopiezoelectric": {
            "T": {"m": heat},
            "piezo": {"m": piezo},
            "expansion": {"m": expansion},
        },
    }[physics]
    temperature = [300.0, 310.0, 300.0, 300.0]
    displacement = [0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0]
    potential = [0, 1, 0, 0]
    conditions = []
    if physics in {"heat", "thermopiezoelectric"}:
        conditions.append(
            {"field": "T", "nodes": [0, 1, 2, 3], "value": temperature, "unit": "K"}
        )
    if physics in {"elasticity", "piezoelectric", "thermopiezoelectric"}:
        conditions.append(
            {"field": "u", "nodes": [0, 1, 2, 3], "value": displacement, "unit": "m"}
        )
    if physics in {"electrostatics", "piezoelectric", "thermopiezoelectric"}:
        conditions.append(
            {"field": "phi", "nodes": [0, 1, 2, 3], "value": potential, "unit": "V"}
        )
    return {
        "schema_version": 1,
        "units": "SI",
        "physics": physics,
        "materials": materials,
        "regions": {field: {"body": "m"} for field in materials},
        "boundary": {"dirichlet": conditions},
        **(
            {"solver": {"electric_thermal_strain": True}}
            if physics == "thermopiezoelectric"
            else {}
        ),
    }


def test_model_and_direct_single_field_solvers_agree() -> None:
    mesh = _mesh()
    nodes = [0, 1, 2, 3]
    temperature = np.array([300.0, 310.0, 300.0, 300.0])
    displacement = np.array([0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0])
    potential = np.array([0, 1, 0, 0])

    heat = Model(mesh=mesh, document=_document("heat"))
    direct_heat = solve_heat(
        mesh,
        heat.assignment.cell_materials,
        [Dirichlet(thermal_dofs(mesh).indices("T", nodes), temperature)],
    )
    model_heat = heat.solve().field_result
    assert_allclose(model_heat.temperature, direct_heat.temperature, atol=1e-12, rtol=0)
    assert_allclose(model_heat.heat_flux, direct_heat.heat_flux, atol=1e-12, rtol=0)
    assert_allclose(model_heat.heat_flux[0], [-20, 0, 0], atol=1e-12, rtol=0)
    assert_allclose(model_heat.reactions, direct_heat.reactions, atol=1e-12, rtol=0)
    assert_allclose(
        model_heat.power_balance, direct_heat.power_balance, atol=1e-12, rtol=0
    )

    elastic = Model(mesh=mesh, document=_document("elasticity"))
    direct_elastic = solve_elasticity(
        mesh,
        elastic.assignment.cell_materials,
        [Dirichlet(mechanical_dofs(mesh).indices("u", nodes), displacement)],
    )
    model_elastic = elastic.solve().field_result
    assert_allclose(
        model_elastic.displacement, direct_elastic.displacement, atol=1e-12, rtol=0
    )
    assert_allclose(model_elastic.stress, direct_elastic.stress, atol=1e-12, rtol=0)
    assert_allclose(model_elastic.stress[0, 0], 0.1, atol=1e-12, rtol=0)
    assert_allclose(
        model_elastic.force_balance, direct_elastic.force_balance, atol=1e-12, rtol=0
    )
    assert_allclose(model_elastic.energy, direct_elastic.energy, atol=1e-12, rtol=0)

    electric = Model(mesh=mesh, document=_document("electrostatics"))
    direct_electric = solve_electrostatics(
        mesh,
        electric.assignment.cell_materials,
        [Dirichlet(electric_dofs(mesh).indices("phi", nodes), potential)],
    )
    model_electric = electric.solve().field_result
    assert_allclose(
        model_electric.potential, direct_electric.potential, atol=1e-12, rtol=0
    )
    assert_allclose(
        model_electric.electric_displacement,
        direct_electric.electric_displacement,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        model_electric.electric_displacement[0], [-4, 0, 0], atol=1e-12, rtol=0
    )
    assert_allclose(
        model_electric.charge_balance,
        direct_electric.charge_balance,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(model_electric.energy, direct_electric.energy, atol=1e-12, rtol=0)


def test_model_and_direct_coupled_solvers_agree() -> None:
    mesh = _mesh()
    nodes = [0, 1, 2, 3]
    displacement = np.array([0, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0])
    potential = np.array([0, 1, 0, 0])
    dofs = piezoelectric_dofs(mesh)
    conditions = [
        Dirichlet(dofs.indices("u", nodes), displacement),
        Dirichlet(dofs.indices("phi", nodes), potential),
    ]
    piezo = Model(mesh=mesh, document=_document("piezoelectric"))
    direct_piezo = solve_piezoelectric(
        mesh, piezo.assignment.cell_materials, conditions
    )
    model_piezo = piezo.solve().field_result
    assert_allclose(
        model_piezo.displacement, direct_piezo.displacement, atol=1e-12, rtol=0
    )
    assert_allclose(model_piezo.potential, direct_piezo.potential, atol=1e-12, rtol=0)
    assert_allclose(model_piezo.stress, direct_piezo.stress, atol=1e-12, rtol=0)
    assert_allclose(model_piezo.stress[0, 0], 2.1, atol=1e-12, rtol=0)
    assert_allclose(
        model_piezo.electric_displacement,
        direct_piezo.electric_displacement,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(model_piezo.electric_displacement[0, 0], -3.98, atol=1e-12, rtol=0)
    assert_allclose(
        model_piezo.charge_balance, direct_piezo.charge_balance, atol=1e-12, rtol=0
    )

    coupled = Model(mesh=mesh, document=_document("thermopiezoelectric"))
    thermal_conditions = [
        Dirichlet(
            thermal_dofs(mesh).indices("T", nodes),
            np.array([300.0, 310.0, 300.0, 300.0]),
        )
    ]
    direct_coupled = solve_thermopiezoelectric(
        mesh,
        coupled.assignments["T"].cell_materials,
        coupled.assignments["piezo"].cell_materials,
        coupled.assignments["expansion"].cell_materials,
        thermal_conditions,
        conditions,
        electric_thermal_strain=True,
    )
    model_coupled = coupled.solve().field_result
    assert_allclose(
        model_coupled.heat.temperature,
        direct_coupled.heat.temperature,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        model_coupled.heat.power_balance,
        direct_coupled.heat.power_balance,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        model_coupled.piezoelectric.stress,
        direct_coupled.piezoelectric.stress,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(model_coupled.piezoelectric.stress[0, 0], 1.85, atol=1e-12, rtol=0)
    assert_allclose(
        model_coupled.piezoelectric.electric_displacement,
        direct_coupled.piezoelectric.electric_displacement,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        model_coupled.piezoelectric.electric_displacement[0, 0],
        -4.03,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        model_coupled.piezoelectric.force_balance,
        direct_coupled.piezoelectric.force_balance,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        model_coupled.piezoelectric.charge_balance,
        direct_coupled.piezoelectric.charge_balance,
        atol=1e-12,
        rtol=0,
    )
