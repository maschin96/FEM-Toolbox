"""Model-level constraint checks before any field assembly or LU solve."""

import numpy as np
import pytest

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.errors import (
    ElectricGaugeError,
    ModelSchemaError,
    ThermalNullspaceError,
    ValidationError,
)
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalMaterial,
)
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.solvers.elasticity import solve_elasticity
from fem_toolbox.solvers.electrostatics import solve_electrostatics
from fem_toolbox.solvers.heat import solve_heat
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric


def _two_components() -> Mesh:
    tetra = np.array([[0, 1, 2, 3], [4, 5, 6, 7]])
    one = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    return Mesh(
        points=np.vstack([one, one + np.array([3, 0, 0])]),
        tetra4=tetra,
        volume_sets={"body": np.array([0, 1])},
    )


def _heat_document() -> dict:
    return {
        "schema_version": 1,
        "units": "SI",
        "physics": "heat",
        "materials": {"T": {"k": {"conductivity": 2.0, "unit": "W/(m K)"}}},
        "regions": {"T": {"body": "k"}},
        "boundary": {
            "dirichlet": [{"field": "T", "nodes": [0], "value": 300.0, "unit": "K"}]
        },
    }


def test_missing_thermal_reference_is_reported_per_component_at_construction() -> None:
    with pytest.raises(ModelSchemaError, match=r"T.*component 1.*nodes"):
        Model(mesh=_two_components(), document=_heat_document())
    with pytest.raises(ThermalNullspaceError, match="component 1"):
        solve_heat(
            _two_components(), ThermalMaterial.isotropic(2), [Dirichlet([0], [300])]
        )


def test_equal_duplicate_dirichlet_is_allowed_but_conflicting_value_is_rejected() -> (
    None
):
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        volume_sets={"body": np.array([0])},
    )
    document = _heat_document()
    document["boundary"]["dirichlet"].append(
        {"field": "T", "nodes": [0], "value": 300.0, "unit": "K"}
    )
    Model(mesh=mesh, document=document).solve()
    document["boundary"]["dirichlet"][1]["value"] = 301.0
    with pytest.raises(ModelSchemaError, match=r"boundary\.dirichlet\[1\].*DOF 0"):
        Model(mesh=mesh, document=document)


def test_missing_electric_gauge_matches_direct_solver_component() -> None:
    mesh = _two_components()
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "electrostatics",
        "materials": {"phi": {"eps": {"permittivity": 2.0, "unit": "F/m"}}},
        "regions": {"phi": {"body": "eps"}},
        "boundary": {
            "dirichlet": [{"field": "phi", "nodes": [0], "value": 0, "unit": "V"}]
        },
    }
    with pytest.raises(
        ModelSchemaError, match=r"phi.*component 1.*nodes \[4, 5, 6, 7\]"
    ):
        Model(mesh=mesh, document=document)
    with pytest.raises(ElectricGaugeError, match="component 1"):
        solve_electrostatics(
            mesh, DielectricMaterial.isotropic(2), [Dirichlet([0], [0])]
        )


def test_missing_mechanical_support_matches_direct_solver_component() -> None:
    mesh = _two_components()
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "elasticity",
        "materials": {"u": {"solid": {"young": 10.0, "poisson": 0.0, "unit": "Pa"}}},
        "regions": {"u": {"body": "solid"}},
        "boundary": {
            "dirichlet": [
                {"field": "u", "nodes": [0, 1, 2, 3], "value": 0, "unit": "m"}
            ]
        },
    }
    with pytest.raises(ModelSchemaError, match=r"u.*component 1.*rank 0/6"):
        Model(mesh=mesh, document=document)
    with pytest.raises(ValidationError, match=r"rigid modes in component 1"):
        solve_elasticity(
            mesh,
            ElasticMaterial.isotropic(10, 0),
            [Dirichlet(np.arange(12), np.zeros(12))],
        )


def _electrode_case() -> tuple[Mesh, dict]:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        triangle3=np.array([[1, 2, 3]]),
        volume_sets={"body": np.array([0])},
        facet_sets={"terminal": np.array([0])},
    )
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "piezoelectric",
        "materials": {
            "piezo": {
                "p": {
                    "stiffness": (np.eye(6) * 10).tolist(),
                    "stress_coupling": np.zeros((3, 6)).tolist(),
                    "permittivity": (np.eye(3) * 2).tolist(),
                    "units": {
                        "stiffness": "Pa",
                        "stress_coupling": "C/m^2",
                        "permittivity": "F/m",
                    },
                }
            }
        },
        "regions": {"piezo": {"body": "p"}},
        "boundary": {
            "dirichlet": [
                {"field": "u", "nodes": [0, 1, 2, 3], "value": 0, "unit": "m"},
                {"field": "phi", "nodes": [0], "value": 0, "unit": "V"},
            ],
            "electrodes": [
                {"name": "terminal", "set": "terminal", "charge": 0, "charge_unit": "C"}
            ],
        },
    }
    return mesh, document


def test_floating_electrode_is_valid_but_dirichlet_collision_is_not() -> None:
    mesh, document = _electrode_case()
    result = Model(mesh=mesh, document=document).solve().field_result
    assert abs(result.electrode_states[0].charge) < 1e-12
    document["boundary"]["dirichlet"].append(
        {"field": "phi", "nodes": [1], "value": 0, "unit": "V"}
    )
    with pytest.raises(ModelSchemaError, match=r"floating electrode.*Dirichlet"):
        Model(mesh=mesh, document=document)


def test_floating_charge_without_reference_does_not_fix_gauge() -> None:
    mesh, document = _electrode_case()
    document["boundary"]["dirichlet"].pop()
    with pytest.raises(
        ModelSchemaError, match=r"phi.*component 0.*missing potential reference"
    ):
        Model(mesh=mesh, document=document)


def test_fixed_electrode_potential_must_match_dirichlet_and_sets_cannot_overlap() -> (
    None
):
    mesh, document = _electrode_case()
    electrode = document["boundary"]["electrodes"][0]
    del electrode["charge"]
    del electrode["charge_unit"]
    electrode["potential"] = 1.0
    electrode["potential_unit"] = "V"
    document["boundary"]["dirichlet"].append(
        {"field": "phi", "nodes": [1], "value": 0, "unit": "V"}
    )
    with pytest.raises(ModelSchemaError, match=r"electrode.*conflicts"):
        Model(mesh=mesh, document=document)
    document["boundary"]["dirichlet"][-1]["value"] = 1.0
    Model(mesh=mesh, document=document).solve()
    document["boundary"]["electrodes"].append(
        {"name": "other", "set": "terminal", "potential": 1.0, "potential_unit": "V"}
    )
    with pytest.raises(ModelSchemaError, match=r"overlaps"):
        Model(mesh=mesh, document=document)


def test_positive_convection_anchors_heat_but_zero_coefficient_does_not() -> None:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        triangle3=np.array([[1, 2, 3]]),
        volume_sets={"body": np.array([0])},
        facet_sets={"cooling": np.array([0])},
    )
    document = _heat_document()
    document["boundary"] = {
        "convection": [
            {
                "set": "cooling",
                "coefficient": 1.0,
                "ambient_temperature": 300.0,
                "unit": "W/(m^2 K)",
                "ambient_unit": "K",
            }
        ]
    }
    result = Model(mesh=mesh, document=document).solve().field_result
    np.testing.assert_allclose(result.temperature, 300, atol=1e-12, rtol=0)
    document["boundary"]["convection"][0]["coefficient"] = 0.0
    with pytest.raises(ModelSchemaError, match=r"T.*component 0.*convection"):
        Model(mesh=mesh, document=document)


def test_overlapping_flux_assignments_report_facet_and_set() -> None:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
        triangle3=np.array([[1, 2, 3]]),
        volume_sets={"body": np.array([0])},
        facet_sets={"outer": np.array([0])},
    )
    document = _heat_document()
    document["boundary"]["flux"] = [
        {"set": "outer", "value": 1, "unit": "W/m^2"},
        {"set": "outer", "value": 2, "unit": "W/m^2"},
    ]
    with pytest.raises(ModelSchemaError, match=r"facet 0.*outer"):
        Model(mesh=mesh, document=document)


def test_surface_charge_on_electrode_is_rejected_by_model_and_solver() -> None:
    mesh, document = _electrode_case()
    document["boundary"]["surface_charge"] = [
        {"set": "terminal", "value": 1.0, "unit": "C/m^2"}
    ]
    with pytest.raises(ModelSchemaError, match=r"terminal.*surface charge"):
        Model(mesh=mesh, document=document)
    dofs = piezoelectric_dofs(mesh)
    material = PiezoelectricMaterial(np.eye(6) * 10, np.zeros((3, 6)), np.eye(3) * 2)
    with pytest.raises(ValidationError, match="surface charge"):
        solve_piezoelectric(
            mesh,
            material,
            [
                Dirichlet(dofs.indices("u", [0, 1, 2, 3]), np.zeros(12)),
                Dirichlet(dofs.indices("phi", [0]), np.zeros(1)),
            ],
            electrodes=[Electrode("terminal", "terminal", charge=0)],
            surface_charges={0: 1.0},
        )
