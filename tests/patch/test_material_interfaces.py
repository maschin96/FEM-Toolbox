"""Independent two-layer 1D flux and traction references on a 3D Tet4 mesh."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
    assign_region_materials,
)
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric


def _layer_mesh(cube_mesh: object) -> Mesh:
    original = cube_mesh(2)
    cell_x = original.points[original.tetra4, 0].mean(axis=1)
    return Mesh(
        points=original.points,
        tetra4=original.tetra4,
        tetra4_ids=np.arange(1001, 1001 + original.n_tetra4),
        volume_sets={
            "left": np.flatnonzero(cell_x < 0.5),
            "right": np.flatnonzero(cell_x > 0.5),
        },
    )


def _piecewise(x: np.ndarray, coefficient: float) -> np.ndarray:
    return np.minimum(x, 0.5) / coefficient + np.maximum(x - 0.5, 0) / (2 * coefficient)


def _scalar_bc(
    mesh: Mesh, field: str, left: float, right: float, dofs: object
) -> list[Dirichlet]:
    left_nodes = np.flatnonzero(mesh.points[:, 0] == 0)
    right_nodes = np.flatnonzero(mesh.points[:, 0] == 1)
    return [
        Dirichlet(dofs.indices(field, left_nodes), np.full(len(left_nodes), left)),
        Dirichlet(dofs.indices(field, right_nodes), np.full(len(right_nodes), right)),
    ]


def _mechanical_bc(mesh: Mesh, dofs: object, extension: float) -> list[Dirichlet]:
    nodes = np.arange(mesh.n_points)
    ids = dofs.indices("u", nodes).reshape(-1, 3)
    return [
        Dirichlet(ids[:, 1:].ravel(), np.zeros(mesh.n_points * 2)),
        Dirichlet(
            ids[mesh.points[:, 0] == 0, 0], np.zeros(np.sum(mesh.points[:, 0] == 0))
        ),
        Dirichlet(
            ids[mesh.points[:, 0] == 1, 0],
            np.full(np.sum(mesh.points[:, 0] == 1), extension),
        ),
    ]


def test_heat_interface_piecewise_solution_and_flux(cube_mesh: object) -> None:
    mesh = _layer_mesh(cube_mesh)
    assignment = assign_region_materials(
        mesh,
        "T",
        {"left": "k2", "right": "k4"},
        {"k2": ThermalMaterial.isotropic(2), "k4": ThermalMaterial.isotropic(4)},
    )
    result = solve_heat(
        mesh,
        assignment.cell_materials,
        _scalar_bc(mesh, "T", 301, 300, thermal_dofs(mesh)),
    )
    q = 1 / (0.5 / 2 + 0.5 / 4)
    expected = 301 - q * _piecewise(mesh.points[:, 0], 2)
    assert_allclose(result.temperature, expected, atol=3e-11, rtol=2e-13)
    assert_allclose(
        result.heat_flux, np.tile([q, 0, 0], (mesh.n_tetra4, 1)), atol=2e-11, rtol=2e-12
    )
    assert abs(result.power_balance) < 1e-10
    assert result.relative_residual < 1e-12


def test_dielectric_interface_piecewise_solution_and_charge(cube_mesh: object) -> None:
    mesh = _layer_mesh(cube_mesh)
    assignment = assign_region_materials(
        mesh,
        "phi",
        {"left": "e2", "right": "e4"},
        {
            "e2": DielectricMaterial.isotropic(2e-9),
            "e4": DielectricMaterial.isotropic(4e-9),
        },
    )
    result = solve_electrostatics(
        mesh,
        assignment.cell_materials,
        _scalar_bc(mesh, "phi", 1, 0, electric_dofs(mesh)),
    )
    charge_density = 1 / (0.5 / 2e-9 + 0.5 / 4e-9)
    expected = 1 - charge_density * _piecewise(mesh.points[:, 0], 2e-9)
    assert_allclose(result.potential, expected, atol=3e-12, rtol=2e-12)
    assert_allclose(
        result.electric_displacement,
        np.tile([charge_density, 0, 0], (mesh.n_tetra4, 1)),
        atol=2e-20,
        rtol=2e-11,
    )
    assert abs(result.charge_balance) < 1e-19
    assert result.relative_residual < 1e-12


def test_elastic_interface_piecewise_extension_and_traction(cube_mesh: object) -> None:
    mesh = _layer_mesh(cube_mesh)
    c2 = np.diag([2e6] * 3 + [1e6] * 3)
    c4 = np.diag([4e6] * 3 + [2e6] * 3)
    assignment = assign_region_materials(
        mesh,
        "u",
        {"left": "c2", "right": "c4"},
        {"c2": ElasticMaterial(c2), "c4": ElasticMaterial(c4)},
    )
    extension = 1e-3
    result = solve_elasticity(
        mesh,
        assignment.cell_materials,
        _mechanical_bc(mesh, mechanical_dofs(mesh), extension),
    )
    traction = extension / (0.5 / 2e6 + 0.5 / 4e6)
    expected = traction * _piecewise(mesh.points[:, 0], 2e6)
    assert_allclose(result.displacement[:, 0], expected, atol=2e-14, rtol=2e-11)
    assert_allclose(result.displacement[:, 1:], 0, atol=2e-14, rtol=0)
    assert_allclose(result.stress[:, 0], traction, atol=2e-7, rtol=2e-11)
    assert_allclose(result.stress[:, 1:], 0, atol=2e-7, rtol=0)
    assert_allclose(result.force_balance, 0, atol=2e-8, rtol=0)
    assert result.relative_residual < 1e-12


def test_uncoupled_piezo_limit_preserves_both_interface_balances(
    cube_mesh: object,
) -> None:
    mesh = _layer_mesh(cube_mesh)
    piezo = assign_region_materials(
        mesh,
        "piezo",
        {"left": "p2", "right": "p4"},
        {
            "p2": PiezoelectricMaterial(
                np.diag([2e6] * 3 + [1e6] * 3), np.zeros((3, 6)), np.eye(3) * 2e-9
            ),
            "p4": PiezoelectricMaterial(
                np.diag([4e6] * 3 + [2e6] * 3), np.zeros((3, 6)), np.eye(3) * 4e-9
            ),
        },
    )
    expansion = assign_region_materials(
        mesh,
        "expansion",
        {"left": "a0", "right": "a0"},
        {"a0": ThermalExpansion(np.zeros(6), 293.15)},
    )
    assert (
        piezo.to_provenance_dict()["catalog"]["p2"]["axes_local_to_global"]
        == np.eye(3).tolist()
    )
    assert (
        expansion.to_provenance_dict()["catalog"]["a0"]["reference_temperature"]
        == 293.15
    )
    dofs = piezoelectric_dofs(mesh)
    extension = 1e-3
    result = solve_piezoelectric(
        mesh,
        piezo.cell_materials,
        _mechanical_bc(mesh, dofs, extension) + _scalar_bc(mesh, "phi", 1, 0, dofs),
        expansion=expansion.cell_materials,
        temperature=np.full(mesh.n_points, 293.15),
        electric_thermal_strain=True,
    )
    traction = extension / (0.5 / 2e6 + 0.5 / 4e6)
    charge_density = 1 / (0.5 / 2e-9 + 0.5 / 4e-9)
    assert_allclose(
        result.displacement[:, 0],
        traction * _piecewise(mesh.points[:, 0], 2e6),
        atol=2e-13,
        rtol=2e-10,
    )
    assert_allclose(
        result.potential,
        1 - charge_density * _piecewise(mesh.points[:, 0], 2e-9),
        atol=2e-11,
        rtol=2e-10,
    )
    assert_allclose(result.stress[:, 0], traction, atol=2e-6, rtol=2e-10)
    assert_allclose(
        result.electric_displacement[:, 0], charge_density, atol=2e-18, rtol=2e-10
    )
    assert_allclose(result.force_balance, 0, atol=2e-7, rtol=0)
    assert abs(result.charge_balance) < 1e-18
