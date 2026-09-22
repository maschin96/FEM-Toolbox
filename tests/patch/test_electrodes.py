"""Independent parallel-plate and piezoelectric electrode limits."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric


def _setup(cube_mesh):
    mesh = cube_mesh(1)
    left = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 0, axis=1))
    right = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 1, axis=1))
    dofs = piezoelectric_dofs(mesh)
    stiffness = np.eye(6) * 10.0
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.5
    material = PiezoelectricMaterial(stiffness, coupling, np.eye(3) * 2.0)
    displacement = np.column_stack(
        [0.1 * mesh.points[:, 0], np.zeros(mesh.n_points), np.zeros(mesh.n_points)]
    )
    bc = Dirichlet(dofs.indices("u", np.arange(mesh.n_points)), displacement.ravel())
    return mesh, material, bc, left, right, dofs


def _state(result, name):
    return next(item for item in result.electrode_states if item.name == name)


def test_fixed_electrodes_equal_manual_dirichlet_and_report_charge(cube_mesh):
    mesh, material, bc, left, right, dofs = _setup(cube_mesh)
    electrodes = [
        Electrode("left", left, potential=0),
        Electrode("right", right, potential=1),
    ]
    result = solve_piezoelectric(mesh, material, [bc], electrodes=electrodes)
    manual = solve_piezoelectric(
        mesh,
        material,
        [
            bc,
            Dirichlet(
                dofs.indices("phi", np.flatnonzero(mesh.points[:, 0] == 0)), np.zeros(4)
            ),
            Dirichlet(
                dofs.indices("phi", np.flatnonzero(mesh.points[:, 0] == 1)), np.ones(4)
            ),
        ],
    )
    assert_allclose(result.potential, manual.potential, atol=1e-13, rtol=0)
    assert_allclose(
        result.electric_displacement, manual.electric_displacement, atol=1e-12, rtol=0
    )
    assert_allclose(_state(result, "right").charge, 1.95, atol=1e-12, rtol=0)
    assert abs(result.charge_balance) < 1e-12
    assert (
        result.provenance.solver_parameters["input_schema"]
        == "piezoelectric-discrete-v2"
    )
    assert (
        manual.provenance.solver_parameters["input_schema"]
        == "piezoelectric-discrete-v1"
    )
    assert result.provenance.input_sha256 != manual.provenance.input_sha256


def test_open_and_short_circuit_have_distinct_analytic_limits(cube_mesh):
    mesh, material, bc, left, right, _ = _setup(cube_mesh)
    common = Electrode("left", left, potential=0)
    open_result = solve_piezoelectric(
        mesh, material, [bc], electrodes=[common, Electrode("right", right, charge=0)]
    )
    short_result = solve_piezoelectric(
        mesh,
        material,
        [bc],
        electrodes=[common, Electrode("right", right, potential=0)],
    )
    assert_allclose(_state(open_result, "right").potential, 0.025, atol=1e-12, rtol=0)
    assert_allclose(_state(open_result, "right").charge, 0, atol=1e-12, rtol=0)
    assert_allclose(
        open_result.potential[mesh.points[:, 0] == 1],
        _state(open_result, "right").potential,
        atol=1e-13,
        rtol=0,
    )
    assert_allclose(_state(short_result, "right").potential, 0, atol=1e-13, rtol=0)
    assert_allclose(_state(short_result, "right").charge, -0.05, atol=1e-12, rtol=0)
    assert_allclose(open_result.electric_displacement, 0, atol=1e-12, rtol=0)
    assert_allclose(short_result.electric_displacement[:, 0], 0.05, atol=1e-12, rtol=0)
    assert abs(open_result.charge_balance) < 1e-12
    assert abs(short_result.charge_balance) < 1e-12


def test_prescribed_floating_charge_matches_capacitance(cube_mesh):
    mesh, material, bc, left, right, _ = _setup(cube_mesh)
    result = solve_piezoelectric(
        mesh,
        material,
        [bc],
        electrodes=[
            Electrode("left", left, potential=0),
            Electrode("right", right, charge=1.95),
        ],
    )
    assert_allclose(_state(result, "right").potential, 1, atol=1e-12, rtol=0)
    assert_allclose(_state(result, "right").charge, 1.95, atol=1e-12, rtol=0)


def test_overlapping_or_unreferenced_electrodes_are_rejected(cube_mesh):
    mesh, material, bc, left, right, _ = _setup(cube_mesh)
    with pytest.raises(ValidationError, match="overlap"):
        solve_piezoelectric(
            mesh,
            material,
            [bc],
            electrodes=[
                Electrode("a", left, potential=0),
                Electrode("b", left, charge=0),
            ],
        )
    with pytest.raises(ValidationError, match="gauge"):
        solve_piezoelectric(
            mesh, material, [bc], electrodes=[Electrode("right", right, charge=0)]
        )


def test_named_facet_sets_and_dirichlet_conflicts(cube_mesh):
    mesh, material, bc, left, right, dofs = _setup(cube_mesh)
    mesh = Mesh(
        points=mesh.points,
        tetra4=mesh.tetra4,
        triangle3=mesh.triangle3,
        facet_sets={"left": left, "right": right},
    )
    result = solve_piezoelectric(
        mesh,
        material,
        [bc],
        electrodes=[
            Electrode("left", "left", potential=0),
            Electrode("right", "right", charge=0),
        ],
    )
    assert_allclose(_state(result, "right").potential, 0.025, atol=1e-12, rtol=0)
    right_node = int(np.flatnonzero(mesh.points[:, 0] == 1)[0])
    with pytest.raises(ValidationError, match="conflicts"):
        solve_piezoelectric(
            mesh,
            material,
            [bc, Dirichlet(dofs.indices("phi", [right_node]), [0])],
            electrodes=[
                Electrode("left", "left", potential=0),
                Electrode("right", "right", charge=0),
            ],
        )
    with pytest.raises(ValidationError, match="conflicting Dirichlet"):
        solve_piezoelectric(
            mesh,
            material,
            [bc, Dirichlet(dofs.indices("phi", [right_node]), [1])],
            electrodes=[Electrode("right", "right", potential=0)],
        )
