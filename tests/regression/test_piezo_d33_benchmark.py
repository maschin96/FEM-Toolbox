"""Analytical d33 bar: voltage drive, force drive, and a simple open circuit."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric

LENGTH = 0.02  # m
WIDTH = 0.01  # m
AREA = WIDTH**2
C33 = 1e9  # Pa
E33 = 5.0  # C/m²
K33 = 1e-8  # F/m


def bar(cube_mesh, n):
    unit = cube_mesh(n)
    return Mesh(
        points=unit.points * [WIDTH, WIDTH, LENGTH],
        tetra4=unit.tetra4,
        triangle3=unit.triangle3,
    )


def material():
    stiffness = np.diag([2e9, 2e9, C33, 5e8, 5e8, 5e8])
    coupling = np.zeros((3, 6))
    coupling[2, 2] = E33
    return PiezoelectricMaterial(stiffness, coupling, np.diag([K33] * 3))


def bar_conditions(mesh, *, voltage_top=None):
    z = mesh.points[:, 2]
    nodes = np.arange(mesh.n_points)
    bottom = np.flatnonzero(z == 0)
    top = np.flatnonzero(z == LENGTH)
    dofs = piezoelectric_dofs(mesh)
    conditions = [
        Dirichlet(dofs.indices("u", nodes, ["x", "y"]), np.zeros(2 * len(nodes))),
        Dirichlet(dofs.indices("u", bottom, ["z"]), np.zeros(len(bottom))),
        Dirichlet(dofs.indices("phi", bottom), np.zeros(len(bottom))),
    ]
    if voltage_top is not None:
        conditions.append(
            Dirichlet(dofs.indices("phi", top), np.full(len(top), voltage_top))
        )
    return conditions, top


def top_tractions(mesh, force):
    top_faces = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 2] == LENGTH, axis=1))
    return {int(i): [0, 0, force / AREA] for i in top_faces}


def test_direct_voltage_and_inverse_short_circuit_reciprocity(cube_mesh):
    mesh = bar(cube_mesh, 2)
    drive_voltage = 100.0
    drive_force = 10.0
    direct_bc, top = bar_conditions(mesh, voltage_top=drive_voltage)
    short_bc, _ = bar_conditions(mesh, voltage_top=0.0)
    direct = solve_piezoelectric(mesh, material(), direct_bc)
    inverse = solve_piezoelectric(
        mesh, material(), short_bc, tractions=top_tractions(mesh, drive_force)
    )
    d33 = E33 / C33
    capacitance_free = AREA / LENGTH * (K33 + E33**2 / C33)
    assert_allclose(direct.displacement[top, 2], -d33 * drive_voltage, rtol=2e-11)
    assert_allclose(
        direct.electrode_charges[top].sum(),
        capacitance_free * drive_voltage,
        rtol=2e-11,
    )
    assert_allclose(
        inverse.displacement[top, 2], drive_force * LENGTH / (AREA * C33), rtol=2e-11
    )
    assert_allclose(
        inverse.electrode_charges[top].sum(), -d33 * drive_force, rtol=2e-11
    )
    assert_allclose(
        direct.displacement[top, 2].mean() / drive_voltage,
        inverse.electrode_charges[top].sum() / drive_force,
        rtol=2e-11,
    )
    assert_allclose(
        direct.total_enthalpy, -0.5 * capacitance_free * drive_voltage**2, rtol=2e-11
    )
    assert_allclose(
        inverse.total_enthalpy, 0.5 * drive_force**2 * LENGTH / (AREA * C33), rtol=2e-11
    )
    assert_allclose(direct.force_balance, 0, atol=1e-12)
    assert_allclose(inverse.force_balance, 0, atol=1e-12)
    assert abs(direct.charge_balance) < 1e-18
    assert abs(inverse.charge_balance) < 1e-18


def test_short_and_open_circuit_have_distinct_force_response(cube_mesh):
    mesh = bar(cube_mesh, 2)
    force = 10.0
    short_bc, top = bar_conditions(mesh, voltage_top=0.0)
    open_bc, _ = bar_conditions(mesh)
    traction = top_tractions(mesh, force)
    short = solve_piezoelectric(mesh, material(), short_bc, tractions=traction)
    open_result = solve_piezoelectric(mesh, material(), open_bc, tractions=traction)
    effective_stiffness = C33 + E33**2 / K33
    assert_allclose(
        short.displacement[top, 2], force * LENGTH / (AREA * C33), rtol=2e-11
    )
    assert_allclose(
        open_result.displacement[top, 2],
        force * LENGTH / (AREA * effective_stiffness),
        rtol=2e-11,
    )
    assert_allclose(
        open_result.potential[top],
        E33 * force * LENGTH / (AREA * K33 * effective_stiffness),
        rtol=2e-11,
    )
    assert_allclose(open_result.electric_displacement[:, 2], 0, atol=1e-17)
    assert_allclose(
        open_result.total_enthalpy,
        0.5 * force**2 * LENGTH / (AREA * effective_stiffness),
        rtol=2e-11,
    )
    assert np.ptp(open_result.potential[top]) < 1e-10
    assert abs(open_result.charge_balance) < 1e-18
    assert open_result.relative_residual_u < 1e-12
    assert open_result.relative_residual_phi < 1e-12
