"""Native PETSc parity and independent patch/oscillator references in SI units."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalMaterial,
)
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from fem_toolbox.solvers.electrostatics import solve_electrostatics
from fem_toolbox.solvers.heat import solve_heat
from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes
from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
from fem_toolbox.solvers.piezo_transient import solve_transient_piezoelectric
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from fem_toolbox.solvers.transient_heat import solve_transient_heat


def test_single_field_patches_and_heat_storage(cube_mesh, physics_backend):
    mesh = cube_mesh(2)
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    x = mesh.points[:, 0]
    thermal = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    for backend in (None, physics_backend):
        heat = solve_heat(
            mesh, thermal, [Dirichlet(boundary, 300 + x[boundary])], backend=backend
        )
        assert_allclose(heat.temperature, 300 + x, rtol=0, atol=1e-10)
        assert_allclose(
            heat.heat_flux,
            np.tile([-2, 0, 0], (len(mesh.tetra4), 1)),
            rtol=0,
            atol=1e-10,
        )
        electric = solve_electrostatics(
            mesh,
            DielectricMaterial.isotropic(2),
            [Dirichlet(boundary, x[boundary])],
            backend=backend,
        )
        assert_allclose(electric.potential, x, rtol=0, atol=1e-13)
        u = 0.01 * mesh.points
        bc = Dirichlet(
            mechanical_dofs(mesh).indices("u", boundary), u[boundary].ravel()
        )
        elastic = solve_elasticity(
            mesh, ElasticMaterial.isotropic(100, 0.25), [bc], backend=backend
        )
        assert_allclose(elastic.displacement, u, rtol=0, atol=1e-13)
        history = solve_transient_heat(
            mesh, thermal, [0, 0.5, 1], 300, source=12, backend=backend
        )
        for frame in history.frames:
            assert_allclose(frame.temperature, 300 + frame.time, rtol=0, atol=1e-10)
        for balance in history.balances:
            assert abs(balance.balance) < 1e-10


@pytest.mark.parametrize("floating", [False, True])
def test_fixed_and_charged_floating_electrodes(cube_mesh, physics_backend, floating):
    mesh = cube_mesh(2)
    x = mesh.points[:, 0]
    left = np.flatnonzero(np.all(x[mesh.triangle3] == 0, axis=1))
    right = np.flatnonzero(np.all(x[mesh.triangle3] == 1, axis=1))
    fixed = piezoelectric_dofs(mesh).indices("u", np.flatnonzero(x == 0))
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.5
    material = PiezoelectricMaterial(np.eye(6) * 10, coupling, np.eye(3) * 2)
    kwargs = dict(
        electrodes=[
            Electrode("left", left, potential=0),
            Electrode("right", right, charge=0.3)
            if floating
            else Electrode("right", right, potential=1),
        ]
    )
    bc = [Dirichlet(fixed, np.zeros(len(fixed)))]
    reference = solve_piezoelectric(mesh, material, bc, **kwargs)
    result = solve_piezoelectric(mesh, material, bc, backend=physics_backend, **kwargs)
    for field in (
        "displacement",
        "potential",
        "stress",
        "electric_displacement",
        "electrode_charges",
        "force_balance",
        "moment_balance",
        "mechanical_energy",
        "coupling_energy",
        "dielectric_enthalpy",
    ):
        assert_allclose(
            getattr(result, field), getattr(reference, field), rtol=1e-11, atol=1e-12
        )
    assert max(result.relative_residual_u, result.relative_residual_phi) < 1e-10
    assert abs(result.charge_balance) < 1e-12
    assert result.provenance.solver_parameters["backend"] == physics_backend.name


def test_modal_newmark_and_complex_harmonic(physics_backend):
    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], tetra4=[[0, 1, 2, 3]]
    )
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 2
    material = PiezoelectricMaterial(
        np.diag([10, 20, 30, 4, 5, 6]), coupling, np.diag([3, 4, 5])
    )
    fixed = np.setdiff1d(np.arange(16), [3, 13])
    bc = [Dirichlet(fixed, np.zeros(len(fixed)))]
    # V=1/6, M11=rho*V/10=0.1; open-circuit stiffness (10+2²/3)*V.
    mass, stiffness = 0.1, (10 + 4 / 3) / 6
    f = np.zeros(16, dtype=complex)
    f[3] = 1 + 2j
    histories = []
    for backend in (None, physics_backend):
        modal = solve_piezoelectric_modes(
            mesh, material, 6, bc, count=1, backend=backend
        )
        assert_allclose(
            modal.angular_frequencies, [np.sqrt(stiffness / mass)], rtol=1e-12
        )
        harmonic = solve_harmonic_piezoelectric(
            mesh, material, 6, bc, 0.1, loads=f, backend=backend
        )
        u = (1 + 2j) / (stiffness - (2 * np.pi * 0.1) ** 2 * mass)
        assert_allclose(harmonic.displacement[1, 0], u, rtol=1e-12)
        assert_allclose(harmonic.potential[1], 2 / 3 * u, rtol=1e-12)
        initial = np.zeros((4, 3))
        initial[1, 0] = 0.1
        transient = solve_transient_piezoelectric(
            mesh,
            material,
            6,
            bc,
            [0, 0.05, 0.1],
            initial_displacement=initial,
            backend=backend,
        )
        histories.append(transient)
        assert_allclose(transient.total_energy, 0.5 * stiffness * 0.1**2, rtol=1e-12)
        for result in (modal, harmonic, transient):
            assert np.max(result.relative_residual_u) < 1e-10
            assert np.max(result.relative_residual_phi) < 1e-10
    for field in (
        "displacement",
        "velocity",
        "acceleration",
        "potential",
        "total_energy",
    ):
        assert_allclose(
            getattr(histories[1], field),
            getattr(histories[0], field),
            rtol=1e-11,
            atol=1e-12,
        )


def test_tet10_heat_stationary_and_transient(cube_mesh, physics_backend):
    from tests.reference.quadratic_mesh import elevate

    mesh = elevate(cube_mesh(2))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    boundary = np.unique(mesh.facets)
    exact = 300 + mesh.points[:, 0] ** 2
    args = [mesh, material, [Dirichlet(boundary, exact[boundary])]]
    reference = solve_heat(*args, source=-4)
    result = solve_heat(*args, source=-4, backend=physics_backend)
    assert_allclose(result.temperature, exact, rtol=0, atol=1e-9)
    assert_allclose(result.temperature, reference.temperature, rtol=0, atol=1e-9)
    assert_allclose(result.heat_flux, reference.heat_flux, rtol=0, atol=1e-9)
    assert abs(result.power_balance) < 1e-9
    history = solve_transient_heat(
        mesh, material, [0, 0.2, 1], 300, source=12, backend=physics_backend
    )
    assert_allclose(history.frames[-1].temperature, 301, rtol=0, atol=1e-9)
    assert max(abs(b.balance) for b in history.balances) < 1e-9


def test_tet10_electrostatics_patch_and_energy(cube_mesh, physics_backend):
    from tests.reference.quadratic_mesh import elevate

    mesh = elevate(cube_mesh(2))
    material = DielectricMaterial.isotropic(2e-11)
    boundary = np.unique(mesh.facets)
    exact = mesh.points[:, 0] ** 2
    args = [mesh, material, [Dirichlet(boundary, exact[boundary])]]
    reference = solve_electrostatics(*args, volume_charge=-4e-11)
    result = solve_electrostatics(*args, volume_charge=-4e-11, backend=physics_backend)
    assert_allclose(result.potential, exact, rtol=0, atol=1e-10)
    assert_allclose(result.electric_field, reference.electric_field, rtol=0, atol=1e-9)
    assert_allclose(
        result.element_energy, reference.element_energy, rtol=1e-9, atol=1e-24
    )
    assert abs(result.energy - 4e-11 / 3) < 1e-21
    assert abs(result.charge_balance) < 1e-21


def test_tet10_elasticity_patch_and_energy(cube_mesh, physics_backend):
    from tests.reference.quadratic_mesh import elevate

    mesh = elevate(cube_mesh(2))
    boundary = np.unique(mesh.facets)
    exact = np.zeros_like(mesh.points)
    exact[:, 0] = 0.01 * mesh.points[:, 0] ** 2
    material = ElasticMaterial.isotropic(100, 0.25)
    bc = [
        Dirichlet(mechanical_dofs(mesh).indices("u", boundary), exact[boundary].ravel())
    ]
    reference = solve_elasticity(mesh, material, bc, body_force=[-2.4, 0, 0])
    result = solve_elasticity(
        mesh, material, bc, body_force=[-2.4, 0, 0], backend=physics_backend
    )
    assert_allclose(result.displacement, exact, rtol=0, atol=1e-11)
    assert_allclose(result.stress, reference.stress, rtol=0, atol=1e-8)
    assert_allclose(
        result.element_energy, reference.element_energy, rtol=1e-9, atol=1e-13
    )
    assert abs(result.energy - 0.008) < 1e-11
    assert_allclose(result.force_balance, 0, rtol=0, atol=1e-9)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=1e-9)


def test_tet10_thermomechanical_transfer_and_energy(cube_mesh, physics_backend):
    from fem_toolbox.materials import ThermalExpansion
    from fem_toolbox.solvers.heat import thermal_dofs
    from fem_toolbox.solvers.thermomechanics import solve_thermomechanics
    from tests.reference.quadratic_mesh import elevate

    mesh = elevate(cube_mesh(2))
    fixed = np.unique(mesh.facets)
    args = (
        mesh,
        ThermalMaterial.isotropic(2),
        ElasticMaterial.isotropic(100, 0.25),
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 1),
        [
            Dirichlet(
                thermal_dofs(mesh).indices("T", fixed), 1 + mesh.points[fixed, 0] ** 2
            )
        ],
        [
            Dirichlet(
                mechanical_dofs(mesh).indices("u", fixed), np.zeros(3 * len(fixed))
            )
        ],
    )
    reference = solve_thermomechanics(*args, source=-4)
    result = solve_thermomechanics(
        *args,
        source=-4,
        thermal_backend=physics_backend,
        mechanical_backend=physics_backend,
    )
    assert_allclose(
        result.heat.temperature, 1 + mesh.points[:, 0] ** 2, rtol=0, atol=1e-10
    )
    assert_allclose(
        result.elasticity.displacement,
        reference.elasticity.displacement,
        rtol=0,
        atol=1e-11,
    )
    assert_allclose(
        result.elasticity.stress, reference.elasticity.stress, rtol=0, atol=1e-8
    )
    assert_allclose(
        result.elasticity.element_energy,
        reference.elasticity.element_energy,
        rtol=1e-9,
        atol=1e-13,
    )
    assert abs(result.heat.power_balance) < 1e-9
    assert_allclose(result.elasticity.force_balance, 0, rtol=0, atol=1e-9)
    assert_allclose(result.elasticity.moment_balance, 0, rtol=0, atol=1e-9)


def test_tet10_quadratic_piezo_and_floating_electrode(cube_mesh, physics_backend):
    from tests.reference.quadratic_mesh import elevate
    from tests.unit.test_tet10_piezoelectric_integration import material, quadratic_case

    mesh = elevate(cube_mesh(2))
    reference = quadratic_case(mesh)
    result = quadratic_case(mesh, physics_backend)
    for field in (
        "displacement",
        "potential",
        "stress",
        "electric_displacement",
        "mechanical_energy",
        "coupling_energy",
        "dielectric_enthalpy",
    ):
        assert_allclose(
            getattr(result, field), getattr(reference, field), rtol=2e-8, atol=2e-10
        )
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    result = solve_piezoelectric(
        mesh,
        material(),
        [
            Dirichlet(
                piezoelectric_dofs(mesh).indices("u", left), np.zeros(3 * len(left))
            )
        ],
        electrodes=[
            Electrode("ground", "left", potential=0),
            Electrode("top", "right", charge=3.04),
        ],
        backend=physics_backend,
    )
    assert_allclose(result.potential, mesh.points[:, 0], rtol=0, atol=2e-9)
    assert_allclose(
        result.displacement[:, 0], -0.02 * mesh.points[:, 0], rtol=0, atol=2e-10
    )
    assert_allclose(result.force_balance, 0, rtol=0, atol=2e-9)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-9)
    assert abs(result.charge_balance) < 2e-9


@pytest.mark.parametrize("electric", [False, True])
def test_tet10_thermopiezo_sequential_backend_parity(
    cube_mesh, physics_backend, electric
):
    from tests.reference.quadratic_mesh import elevate
    from tests.unit.test_tet10_thermopiezoelectric_integration import coupled_case

    mesh = elevate(cube_mesh(2))
    reference = coupled_case(mesh, electric)
    result = coupled_case(mesh, electric, physics_backend, physics_backend)
    assert_allclose(
        result.heat.temperature, reference.heat.temperature, rtol=0, atol=2e-9
    )
    for field in (
        "displacement",
        "potential",
        "stress",
        "electric_displacement",
        "mechanical_energy",
        "coupling_energy",
        "dielectric_enthalpy",
    ):
        assert_allclose(
            getattr(result.piezoelectric, field),
            getattr(reference.piezoelectric, field),
            rtol=2e-8,
            atol=2e-10,
        )
    assert abs(result.heat.power_balance) < 2e-9
    assert abs(result.piezoelectric.charge_balance) < 2e-9
    assert_allclose(result.piezoelectric.force_balance, 0, rtol=0, atol=2e-9)
    assert_allclose(result.piezoelectric.moment_balance, 0, rtol=0, atol=2e-9)
