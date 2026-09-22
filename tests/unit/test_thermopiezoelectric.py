"""Independent one-tetra checks for stationary T -> (u, phi) transfer."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import (
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
)
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.heat import solve_heat
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from fem_toolbox.solvers.thermopiezoelectric import solve_thermopiezoelectric


def test_clamped_thermal_open_circuit_voltage_from_constitutive_law() -> None:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
    )
    e = np.zeros((3, 6))
    e[0, 0] = 2.0
    material = PiezoelectricMaterial(np.eye(6) * 10, e, np.eye(3) * 4)
    dofs = piezoelectric_dofs(mesh)
    clamped = Dirichlet(dofs.indices("u", [0, 1, 2, 3]), np.zeros(12))
    grounded = Dirichlet(dofs.indices("phi", [0, 2, 3]), np.zeros(3))
    result = solve_thermopiezoelectric(
        mesh,
        ThermalMaterial.isotropic(2),
        material,
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
        [Dirichlet([0, 1, 2, 3], [310, 310, 310, 310])],
        [clamped, grounded],
        electric_thermal_strain=True,
    )
    assert_allclose(result.heat.temperature, 310, atol=1e-12, rtol=0)
    # D_x=e_xx*(epsilon_x-alpha*dT)+kappa*E_x=0; E_x=-phi(1).
    assert_allclose(result.piezoelectric.potential[1], -0.05, atol=1e-12, rtol=0)
    assert_allclose(result.piezoelectric.electric_displacement, 0, atol=1e-12, rtol=0)
    assert abs(result.heat.power_balance) < 1e-12
    assert np.linalg.norm(result.piezoelectric.force_balance) < 1e-12
    assert abs(result.piezoelectric.charge_balance) < 1e-12
    without_electric_thermal_term = solve_thermopiezoelectric(
        mesh,
        ThermalMaterial.isotropic(2),
        material,
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
        [Dirichlet([0, 1, 2, 3], [310, 310, 310, 310])],
        [clamped, grounded],
        electric_thermal_strain=False,
    )
    assert_allclose(
        without_electric_thermal_term.piezoelectric.potential, 0, atol=1e-12, rtol=0
    )


def test_thermal_and_electric_loads_superpose_with_fixed_fields() -> None:
    mesh = Mesh(
        points=np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float),
        tetra4=np.array([[0, 1, 2, 3]]),
    )
    e = np.zeros((3, 6))
    e[0, 0] = 2
    material = PiezoelectricMaterial(np.eye(6) * 10, e, np.eye(3) * 4)
    dofs = piezoelectric_dofs(mesh)
    u = Dirichlet(dofs.indices("u", [0, 1, 2, 3]), np.zeros(12))

    def run(temperature: float, voltage: float):
        phi = Dirichlet(dofs.indices("phi", [0, 1, 2, 3]), [0, voltage, 0, 0])
        return solve_thermopiezoelectric(
            mesh,
            ThermalMaterial.isotropic(2),
            material,
            ThermalExpansion([0.01, 0, 0, 0, 0, 0], 300),
            [Dirichlet([0, 1, 2, 3], [temperature] * 4)],
            [u, phi],
            electric_thermal_strain=True,
        )

    reference = run(300, 0).piezoelectric
    thermal = run(310, 0).piezoelectric
    electric = run(300, 1).piezoelectric
    combined = run(310, 1).piezoelectric
    assert_allclose(
        combined.stress - reference.stress,
        thermal.stress - reference.stress + electric.stress - reference.stress,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        combined.electric_displacement - reference.electric_displacement,
        thermal.electric_displacement
        - reference.electric_displacement
        + electric.electric_displacement
        - reference.electric_displacement,
        atol=1e-12,
        rtol=0,
    )
    # Independent constitutive values: sigma_x=-C_xx*alpha*dT+e*voltage;
    # D_x=-e*alpha*dT-kappa*voltage for E_x=-voltage.
    assert_allclose(combined.stress[0, 0], 1.0, atol=1e-12, rtol=0)
    assert_allclose(combined.electric_displacement[0, 0], -4.2, atol=1e-12, rtol=0)
    assert_allclose(
        reference.stress,
        solve_piezoelectric(
            mesh,
            material,
            [u, Dirichlet(dofs.indices("phi", [0, 1, 2, 3]), np.zeros(4))],
        ).stress,
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        run(300, 0).heat.temperature,
        solve_heat(
            mesh,
            ThermalMaterial.isotropic(2),
            [Dirichlet([0, 1, 2, 3], [300] * 4)],
        ).temperature,
        atol=1e-12,
        rtol=0,
    )
