"""Transient heat conservation, timed loads and independently solvable cases."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.solvers.transient_heat import (
    TemperatureBC,
    TransientConvection,
    solve_transient_heat,
)
from fem_toolbox.time import Ramp, Step


def material():
    return ThermalMaterial.isotropic(2.0, density=3.0, specific_heat=4.0)


def test_isolated_heating_step_and_probes(cube_mesh):
    mesh = cube_mesh(1)
    result = solve_transient_heat(
        mesh,
        material(),
        [0.0, 1.0],
        300.0,
        source=Step(0.25, 0.0, 12.0),
        probes={"center": [0.5, 0.5, 0.5], "corner": 0},
    )
    np.testing.assert_array_equal(result.times, [0.0, 0.25, 1.0])
    np.testing.assert_allclose(result.frames[-1].temperature, 300.75, atol=1e-11)
    np.testing.assert_allclose(
        result.probe_temperature[:, 0], [300.0, 300.0, 300.75], atol=1e-11
    )
    assert max(abs(b.balance) for b in result.balances) < 1e-9
    assert result.balances[-1].storage_power == pytest.approx(12.0, abs=1e-9)


def test_moving_dirichlet_storage_and_no_initial_projection(cube_mesh):
    mesh = cube_mesh(1)
    bcs = [TemperatureBC(np.arange(mesh.n_points), Ramp(0, 1, 300, 302))]
    result = solve_transient_heat(mesh, material(), [0, 0.5, 1], 300.0, conditions=bcs)
    np.testing.assert_allclose(result.frames[-1].temperature, 302.0)
    assert result.balances[-1].dirichlet_inflow == pytest.approx(24.0, abs=1e-9)
    with pytest.raises(ValidationError, match="initial"):
        solve_transient_heat(mesh, material(), [0, 1], 290.0, conditions=bcs)


def test_time_varying_flux_and_convection(cube_mesh):
    mesh = cube_mesh(1)
    # All temperatures prescribed: evaluate independent boundary power exactly.
    result = solve_transient_heat(
        mesh,
        material(),
        [0.0, 1.0],
        300.0,
        conditions=[TemperatureBC(np.arange(mesh.n_points), 300.0)],
        outward_fluxes={i: Ramp(0, 1, 0, 2) for i in range(mesh.n_triangle3)},
        convection={
            i: TransientConvection(Ramp(0, 1, 0, 3), 290.0)
            for i in range(mesh.n_triangle3)
        },
    )
    b = result.balances[-1]
    assert b.prescribed_outflow == pytest.approx(12.0)  # six m²
    assert b.convection_outflow == pytest.approx(180.0)
    assert b.dirichlet_inflow == pytest.approx(192.0)


def test_multimaterial_uniform_heating(cube_mesh):
    mesh = cube_mesh(1)
    mats = [
        ThermalMaterial.isotropic(2, density=2, specific_heat=3)
        if i % 2
        else material()
        for i in range(mesh.n_tetra4)
    ]
    result = solve_transient_heat(
        mesh,
        mats,
        [0, 1],
        300.0,
        conditions=[TemperatureBC(np.arange(mesh.n_points), Ramp(0, 1, 300, 301))],
    )
    assert result.balances[-1].storage_power == pytest.approx(9.0, abs=1e-9)


def test_invalid_storage_probe_and_memory(cube_mesh):
    mesh = cube_mesh(1)
    with pytest.raises(ValidationError, match="storage"):
        solve_transient_heat(mesh, ThermalMaterial.isotropic(1.0), [0, 1], 300.0)
    with pytest.raises(ValidationError, match="outside"):
        solve_transient_heat(mesh, material(), [0, 1], 300.0, probes={"bad": [2, 2, 2]})
    with pytest.raises(ValidationError, match="budget"):
        solve_transient_heat(mesh, material(), [0, 1], 300.0, max_history_bytes=1)


def test_spatial_time_callable_step_does_not_heat_before_event(cube_mesh):
    result = solve_transient_heat(
        cube_mesh(1),
        material(),
        [0, 0.25, 1],
        300.0,
        source=lambda x, t: 12.0 if t >= 0.25 else 0.0,
    )
    np.testing.assert_allclose(result.frames[1].temperature, 300.0, rtol=0, atol=1e-11)
    np.testing.assert_allclose(
        result.frames[-1].temperature, 300.75, rtol=0, atol=1e-11
    )


def test_insulated_diffusion_conserves_heat_and_dissipates_variance():
    from fem_toolbox.mesh import Mesh

    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], tetra4=[[0, 1, 2, 3]]
    )
    result = solve_transient_heat(
        mesh,
        ThermalMaterial.isotropic(1, density=1, specific_heat=1),
        [0, 0.01, 0.1, 1, 10],
        [300, 310, 320, 330],
    )
    # Analytic integral Ni Nj on unit reference tetra: V=1/6.
    capacity = (np.ones((4, 4)) + np.eye(4)) / 120
    energies = []
    for frame in result.frames:
        assert np.sum(capacity @ frame.temperature) == pytest.approx(315 / 6, abs=1e-10)
        deviation = frame.temperature - 315
        energies.append(deviation @ capacity @ deviation)
    assert np.all(np.diff(energies) < 0)
    assert max(abs(b.stored_energy_change) for b in result.balances) < 1e-10
