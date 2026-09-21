"""Independent thermal patch, power balance and nullspace references."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ThermalNullspaceError, ValidationError
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.heat import Convection
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs


def faces_on(mesh, axis, value):
    return np.flatnonzero(np.all(mesh.points[mesh.triangle3, axis] == value, axis=1))


def test_anisotropic_affine_patch(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    mat = ThermalMaterial([[2, 0.5, 0], [0.5, 3, 0.2], [0, 0.2, 4]])
    gradient = np.array([1.0, 2, -3])
    exact = 300 + mesh.points @ gradient
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    dofs = thermal_dofs(mesh)
    result = solve_heat(
        mesh, mat, [Dirichlet(dofs.indices("T", boundary), exact[boundary])]
    )
    assert_allclose(result.temperature, exact, atol=3e-12, rtol=0)
    assert_allclose(
        result.gradient, np.tile(gradient, (mesh.n_tetra4, 1)), atol=1e-11, rtol=0
    )
    assert_allclose(
        result.heat_flux,
        np.tile([-3, -5.9, 11.6], (mesh.n_tetra4, 1)),
        atol=5e-11,
        rtol=0,
    )
    assert abs(result.power_balance) < 1e-10
    assert result.relative_residual < 1e-12
    with pytest.raises(ValueError):
        result.temperature[0] = 0


def test_prescribed_outward_flux_sign(cube_mesh):
    mesh = cube_mesh(2)
    left = np.flatnonzero(mesh.points[:, 0] == 0)
    right = faces_on(mesh, 0, 1)
    result = solve_heat(
        mesh,
        ThermalMaterial.isotropic(2),
        [Dirichlet(thermal_dofs(mesh).indices("T", left), np.full(len(left), 400.0))],
        outward_fluxes={int(f): 20.0 for f in right},
    )
    assert_allclose(
        result.temperature, 400 - 10 * mesh.points[:, 0], atol=3e-12, rtol=0
    )
    assert_allclose(
        result.heat_flux, np.tile([20, 0, 0], (mesh.n_tetra4, 1)), atol=2e-11, rtol=0
    )
    assert_allclose(
        [result.prescribed_outflow, result.dirichlet_inflow],
        [20, 20],
        atol=1e-10,
        rtol=0,
    )
    assert abs(result.power_balance) < 1e-10


@pytest.mark.parametrize("lengths", [(1.0, 1.0, 1.0), (2.0, 3.0, 4.0)])
def test_convection_quader(cube_mesh, lengths):
    original = cube_mesh(2)
    mesh = Mesh(
        points=original.points * np.array(lengths),
        tetra4=original.tetra4,
        triangle3=original.triangle3,
    )
    left = np.flatnonzero(mesh.points[:, 0] == 0)
    right = faces_on(mesh, 0, lengths[0])
    result = solve_heat(
        mesh,
        ThermalMaterial.isotropic(2),
        [Dirichlet(thermal_dofs(mesh).indices("T", left), np.full(len(left), 400.0))],
        convection={int(f): Convection(4, 300) for f in right},
    )
    q = 100 / (lengths[0] / 2 + 1 / 4)
    assert_allclose(
        result.temperature, 400 - q / 2 * mesh.points[:, 0], atol=3e-12, rtol=0
    )
    assert_allclose(
        [result.convection_outflow, result.dirichlet_inflow],
        [q * lengths[1] * lengths[2]] * 2,
        atol=1e-10,
        rtol=0,
    )
    assert abs(result.power_balance) < 1e-10


def test_convection_removes_nullspace_and_constant_temperature(cube_mesh):
    mesh = cube_mesh(2)
    result = solve_heat(
        mesh,
        ThermalMaterial.isotropic(2),
        [],
        convection={int(f): Convection(4, 300) for f in faces_on(mesh, 0, 1)},
    )
    assert_allclose(result.temperature, 300, atol=3e-12, rtol=0)
    assert_allclose(result.heat_flux, 0, atol=2e-11, rtol=0)
    assert abs(result.convection_outflow) < 1e-10


def test_source_reactions_all_fixed_and_provenance(cube_mesh):
    mesh = cube_mesh(1)
    bc = [Dirichlet(np.arange(mesh.n_points), np.full(mesh.n_points, 300.0))]
    first = solve_heat(mesh, ThermalMaterial.isotropic(2), bc, source=6.0)
    again = solve_heat(mesh, ThermalMaterial.isotropic(2), bc, source=6.0)
    changed = solve_heat(mesh, ThermalMaterial.isotropic(2), bc, source=3.0)
    assert_allclose(
        [first.source_power, first.dirichlet_inflow], [6, -6], atol=1e-11, rtol=0
    )
    assert abs(first.power_balance) < 1e-11
    assert first.provenance.to_dict() == again.provenance.to_dict()
    assert first.provenance.input_sha256 != changed.provenance.input_sha256


@pytest.mark.parametrize("source,balanced", [(0.0, True), (6.0, True), (1.0, False)])
def test_neumann_compatibility_is_diagnosed(cube_mesh, source, balanced):
    mesh = cube_mesh(2)
    fluxes = {int(f): 6.0 for f in faces_on(mesh, 0, 1)} if source == 6 else {}
    with pytest.raises(ThermalNullspaceError) as failure:
        solve_heat(
            mesh, ThermalMaterial.isotropic(2), [], source=source, outward_fluxes=fluxes
        )
    assert failure.value.compatible is balanced
    assert failure.value.component == 0
    assert_allclose(failure.value.net_power, 0 if balanced else 1, atol=1e-13)
    with pytest.raises(ThermalNullspaceError):
        solve_heat(
            mesh, ThermalMaterial.isotropic(2), [], convection={0: Convection(0, 300)}
        )


def test_disconnected_components_checked_individually(cube_mesh):
    original = cube_mesh(1)
    mesh = Mesh(
        points=np.vstack([original.points, original.points + np.array([2, 0, 0])]),
        tetra4=np.vstack([original.tetra4, original.tetra4 + original.n_points]),
    )
    with pytest.raises(ThermalNullspaceError) as failure:
        solve_heat(
            mesh,
            ThermalMaterial.isotropic(2),
            [],
            source=lambda x: 1.0 if x[0] < 1.5 else -1.0,
        )
    assert (
        not failure.value.compatible
    )  # Global cancellation cannot anchor either cube.
    with pytest.raises(ThermalNullspaceError) as failure:
        solve_heat(mesh, ThermalMaterial.isotropic(2), [Dirichlet([0], [300])])
    assert failure.value.component == 1


def test_invalid_boundary_and_inactive_nodes(cube_mesh):
    original = cube_mesh(1)
    mesh = Mesh(
        points=np.vstack([original.points, [5, 5, 5]]),
        tetra4=original.tetra4,
        triangle3=[[0, 4, 7]],
    )  # Internal face.
    assert thermal_dofs(mesh).size == 8
    with pytest.raises(ValidationError, match="exterior"):
        solve_heat(mesh, ThermalMaterial.isotropic(2), [], outward_fluxes={0: 1.0})
    with pytest.raises(ValidationError):
        solve_heat(mesh, ThermalMaterial.isotropic(2), [Dirichlet([8], [300])])
    with pytest.raises(ValidationError):
        solve_heat(mesh, [ThermalMaterial.isotropic(2)], [])


def test_two_material_layers(cube_mesh):
    mesh = cube_mesh(2)
    centroids = mesh.points[mesh.tetra4].mean(axis=1)
    materials = [
        ThermalMaterial.isotropic(2 if x < 0.5 else 4) for x in centroids[:, 0]
    ]
    left = np.flatnonzero(mesh.points[:, 0] == 0)
    right = np.flatnonzero(mesh.points[:, 0] == 1)
    dofs = thermal_dofs(mesh)
    result = solve_heat(
        mesh,
        materials,
        [
            Dirichlet(dofs.indices("T", left), np.full(len(left), 400.0)),
            Dirichlet(dofs.indices("T", right), np.full(len(right), 300.0)),
        ],
    )
    x = mesh.points[:, 0]
    expected = 400 - (800 / 3) * np.where(x <= 0.5, x / 2, 0.25 + (x - 0.5) / 4)
    assert_allclose(result.temperature, expected, atol=3e-12, rtol=0)
    assert_allclose(
        result.heat_flux,
        np.tile([800 / 3, 0, 0], (mesh.n_tetra4, 1)),
        atol=5e-11,
        rtol=0,
    )
    assert_allclose(
        [result.reactions[left].sum(), result.reactions[right].sum()],
        [800 / 3, -800 / 3],
        atol=1e-10,
        rtol=0,
    )
