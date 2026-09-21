import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ElectricGaugeError, ValidationError
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics


def test_anisotropic_affine_patch_and_gauge_shift(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    material = DielectricMaterial(
        1e-11 * np.array([[2, 0.5, 0], [0.5, 3, 0.2], [0, 0.2, 4]])
    )
    gradient = np.array([1.0, 2, -3])
    exact = mesh.points @ gradient
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    dofs = electric_dofs(mesh)
    result = solve_electrostatics(
        mesh, material, [Dirichlet(dofs.indices("phi", boundary), exact[boundary])]
    )
    shifted = solve_electrostatics(
        mesh,
        material,
        [Dirichlet(dofs.indices("phi", boundary), exact[boundary] + 100)],
    )
    assert_allclose(result.potential, exact, atol=2e-14, rtol=0)
    assert_allclose(
        result.electric_field,
        np.tile(-gradient, (mesh.n_tetra4, 1)),
        atol=3e-14,
        rtol=0,
    )
    assert_allclose(
        result.electric_displacement,
        np.tile(np.array([-3, -5.9, 11.6]) * 1e-11, (mesh.n_tetra4, 1)),
        atol=2e-24,
        rtol=0,
    )
    assert_allclose(shifted.electric_field, result.electric_field, atol=3e-12, rtol=0)
    assert_allclose(
        result.energy,
        0.5 * gradient @ material.permittivity @ gradient,
        atol=2e-24,
        rtol=0,
    )
    assert abs(result.charge_balance) < 1e-23
    assert_allclose(
        2 * result.energy,
        np.dot(result.potential, result.electrode_charges),
        atol=2e-24,
        rtol=0,
    )
    with pytest.raises(ValueError):
        result.potential[0] = 0


def test_constant_potential(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    r = solve_electrostatics(
        mesh, DielectricMaterial.isotropic(2e-11), [Dirichlet([0], [25.0])]
    )
    assert_allclose(r.potential, 25, atol=3e-12, rtol=0)
    assert_allclose(r.electric_field, 0, atol=3e-12, rtol=0)
    assert r.energy < 1e-33


def test_parallel_plate_capacitance(cube_mesh):
    original = cube_mesh(2)
    lengths = np.array([0.002, 0.03, 0.04])
    mesh = Mesh(points=original.points * lengths, tetra4=original.tetra4)
    low = np.flatnonzero(mesh.points[:, 0] == 0)
    high = np.flatnonzero(mesh.points[:, 0] == lengths[0])
    dofs = electric_dofs(mesh)
    r = solve_electrostatics(
        mesh,
        DielectricMaterial.isotropic(2e-11),
        [
            Dirichlet(dofs.indices("phi", low), np.zeros(len(low))),
            Dirichlet(dofs.indices("phi", high), np.full(len(high), 100.0)),
        ],
    )
    capacitance = 2e-11 * 0.03 * 0.04 / 0.002
    assert_allclose(r.potential, 100 * mesh.points[:, 0] / 0.002, atol=3e-12, rtol=0)
    assert_allclose(
        r.electric_field, np.tile([-50000, 0, 0], (mesh.n_tetra4, 1)), atol=3e-9, rtol=0
    )
    assert_allclose(
        r.electrode_charges[high].sum() / 100, capacitance, atol=2e-24, rtol=0
    )
    assert_allclose(
        r.electrode_charges[low].sum(), -capacitance * 100, atol=2e-22, rtol=0
    )
    assert_allclose(r.energy, 0.5 * capacitance * 100**2, atol=2e-20, rtol=0)
    assert abs(r.charge_balance) < 2e-22


def test_single_surface_charge_sign(cube_mesh):
    mesh = cube_mesh(2)
    left = np.flatnonzero(mesh.points[:, 0] == 0)
    right = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 1, axis=1))
    r = solve_electrostatics(
        mesh,
        DielectricMaterial.isotropic(2e-11),
        [Dirichlet(electric_dofs(mesh).indices("phi", left), np.zeros(len(left)))],
        surface_charges={int(f): 2e-11 for f in right},
    )
    assert_allclose(r.potential, -mesh.points[:, 0], atol=2e-14, rtol=0)
    assert_allclose(
        r.electric_displacement,
        np.tile([2e-11, 0, 0], (mesh.n_tetra4, 1)),
        atol=2e-24,
        rtol=0,
    )
    assert_allclose(
        [r.prescribed_surface_charge, r.dirichlet_outward_charge],
        [2e-11, -2e-11],
        atol=2e-24,
        rtol=0,
    )
    assert abs(r.charge_balance) < 2e-24


@pytest.mark.parametrize("rho,compatible", [(0.0, True), (6e-11, True), (1e-11, False)])
def test_gauge_compatibility(cube_mesh, rho, compatible):
    mesh = cube_mesh(2)
    right = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 1, axis=1))
    surfaces = {int(f): 6e-11 for f in right} if rho == 6e-11 else {}
    with pytest.raises(ElectricGaugeError) as failure:
        solve_electrostatics(
            mesh,
            DielectricMaterial.isotropic(2e-11),
            [],
            volume_charge=rho,
            surface_charges=surfaces,
        )
    assert failure.value.compatible is compatible
    assert_allclose(
        failure.value.net_charge, 0 if compatible else 1e-11, atol=1e-25, rtol=0
    )


def test_disconnected_and_invalid_inputs(cube_mesh):
    base = cube_mesh(1)
    mesh = Mesh(
        points=np.vstack([base.points, base.points + np.array([2, 0, 0])]),
        tetra4=np.vstack([base.tetra4, base.tetra4 + 8]),
    )
    with pytest.raises(ElectricGaugeError) as failure:
        solve_electrostatics(
            mesh, DielectricMaterial.isotropic(2e-11), [Dirichlet([0], [0])]
        )
    assert failure.value.component == 1
    with pytest.raises(ValidationError):
        solve_electrostatics(base, [DielectricMaterial.isotropic(2e-11)], [])
    with pytest.raises(ValidationError):
        solve_electrostatics(
            base, DielectricMaterial.isotropic(2e-11), [Dirichlet([0, 0], [0, 1])]
        )
    with pytest.raises(ValidationError):
        solve_electrostatics(
            base, DielectricMaterial.isotropic(2e-11), [], surface_charges={999: 0}
        )


def test_volume_charge_reactions_and_provenance(cube_mesh):
    mesh = cube_mesh(1)
    material = DielectricMaterial.isotropic(2e-11)
    bc = [Dirichlet(np.arange(mesh.n_points), np.zeros(mesh.n_points))]
    r = solve_electrostatics(mesh, material, bc, volume_charge=6e-11)
    again = solve_electrostatics(mesh, material, bc, volume_charge=6e-11)
    changed = solve_electrostatics(mesh, material, bc, volume_charge=3e-11)
    assert_allclose(
        [r.volume_charge, r.dirichlet_outward_charge, r.electrode_charges.sum()],
        [6e-11, 6e-11, -6e-11],
        atol=1e-25,
        rtol=0,
    )
    assert r.provenance.to_dict() == again.provenance.to_dict()
    assert r.provenance.input_sha256 != changed.provenance.input_sha256


def test_series_dielectric_layers(cube_mesh):
    mesh = cube_mesh(2)
    centroids = mesh.points[mesh.tetra4].mean(axis=1)
    materials = [
        DielectricMaterial.isotropic(2e-11 if x < 0.5 else 4e-11)
        for x in centroids[:, 0]
    ]
    left = np.flatnonzero(mesh.points[:, 0] == 0)
    right = np.flatnonzero(mesh.points[:, 0] == 1)
    dofs = electric_dofs(mesh)
    r = solve_electrostatics(
        mesh,
        materials,
        [
            Dirichlet(dofs.indices("phi", left), np.ones(len(left))),
            Dirichlet(dofs.indices("phi", right), np.zeros(len(right))),
        ],
    )
    x = mesh.points[:, 0]
    assert_allclose(
        r.potential,
        1 - (8 / 3) * np.where(x <= 0.5, x / 2, 0.25 + (x - 0.5) / 4),
        atol=2e-14,
        rtol=0,
    )
    assert_allclose(
        r.electric_displacement,
        np.tile([8e-11 / 3, 0, 0], (mesh.n_tetra4, 1)),
        atol=2e-24,
        rtol=0,
    )
    assert_allclose(r.electrode_charges[left].sum(), 8e-11 / 3, atol=2e-24, rtol=0)


def test_active_node_mapping_and_interior_surface_rejection(cube_mesh):
    base = cube_mesh(1)
    mesh = Mesh(points=np.vstack([[5, 5, 5], base.points]), tetra4=base.tetra4 + 1)
    dofs = electric_dofs(mesh)
    assert dofs.size == 8
    r = solve_electrostatics(
        mesh,
        DielectricMaterial.isotropic(2e-11),
        [Dirichlet(dofs.indices("phi", np.arange(1, 9)), base.points[:, 0])],
    )
    assert_allclose(r.node_ids, np.arange(1, 9), atol=0)
    assert_allclose(
        r.electric_field, np.tile([-1, 0, 0], (mesh.n_tetra4, 1)), atol=1e-14, rtol=0
    )
    internal = Mesh(points=base.points, tetra4=base.tetra4, triangle3=[[0, 4, 7]])
    with pytest.raises(ValidationError, match="exterior"):
        solve_electrostatics(
            internal,
            DielectricMaterial.isotropic(2e-11),
            [],
            surface_charges={0: 1e-11},
        )
