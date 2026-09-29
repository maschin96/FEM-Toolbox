"""Independent residual identities and deterministic bulk marking."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.heat import Convection
from fem_toolbox.post import heat_indicator, mark_bulk


def test_affine_patch_and_material_flux_continuity(cube_mesh):
    mesh = cube_mesh(2)
    left = mesh.points[mesh.tetra4].mean(axis=1)[:, 0] < 0.5
    materials = [ThermalMaterial.isotropic(1 if side else 8) for side in left]
    x = mesh.points[:, 0]
    temperature = np.where(x <= 0.5, x, 0.5 + (x - 0.5) / 8)
    estimate = heat_indicator(
        mesh, temperature, materials, dirichlet_facets=np.arange(mesh.n_triangle3)
    )
    # Gradient jumps eightfold; physical normal heat flux remains continuous.
    assert estimate.total < 1e-13
    assert mark_bulk(estimate.cell).size == 0


def test_single_tetrahedron_volume_and_boundary_terms():
    mesh = Mesh(
        points=np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        ),
        tetra4=np.array([[0, 1, 2, 3]]),
        triangle3=np.array([[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]),
    )
    material = ThermalMaterial.isotropic(2)
    # h_K²=2, V=1/6, Q=3 -> 2/2 * 9/6 = 1.5 W K.
    estimate = heat_indicator(
        mesh, np.zeros(4), material, source=3.0, dirichlet_facets=[0, 1, 2, 3]
    )
    np.testing.assert_allclose(estimate.volume_squared, [1.5], rtol=1e-14)
    np.testing.assert_allclose(estimate.cell, [np.sqrt(1.5)], rtol=1e-14)
    # z=0 face: q_out=2 for T=z; impose 1 -> mismatch²=1,
    # contribution h_F/k * A = sqrt(2)/4.
    estimate = heat_indicator(
        mesh,
        mesh.points[:, 2],
        material,
        dirichlet_facets=[1, 2, 3],
        outward_fluxes={0: 1.0},
    )
    np.testing.assert_allclose(estimate.boundary_squared, [np.sqrt(2) / 4], rtol=1e-14)
    # Robin + prescribed flux is additive, just as in solve_heat.
    estimate = heat_indicator(
        mesh,
        mesh.points[:, 2],
        material,
        dirichlet_facets=[1, 2, 3],
        outward_fluxes={0: 1.0},
        convection={0: Convection(1.0, -1.0)},
    )
    assert estimate.total < 1e-14
    with pytest.raises(ValueError):
        estimate.cell[0] = 1


def test_interior_jump_is_counted_and_split_equally():
    mesh = Mesh(
        points=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
            ]
        ),
        tetra4=np.array([[0, 1, 2, 3], [0, 2, 1, 4]]),
        triangle3=np.array(
            [[0, 1, 3], [0, 2, 3], [1, 2, 3], [0, 1, 4], [0, 2, 4], [1, 2, 4]]
        ),
    )
    estimate = heat_indicator(
        mesh,
        abs(mesh.points[:, 2]),
        ThermalMaterial.isotropic(1),
        dirichlet_facets=np.arange(6),
    )
    # Shared z=0 area=1/2, diameter=sqrt(2), normal jump=2.
    np.testing.assert_allclose(estimate.interior_squared, [np.sqrt(2)] * 2, rtol=1e-14)


def test_bulk_marking_minimal_stable_and_scale_invariant():
    np.testing.assert_array_equal(mark_bulk([1, 2, 2, 0], fraction=0.5), [1, 2])
    np.testing.assert_array_equal(mark_bulk([1, 2, 2, 0], fraction=1), [0, 1, 2])
    assert mark_bulk([0, 0]).size == 0
    np.testing.assert_array_equal(
        mark_bulk([1e200, 2e200, 2e200], fraction=0.5), [1, 2]
    )
    for values, fraction in [([-1], 0.5), ([np.nan], 0.5), ([1], 0), ([1], 1.1)]:
        with pytest.raises(ValidationError):
            mark_bulk(values, fraction=fraction)


def test_invalid_input_is_rejected(cube_mesh):
    mesh = cube_mesh(1)
    material = ThermalMaterial.isotropic(1)
    for kwargs in [
        dict(dirichlet_facets=[0], outward_fluxes={0: 1}),
        dict(dirichlet_facets=[999]),
        dict(source=np.nan),
        dict(convection={0: 3}),
        dict(outward_fluxes={999: 0}),
    ]:
        with pytest.raises(ValidationError):
            heat_indicator(mesh, np.zeros(mesh.n_points), material, **kwargs)
    with pytest.raises(ValidationError):
        heat_indicator(mesh, np.zeros(1), material)
    with pytest.raises(ValidationError):
        heat_indicator(mesh, np.zeros(mesh.n_points), [])
    incomplete = Mesh(points=mesh.points, tetra4=mesh.tetra4)
    with pytest.raises(ValidationError):
        heat_indicator(incomplete, np.zeros(mesh.n_points), material)


def test_distorted_affine_patch_and_cell_order(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    material = ThermalMaterial([[2, 0.2, 0.1], [0.2, 3, 0.1], [0.1, 0.1, 4]])
    temperature = mesh.points @ np.array([1.0, 2.0, -3.0])
    estimate = heat_indicator(
        mesh, temperature, material, dirichlet_facets=np.arange(mesh.n_triangle3)
    )
    assert estimate.total < 2e-13
    # External IDs and local orientation must not become array indices.
    reversed_mesh = Mesh(
        points=mesh.points,
        tetra4=mesh.tetra4[::-1, ::-1],
        tetra4_ids=np.arange(mesh.n_tetra4) + 1000,
        triangle3=mesh.triangle3[:, ::-1],
    )
    original = heat_indicator(
        mesh, temperature**2, material, dirichlet_facets=np.arange(mesh.n_triangle3)
    )
    reversed_result = heat_indicator(
        reversed_mesh,
        temperature**2,
        material,
        dirichlet_facets=np.arange(mesh.n_triangle3),
    )
    np.testing.assert_allclose(
        reversed_result.cell[::-1], original.cell, rtol=1e-13, atol=1e-13
    )


def test_material_interface_solver_patch(cube_mesh):
    from fem_toolbox.boundary import Dirichlet
    from fem_toolbox.solvers.heat import solve_heat, thermal_dofs

    mesh = cube_mesh(4)
    x = mesh.points[:, 0]
    material = [
        ThermalMaterial.isotropic(1 if p < 0.5 else 8)
        for p in mesh.points[mesh.tetra4].mean(axis=1)[:, 0]
    ]
    nodes = np.flatnonzero((x == 0) | (x == 1))
    result = solve_heat(
        mesh,
        material,
        [
            Dirichlet(
                thermal_dofs(mesh).indices("T", nodes),
                np.where(x[nodes] == 0, 0.0, 0.5625),
            )
        ],
    )
    dirichlet = np.flatnonzero(
        np.all(x[mesh.triangle3] == 0, axis=1) | np.all(x[mesh.triangle3] == 1, axis=1)
    )
    estimate = heat_indicator(
        mesh, result.temperature, material, dirichlet_facets=dirichlet
    )
    np.testing.assert_allclose(
        result.temperature,
        np.where(x <= 0.5, x, 0.5 + (x - 0.5) / 8),
        rtol=0,
        atol=2e-14,
    )
    assert estimate.total < 1e-12
    assert abs(result.power_balance) < 1e-12


def test_quadratic_mesh_is_explicitly_rejected():
    from fem_toolbox.mesh.quadratic import TET10_EDGES

    corners = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )
    points = np.vstack(
        (corners, [(corners[a] + corners[b]) / 2 for a, b in TET10_EDGES])
    )
    mesh = Mesh(points=points, tetra10=np.arange(10)[None, :])
    with pytest.raises(ValidationError, match="Tet4"):
        heat_indicator(mesh, np.zeros(10), ThermalMaterial.isotropic(1))
