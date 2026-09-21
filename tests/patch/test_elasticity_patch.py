import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ElasticMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity


def test_affine_patch_on_distorted_mesh(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    c = np.diag([11.0, 12, 13, 4, 5, 6])
    c[0, 1] = c[1, 0] = 2
    material = ElasticMaterial(c)
    a = np.array([[0.1, 0.2, 0.3], [-0.1, 0.15, 0.2], [0.12, 0.13, -0.1]])
    expected_u = mesh.points @ a.T + [0.2, -0.1, 0.3]
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    dofs = mechanical_dofs(mesh)
    result = solve_elasticity(
        mesh,
        material,
        [Dirichlet(dofs.indices("u", boundary), expected_u[boundary].ravel())],
    )
    expected_e = np.array([0.1, 0.15, -0.1, 0.33, 0.42, 0.1])
    assert_allclose(result.displacement, expected_u, rtol=0, atol=3e-14)
    assert_allclose(
        result.strain, np.tile(expected_e, (mesh.n_tetra4, 1)), rtol=0, atol=3e-14
    )
    assert_allclose(
        result.stress, np.tile(c @ expected_e, (mesh.n_tetra4, 1)), rtol=0, atol=5e-13
    )
    assert_allclose(result.force_balance, 0, atol=2e-13)
    assert_allclose(result.moment_balance, 0, atol=2e-13)
    assert_allclose(result.energy, 0.5 * expected_e @ c @ expected_e, atol=2e-14)
    assert_allclose(
        2 * result.energy,
        np.sum(result.displacement * (result.applied_loads + result.reactions)),
        atol=2e-13,
    )
    with pytest.raises(ValueError):
        result.displacement[0, 0] = 1


def test_tension_bar_and_point_load_superposition(cube_mesh):
    mesh = cube_mesh(2)
    dofs = mechanical_dofs(mesh)
    conditions = [
        Dirichlet(
            dofs.indices("u", np.flatnonzero(mesh.points[:, axis] == 0), [component]),
            np.zeros(np.count_nonzero(mesh.points[:, axis] == 0)),
        )
        for axis, component in enumerate(("x", "y", "z"))
    ]
    right = np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == 1, axis=1))
    tractions = {int(face): [10.0, 0, 0] for face in right}
    result = solve_elasticity(
        mesh, ElasticMaterial.isotropic(1000, 0.25), conditions, tractions=tractions
    )
    assert_allclose(
        result.displacement, mesh.points * [0.01, -0.0025, -0.0025], atol=2e-16
    )
    assert_allclose(
        result.stress, np.tile([10, 0, 0, 0, 0, 0], (mesh.n_tetra4, 1)), atol=2e-12
    )
    assert_allclose(result.reactions.sum(axis=0), [-10, 0, 0], atol=2e-12)
    assert_allclose(result.energy, 0.05, atol=2e-15)
    assert_allclose(result.force_balance, 0, atol=2e-12)
    assert_allclose(result.moment_balance, 0, atol=2e-12)
    point_loads = {
        int(node): f
        for node, f in zip(result.node_ids, result.applied_loads, strict=True)
        if np.any(f)
    }
    point_result = solve_elasticity(
        mesh, ElasticMaterial.isotropic(1000, 0.25), conditions, point_loads=point_loads
    )
    assert_allclose(point_result.displacement, result.displacement, atol=2e-16)


def test_missing_support_and_disconnected_body(cube_mesh):
    mesh = cube_mesh(1)
    material = ElasticMaterial.isotropic(1000, 0.25)
    with pytest.raises(ValidationError, match="rigid"):
        solve_elasticity(mesh, material, [])
    with pytest.raises(ValidationError, match="rigid"):
        solve_elasticity(mesh, material, [Dirichlet([0, 1, 2], [0, 0, 0])])
    split = Mesh(
        points=np.vstack([mesh.points, mesh.points + np.array([2, 0, 0])]),
        tetra4=np.vstack([mesh.tetra4, mesh.tetra4 + mesh.n_points]),
    )
    with pytest.raises(ValidationError, match="rigid"):
        solve_elasticity(
            split,
            material,
            [Dirichlet(np.arange(3 * mesh.n_points), np.zeros(3 * mesh.n_points))],
        )


def test_unused_nodes_and_invalid_load_selection(cube_mesh):
    original = cube_mesh(1)
    mesh = Mesh(
        points=np.vstack([original.points, [5, 5, 5]]),
        tetra4=original.tetra4,
        triangle3=original.triangle3,
    )
    dofs = mechanical_dofs(mesh)
    assert dofs.size == 3 * original.n_points
    bc = [Dirichlet(np.arange(dofs.size), np.zeros(dofs.size))]
    material = ElasticMaterial.isotropic(1000, 0.25)
    result = solve_elasticity(mesh, material, bc, body_force=[1, 2, 3])
    assert_allclose(result.reactions.sum(axis=0), [-1, -2, -3], atol=2e-15)
    for kwargs in [{"point_loads": {8: [1, 0, 0]}}, {"tractions": {999: [1, 0, 0]}}]:
        with pytest.raises(ValidationError):
            solve_elasticity(mesh, material, bc, **kwargs)
    with pytest.raises(ValidationError):
        solve_elasticity(mesh, [material], bc)


def test_prescribed_rigid_motion_has_no_stress(cube_mesh):
    mesh = cube_mesh(2, distort=True)
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    expected = np.cross([0.02, 0.01, -0.03], mesh.points) + np.array([0.1, -0.2, 0.3])
    dofs = mechanical_dofs(mesh)
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(1000, 0.25),
        [Dirichlet(dofs.indices("u", boundary), expected[boundary].ravel())],
    )
    assert_allclose(result.displacement, expected, atol=2e-15, rtol=0)
    assert_allclose(result.stress, 0, atol=2e-12, rtol=0)
    assert result.energy < 1e-27


def test_nonuniform_materials_and_zero_solution(cube_mesh):
    mesh = cube_mesh(2)
    dofs = mechanical_dofs(mesh)
    nodes = np.flatnonzero(mesh.points[:, 0] == 0)
    materials = [
        ElasticMaterial.isotropic(100 * (1 + i), 0.25) for i in range(mesh.n_tetra4)
    ]
    result = solve_elasticity(
        mesh, materials, [Dirichlet(dofs.indices("u", nodes), np.zeros(3 * len(nodes)))]
    )
    assert_allclose(result.displacement, 0, atol=0)
    assert result.energy == 0


def test_internal_triangle_cannot_receive_boundary_traction(cube_mesh):
    original = cube_mesh(1)
    # Tet 0 and 1 share the (0,4,7) face of the unit cube.
    mesh = Mesh(points=original.points, tetra4=original.tetra4, triangle3=[[0, 4, 7]])
    with pytest.raises(ValidationError, match="exterior"):
        solve_elasticity(
            mesh, ElasticMaterial.isotropic(100, 0.25), [], tractions={0: [1, 0, 0]}
        )


def test_point_connected_mechanism_is_rejected():
    mesh = Mesh(
        points=[
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [-1, 0, 0],
            [0, -1, 0],
            [0, 0, -1],
        ],
        tetra4=[[0, 1, 2, 3], [0, 4, 5, 6]],
    )
    with pytest.raises(ValidationError, match="singular"):
        solve_elasticity(
            mesh,
            ElasticMaterial.isotropic(100, 0.25),
            [Dirichlet(np.arange(12), np.zeros(12))],
        )


def test_deterministic_run_provenance(cube_mesh):
    mesh = cube_mesh(1)
    material = ElasticMaterial.isotropic(100, 0.25)
    bc = [Dirichlet(np.arange(24), np.zeros(24))]
    first = solve_elasticity(mesh, material, bc)
    again = solve_elasticity(mesh, material, bc)
    changed = solve_elasticity(mesh, material, bc, body_force=[1, 0, 0])
    assert first.provenance.to_dict() == again.provenance.to_dict()
    assert first.provenance.input_sha256 != changed.provenance.input_sha256
