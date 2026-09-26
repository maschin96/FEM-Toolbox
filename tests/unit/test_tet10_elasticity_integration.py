"""Independent quadratic patches and pure bending through the production solver."""

import json
from dataclasses import replace

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ValidationError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion
from fem_toolbox.model import Model
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
from tests.reference.quadratic_mesh import elevate


def mechanical_document():
    return {
        "schema_version": 1,
        "units": "SI",
        "physics": "elasticity",
        "materials": {"u": {"solid": {"young": 100.0, "poisson": 0.0, "unit": "Pa"}}},
        "regions": {"u": {"body": "solid"}},
        "loads": {"body_force": {"value": [-2, 0, 0], "unit": "N/m^3"}},
        "boundary": {
            "dirichlet": [
                {"field": "u", "set": "left", "value": [0, 0, 0], "unit": "m"}
            ],
            "traction": [{"set": "right", "value": [2, 0, 0], "unit": "N/m^2"}],
        },
    }


def test_quadratic_displacement_model_and_integral_export(cube_mesh, tmp_path):
    mesh = elevate(cube_mesh(2))
    model = Model(mesh=mesh, document=mechanical_document())
    solved = model.solve()
    result = solved.field_result
    expected = np.zeros_like(mesh.points)
    expected[:, 0] = 0.01 * mesh.points[:, 0] ** 2
    assert_allclose(result.displacement, expected, rtol=0, atol=2e-14)
    centers = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    strain = np.zeros((mesh.n_cells, 6))
    strain[:, 0] = 0.02 * centers[:, 0]
    assert_allclose(result.strain, strain, rtol=0, atol=2e-13)
    assert_allclose(result.stress, 100 * strain, rtol=0, atol=2e-11)
    assert_allclose(result.field_points, centers, rtol=0, atol=1e-15)
    assert result.energy == pytest.approx(0.02 / 3, rel=0, abs=2e-15)
    assert_allclose(result.force_balance, 0, rtol=0, atol=2e-12)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-12)
    replay = Model(mesh=mesh, document=model.to_document()).solve().field_result
    assert_allclose(replay.displacement, result.displacement, rtol=0, atol=0)
    vtu, manifest = export_result_bundle(tmp_path / "solid", mesh, solved)
    data = json.loads(manifest.read_text())
    assert data["bundle_schema_version"] == 2
    fields = data["fields"]["tetra10"]
    assert fields["stress_Pa_at_centroid"]["components"] == [
        "xx",
        "yy",
        "zz",
        "yz",
        "xz",
        "xy",
    ]
    assert fields["strain_1_at_centroid"]["location"] == "cell_centroid"
    assert fields["mechanical_energy_J"]["kind"] == "cell_integral"
    assert_allclose(
        meshio.read(vtu).cell_data_dict["stress_Pa_at_centroid"]["tetra10"],
        result.stress,
    )
    with pytest.raises(ValueError):
        result.field_points[0, 0] = 1
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(tmp_path / "bad", mesh, replace(result, field_points=None))


def test_pure_bending_quadratic_exactness_and_balance(cube_mesh):
    mesh = elevate(cube_mesh(1))
    curvature, young = 0.01, 100.0
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(young, 0),
        [Dirichlet(mechanical_dofs(mesh).indices("u", left), np.zeros(3 * len(left)))],
        tractions={
            int(f): lambda p: [-young * curvature * (p[1] - 0.5), 0, 0]
            for f in mesh.facet_sets["right"]
        },
    )
    expected = np.zeros_like(mesh.points)
    expected[:, 0] = -curvature * mesh.points[:, 0] * (mesh.points[:, 1] - 0.5)
    expected[:, 1] = curvature * mesh.points[:, 0] ** 2 / 2
    assert_allclose(result.displacement, expected, rtol=0, atol=2e-14)
    assert result.energy == pytest.approx(young * curvature**2 / 24, rel=0, abs=2e-16)
    assert_allclose(result.force_balance, 0, rtol=0, atol=2e-13)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-13)
    moment = np.cross(mesh.points, result.applied_loads).sum(axis=0)
    assert_allclose(moment, [0, 0, young * curvature / 12], rtol=0, atol=2e-14)


def test_tet10_rigid_modes_and_explicit_thermal_boundary(cube_mesh):
    mesh = elevate(cube_mesh(1))
    material = ElasticMaterial.isotropic(100, 0.25)
    with pytest.raises(ValidationError, match="rigid"):
        solve_elasticity(mesh, material, [])
    raw = mechanical_document()
    raw["boundary"].pop("dirichlet")
    with pytest.raises(ValidationError, match="rigid"):
        Model(mesh=mesh, document=raw)
    with pytest.raises(ValidationError, match="Tet10 thermal expansion"):
        solve_elasticity(
            mesh,
            material,
            [],
            temperature=np.full(mesh.n_points, 300),
            expansion=ThermalExpansion([1e-5, 1e-5, 1e-5, 0, 0, 0], 293),
        )


@pytest.mark.parametrize("reverse", [False, True])
def test_anisotropic_affine_patch_and_orientation(cube_mesh, reverse):
    from fem_toolbox.mesh import Mesh

    original = elevate(cube_mesh(2))
    mesh = Mesh(
        points=np.vstack([original.points, [2, 2, 2]]),
        tetra10=original.cells[:, [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]]
        if reverse
        else original.cells,
        triangle6=original.facets[:, [0, 2, 1, 5, 4, 3]]
        if reverse
        else original.facets,
    )
    a = np.arange(36).reshape(6, 6) / 10
    tensor = a.T @ a + 100 * np.eye(6)
    h = 0.001 * np.array([[1, 2, 3], [-4, 5, 6], [7, -8, 9]])
    expected = mesh.points @ h.T + [0.3, -0.2, 0.1]
    boundary = np.unique(mesh.facets)
    result = solve_elasticity(
        mesh,
        ElasticMaterial(tensor),
        [
            Dirichlet(
                mechanical_dofs(mesh).indices("u", boundary), expected[boundary].ravel()
            )
        ],
    )
    strain = 0.001 * np.array([1, 5, 9, -2, 10, -2])
    assert len(result.node_ids) == mesh.n_points - 1
    assert_allclose(result.displacement, expected[result.node_ids], rtol=0, atol=2e-14)
    assert_allclose(
        result.strain, np.tile(strain, (mesh.n_cells, 1)), rtol=0, atol=2e-13
    )
    assert_allclose(
        result.stress, np.tile(tensor @ strain, (mesh.n_cells, 1)), rtol=0, atol=1e-10
    )
    assert result.energy == pytest.approx(
        0.5 * strain @ tensor @ strain, rel=0, abs=2e-14
    )


def test_tet10_point_load_and_curvature_validation(cube_mesh):
    from fem_toolbox.mesh import Mesh

    mesh = elevate(cube_mesh(1))
    material = ElasticMaterial.isotropic(100, 0.25)
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    bc = [Dirichlet(mechanical_dofs(mesh).indices("u", left), np.zeros(3 * len(left)))]
    tip = int(np.flatnonzero(np.all(mesh.points == [1, 1, 1], axis=1))[0])
    result = solve_elasticity(mesh, material, bc, point_loads={tip: [0, -1, 0]})
    assert_allclose(result.force_balance, 0, rtol=0, atol=2e-12)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-12)
    assert result.energy > 0
    assert_allclose(
        result.energy, -0.5 * result.displacement[tip, 1], rtol=0, atol=2e-14
    )
    points = mesh.points.copy()
    points[mesh.cells[0, 4], 2] += 0.01
    with pytest.raises(ValidationError):
        solve_elasticity(Mesh(points=points, tetra10=mesh.cells), material, [])


def test_tet10_rigid_motion_has_no_strain_energy(cube_mesh):
    mesh = elevate(cube_mesh(2))
    exact = np.cross([0.01, -0.02, 0.03], mesh.points) + np.array([0.1, -0.2, 0.3])
    boundary = np.unique(mesh.facets)
    result = solve_elasticity(
        mesh,
        ElasticMaterial.isotropic(100, 0.25),
        [
            Dirichlet(
                mechanical_dofs(mesh).indices("u", boundary), exact[boundary].ravel()
            )
        ],
    )
    assert_allclose(result.displacement, exact, rtol=0, atol=2e-14)
    assert_allclose(result.stress, 0, rtol=0, atol=1e-10)
    assert result.energy < 1e-24


def test_tet10_mechanics_gmsh_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_elasticity/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "example")], check=True
    )
    data = json.loads((tmp_path / "example" / "solid.json").read_text())
    assert data["physics"] == "elasticity"
    assert data["bundle_schema_version"] == 2


def test_tet10_two_material_model_traction_continuity(cube_mesh):
    from fem_toolbox.mesh import Mesh

    original = elevate(cube_mesh(2))
    centers = original.points[original.cells[:, :4]].mean(axis=1)
    mesh = Mesh(
        points=original.points,
        tetra10=original.cells,
        triangle6=original.facets,
        facet_sets=original.facet_sets,
        volume_sets={
            "soft": np.flatnonzero(centers[:, 0] < 0.5),
            "stiff": np.flatnonzero(centers[:, 0] > 0.5),
        },
    )
    raw = mechanical_document()
    raw["materials"]["u"]["stiff"] = {"young": 800, "poisson": 0, "unit": "Pa"}
    raw["regions"]["u"] = {"soft": "solid", "stiff": "stiff"}
    raw.pop("loads")
    raw["boundary"]["traction"][0]["value"] = [1, 0, 0]
    result = Model(mesh=mesh, document=raw).solve().field_result
    expected = np.zeros_like(mesh.points)
    x = mesh.points[:, 0]
    expected[:, 0] = np.where(x <= 0.5, x / 100, 0.005 + (x - 0.5) / 800)
    assert_allclose(result.displacement, expected, rtol=0, atol=2e-14)
    assert_allclose(
        result.stress,
        np.tile([1, 0, 0, 0, 0, 0], (mesh.n_cells, 1)),
        rtol=0,
        atol=2e-11,
    )
    assert result.energy == pytest.approx(0.0028125, rel=0, abs=2e-15)
