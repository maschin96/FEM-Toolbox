"""Analytical quadratic electric fields, integrated energy and Model export."""

import json
from dataclasses import replace

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ElectricGaugeError, ValidationError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.model import Model
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics
from tests.reference.quadratic_mesh import elevate


def electric_document():
    return {
        "schema_version": 1,
        "units": "SI",
        "physics": "electrostatics",
        "materials": {"phi": {"dielectric": {"permittivity": 2e-11, "unit": "F/m"}}},
        "regions": {"phi": {"body": "dielectric"}},
        "loads": {"volume_charge": {"value": -4e-11, "unit": "C/m^3"}},
        "boundary": {
            "dirichlet": [{"field": "phi", "set": "left", "value": 0, "unit": "V"}],
            "surface_charge": [{"set": "right", "value": -4e-11, "unit": "C/m^2"}],
        },
    }


def test_quadratic_electric_patch_integrates_energy(cube_mesh):
    mesh = elevate(cube_mesh(2))
    tensor = 1e-11 * np.array([[2, 0.3, 0.2], [0.3, 3, 0.1], [0.2, 0.1, 4]])
    coefficients = np.array([1, 0.5, 0.25])
    exact = mesh.points**2 @ coefficients
    boundary = np.unique(mesh.facets)
    result = solve_electrostatics(
        mesh,
        DielectricMaterial(tensor),
        [Dirichlet(electric_dofs(mesh).indices("phi", boundary), exact[boundary])],
        volume_charge=-9e-11,
    )
    centers = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    expected_e = -2 * centers * coefficients
    assert_allclose(result.potential, exact, rtol=0, atol=2e-13)
    assert_allclose(result.electric_field, expected_e, rtol=0, atol=2e-12)
    assert_allclose(
        result.electric_displacement, expected_e @ tensor, rtol=0, atol=1e-22
    )
    assert_allclose(result.field_points, centers, rtol=0, atol=1e-15)
    # Exact integral of 1/2 grad(phi).kappa.grad(phi) over the unit cube:
    # int x_i²=1/3, int x_i*x_j=1/4 for i!=j.
    expected_energy = 2 / 3 * np.dot(np.diag(tensor), coefficients**2)
    expected_energy += sum(
        tensor[i, j] * coefficients[i] * coefficients[j]
        for i in range(3)
        for j in range(i + 1, 3)
    )
    assert result.energy == pytest.approx(expected_energy, rel=0, abs=2e-24)
    assert abs(result.charge_balance) < 2e-23
    assert result.volume_charge == pytest.approx(-9e-11, rel=0, abs=2e-24)
    assert_allclose(result.electrode_charges, -result.reactions, rtol=0, atol=0)
    with pytest.raises(ValueError):
        result.field_points[0, 0] = 0


def test_triangle6_model_and_export_distinguish_energy_from_samples(
    cube_mesh, tmp_path
):
    mesh = elevate(cube_mesh(1))
    model = Model(mesh=mesh, document=electric_document())
    solved = model.solve()
    result = solved.field_result
    assert_allclose(result.potential, mesh.points[:, 0] ** 2, rtol=0, atol=2e-13)
    assert result.energy == pytest.approx(4e-11 / 3, rel=0, abs=2e-24)
    assert result.prescribed_surface_charge == pytest.approx(-4e-11, rel=0, abs=2e-24)
    assert abs(result.dirichlet_outward_charge) < 2e-23
    assert abs(result.charge_balance) < 2e-23
    replay = Model(mesh=mesh, document=model.to_document()).solve().field_result
    assert_allclose(replay.potential, result.potential, rtol=0, atol=0)
    vtu, manifest = export_result_bundle(tmp_path / "electric", mesh, solved)
    data = json.loads(manifest.read_text())
    assert data["bundle_schema_version"] == 2
    fields = data["fields"]["tetra10"]
    assert fields["electric_field_V_m_at_centroid"]["location"] == "cell_centroid"
    assert fields["electric_energy_J"]["kind"] == "cell_integral"
    assert fields["electric_energy_J"]["location"] == "cell_volume"
    exported = meshio.read(vtu)
    assert_allclose(
        exported.cell_data_dict["electric_energy_J"]["tetra10"], result.element_energy
    )
    assert "electric_energy_J_at_centroid" not in exported.cell_data_dict
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(tmp_path / "bad", mesh, replace(result, field_points=None))


def test_tet10_electric_gauge_and_material_interface(cube_mesh):
    mesh = elevate(cube_mesh(2))
    with pytest.raises(ElectricGaugeError):
        solve_electrostatics(mesh, DielectricMaterial.isotropic(2e-11), [])
    raw = electric_document()
    raw["boundary"].pop("dirichlet")
    with pytest.raises(ValidationError, match=r"reference|gauge"):
        Model(mesh=mesh, document=raw)
    centers = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    materials = [
        DielectricMaterial.isotropic(1e-11 if x < 0.5 else 8e-11) for x in centers[:, 0]
    ]
    x = mesh.points[:, 0]
    exact = np.where(x <= 0.5, x, 0.5 + (x - 0.5) / 8)
    fixed = np.flatnonzero((x == 0) | (x == 1))
    result = solve_electrostatics(mesh, materials, [Dirichlet(fixed, exact[fixed])])
    assert_allclose(result.potential, exact, rtol=0, atol=3e-13)
    assert_allclose(
        result.electric_displacement,
        np.tile([-1e-11, 0, 0], (mesh.n_cells, 1)),
        rtol=0,
        atol=2e-23,
    )
    assert result.energy == pytest.approx(2.8125e-12, rel=0, abs=2e-24)
    assert abs(result.charge_balance) < 2e-23


@pytest.mark.parametrize("reverse", [False, True])
def test_tet10_capacitor_charge_signs_and_orientation(cube_mesh, reverse):
    from fem_toolbox.mesh import Mesh

    original = elevate(cube_mesh(1))
    mesh = Mesh(
        points=np.vstack([original.points, [2, 2, 2]]),  # inactive node
        tetra10=original.cells[:, [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]]
        if reverse
        else original.cells,
        triangle6=original.facets[:, [0, 2, 1, 5, 4, 3]]
        if reverse
        else original.facets,
        facet_sets=original.facet_sets,
    )
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    result = solve_electrostatics(
        mesh,
        DielectricMaterial.isotropic(2e-11),
        [Dirichlet(electric_dofs(mesh).indices("phi", left), np.zeros(len(left)))],
        surface_charges={int(f): -2e-11 for f in mesh.facet_sets["right"]},
    )
    assert len(result.potential) == mesh.n_points - 1
    assert_allclose(
        result.potential, mesh.points[result.node_ids, 0], rtol=0, atol=2e-13
    )
    assert result.energy == pytest.approx(1e-11, rel=0, abs=2e-24)
    assert result.dirichlet_outward_charge == pytest.approx(2e-11, rel=0, abs=2e-24)
    assert result.electrode_charges.sum() == pytest.approx(-2e-11, rel=0, abs=2e-24)


def test_tet10_electric_rejects_curvature_and_bad_export_points(cube_mesh, tmp_path):
    from fem_toolbox.mesh import Mesh

    mesh = elevate(cube_mesh(1))
    points = mesh.points.copy()
    points[mesh.cells[0, 4], 2] += 0.01
    curved = Mesh(points=points, tetra10=mesh.cells, triangle6=mesh.facets)
    with pytest.raises(ValidationError):
        solve_electrostatics(curved, DielectricMaterial.isotropic(2e-11), [])
    result = Model(mesh=mesh, document=electric_document()).solve().field_result
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(
            tmp_path / "wrong",
            mesh,
            replace(result, field_points=result.field_points + 1),
        )


def test_tet10_electric_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_electrostatics/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "example")], check=True
    )
    data = json.loads((tmp_path / "example" / "electric.json").read_text())
    assert data["physics"] == "electrostatics"
    assert data["bundle_schema_version"] == 2
