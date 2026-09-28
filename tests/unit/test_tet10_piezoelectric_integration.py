"""Independent coupled polynomial and conducting-boundary benchmarks."""

import json
from dataclasses import replace

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.errors import ValidationError
from fem_toolbox.io import export_result_bundle
from fem_toolbox.materials import ElasticMaterial, PiezoelectricMaterial
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from tests.reference.quadratic_mesh import elevate


def material():
    e = np.zeros((3, 6))
    e[0, 0] = 2
    return PiezoelectricMaterial(
        ElasticMaterial.isotropic(100, 0).stiffness, e, 3 * np.eye(3)
    )


def quadratic_case(mesh, backend=None):
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    dofs = piezoelectric_dofs(mesh)
    return solve_piezoelectric(
        mesh,
        material(),
        [
            Dirichlet(dofs.indices("u", left), np.zeros(3 * len(left))),
            Dirichlet(dofs.indices("phi", left), np.zeros(len(left))),
        ],
        body_force=[-6.0, 0, 0],
        volume_charge=-5.96,
        tractions={int(f): [6.0, 0, 0] for f in mesh.facet_sets["right"]},
        surface_charges={int(f): -5.96 for f in mesh.facet_sets["right"]},
        backend=backend,
    )


def test_quadratic_coupled_fields_energy_balance_and_export(cube_mesh, tmp_path):
    mesh = elevate(cube_mesh(2))
    result = quadratic_case(mesh)
    u = np.zeros_like(mesh.points)
    u[:, 0] = 0.01 * mesh.points[:, 0] ** 2
    assert_allclose(result.displacement, u, rtol=0, atol=2e-13)
    assert_allclose(result.potential, mesh.points[:, 0] ** 2, rtol=0, atol=2e-12)
    centers = mesh.points[mesh.cells[:, :4]].mean(axis=1)
    assert_allclose(result.field_points, centers, rtol=0, atol=1e-15)
    stress = np.zeros((mesh.n_cells, 6))
    stress[:, 0] = 6.0 * centers[:, 0]
    assert_allclose(result.stress, stress, rtol=0, atol=2e-11)
    electric = np.zeros((mesh.n_cells, 3))
    electric[:, 0] = -2 * centers[:, 0]
    assert_allclose(result.electric_field, electric, rtol=0, atol=2e-11)
    assert_allclose(result.electric_displacement, 2.98 * electric, rtol=0, atol=2e-11)
    assert_allclose(
        [
            result.mechanical_energy.sum(),
            result.coupling_energy.sum(),
            result.dielectric_enthalpy.sum(),
        ],
        [0.02 / 3, 0.08 / 3, -2],
        rtol=0,
        atol=2e-12,
    )
    assert_allclose(result.force_balance, 0, rtol=0, atol=2e-11)
    assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-11)
    assert abs(result.charge_balance) < 2e-11
    vtu, report = export_result_bundle(tmp_path / "piezo", mesh, result)
    loaded = meshio.read(vtu)
    assert_allclose(
        loaded.cell_data_dict["stress_Pa_at_centroid"]["tetra10"], result.stress
    )
    assert_allclose(
        loaded.cell_data_dict["coupling_energy_J"]["tetra10"], result.coupling_energy
    )
    data = json.loads(report.read_text())
    assert data["bundle_schema_version"] == 2
    fields = data["fields"]["tetra10"]
    for name in ("mechanical_energy_J", "coupling_energy_J", "dielectric_enthalpy_J"):
        assert fields[name]["kind"] == "cell_integral"
    assert fields["stress_Pa_at_centroid"]["location"] == "cell_centroid"
    with pytest.raises(ValidationError, match="centroid"):
        export_result_bundle(tmp_path / "bad", mesh, replace(result, field_points=None))


@pytest.mark.parametrize("floating", [False, True])
def test_fixed_and_floating_electrodes_include_midside_nodes(cube_mesh, floating):
    mesh = elevate(cube_mesh(2))
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
            Electrode("top", "right", charge=3.04)
            if floating
            else Electrode("top", "right", potential=1),
        ],
    )
    expected = np.zeros_like(mesh.points)
    expected[:, 0] = -0.02 * mesh.points[:, 0]
    assert_allclose(result.displacement, expected, rtol=0, atol=2e-13)
    assert_allclose(result.potential, mesh.points[:, 0], rtol=0, atol=2e-12)
    assert_allclose(result.stress, 0, rtol=0, atol=2e-11)
    assert result.electrode_states[1].charge == pytest.approx(3.04, abs=2e-11)
    assert result.total_enthalpy == pytest.approx(-1.52, abs=2e-12)


def document():
    mat = material()
    return {
        "schema_version": 1,
        "units": "SI",
        "physics": "piezoelectric",
        "materials": {
            "piezo": {
                "p": {
                    "stiffness": mat.stiffness.tolist(),
                    "stress_coupling": mat.stress_coupling.tolist(),
                    "permittivity": mat.permittivity.tolist(),
                    "units": {
                        "stiffness": "Pa",
                        "stress_coupling": "C/m^2",
                        "permittivity": "F/m",
                    },
                }
            }
        },
        "regions": {"piezo": {"body": "p"}},
        "boundary": {
            "dirichlet": [
                {"field": "u", "set": "left", "value": [0, 0, 0], "unit": "m"}
            ],
            "electrodes": [
                {
                    "name": "ground",
                    "set": "left",
                    "potential": 0,
                    "potential_unit": "V",
                },
                {"name": "top", "set": "right", "charge": 3.04, "charge_unit": "C"},
            ],
        },
    }


def test_model_roundtrip_and_preflight(cube_mesh, tmp_path):
    from fem_toolbox.model import Model

    mesh = elevate(cube_mesh(2))
    model = Model(mesh=mesh, document=document())
    result = model.solve()
    replay = Model(mesh=mesh, document=model.to_document()).solve()
    assert_allclose(
        result.field_result.potential, mesh.points[:, 0], rtol=0, atol=2e-12
    )
    assert_allclose(
        result.field_result.displacement,
        replay.field_result.displacement,
        rtol=0,
        atol=0,
    )
    export_result_bundle(tmp_path / "model", mesh, result)
    missing = document()
    missing["boundary"]["dirichlet"] = []
    with pytest.raises(ValidationError, match="rigid"):
        Model(mesh=mesh, document=missing)
    missing = document()
    missing["boundary"]["electrodes"].pop(0)
    with pytest.raises(ValidationError, match="reference"):
        Model(mesh=mesh, document=missing)


def test_gmsh_toml_example(tmp_path):
    import subprocess
    import sys
    from pathlib import Path

    script = Path(__file__).parents[2] / "examples/tet10_piezoelectric/run.py"
    subprocess.run(
        [sys.executable, str(script), "--output", str(tmp_path / "example")], check=True
    )
    report = json.loads((tmp_path / "example/piezo.json").read_text())
    assert report["physics"] == "piezoelectric"
    assert report["bundle_schema_version"] == 2


@pytest.mark.parametrize("reverse", [False, True])
def test_anisotropic_affine_patch_with_unused_node(cube_mesh, reverse):
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
    c = a.T @ a + 100 * np.eye(6)
    e = np.arange(18).reshape(3, 6) / 10
    k = np.array([[3, 0.1, 0.2], [0.1, 4, 0.3], [0.2, 0.3, 5]])
    mat = PiezoelectricMaterial(c, e, k)
    h = 0.001 * np.array([[1, 2, 3], [-4, 5, 6], [7, -8, 9]])
    u = mesh.points @ h.T + [0.3, -0.2, 0.1]
    phi = mesh.points @ [0.5, -0.2, 0.3] + 4
    fixed = np.unique(mesh.facets)
    dofs = piezoelectric_dofs(mesh)
    result = solve_piezoelectric(
        mesh,
        [mat] * mesh.n_cells,
        [
            Dirichlet(dofs.indices("u", fixed), u[fixed].ravel()),
            Dirichlet(dofs.indices("phi", fixed), phi[fixed]),
        ],
    )
    strain = 0.001 * np.array([1, 5, 9, -2, 10, -2])
    field = np.array([-0.5, 0.2, -0.3])
    assert len(result.node_ids) == mesh.n_points - 1
    assert_allclose(result.displacement, u[result.node_ids], rtol=0, atol=2e-13)
    assert_allclose(result.potential, phi[result.node_ids], rtol=0, atol=2e-12)
    assert_allclose(
        result.stress,
        np.tile(c @ strain - e.T @ field, (mesh.n_cells, 1)),
        rtol=0,
        atol=2e-10,
    )
    assert_allclose(
        result.electric_displacement,
        np.tile(e @ strain + k @ field, (mesh.n_cells, 1)),
        rtol=0,
        atol=2e-11,
    )
    assert result.total_enthalpy == pytest.approx(
        0.5 * strain @ c @ strain - strain @ e.T @ field - 0.5 * field @ k @ field,
        abs=2e-12,
    )
    with pytest.raises(ValueError):
        result.field_points[0, 0] = 1


def test_point_load_reciprocity_and_invalid_inputs(cube_mesh):
    from fem_toolbox.mesh import Mesh

    mesh = elevate(cube_mesh(1))
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    right = np.unique(mesh.facets[mesh.facet_sets["right"]])
    dofs = piezoelectric_dofs(mesh)
    bc = [
        Dirichlet(dofs.indices("u", left), np.zeros(3 * len(left))),
        Dirichlet(dofs.indices("phi", left), np.zeros(len(left))),
    ]
    first, second = int(right[0]), int(right[-1])
    a = solve_piezoelectric(mesh, material(), bc, point_loads={first: [1, 0, 0]})
    b = solve_piezoelectric(mesh, material(), bc, point_loads={second: [1, 0, 0]})
    assert_allclose(
        a.displacement[second, 0], b.displacement[first, 0], rtol=0, atol=2e-13
    )
    assert_allclose(a.force_balance, 0, rtol=0, atol=2e-11)
    assert_allclose(a.moment_balance, 0, rtol=0, atol=2e-11)
    with pytest.raises(ValidationError, match="gauge"):
        solve_piezoelectric(mesh, material(), bc[:1])
    with pytest.raises(ValidationError, match="per cell"):
        solve_piezoelectric(mesh, [material()], bc)
    with pytest.raises(ValidationError, match="conflicts"):
        solve_piezoelectric(
            mesh,
            material(),
            bc,
            electrodes=[Electrode("r", "right", potential=1)],
            surface_charges={int(mesh.facet_sets["right"][0]): 1},
        )
    points = mesh.points.copy()
    points[mesh.cells[0, 4], 2] += 0.01
    with pytest.raises(ValidationError):
        solve_piezoelectric(Mesh(points=points, tetra10=mesh.cells), material(), [])


def test_two_material_model_interface_continuity(cube_mesh):
    from fem_toolbox.mesh import Mesh
    from fem_toolbox.model import Model

    original = elevate(cube_mesh(2))
    center = original.points[original.cells[:, :4]].mean(axis=1)
    mesh = Mesh(
        points=original.points,
        tetra10=original.cells,
        triangle6=original.facets,
        facet_sets=original.facet_sets,
        volume_sets={
            "a": np.flatnonzero(center[:, 0] < 0.5),
            "b": np.flatnonzero(center[:, 0] > 0.5),
        },
    )
    raw = document()
    first = raw["materials"]["piezo"]["p"]
    raw["materials"]["piezo"]["q"] = {
        **{
            key: (2 * np.array(first[key])).tolist()
            for key in ("stiffness", "stress_coupling", "permittivity")
        },
        "units": first["units"],
    }
    raw["regions"]["piezo"] = {"a": "p", "b": "q"}
    raw["solver"] = {"method": "direct", "block_scaling": False}
    result = Model(mesh=mesh, document=raw).solve().field_result
    x = mesh.points[:, 0]
    phi = np.where(x <= 0.5, x, 0.5 + 0.5 * (x - 0.5))
    assert_allclose(result.potential, phi, rtol=0, atol=2e-12)
    assert_allclose(result.displacement[:, 0], -0.02 * phi, rtol=0, atol=2e-13)
    assert_allclose(result.stress, 0, rtol=0, atol=2e-11)
    assert_allclose(
        result.electric_displacement,
        np.tile([-3.04, 0, 0], (mesh.n_cells, 1)),
        rtol=0,
        atol=2e-11,
    )
