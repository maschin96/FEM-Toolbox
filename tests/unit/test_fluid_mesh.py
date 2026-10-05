"""Independent geometry, topology, wall and SI preflight expectations."""

import json

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain, FluidMesh2D
from fem_toolbox.errors import ValidationError


def channel_document():
    return {
        "fluid_mesh_schema_version": 1,
        "length_unit": "m",
        "points": [
            [0, 0],
            [1, 0],
            [2, 0],
            [0, 1],
            [1, 1],
            [2, 1],
            [0, 2],
            [1, 2],
            [2, 2],
        ],
        "quad4": [[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]],
        "line2": [[0, 1], [1, 2], [2, 5], [5, 8], [8, 7], [7, 6], [6, 3], [3, 0]],
        "boundary_sets": {
            "bottom": [0, 1],
            "outlet": [2, 3],
            "top": [4, 5],
            "inlet": [6, 7],
        },
        "volume_sets": {"fluid": [0, 1, 2, 3]},
    }


def problem_document():
    return {
        "fluid_domain_schema_version": 1,
        "mesh": "mesh.json",
        "width_m": 1.0,
        "fluid": {"density_kg_m3": 1.0, "viscosity_Pa_s": 1.0},
        "boundaries": {
            "bottom": {"type": "wall"},
            "top": {"type": "wall", "fsi": True},
            "inlet": {
                "type": "inlet",
                "velocity_m_s": [1.0, 0.0],
                "nu_tilde_m2_s": 3.0,
            },
            "outlet": {"type": "pressure_outlet", "pressure_Pa": 0.0},
        },
    }


def test_channel_geometry_groups_and_owned_arrays():
    raw = channel_document()
    mesh = FluidMesh2D.from_document(raw)
    np.testing.assert_allclose(mesh.cell_areas_m2, 1, rtol=0, atol=1e-14)
    np.testing.assert_allclose(
        mesh.cell_centres_m,
        [[0.5, 0.5], [1.5, 0.5], [0.5, 1.5], [1.5, 1.5]],
        rtol=0,
        atol=1e-14,
    )
    assert mesh.volume_sets["fluid"].tolist() == [0, 1, 2, 3]
    assert mesh.boundary_sets["inlet"].tolist() == [6, 7]
    np.testing.assert_allclose(mesh.boundary_normals[0], [0, -1], atol=1e-14)
    assert mesh.face_neighbours.tolist().count(-1) == 8
    assert mesh.max_nonorthogonality_deg == pytest.approx(0, abs=1e-12)
    assert mesh.max_skewness == pytest.approx(0, abs=1e-12)
    raw["points"][0][0] = 10
    assert mesh.points[0, 0] == 0
    with pytest.raises(ValueError):
        mesh.points.setflags(write=True)
    with pytest.raises(TypeError):
        mesh.boundary_sets["new"] = np.array([0])
    assert (
        FluidMesh2D.from_document(mesh.to_document()).to_document()
        == mesh.to_document()
    )


@pytest.mark.parametrize(
    "damage",
    [
        "clockwise",
        "degenerate",
        "concave",
        "missing_boundary",
        "inner_boundary",
        "duplicate_cell",
        "overlap",
        "float_index",
        "group_index",
        "unit",
        "nan",
        "unknown",
    ],
)
def test_invalid_mesh_diagnosed(damage):
    d = channel_document()
    if damage == "clockwise":
        d["quad4"][0].reverse()
    elif damage == "degenerate":
        d["points"][4] = d["points"][1]
    elif damage == "concave":
        d["points"][4] = [0.2, 0.2]
    elif damage == "missing_boundary":
        d["line2"].pop()
    elif damage == "inner_boundary":
        d["line2"][0] = [1, 4]
    elif damage == "duplicate_cell":
        d["quad4"].append(d["quad4"][0])
    elif damage == "overlap":
        d["points"] += [[0.2, 0.2], [0.8, 0.2], [0.8, 0.8], [0.2, 0.8]]
        d["quad4"].append([9, 10, 11, 12])
        d["line2"] += [[9, 10], [10, 11], [11, 12], [12, 9]]
    elif damage == "float_index":
        d["quad4"][0][0] = 0.5
    elif damage == "group_index":
        d["boundary_sets"]["inlet"] = [100]
    elif damage == "unit":
        d["length_unit"] = "mm"
    elif damage == "nan":
        d["points"][0][0] = float("nan")
    else:
        d["unexpected"] = 0
    with pytest.raises(ValidationError):
        FluidMesh2D.from_document(d)


def test_wall_distance_is_to_segment_not_vertex_and_y_plus():
    mesh = FluidMesh2D.from_document(channel_document())
    wall = mesh.wall_quality(["bottom", "top"])
    np.testing.assert_allclose(wall.wall_distance_m, 0.5, rtol=0, atol=1e-14)
    np.testing.assert_allclose(wall.first_cell_distance_m, 0.5, rtol=0, atol=1e-14)
    np.testing.assert_allclose(wall.layer_growth, 1, rtol=0, atol=1e-14)
    assert wall.boundary_indices.tolist() == [0, 1, 4, 5]
    np.testing.assert_allclose(
        wall.y_plus(density_kg_m3=1, viscosity_Pa_s=1, wall_shear_Pa=[4] * 4),
        1,
        rtol=0,
        atol=1e-14,
    )
    with pytest.raises(ValidationError, match="y_plus"):
        wall.y_plus(density_kg_m3=1, viscosity_Pa_s=1, wall_shear_Pa=[9] * 4)
    with pytest.raises(ValidationError):
        mesh.wall_quality([])


def test_nonuniform_layers_growth_and_anisotropy():
    d = channel_document()
    for point in d["points"]:
        point[0] *= 100
        point[1] *= 0.01
    mesh = FluidMesh2D.from_document(d)
    np.testing.assert_allclose(
        mesh.wall_quality(["bottom", "top"]).wall_distance_m, 0.005, rtol=0, atol=1e-14
    )
    d = channel_document()
    for point in d["points"]:
        if point[1] == 2:
            point[1] = 3
    mesh = FluidMesh2D.from_document(d)
    with pytest.raises(ValidationError, match="growth"):
        mesh.wall_quality(["bottom", "top"])


def test_file_import_and_domain_preflight(tmp_path):
    (tmp_path / "mesh.json").write_text(json.dumps(channel_document()))
    (tmp_path / "domain.json").write_text(json.dumps(problem_document()))
    domain = FluidDomain.from_file(tmp_path / "domain.json")
    assert domain.inlet_mass_flow_kg_s == pytest.approx(2, abs=1e-14)
    assert domain.wall_quality.boundary_indices.tolist() == [0, 1, 4, 5]
    assert domain.boundaries["top"]["fsi"] is True
    with pytest.raises(TypeError):
        domain.boundaries["top"]["fsi"] = False


@pytest.mark.parametrize(
    "damage",
    [
        "missing_wall",
        "duplicate_bc",
        "negative_density",
        "missing_sa",
        "outward_inlet",
        "gauge_and_outlet",
        "unknown_bc",
        "invalid_fsi",
    ],
)
def test_domain_rejects_invalid_bc_and_si(tmp_path, damage):
    mesh = FluidMesh2D.from_document(channel_document())
    d = problem_document()
    if damage == "missing_wall":
        del d["boundaries"]["bottom"]
    elif damage == "duplicate_bc":
        raw = channel_document()
        raw["boundary_sets"]["other"] = [0]
        mesh = FluidMesh2D.from_document(raw)
        d["boundaries"]["other"] = {"type": "wall"}
    elif damage == "negative_density":
        d["fluid"]["density_kg_m3"] = -1
    elif damage == "missing_sa":
        del d["boundaries"]["inlet"]["nu_tilde_m2_s"]
    elif damage == "outward_inlet":
        d["boundaries"]["inlet"]["velocity_m_s"] = [-1, 0]
    elif damage == "gauge_and_outlet":
        d["pressure_reference_cell"] = 0
    elif damage == "unknown_bc":
        d["boundaries"]["bottom"]["type"] = "magic"
    else:
        d["boundaries"]["inlet"]["fsi"] = True
    with pytest.raises(ValidationError):
        FluidDomain.from_document(d, mesh=mesh)


def test_explicit_width_required_and_mass_scales_with_width():
    d = problem_document()
    mesh = FluidMesh2D.from_document(channel_document())
    d["width_m"] = 4
    assert FluidDomain.from_document(
        d, mesh=mesh
    ).inlet_mass_flow_kg_s == pytest.approx(8, rel=0, abs=1e-14)
    del d["width_m"]
    with pytest.raises(ValidationError, match="width_m"):
        FluidDomain.from_document(d, mesh=mesh)


def test_partial_face_contacts_are_nonconforming():
    d = {
        "fluid_mesh_schema_version": 1,
        "length_unit": "m",
        "points": [
            [0, 0],
            [1, 0],
            [2, 0],
            [0, 1],
            [1, 1],
            [2, 1],
            [0, 2],
            [1, 2],
            [2, 2],
            [1.5, 1],
            [1.5, 2],
        ],
        "quad4": [
            [0, 1, 4, 3],
            [1, 2, 5, 4],
            [3, 4, 7, 6],
            [4, 9, 10, 7],
            [9, 5, 8, 10],
        ],
        "line2": [
            [0, 1],
            [1, 2],
            [2, 5],
            [5, 8],
            [8, 10],
            [10, 7],
            [7, 6],
            [6, 3],
            [3, 0],
            [4, 5],
            [4, 9],
            [9, 5],
        ],
        "boundary_sets": {"all": list(range(12))},
    }
    with pytest.raises(ValidationError, match="conform"):
        FluidMesh2D.from_document(d)


@pytest.mark.parametrize("field", ["points", "quad4"])
def test_mixed_booleans_are_not_numeric_inputs(field):
    d = channel_document()
    d[field][0][0] = False
    with pytest.raises(ValidationError):
        FluidMesh2D.from_document(d)


def periodic_document():
    d = problem_document()
    d["pressure_reference_cell"] = 0
    d["boundaries"]["inlet"] = {
        "type": "periodic",
        "partner": "outlet",
        "translation_m": [2, 0],
        "pressure_drop_Pa": 1,
    }
    d["boundaries"]["outlet"] = {
        "type": "periodic",
        "partner": "inlet",
        "translation_m": [-2, 0],
        "pressure_drop_Pa": -1,
    }
    return d


def test_periodic_pairs_geometry_gauge_and_symmetry():
    mesh = FluidMesh2D.from_document(channel_document())
    d = periodic_document()
    domain = FluidDomain.from_document(d, mesh=mesh)
    assert dict(domain.periodic_partners) == {6: 3, 7: 2, 2: 7, 3: 6}
    assert domain.inlet_mass_flow_kg_s == 0
    d["boundaries"]["top"] = {"type": "symmetry"}
    domain = FluidDomain.from_document(d, mesh=mesh)
    assert domain.wall_quality.boundary_indices.tolist() == [0, 1]


@pytest.mark.parametrize(
    "damage", ["translation", "pressure", "reciprocal", "gauge", "congruence"]
)
def test_periodic_mismatch_rejected(damage):
    d = periodic_document()
    raw = channel_document()
    if damage == "translation":
        d["boundaries"]["outlet"]["translation_m"] = [-3, 0]
    elif damage == "pressure":
        d["boundaries"]["outlet"]["pressure_drop_Pa"] = 1
    elif damage == "reciprocal":
        d["boundaries"]["outlet"]["partner"] = "top"
    elif damage == "gauge":
        del d["pressure_reference_cell"]
    else:
        raw["points"][5][0] = 2.1
    with pytest.raises(ValidationError):
        FluidDomain.from_document(d, mesh=FluidMesh2D.from_document(raw))


def test_sheared_cells_have_analytic_normals_angles_and_distances():
    d = channel_document()
    for p in d["points"]:
        p[0] += 0.3 * p[1]
    mesh = FluidMesh2D.from_document(d)
    np.testing.assert_allclose(mesh.cell_areas_m2, 1, rtol=0, atol=1e-14)
    assert mesh.max_nonorthogonality_deg == pytest.approx(
        np.degrees(np.arctan(0.3)), rel=0, abs=1e-12
    )
    assert mesh.max_skewness < 1e-14
    np.testing.assert_allclose(
        mesh.wall_quality(["bottom", "top"]).wall_distance_m, 0.5, rtol=0, atol=1e-14
    )
    for p in d["points"]:
        p[0] += 2 * p[1]
    with pytest.raises(ValidationError, match="quality"):
        FluidMesh2D.from_document(d)


def test_closest_wall_can_be_segment_endpoint():
    d = channel_document()
    d["boundary_sets"]["short"] = [0]
    mesh = FluidMesh2D.from_document(d)
    np.testing.assert_allclose(
        mesh.wall_quality(["short"]).wall_distance_m,
        [0.5, np.sqrt(0.5), 1.5, np.sqrt(2.5)],
        rtol=0,
        atol=1e-14,
    )


@pytest.mark.parametrize(
    "damage",
    [
        "schema",
        "shape",
        "duplicate_face",
        "nonmanifold",
        "empty_group",
        "negative_index",
        "infinite",
    ],
)
def test_additional_import_errors(damage):
    d = channel_document()
    if damage == "schema":
        d["fluid_mesh_schema_version"] = True
    elif damage == "shape":
        d["points"][0] = [0, 0, 0]
    elif damage == "duplicate_face":
        d["line2"].append(d["line2"][0])
    elif damage == "nonmanifold":
        d["quad4"].append([1, 2, 8, 7])
    elif damage == "empty_group":
        d["boundary_sets"]["bottom"] = []
    elif damage == "negative_index":
        d["quad4"][0][0] = -1
    else:
        d["points"][0][0] = float("inf")
    with pytest.raises(ValidationError):
        FluidMesh2D.from_document(d)


@pytest.mark.parametrize(
    "damage",
    [
        "width",
        "viscosity",
        "velocity_shape",
        "sa_negative",
        "fsi_type",
        "no_wall",
        "unknown",
        "mesh_escape",
    ],
)
def test_additional_domain_errors(damage):
    d = problem_document()
    if damage == "width":
        d["width_m"] = 0
    elif damage == "viscosity":
        d["fluid"]["viscosity_Pa_s"] = float("nan")
    elif damage == "velocity_shape":
        d["boundaries"]["inlet"]["velocity_m_s"] = [[1, 0]]
    elif damage == "sa_negative":
        d["boundaries"]["inlet"]["nu_tilde_m2_s"] = [-1, 0]
    elif damage == "fsi_type":
        d["boundaries"]["top"]["fsi"] = 1
    elif damage == "no_wall":
        d["boundaries"]["top"] = {"type": "symmetry"}
        d["boundaries"]["bottom"] = {"type": "symmetry"}
    elif damage == "unknown":
        d["unexpected"] = 1
    else:
        d["mesh"] = "../mesh.json"
    with pytest.raises(ValidationError):
        FluidDomain.from_document(d, mesh=FluidMesh2D.from_document(channel_document()))


def test_unreadable_json_and_escaping_symlink(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("{")
    with pytest.raises(ValidationError):
        FluidMesh2D.from_file(path)
    with pytest.raises(ValidationError):
        FluidDomain.from_file(path)
    with pytest.raises(ValidationError):
        FluidMesh2D.from_file(tmp_path / "missing.json")
    inner = tmp_path / "input"
    inner.mkdir()
    (tmp_path / "outside.json").write_text(json.dumps(channel_document()))
    (inner / "mesh.json").symlink_to(tmp_path / "outside.json")
    (inner / "domain.json").write_text(json.dumps(problem_document()))
    with pytest.raises(ValidationError, match="symlink"):
        FluidDomain.from_file(inner / "domain.json")


def test_anisotropic_aspect_ratio_is_reported_without_false_rejection():
    d = channel_document()
    for p in d["points"]:
        p[0] *= 100
        p[1] *= 0.01
    mesh = FluidMesh2D.from_document(d)
    np.testing.assert_allclose(mesh.cell_aspect_ratio, 10000, rtol=0, atol=1e-10)


def test_oriented_faces_close_each_control_volume():
    mesh = FluidMesh2D.from_document(channel_document())
    for cell, area in enumerate(mesh.cell_areas_m2):
        signs = (mesh.face_owners == cell).astype(float) - (
            mesh.face_neighbours == cell
        ).astype(float)
        vector = signs[:, None] * mesh.face_normals * mesh.face_lengths_m[:, None]
        np.testing.assert_allclose(vector.sum(axis=0), [0, 0], rtol=0, atol=1e-14)
        boundary_area = np.sum(mesh.face_centres_m * vector) / 2
        assert boundary_area == pytest.approx(area, rel=0, abs=1e-14)


def test_domain_detaches_inlet_arrays_and_rejects_wrong_mesh_type():
    d = problem_document()
    domain = FluidDomain.from_document(
        d, mesh=FluidMesh2D.from_document(channel_document())
    )
    d["boundaries"]["inlet"]["velocity_m_s"][0] = 99
    assert domain.boundaries["inlet"]["velocity_m_s"][0, 0] == 1
    with pytest.raises(ValueError):
        domain.boundaries["inlet"]["velocity_m_s"].setflags(write=True)
    with pytest.raises(ValidationError):
        FluidDomain.from_document(d, mesh=None)
