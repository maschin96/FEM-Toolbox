"""Independent harmonic patches and reference-geometry mesh quality failures."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain
from fem_toolbox.cfd.mesh_motion import MeshMotionError, move_fluid_domain


def channel(mesh, *, periodic=False, closed=False):
    bc = {"bottom": {"type": "wall"}, "top": {"type": "wall", "fsi": True}}
    if periodic:
        bc.update(
            left={
                "type": "periodic",
                "partner": "right",
                "translation_m": [1.0, 0.0],
                "pressure_drop_Pa": 0.0,
            },
            right={
                "type": "periodic",
                "partner": "left",
                "translation_m": [-1.0, 0.0],
                "pressure_drop_Pa": 0.0,
            },
        )
    elif closed:
        bc.update(left={"type": "wall"}, right={"type": "wall"})
    else:
        bc.update(
            left={"type": "inlet", "velocity_m_s": [1.0, 0.0], "nu_tilde_m2_s": 0.1},
            right={"type": "pressure_outlet", "pressure_Pa": 0.0},
        )
    document = {
        "fluid_domain_schema_version": 1,
        "mesh": "fluid.json",
        "width_m": 0.7,
        "fluid": {"density_kg_m3": 1.0, "viscosity_Pa_s": 1.0},
        "boundaries": bc,
    }
    if periodic or closed:
        document["pressure_reference_cell"] = 0
    return FluidDomain.from_document(document, mesh=mesh)


def top_displacement(domain, dy=0.001, dx=0.0):
    return {"top": np.tile([dx, dy], (len(domain.mesh.boundary_sets["top"]), 2, 1))}


@pytest.mark.parametrize("shear", [0.0, 0.2])
def test_harmonic_affine_patch_absolute_coordinates_and_fixed_flow(quad_mesh, shear):
    domain = channel(quad_mesh(5, 6, shear=shear))
    imposed = top_displacement(domain, dy=0.001, dx=0.0004)
    result = move_fluid_domain(domain, imposed, length_scale_m=1.0)
    expected = domain.mesh.points[:, 1, None] * [0.0004, 0.001]
    np.testing.assert_allclose(result.displacement_m, expected, atol=1e-13)
    np.testing.assert_allclose(
        result.domain.mesh.points, domain.mesh.points + expected, atol=1e-13
    )
    assert result.minimum_area_ratio == pytest.approx(1.001)
    assert result.domain.inlet_mass_flow_kg_s == pytest.approx(
        domain.inlet_mass_flow_kg_s, rel=1e-13
    )
    np.testing.assert_array_equal(
        result.domain.boundaries["left"]["nu_tilde_m2_s"],
        domain.boundaries["left"]["nu_tilde_m2_s"],
    )
    np.testing.assert_array_equal(result.domain.mesh.quad4, domain.mesh.quad4)
    np.testing.assert_array_equal(result.domain.mesh.line2, domain.mesh.line2)
    assert result.maximum_relative_laplace_residual < 1e-12
    assert result.reference_mesh_sha256 != result.mesh_sha256
    assert len(result.mesh_sha256) == 64
    repeat = move_fluid_domain(domain, imposed, length_scale_m=1.0)
    np.testing.assert_array_equal(result.domain.mesh.points, repeat.domain.mesh.points)
    assert repeat.mesh_sha256 == result.mesh_sha256
    with pytest.raises(ValueError):
        result.displacement_m.setflags(write=True)


def test_periodic_displacements_remain_congruent(quad_mesh):
    domain = channel(quad_mesh(4, 6), periodic=True)
    result = move_fluid_domain(domain, top_displacement(domain), length_scale_m=1)
    np.testing.assert_allclose(
        result.displacement_m[:, 1], 0.001 * domain.mesh.points[:, 1], atol=1e-13
    )
    assert len(result.domain.periodic_partners) == len(domain.periodic_partners)


def test_inconsistent_periodic_corners_are_rejected(quad_mesh):
    domain = channel(quad_mesh(4, 4), periodic=True)
    imposed = top_displacement(domain)
    x = domain.mesh.points[domain.mesh.line2[domain.mesh.boundary_sets["top"]], 0]
    imposed["top"][:, :, 1] = 0.001 * x
    with pytest.raises(MeshMotionError, match="periodic"):
        move_fluid_domain(domain, imposed, length_scale_m=1)


@pytest.mark.parametrize(
    "dy,stage", [(0.03, "displacement"), (-0.002, "geometry"), (-0.0006, "area")]
)
def test_displacement_inversion_and_minimum_area_are_separate_failures(
    quad_mesh, dy, stage
):
    domain = channel(quad_mesh(4, 4, size_m=(1, 0.001)))
    with pytest.raises(MeshMotionError) as caught:
        move_fluid_domain(domain, top_displacement(domain, dy), length_scale_m=1)
    assert caught.value.stage == stage
    assert caught.value.displacement_m.shape == domain.mesh.points.shape
    assert len(caught.value.reference_mesh_sha256) == 64


def test_fixed_wall_corner_conflict_and_missing_interface_are_errors(quad_mesh):
    domain = channel(quad_mesh(4, 4), closed=True)
    with pytest.raises(MeshMotionError, match="corner"):
        move_fluid_domain(domain, top_displacement(domain), length_scale_m=1)
    with pytest.raises(MeshMotionError, match="interface"):
        move_fluid_domain(domain, {}, length_scale_m=1)


def test_resource_budget_preserves_reference_identity(quad_mesh):
    domain = channel(quad_mesh(4, 4))
    with pytest.raises(MeshMotionError, match="budget"):
        move_fluid_domain(
            domain, top_displacement(domain), length_scale_m=1, max_nodes=4
        )


def test_warped_patch_retains_volume_groups_and_recomputes_wall_distance(quad_mesh):
    from fem_toolbox.cfd import FluidMesh2D

    original = quad_mesh(5, 6, warp=0.02)
    mesh = FluidMesh2D(
        original.points,
        original.quad4,
        original.line2,
        original.boundary_sets,
        volume_sets={"upstream": np.flatnonzero(original.cell_centres_m[:, 0] < 0.5)},
    )
    domain = channel(mesh)
    result = move_fluid_domain(domain, top_displacement(domain), length_scale_m=1)
    np.testing.assert_allclose(
        result.displacement_m[:, 1], 0.001 * mesh.points[:, 1], atol=1e-13
    )
    np.testing.assert_array_equal(
        result.domain.mesh.volume_sets["upstream"], mesh.volume_sets["upstream"]
    )
    np.testing.assert_allclose(
        result.domain.wall_quality.wall_distance_m,
        1.001 * domain.wall_quality.wall_distance_m,
        atol=1e-13,
    )


def test_nonmatching_structural_trace_moves_fluid_and_rebuilds_current_overlap(
    cube_mesh, quad_mesh
):
    from fem_toolbox.coupling import build_domain_interface
    from fem_toolbox.mesh import Mesh

    domain = channel(quad_mesh(5, 6))
    original = cube_mesh(2)
    points = original.points * [1, 1, 0.7] + [0, 1, 0]
    rows = np.flatnonzero(np.all(points[original.triangle3, 1] == 1, axis=1))
    solid = Mesh(
        points=points,
        tetra4=original.tetra4,
        triangle3=original.triangle3,
        facet_sets={"wet": rows},
    )
    mapping = build_domain_interface(
        domain, wall_name="top", structure=solid, facet_set="wet", length_scale_m=1
    )
    u = np.zeros_like(points)
    u[:, 1] = 0.001
    endpoints = mapping.segment_endpoint_displacements(u)
    moved = move_fluid_domain(domain, {"top": endpoints}, length_scale_m=1)
    current_solid = Mesh(
        points=points + u,
        tetra4=solid.tetra4,
        triangle3=solid.triangle3,
        facet_sets=solid.facet_sets,
    )
    current = build_domain_interface(
        moved.domain,
        wall_name="top",
        structure=current_solid,
        facet_set="wet",
        length_scale_m=1,
    )
    assert current.maximum_relative_coverage_error < 1e-10
    np.testing.assert_allclose(current.quadrature_points_m[:, 1], 1.001, atol=1e-13)


def test_wall_layer_quality_can_fail_even_with_positive_sufficient_cell_areas(
    quad_mesh,
):
    domain = channel(quad_mesh(8, 8))
    imposed = top_displacement(domain, dy=0.0)
    x = domain.mesh.points[domain.mesh.line2[domain.mesh.boundary_sets["top"]], 0]
    imposed["top"][:, :, 1] = -0.18 * np.sin(np.pi * x)
    with pytest.raises(MeshMotionError) as caught:
        move_fluid_domain(domain, imposed, length_scale_m=10)
    assert caught.value.stage == "quality"
    assert "growth" in str(caught.value)


def test_settings_identity_and_explicit_unscaled_inlet_profile(quad_mesh):
    domain = channel(quad_mesh(4, 4))
    imposed = top_displacement(domain)
    fixed = move_fluid_domain(domain, imposed, length_scale_m=1)
    unscaled = move_fluid_domain(
        domain, imposed, length_scale_m=1, preserve_inlet_flow=False
    )
    assert fixed.mesh_sha256 == unscaled.mesh_sha256
    assert fixed.input_sha256 != unscaled.input_sha256
    assert unscaled.domain.inlet_mass_flow_kg_s == pytest.approx(
        1.001 * domain.inlet_mass_flow_kg_s
    )
    assert unscaled.settings["preserve_inlet_flow"] is False
    with pytest.raises(TypeError):
        unscaled.settings["preserve_inlet_flow"] = True
