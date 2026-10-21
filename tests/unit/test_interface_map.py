"""Analytic overlap integrals, adjoint work and explicit unsupported geometry."""

import numpy as np
import pytest

from fem_toolbox.coupling import build_interface_map
from fem_toolbox.errors import ValidationError
from fem_toolbox.mesh import Mesh


def interface(cube_mesh, *, width=4.0, x_end=3.0, normal=(0, 1), shift=0):
    base = cube_mesh(2)
    points = base.points * np.array([3.0, 1.0, 4.0]) + [0, 1, 0]
    points[:, 0] += shift
    faces = np.flatnonzero(np.all(np.isclose(points[base.triangle3, 1], 1), axis=1))
    solid = Mesh(
        points=points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        facet_sets={"interface": faces},
    )
    segments = np.array([[[0, 1], [0.7, 1]], [[0.7, 1], [x_end, 1]]], float)
    return build_interface_map(
        segments_m=segments,
        fluid_normals=np.tile(normal, (2, 1)),
        width_m=width,
        structure=solid,
        facet_set="interface",
        length_scale_m=3,
    ), solid


def test_constant_pressure_and_shear_force_all_moments_and_gauge(cube_mesh):
    mapping, solid = interface(cube_mesh)
    # n_f points upward. -sigma*n_f: drag +3 x and pressure +2 y.
    loads = mapping.loads_from_fluid_traction(
        np.tile([-3.0, -7.0], (2, 1)), external_pressure_Pa=5
    )
    np.testing.assert_allclose(
        loads.nodal_forces_N.sum(axis=0), [36, 24, 0], atol=1e-12
    )
    np.testing.assert_allclose(
        np.cross(solid.points, loads.nodal_forces_N).sum(axis=0),
        [-48, 72, 0],
        atol=1e-11,
    )
    np.testing.assert_allclose(loads.fluid_force_N, [36, 84, 0], atol=1e-12)
    np.testing.assert_allclose(loads.external_pressure_force_N, [0, -60, 0], atol=1e-12)
    shifted = mapping.loads_from_fluid_traction(
        np.tile([-3.0, -10.0], (2, 1)), external_pressure_Pa=8
    )
    np.testing.assert_allclose(shifted.nodal_forces_N, loads.nodal_forces_N, atol=1e-12)
    assert loads.audit.maximum_relative_error < 1e-10
    with pytest.raises(ValueError):
        loads.nodal_forces_N.setflags(write=True)


def test_affine_traction_exact_integrals_and_arbitrary_virtual_work(cube_mesh):
    mapping, solid = interface(cube_mesh)
    q = mapping.quadrature_points_m
    traction = np.column_stack([q[:, 0], 2 + q[:, 2], np.zeros(len(q))])
    force = traction * mapping.quadrature_areas_m2[:, None]
    nodal = mapping.transfer_forces(force)
    # Independent rectangular integrals: int x=18; int(2+z)=48.
    np.testing.assert_allclose(nodal.sum(axis=0), [18, 48, 0], atol=1e-11)
    # Mx=-int z(2+z)=-112, My=int zx=36, Mz=int[x(2+z)-y*x]=54.
    np.testing.assert_allclose(
        np.cross(solid.points, nodal).sum(axis=0), [-112, 36, 54], atol=1e-10
    )
    rng = np.random.default_rng(89)
    virtual = rng.normal(size=solid.points.shape)
    assert np.sum(virtual * nodal) == pytest.approx(
        np.sum(mapping.interpolate(virtual) * force), abs=1e-12
    )
    assert mapping.audit(force, virtual).maximum_relative_error < 1e-10
    affine = solid.points @ np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).T + [1, 2, 3]
    np.testing.assert_allclose(
        mapping.interpolate(affine),
        q @ np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).T + [1, 2, 3],
        atol=1e-12,
    )
    matrix = mapping.interpolation_matrix
    matrix.data[:] = 0
    np.testing.assert_allclose(
        mapping.interpolate(np.ones_like(solid.points)), 1, atol=1e-14
    )


@pytest.mark.parametrize(
    "change",
    [
        {"width": 3.9},
        {"width": 4.1},
        {"x_end": 2.9},
        {"normal": (0, -1)},
        {"shift": 0.01},
    ],
)
def test_gaps_wrong_width_and_orientation_are_errors(cube_mesh, change):
    with pytest.raises(ValidationError):
        interface(cube_mesh, **change)


def test_planar_trace_and_span_dependence_diagnosis(cube_mesh):
    mapping, solid = interface(cube_mesh)
    u = np.zeros_like(solid.points)
    u[:, 1] = 0.001 + 0.002 * solid.points[:, 0]
    endpoints = mapping.segment_endpoint_displacements(u)
    np.testing.assert_allclose(
        endpoints[:, :, 1], [[0.001, 0.0024], [0.0024, 0.007]], atol=1e-14
    )
    u[:, 1] += solid.points[:, 2] * 1e-5
    with pytest.raises(ValidationError, match="span"):
        mapping.segment_endpoint_displacements(u)
    u[:] = 0
    u[:, 2] = 1e-3
    with pytest.raises(ValidationError, match="plane"):
        mapping.segment_endpoint_displacements(u)


def test_balanced_holes_and_double_fluid_coverage_cannot_cancel(cube_mesh):
    base = cube_mesh(1)
    points = base.points * [3, 1, 4] + [0, 1, 0]
    selected = np.flatnonzero(np.all(points[base.triangle3, 1] == 1, axis=1))
    solid = Mesh(
        points=points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        facet_sets={"interface": selected},
    )
    segments = np.array(
        [
            [[0, 1], [0.75, 1]],
            [[0, 1], [0.75, 1]],
            [[2.25, 1], [3, 1]],
            [[2.25, 1], [3, 1]],
        ]
    )
    with pytest.raises(ValidationError, match="double"):
        build_interface_map(
            segments_m=segments,
            fluid_normals=np.tile([0, 1], (4, 1)),
            width_m=4,
            structure=solid,
            facet_set="interface",
            length_scale_m=3,
        )


@pytest.mark.parametrize("use_fluid_solver", [False, True])
def test_one_way_pressure_uses_existing_solid_solver_and_balances_reactions(
    cube_mesh, quad_mesh, use_fluid_solver
):
    from fem_toolbox.coupling import solve_one_way
    from fem_toolbox.model import Model

    mapping, original = interface(cube_mesh)
    solid = Mesh(
        points=original.points,
        tetra4=original.tetra4,
        triangle3=original.triangle3,
        facet_sets=original.facet_sets,
        volume_sets={"solid": np.arange(original.n_cells)},
    )
    nodes = np.arange(len(solid.points)).tolist()
    top = np.flatnonzero(solid.points[:, 1] == 2).tolist()
    doc = {
        "schema_version": 1,
        "units": "SI",
        "physics": "elasticity",
        "mesh": {"path": "solid.msh", "length_scale_to_m": 1.0},
        "materials": {"u": {"solid": {"young": 1000.0, "poisson": 0.0, "unit": "Pa"}}},
        "regions": {"u": {"solid": "solid"}},
        "boundary": {
            "dirichlet": [
                {
                    "field": "u",
                    "nodes": nodes,
                    "components": ["x", "z"],
                    "value": 0.0,
                    "unit": "m",
                },
                {
                    "field": "u",
                    "nodes": top,
                    "components": ["y"],
                    "value": 0.0,
                    "unit": "m",
                },
            ]
        },
    }
    model = Model(mesh=solid, document=doc)
    traction = np.tile([0.0, -2.0], (2, 1))
    external = 0.0
    if use_fluid_solver:
        from fem_toolbox.cfd import FluidDomain, solve_quad_flow
        from fem_toolbox.coupling import build_domain_interface

        fluid = quad_mesh(4, 4, size_m=(3, 1))
        domain = FluidDomain.from_document(
            {
                "fluid_domain_schema_version": 1,
                "mesh": "fluid.json",
                "width_m": 4.0,
                "fluid": {"density_kg_m3": 1.0, "viscosity_Pa_s": 1.0},
                "boundaries": {
                    name: {"type": "wall", "fsi": name == "top"}
                    for name in fluid.boundary_sets
                },
                "pressure_reference_cell": 0,
            },
            mesh=fluid,
        )
        flow = solve_quad_flow(
            domain, velocity_scale_m_s=1, length_scale_m=3, pressure_reference_Pa=7
        )
        mapping = build_domain_interface(
            domain,
            wall_name="top",
            structure=solid,
            facet_set="interface",
            length_scale_m=3,
        )
        lookup = {
            int(row): i for i, row in enumerate(domain.wall_quality.boundary_indices)
        }
        traction = flow.wall_traction_Pa[
            [lookup[int(row)] for row in fluid.boundary_sets["top"]]
        ]
        external = 5.0
    result = solve_one_way(mapping, model, traction, external_pressure_Pa=external)
    expected = np.zeros_like(solid.points)
    expected[:, 1] = 0.002 * (2 - solid.points[:, 1])
    np.testing.assert_allclose(result.structure.displacement, expected, atol=1e-13)
    np.testing.assert_allclose(
        result.structure.reactions.sum(axis=0), [0, -24, 0], atol=1e-11
    )
    np.testing.assert_allclose(result.structure.force_balance, 0, atol=1e-11)
    np.testing.assert_allclose(result.structure.moment_balance, 0, atol=1e-11)
    assert result.structure.relative_residual < 1e-12
    assert result.maximum_displacement_ratio < 0.02
    assert result.maximum_displacement_gradient < 0.02


def test_quadrature_resource_budget_is_an_explicit_error(cube_mesh):
    mapping, solid = interface(cube_mesh)
    with pytest.raises(ValidationError, match="budget"):
        build_interface_map(
            segments_m=mapping.segments_m,
            fluid_normals=mapping.fluid_normals,
            width_m=4,
            structure=solid,
            facet_set="interface",
            length_scale_m=3,
            max_quadrature_points=2,
        )


@pytest.mark.parametrize("bad", [True, 0, -1, np.nan, np.inf, "4"])
def test_invalid_physical_width_does_not_enter_geometry(cube_mesh, bad):
    with pytest.raises(ValidationError):
        interface(cube_mesh, width=bad)


def test_overlapping_separate_solid_bodies_are_double_coverage(cube_mesh):
    mapping, solid = interface(cube_mesh)
    count = len(solid.points)
    points = np.vstack([solid.points, solid.points])
    cells = np.vstack([solid.tetra4, solid.tetra4 + count])
    facets = np.vstack([solid.triangle3, solid.triangle3 + count])
    selected = np.r_[
        solid.facet_sets["interface"],
        solid.facet_sets["interface"] + len(solid.triangle3),
    ]
    doubled = Mesh(
        points=points,
        tetra4=cells,
        triangle3=facets,
        facet_sets={"interface": selected},
    )
    with pytest.raises(ValidationError, match="double"):
        build_interface_map(
            segments_m=mapping.segments_m,
            fluid_normals=mapping.fluid_normals,
            width_m=4,
            structure=doubled,
            facet_set="interface",
            length_scale_m=3,
        )


def test_rigid_rotation_and_reversed_fluid_endpoints_preserve_physical_load(cube_mesh):
    base, solid = interface(cube_mesh)
    angle = 0.37
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
    )
    rotation3 = np.eye(3)
    rotation3[:2, :2] = rotation
    points = solid.points @ rotation3.T
    rotated = Mesh(
        points=points,
        tetra4=solid.tetra4,
        triangle3=solid.triangle3,
        facet_sets=solid.facet_sets,
    )
    mapping = build_interface_map(
        segments_m=base.segments_m[:, ::-1] @ rotation.T,
        fluid_normals=base.fluid_normals @ rotation.T,
        width_m=4,
        structure=rotated,
        facet_set="interface",
        length_scale_m=3,
    )
    original = base.loads_from_fluid_traction(np.tile([-3.0, -2.0], (2, 1)))
    transformed = mapping.loads_from_fluid_traction(
        np.tile([-3.0, -2.0], (2, 1)) @ rotation.T
    )
    np.testing.assert_allclose(
        transformed.nodal_forces_N, original.nodal_forces_N @ rotation3.T, atol=1e-12
    )
