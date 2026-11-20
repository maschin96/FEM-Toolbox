"""Stationary partitioned FSI: analytic fixed points, history and stop kinds (#91)."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain
from fem_toolbox.coupling import (
    FsiCouplingError,
    StationaryFsiProblem,
    solve_stationary_fsi,
)
from fem_toolbox.coupling.fsi import _diverging
from fem_toolbox.errors import ValidationError
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model


def structure(cube_mesh, *, length=3.0, width=4.0, young=1000.0):
    """Tet4 block y in [1,2] over the fluid lid; x/z guided, top y fixed."""
    base = cube_mesh(2)
    points = base.points * np.array([length, 1.0, width]) + [0, 1, 0]
    wet = np.flatnonzero(np.all(np.isclose(points[base.triangle3, 1], 1), axis=1))
    mesh = Mesh(
        points=points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        facet_sets={"wet": wet},
        volume_sets={"solid": np.arange(base.n_cells)},
    )
    top = np.flatnonzero(points[:, 1] == 2).tolist()
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "elasticity",
        "mesh": {"path": "solid.msh", "length_scale_to_m": 1.0},
        "materials": {"u": {"solid": {"young": young, "poisson": 0.0, "unit": "Pa"}}},
        "regions": {"u": {"solid": "solid"}},
        "boundary": {
            "dirichlet": [
                {
                    "field": "u",
                    "nodes": list(range(len(points))),
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
    return Model(mesh=mesh, document=document)


def domain(mesh, sides, *, width=4.0, mu=1.0, reference=False):
    document = {
        "fluid_domain_schema_version": 1,
        "mesh": "fluid.json",
        "width_m": width,
        "fluid": {"density_kg_m3": 1.0, "viscosity_Pa_s": mu},
        "boundaries": {
            "bottom": {"type": "wall"},
            "top": {"type": "wall", "fsi": True},
            **sides,
        },
    }
    if reference:
        document["pressure_reference_cell"] = 0
    return FluidDomain.from_document(document, mesh=mesh)


def periodic(length, drop=0.0):
    return {
        "left": {
            "type": "periodic",
            "partner": "right",
            "translation_m": [length, 0.0],
            "pressure_drop_Pa": drop,
        },
        "right": {
            "type": "periodic",
            "partner": "left",
            "translation_m": [-length, 0.0],
            "pressure_drop_Pa": -drop,
        },
    }


CHANNEL = {
    "left": {"type": "inlet", "velocity_m_s": [1.0, 0.0], "nu_tilde_m2_s": 0.0},
    "right": {"type": "pressure_outlet", "pressure_Pa": 0.0},
}


@pytest.fixture
def resting(quad_mesh, cube_mesh):
    """Fluid at rest at 7 Pa, 5 Pa outside: d_y=0.002(2-y) m independently."""
    fluid = domain(quad_mesh(4, 4, size_m=(3, 1)), periodic(3.0), reference=True)
    return StationaryFsiProblem(
        fluid_domain=fluid,
        wall_name="top",
        structure=structure(cube_mesh),
        facet_set="wet",
        length_scale_m=3.0,
        velocity_scale_m_s=1.0,
        fluid_options={"pressure_reference_Pa": 7.0},
        external_pressure_Pa=5.0,
    )


def test_geometry_independent_load_converges_geometrically_to_analytic_state(
    resting,
):
    result = solve_stationary_fsi(resting)
    points = resting.structure.elasticity_inputs()["mesh"].points
    expected = np.zeros_like(points)
    expected[:, 1] = 0.002 * (2 - points[:, 1])
    assert result.converged
    np.testing.assert_allclose(result.structure_displacement_m, expected, atol=1e-13)
    np.testing.assert_allclose(result.geometry_displacement_m, expected, atol=3e-8)
    # Constant load: r_k = (1-omega)^k d* exactly, and smallest k with r<=1e-8 Lref.
    residuals = np.array([row.displacement_residual_m for row in result.history])
    k = np.arange(len(residuals))
    np.testing.assert_allclose(residuals, 0.002 * 0.8**k, rtol=1e-9, atol=1e-15)
    assert result.iterations == int(np.ceil(np.log(1.5e-5) / np.log(0.8))) + 1
    assert residuals[-1] <= 3e-8 < residuals[-2]
    assert result.history[0].load_change_N == np.inf
    final = result.history[-1]
    assert final.load_change_N <= result.settings["load_gate_N"] == 1.2e-5
    np.testing.assert_allclose(final.fluid_force_N, [0, 84, 0], atol=1e-9)
    np.testing.assert_allclose(final.external_pressure_force_N, [0, -60, 0], atol=1e-9)
    np.testing.assert_allclose(final.structure_reaction_N, [0, -24, 0], atol=1e-9)
    assert final.transfer_relative_error <= 1e-10
    np.testing.assert_allclose(
        result.mesh_motion.domain.mesh.points[
            result.mesh_motion.domain.mesh.line2[
                result.mesh_motion.domain.mesh.boundary_sets["top"]
            ]
        ][..., 1],
        1 + expected[0, 1] * (1 - 0.8 ** (result.iterations - 1)),
        atol=1e-13,
    )
    assert result.settings["fsi_schema_version"] == 1
    with pytest.raises(ValueError):
        result.geometry_displacement_m.setflags(write=True)
    with pytest.raises(TypeError):
        result.settings["relaxation"] = 1.0  # type: ignore[index]


def test_hydrostatic_lid_load_feeds_back_through_moved_geometry(quad_mesh, cube_mesh):
    """Gravity: lifting the lid by d lowers its pressure, an exact linear loop.

    Harmonic motion of a uniform lid is affine, so the reference cell centre
    y_c rises by d*y_c/H0 and p_lid = p_ref - rho g (H0 - y_c) - rho g d (1 - y_c/H0).
    With d = (p_lid - p_ext) h/E the fixed point is independent of the solver.
    """
    rho_g, p_ref, p_ext, young, y_c = 100.0, 100.0, 5.0, 1000.0, 0.125
    fluid = domain(quad_mesh(4, 4, size_m=(3, 1)), periodic(3.0), reference=True)
    problem = StationaryFsiProblem(
        fluid_domain=fluid,
        wall_name="top",
        structure=structure(cube_mesh, young=young),
        facet_set="wet",
        length_scale_m=3.0,
        velocity_scale_m_s=1.0,
        fluid_options={"pressure_reference_Pa": p_ref, "force_N_m3": (0.0, -rho_g)},
        external_pressure_Pa=p_ext,
    )
    result = solve_stationary_fsi(problem)
    one_way = (p_ref - rho_g * (1 - y_c) - p_ext) / young
    coupled = one_way * young / (young + rho_g * (1 - y_c))
    lid = result.structure_displacement_m[result.interface._triangles.ravel(), 1]
    first, final = result.history[0], result.history[-1]
    assert result.converged
    # Stop gate |r|<=1e-8 Lref bounds |g(d^k)-d*| by |g'| |r|/(1+|g'|)=2.4e-9 m.
    np.testing.assert_allclose(lid, coupled, rtol=0, atol=2.5e-9)
    np.testing.assert_allclose(first.displacement_residual_m, one_way, rtol=1e-9)
    assert coupled < 0.93 * one_way
    # Lid force follows the moved geometry, not the reference pressure.
    np.testing.assert_allclose(
        final.fluid_force_N[1] + final.external_pressure_force_N[1],
        young * lid.mean() * 12.0,
        rtol=1e-10,
    )
    # Linear map: constant contraction 1 - omega (1 + rho g (1-y_c)/E).
    residuals = np.array([row.displacement_residual_m for row in result.history])
    np.testing.assert_allclose(
        residuals[1:] / residuals[:-1], 1 - 0.2 * (1 + 0.0875), rtol=1e-6
    )
    assert final.load_change_N <= result.settings["load_gate_N"]
    np.testing.assert_allclose(
        result.structure_displacement_m, result.geometry_displacement_m, atol=3e-8
    )


def test_identical_inputs_reproduce_history_and_identity(resting):
    short = StationaryFsiProblem(
        **{
            **{name: getattr(resting, name) for name in resting.__slots__},
            "relaxation": 0.9,
        }
    )
    a, b = solve_stationary_fsi(short), solve_stationary_fsi(short)
    assert a.input_sha256 == b.input_sha256
    assert [r.displacement_residual_m for r in a.history] == [
        r.displacement_residual_m for r in b.history
    ]
    np.testing.assert_array_equal(a.geometry_displacement_m, b.geometry_displacement_m)
    other = StationaryFsiProblem(
        **{
            **{name: getattr(resting, name) for name in resting.__slots__},
            "relaxation": 0.8,
        }
    )
    assert solve_stationary_fsi(other).input_sha256 != a.input_sha256


def test_budget_exhaustion_keeps_complete_history_and_last_iterate(resting):
    problem = StationaryFsiProblem(
        **{
            **{name: getattr(resting, name) for name in resting.__slots__},
            "max_coupling_iterations": 3,
        }
    )
    with pytest.raises(FsiCouplingError) as caught:
        solve_stationary_fsi(problem)
    error = caught.value
    assert error.kind == "budget" and len(error.history) == 3
    # d^3 = (1-0.8^3) d*, the last relaxed iterate, never declared converged.
    np.testing.assert_allclose(error.displacement_m[:, 1].max(), 0.002 * (1 - 0.8**3))
    with pytest.raises(ValueError):
        error.displacement_m.setflags(write=True)


def test_structure_gate_fluid_failure_and_fixed_wall_conflict_are_classified(
    resting, quad_mesh
):
    fields = {name: getattr(resting, name) for name in resting.__slots__}
    with pytest.raises(FsiCouplingError) as caught:
        solve_stationary_fsi(
            StationaryFsiProblem(**{**fields, "external_pressure_Pa": -500.0})
        )
    assert caught.value.kind == "structure" and caught.value.history == ()
    with pytest.raises(FsiCouplingError) as caught:
        solve_stationary_fsi(
            StationaryFsiProblem(
                **{
                    **fields,
                    "fluid_domain": domain(
                        quad_mesh(4, 4, size_m=(3, 1)), CHANNEL, width=4.0
                    ),
                    "fluid_options": {"max_iterations": 1},
                }
            )
        )
    assert caught.value.kind == "fluid"
    closed = domain(
        quad_mesh(4, 4, size_m=(3, 1)),
        {"left": {"type": "wall"}, "right": {"type": "wall"}},
        reference=True,
    )
    with pytest.raises(FsiCouplingError) as caught:
        solve_stationary_fsi(StationaryFsiProblem(**{**fields, "fluid_domain": closed}))
    assert caught.value.kind == "mesh" and len(caught.value.history) == 1


def test_divergence_requires_consecutive_unrelaxed_growth():
    assert not _diverging([1.0, 2.0, 3.0], 3)
    assert _diverging([1.0, 2.0, 3.0, 4.0], 3)
    assert not _diverging([1.0, 2.0, 1.5, 4.0, 5.0], 3)


@pytest.mark.parametrize(
    ("change", "match"),
    [
        ({"relaxation": 0.0}, "relaxation"),
        ({"relaxation": 1.5}, "relaxation"),
        ({"max_coupling_iterations": 0}, "max_coupling_iterations"),
        ({"length_scale_m": -1.0}, "scales"),
        ({"fluid_model": "k-epsilon"}, "fluid_model"),
        ({"fluid_model": "sa-noft2"}, "initial_nu_tilde"),
        ({"initial_nu_tilde_m2_s": 0.1}, "only valid"),
        ({"fluid_options": {"length_scale_m": 1.0}}, "reserved"),
        ({"wall_name": "bottom"}, "FSI-marked"),
        ({"facet_set": "missing"}, "facet set"),
        ({"preserve_inlet_flow": 1}, "boolean"),
    ],
)
def test_invalid_inputs_are_rejected_before_any_solve(resting, change, match):
    fields = {name: getattr(resting, name) for name in resting.__slots__}
    with pytest.raises(ValidationError, match=match):
        solve_stationary_fsi(StationaryFsiProblem(**{**fields, **change}))


def test_wall_resolved_sa_path_couples_through_the_same_loop(quad_mesh, cube_mesh):
    fluid = domain(
        quad_mesh(6, 8), periodic(1.0, 0.01), width=0.7, mu=0.1, reference=True
    )
    problem = StationaryFsiProblem(
        fluid_domain=fluid,
        wall_name="top",
        structure=structure(cube_mesh, length=1.0, width=0.7, young=100.0),
        facet_set="wet",
        length_scale_m=1.0,
        velocity_scale_m_s=1.0,
        fluid_model="sa-noft2",
        initial_nu_tilde_m2_s=0.3,
        fluid_options={"pressure_reference_Pa": 0.5},
        relaxation=1.0,
    )
    result = solve_stationary_fsi(problem)
    assert result.converged
    assert np.max(result.fluid.wall_y_plus) <= 1
    # Uniform 0.5 Pa periodic pressure: d_y = 0.5/100 (2-y) m.
    lift = result.structure_displacement_m[:, 1]
    points = problem.structure.elasticity_inputs()["mesh"].points
    np.testing.assert_allclose(lift, 0.005 * (2 - points[:, 1]), atol=1e-9)
