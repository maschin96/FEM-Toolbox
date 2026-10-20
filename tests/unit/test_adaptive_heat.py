"""Adaptive control decisions, diagnostics and boundary inheritance."""

import numpy as np
import pytest

from fem_toolbox.adaptive import AdaptiveOptions, HeatProblem, solve_adaptive_heat
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial


def problem(mesh, source=1.0):
    return HeatProblem(
        ThermalMaterial.isotropic(1),
        {i: 300.0 for i in range(mesh.n_triangle3)},
        source=source,
    )


def test_tolerance_budget_stagnation_and_repeat(cube_mesh):
    mesh = cube_mesh(2)
    p = problem(mesh)
    done = solve_adaptive_heat(mesh, p, AdaptiveOptions(tolerance=100))
    assert done.status == "tolerance"
    assert len(done.stages) == 1
    budget = solve_adaptive_heat(mesh, p, AdaptiveOptions(tolerance=0, max_dofs=27))
    assert budget.status == "dof_budget"
    assert len(budget.stages) == 1
    options = AdaptiveOptions(tolerance=0, max_steps=3)
    run = solve_adaptive_heat(mesh, p, options)
    assert run.status == "step_budget"
    assert len(run.stages) == 3
    assert all(abs(s.result.power_balance) < 1e-9 for s in run.stages)
    assert run.to_dict() == solve_adaptive_heat(mesh, p, options).to_dict()
    stagnant = solve_adaptive_heat(
        mesh,
        p,
        AdaptiveOptions(
            tolerance=0, max_steps=4, min_reduction=0.99, stagnation_window=1
        ),
    )
    assert stagnant.status == "stagnation"
    assert len(stagnant.stages) == 2


def test_errors_are_diagnosed_without_losing_completed_stages(cube_mesh):
    mesh = cube_mesh(1)
    invalid = solve_adaptive_heat(mesh, HeatProblem(ThermalMaterial.isotropic(1), {}))
    assert invalid.status == "error"
    assert invalid.stages == ()
    assert "nullspace" in invalid.diagnostic.lower()

    def broken(*args, **kwargs):
        raise ValidationError("adapter failure", field="remesh")

    run = solve_adaptive_heat(
        mesh, problem(mesh), AdaptiveOptions(tolerance=0), remesher=broken
    )
    assert run.status == "error"
    assert len(run.stages) == 1
    assert "adapter failure" in run.diagnostic


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tolerance": -1},
        {"max_steps": 0},
        {"fraction": 0},
        {"max_dofs": True},
        {"min_reduction": float("nan")},
    ],
)
def test_invalid_options(kwargs):
    with pytest.raises(ValidationError):
        AdaptiveOptions(**kwargs)


def test_layered_material_and_boundary_inheritance(cube_mesh):
    from fem_toolbox.mesh import Mesh
    from fem_toolbox.mesh.refinement import refine_tet4

    base = cube_mesh(2)
    tags = (base.points[base.tetra4].mean(axis=1)[:, 0] >= 0.5).astype(int)
    mesh = Mesh(
        points=base.points,
        tetra4=base.tetra4,
        triangle3=base.triangle3,
        cell_tags=tags,
        volume_sets={"left": np.flatnonzero(tags == 0)},
    )

    def exact(x):
        return 300 + min(x[0], 0.5) + max(x[0] - 0.5, 0) / 2

    p = HeatProblem(
        tuple(ThermalMaterial.isotropic(1 + int(t)) for t in tags),
        {i: exact for i in range(mesh.n_triangle3)},
    )
    refinement = refine_tet4(mesh, [0, 25])
    run = solve_adaptive_heat(
        refinement.mesh, p.refined(refinement), AdaptiveOptions(tolerance=1e-8)
    )
    assert run.status == "tolerance"
    stage = run.stages[0]
    np.testing.assert_allclose(
        stage.result.temperature,
        [exact(x) for x in stage.mesh.points[stage.result.node_ids]],
        atol=1e-10,
    )
    np.testing.assert_allclose(stage.result.heat_flux[:, 0], -1, atol=1e-10)


def test_flux_robin_and_dirichlet_inheritance(cube_mesh):
    from fem_toolbox.mesh.refinement import refine_tet4
    from fem_toolbox.physics.heat import Convection

    mesh = cube_mesh(2)
    centers = mesh.points[mesh.triangle3].mean(axis=1)
    fixed = {int(i): 300.0 for i in np.flatnonzero(centers[:, 0] == 0)}
    flux = {int(i): -1.0 for i in np.flatnonzero(centers[:, 0] == 1)}
    robin = {
        int(i): Convection(1.0, lambda x: 300 + x[0])
        for i in np.flatnonzero(centers[:, 1] == 1)
    }
    p = HeatProblem(
        ThermalMaterial.isotropic(1), fixed, outward_fluxes=flux, convection=robin
    )
    refinement = refine_tet4(mesh, np.arange(mesh.n_tetra4))
    run = solve_adaptive_heat(
        refinement.mesh, p.refined(refinement), AdaptiveOptions(tolerance=1e-8)
    )
    assert run.status == "tolerance"
    stage = run.stages[0]
    np.testing.assert_allclose(
        stage.result.temperature,
        300 + stage.mesh.points[stage.result.node_ids, 0],
        atol=1e-10,
    )
    assert abs(stage.result.power_balance) < 1e-9


def test_invalid_remesher_parent_maps_are_rejected(cube_mesh):
    from fem_toolbox.mesh.refinement import Refinement, refine_tet4

    mesh = cube_mesh(1)

    def wrong(mesh, marked, **kwargs):
        valid = refine_tet4(mesh, marked, **kwargs)
        return Refinement(
            valid.mesh, valid.cell_parents + mesh.n_tetra4, valid.facet_parents
        )

    run = solve_adaptive_heat(
        mesh, problem(mesh), AdaptiveOptions(tolerance=0), remesher=wrong
    )
    assert run.status == "error"
    assert len(run.stages) == 1
    assert "parent map" in run.diagnostic
