"""Independent manufactured error and a frozen localized error target."""

from pathlib import Path
from runpy import run_path

import numpy as np

from fem_toolbox.adaptive import AdaptiveOptions, HeatProblem, solve_adaptive_heat
from fem_toolbox.materials import ThermalMaterial

reference = run_path(
    str(Path(__file__).resolve().parents[2] / "examples/adaptive_heat/run_adaptive.py")
)


def test_localized_adaptive_uses_fewer_dofs_at_fixed_error():
    # Frozen after pilot exploration, before this regression; not a universal claim.
    target = 1.5  # sqrt(W K), independently integrated true energy error
    mesh = reference["cube"](4)
    exact, source, gradient = reference["functions"](True)
    problem = HeatProblem(
        ThermalMaterial.isotropic(1),
        {i: exact for i in range(mesh.n_triangle3)},
        source,
    )
    counts = {}
    for strategy, fraction in (("adaptive", 0.5), ("uniform", 1.0)):
        run = solve_adaptive_heat(
            mesh,
            problem,
            AdaptiveOptions(
                tolerance=0,
                max_steps=6,
                max_dofs=10000,
                fraction=fraction,
                stagnation_window=6,
            ),
        )
        assert run.status == "step_budget"
        errors = [
            reference["energy_error"](s.mesh, s.result, gradient, order=7)
            for s in run.stages
        ]
        assert np.all(np.diff(errors) < 0)
        first = next(i for i, error in enumerate(errors) if error <= target)
        counts[strategy] = len(run.stages[first].result.node_ids)
        for stage in run.stages:
            assert abs(stage.result.power_balance) < 1e-8
            assert len(stage.mesh.facet_sets["boundary"]) == stage.mesh.n_triangle3
            assert (
                sum(len(s) for s in stage.mesh.volume_sets.values())
                == stage.mesh.n_tetra4
            )
        assert len(run.stages) >= 4
    assert counts["adaptive"] < counts["uniform"]
