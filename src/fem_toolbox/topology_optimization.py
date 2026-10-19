"""Bounded, deterministic SLSQP driver for filtered SIMP compliance."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize  # type: ignore[import-untyped]

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError
from fem_toolbox.topology import DensityEvaluation, TopologyProblem


@dataclass(frozen=True)
class OptimizationOptions:
    """Iteration/design-evaluation budgets and dimensionless SLSQP tolerance."""

    max_iterations: int = 100
    max_evaluations: int = 300
    tolerance: float = 1e-7

    def __post_init__(self) -> None:
        for value in (self.max_iterations, self.max_evaluations):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValidationError(
                    "positive integer optimization budget required", field="options"
                )
        if not np.isfinite(self.tolerance) or not 0 < self.tolerance < 1:
            raise ValidationError("finite tolerance in (0,1) required", field="options")


@dataclass(frozen=True)
class OptimizationResult:
    """Final feasible state, explicit stop status and accepted-step diagnostics.

    A budget stop retains the best feasible state without claiming convergence.
    Each history row contains SI compliance, fraction, residual and step change.
    """

    status: str
    message: str
    final: DensityEvaluation
    initial_objective: float
    volume_limit: float
    evaluations: int
    history: tuple[dict[str, Any], ...]
    options: OptimizationOptions


class _EvaluationLimit(Exception):
    pass


def optimize_topology(
    problem: TopologyProblem,
    volume_fraction: float,
    *,
    initial: ArrayLike | None = None,
    options: OptimizationOptions | None = None,
) -> OptimizationResult:
    """Minimize weighted compliance with exact filter gradients and SI states.

    Invalid configurations raise ValidationError. Numerical failures during the
    search return solver_failed with the last best feasible state and message.
    The initial solve itself must succeed; no invalid result is fabricated.
    """
    options = options or OptimizationOptions()
    feasible = problem.initial_density(volume_fraction)
    x = feasible if initial is None else real_array(initial, (problem.mesh.n_cells,))
    start = problem.evaluate(x)
    if start.volume_fraction > volume_fraction + 1e-10:
        raise ValidationError("initial design violates volume budget", field="initial")
    cached = start
    best = start
    evaluations = 1
    history: list[dict[str, Any]] = []
    previous = start

    def record(state: DensityEvaluation) -> None:
        nonlocal previous
        history.append(
            {
                "iteration": len(history),
                "compliance_J": state.objective,
                "volume_fraction": state.volume_fraction,
                "volume_violation": max(0.0, state.volume_fraction - volume_fraction),
                "maximum_relative_residual": max(
                    s.relative_residual for s in state.states
                ),
                "maximum_force_balance_N": max(
                    float(np.linalg.norm(s.force_balance)) for s in state.states
                ),
                "maximum_density_change": float(
                    np.max(np.abs(state.density - previous.density))
                ),
                "relative_objective_change": abs(state.objective - previous.objective)
                / start.objective,
                "state_status": "succeeded",
            }
        )
        previous = state

    record(start)

    def evaluate(z: NDArray[np.float64]) -> DensityEvaluation:
        nonlocal cached, best, evaluations
        trial = x.copy()
        trial[problem.free] = np.clip(z, 0, 1)
        if np.array_equal(trial, cached.density):
            return cached
        if evaluations >= options.max_evaluations:
            raise _EvaluationLimit("design evaluation budget exhausted")
        evaluations += 1
        cached = problem.evaluate(trial)
        if (
            cached.volume_fraction <= volume_fraction + 1e-8
            and cached.objective < best.objective
        ):
            best = cached
        return cached

    def objective(z: NDArray[np.float64]) -> tuple[float, NDArray[np.float64]]:
        state = evaluate(z)
        return state.objective / start.objective, state.objective_gradient[
            problem.free
        ] / start.objective

    def constraint(z: NDArray[np.float64]) -> float:
        # The density filter and fixed-region overrides are affine.
        return float(
            volume_fraction
            - start.volume_fraction
            - problem.volume_gradient[problem.free] @ (z - x[problem.free])
        )

    status = "converged"
    message = "no free design variables"
    if len(problem.free):
        try:
            solved = minimize(
                objective,
                x[problem.free],
                jac=True,
                method="SLSQP",
                bounds=[(0.0, 1.0)] * len(problem.free),
                constraints=[
                    {
                        "type": "ineq",
                        "fun": constraint,
                        "jac": lambda z: -problem.volume_gradient[problem.free],
                    }
                ],
                callback=lambda z: record(evaluate(z)),
                options={"maxiter": options.max_iterations, "ftol": options.tolerance},
            )
            message = str(solved.message)
            state = evaluate(solved.x)
            if solved.success and state.volume_fraction <= volume_fraction + 1e-8:
                status = "converged"
                best = state
            else:
                status = "iteration_limit" if solved.status == 9 else "solver_failed"
        except _EvaluationLimit as error:
            status = "evaluation_limit"
            message = str(error)
        except ValidationError as error:
            status = "solver_failed"
            message = str(error)
    if not np.array_equal(previous.density, best.density):
        record(best)
    return OptimizationResult(
        status,
        message,
        best,
        start.objective,
        float(volume_fraction),
        evaluations,
        tuple(history),
        options,
    )
