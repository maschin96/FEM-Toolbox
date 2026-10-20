"""Reproducible adaptive/uniform manufactured-solution comparison.

Run: python examples/adaptive_heat/run_adaptive.py --output results/adaptive_heat
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import permutations, product
from pathlib import Path
from time import perf_counter

import numpy as np

from fem_toolbox.adaptive import AdaptiveOptions, HeatProblem, solve_adaptive_heat
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh


def cube(n: int) -> Mesh:
    """Freudenthal triangulation, independent of production remeshing."""
    points = np.array(list(product(np.linspace(0, 1, n + 1), repeat=3)))
    cells = []
    for base in product(range(n), repeat=3):
        for order in permutations(range(3)):
            position = np.array(base)
            nodes = []
            for axis in (None, *order):
                if axis is not None:
                    position = position.copy()
                    position[axis] += 1
                nodes.append(
                    int((position[0] * (n + 1) + position[1]) * (n + 1) + position[2])
                )
            cells.append(nodes)
    cells_array = np.array(cells)
    faces = Counter(
        tuple(sorted(np.delete(c, i))) for c in cells_array for i in range(4)
    )
    facets = np.array([f for f, count in faces.items() if count == 1])
    tags = (points[cells_array].mean(axis=1)[:, 0] >= 0.5).astype(int)
    return Mesh(
        points=points,
        tetra4=cells_array,
        triangle3=facets,
        cell_tags=tags,
        facet_tags=np.ones(len(facets), dtype=int),
        volume_sets={
            "left": np.flatnonzero(tags == 0),
            "right": np.flatnonzero(tags == 1),
        },
        facet_sets={"boundary": np.arange(len(facets))},
    )


def energy_error(mesh, result, gradient, order=5):
    """Independent tensor Gauss/Duffy integration of |grad(T)-grad(T_h)|²."""
    roots, weights = np.polynomial.legendre.leggauss(order)
    roots, weights = (roots + 1) / 2, weights / 2
    bary, quadrature = [], []
    for a, b, c in product(range(order), repeat=3):
        r, s, t = roots[[a, b, c]]
        xyz = [r, (1 - r) * s, (1 - r) * (1 - s) * t]
        bary.append([1 - sum(xyz), *xyz])
        quadrature.append(
            6 * weights[a] * weights[b] * weights[c] * (1 - r) ** 2 * (1 - s)
        )
    bary, quadrature = np.array(bary), np.array(quadrature)
    squared = 0.0
    for i, cell in enumerate(mesh.tetra4):
        g = Tet4Geometry(mesh.points[cell])
        difference = result.gradient[i] - gradient(bary @ mesh.points[cell])
        squared += g.volume * float(quadrature @ np.sum(difference**2, axis=1))
    return np.sqrt(squared)


def functions(localized):
    if localized:

        def exact(x):
            return float(300 + np.exp(-80 * (x[0] - 0.5) ** 2))

        def source(x):
            return float(
                (160 - 25600 * (x[0] - 0.5) ** 2) * np.exp(-80 * (x[0] - 0.5) ** 2)
            )

        def gradient(x):
            g = np.zeros_like(x)
            g[:, 0] = -160 * (x[:, 0] - 0.5) * np.exp(-80 * (x[:, 0] - 0.5) ** 2)
            return g
    else:

        def exact(x):
            return float(300 + np.sum(x * (1 - x)))

        def source(x):
            return 6.0

        def gradient(x):
            return 1 - 2 * x

    return exact, source, gradient


def compare(output: Path, *, steps=7, n=4):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for case, localized in (("smooth", False), ("localized", True)):
        exact, source, gradient = functions(localized)
        mesh = cube(n)
        problem = HeatProblem(
            ThermalMaterial.isotropic(1),
            {i: exact for i in range(mesh.n_triangle3)},
            source,
        )
        for strategy, fraction in (("adaptive", 0.5), ("uniform", 1.0)):
            start = perf_counter()
            run = solve_adaptive_heat(
                mesh,
                problem,
                AdaptiveOptions(
                    tolerance=0,
                    max_steps=steps,
                    max_dofs=50000,
                    fraction=fraction,
                    stagnation_window=steps,
                ),
            )
            elapsed = perf_counter() - start
            run.export(output / f"{case}-{strategy}.json")
            if run.status != "step_budget":
                raise RuntimeError(run.diagnostic)
            for level, stage in enumerate(run.stages):
                error = energy_error(stage.mesh, stage.result, gradient, order=7)
                lower = energy_error(stage.mesh, stage.result, gradient, order=6)
                if not np.isclose(error, lower, rtol=2e-4):
                    raise RuntimeError("independent error quadrature unresolved")
                rows.append(
                    {
                        "case": case,
                        "strategy": strategy,
                        "level": level,
                        "dofs": len(stage.result.node_ids),
                        "cells": stage.mesh.n_tetra4,
                        "error": float(error),
                        "indicator": stage.indicator.total,
                        "effectivity": stage.indicator.total / error,
                        "power_balance": stage.result.power_balance,
                        "mesh_sha256": stage.to_dict()["mesh_sha256"],
                        "solve_refine_seconds_total": elapsed,
                    }
                )
                print(
                    case,
                    strategy,
                    level,
                    len(stage.result.node_ids),
                    f"error={error:.6g}",
                    flush=True,
                )
    (output / "comparison.json").write_text(json.dumps(rows, indent=2) + "\n")
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("results/adaptive_heat"))
    parser.add_argument("--steps", type=int, default=7)
    parser.add_argument("--n", type=int, default=4)
    args = parser.parse_args()
    compare(args.output, steps=args.steps, n=args.n)
