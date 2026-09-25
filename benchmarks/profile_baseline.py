"""Reproducible developer baseline; run from repository root with PYTHONPATH=src.

Each sample runs in a fresh process. No timing threshold is a correctness gate.
Instrumentation patches only this single-threaded benchmark's solver bindings.
"""

from __future__ import annotations

import argparse
import cProfile
import json
import os
import platform
import pstats
import statistics
import subprocess
import sys
import tracemalloc
from collections import Counter, defaultdict
from contextlib import ExitStack
from importlib.metadata import version
from itertools import permutations, product
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import numpy as np

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.io import write_vtu
from fem_toolbox.io.heat_series import HeatSeriesWriter
from fem_toolbox.materials import PiezoelectricMaterial, ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers import heat, piezoelectric, transient_heat

CASES = ("heat", "electrode", "transient")
THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
)


def cube(n):
    """Unit cube, six conforming Tet4 per Cartesian cell; complete exterior."""
    points = np.array(list(product(np.linspace(0, 1, n + 1), repeat=3)))

    def index(p):
        return (p[0] * (n + 1) + p[1]) * (n + 1) + p[2]

    cells = []
    for base in product(range(n), repeat=3):
        for order in permutations(range(3)):
            p = list(base)
            cell = [index(p)]
            for axis in order:
                p[axis] += 1
                cell.append(index(p))
            cells.append(cell)
    cells = np.array(cells)
    faces = Counter(tuple(sorted(np.delete(c, i))) for c in cells for i in range(4))
    return Mesh(
        points=points,
        tetra4=cells,
        triangle3=np.array([f for f, count in faces.items() if count == 1]),
    )


def problem(case, mesh, steps):
    """Return unchanged production solve, primary field, and analytic expectation."""
    x = mesh.points[:, 0]
    left = np.flatnonzero(x == 0)
    right = np.flatnonzero(x == 1)
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    if case == "heat":
        conditions = [
            Dirichlet(left, np.full(len(left), 300.0)),
            Dirichlet(right, np.full(len(right), 301.0)),
        ]
        return (
            heat,
            lambda: heat.solve_heat(mesh, material, conditions),
            lambda r: r.temperature,
            300 + x,
        )
    if case == "transient":
        return (
            transient_heat,
            lambda: transient_heat.solve_transient_heat(
                mesh, material, np.linspace(0, 1, steps + 1), 300.0, source=12.0
            ),
            lambda r: r.frames[-1].temperature,
            np.full(mesh.n_points, 301.0),
        )
    stiffness = np.eye(6) * 10
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.5
    piezo = PiezoelectricMaterial(stiffness, coupling, np.eye(3) * 2)
    dofs = piezoelectric.piezoelectric_dofs(mesh)
    nodes = np.arange(mesh.n_points)
    conditions = [
        Dirichlet(dofs.indices("u", nodes, ["y", "z"]), np.zeros(2 * len(nodes))),
        Dirichlet(dofs.indices("u", left, ["x"]), np.zeros(len(left))),
    ]
    left_faces = np.flatnonzero(np.all(x[mesh.triangle3] == 0, axis=1))
    right_faces = np.flatnonzero(np.all(x[mesh.triangle3] == 1, axis=1))
    electrodes = [
        Electrode("ground", left_faces, potential=0),
        Electrode("floating", right_faces, charge=0),
    ]
    return (
        piezoelectric,
        lambda: piezoelectric.solve_piezoelectric(
            mesh,
            piezo,
            conditions,
            electrodes=electrodes,
            tractions={int(i): [2, 0, 0] for i in right_faces},
        ),
        lambda r: r.potential,
        x / 20.25,
    )


class Measurements:
    """Inclusive function-boundary timings, excluding nested measured calls."""

    def __init__(self):
        self.seconds = defaultdict(float)
        self.depth = 0
        self.systems = []

    def wrap(self, name, function):
        def measured(*args, **kwargs):
            outer = self.depth == 0
            self.depth += 1
            start = perf_counter()
            try:
                return function(*args, **kwargs)
            finally:
                elapsed = perf_counter() - start
                self.depth -= 1
                if outer:
                    self.seconds[name] += elapsed

        return measured

    def factor(self, function):
        def measured(matrix, *args, **kwargs):
            factor = self.wrap("factorization", function)(matrix, *args, **kwargs)
            self.systems.append({"dofs": matrix.shape[0], "nnz": matrix.nnz})
            owner = self

            class Factor:
                def solve(self, *args, **kwargs):
                    return owner.wrap("solve", factor.solve)(*args, **kwargs)

                def __getattr__(self, name):
                    return getattr(factor, name)

            return Factor()

        return measured

    def install(self, stack, module):
        groups = {
            "geometry": ["Tet4Geometry"],
            "assembly": [
                "assemble",
                "conductance_matrix",
                "capacity_matrix",
                "piezoelectric_matrix",
                "source_load",
                "traction_load",
            ],
            "boundary": ["eliminate_dirichlet", "_electrode_reduction"],
        }
        for phase, names in groups.items():
            for name in names:
                if hasattr(module, name):
                    stack.enter_context(
                        patch.object(
                            module, name, self.wrap(phase, getattr(module, name))
                        )
                    )
        stack.enter_context(patch.object(module, "splu", self.factor(module.splu)))


def environment():
    cpu = platform.processor() or platform.machine()
    cpu_probe_error = None
    if sys.platform == "darwin":
        probe = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"], capture_output=True, text=True
        )
        if probe.returncode == 0:
            cpu = probe.stdout.strip()
        else:
            cpu_probe_error = probe.stderr.strip()
    return {
        "platform": platform.platform(),
        "cpu": cpu,
        "cpu_probe_error": cpu_probe_error,
        "logical_cpus": os.cpu_count(),
        "python": platform.python_version(),
        "versions": {n: version(n) for n in ("numpy", "scipy", "meshio")},
        "thread_environment": {n: os.environ.get(n) for n in THREAD_VARS},
        "numpy_configuration": np.show_config(mode="dicts"),
        "git_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "git_dirty": bool(subprocess.check_output(["git", "status", "--porcelain"])),
    }


def verify(case, mesh, result):
    """Independent affine/uniform solutions and SI balance checks, not timings."""
    x = mesh.points[:, 0]
    if case == "heat":
        np.testing.assert_allclose(
            result.heat_flux, np.tile([-2, 0, 0], (mesh.n_tetra4, 1)), rtol=0, atol=1e-9
        )
        balance = abs(result.power_balance)
        residual = result.relative_residual
    elif case == "electrode":
        # D_x=.5*epsilon_xx - 2*phi'=0; sigma_xx=10*epsilon_xx+.5*phi'=2.
        displacement = np.column_stack(
            [2 * x / 10.125, np.zeros_like(x), np.zeros_like(x)]
        )
        np.testing.assert_allclose(result.displacement, displacement, rtol=0, atol=1e-9)
        np.testing.assert_allclose(result.electric_displacement, 0, rtol=0, atol=1e-9)
        balance = max(
            float(np.max(np.abs(result.force_balance))), abs(result.charge_balance)
        )
        residual = max(result.relative_residual_u, result.relative_residual_phi)
    else:
        for frame in result.frames:
            np.testing.assert_allclose(
                frame.temperature, 300 + frame.time, rtol=0, atol=1e-9
            )
        balance = max(abs(b.balance) for b in result.balances)
        residual = max(b.relative_residual for b in result.balances)
    if balance > 1e-9 or residual > 1e-9:
        raise RuntimeError("benchmark balance/residual gate failed")
    return {"balance_error": balance, "relative_residual": residual}


def export(case, output, mesh, result, field):
    write_vtu(output / "solution.vtu", mesh, point_data={"field": field(result)})
    if case == "transient":
        writer = HeatSeriesWriter(output / "series", mesh)
        for frame in result.frames:
            writer.write(frame)
        writer.finish(result)


def worker(case, n, steps, output):
    output.mkdir(parents=True, exist_ok=False)
    mesh = cube(n)
    module, solve, field, expected = problem(case, mesh, steps)
    start = perf_counter()
    reference = solve()
    wall = perf_counter() - start
    error = float(np.max(np.abs(field(reference) - expected)))
    # Absolute 1e-9 K/V: O(1) manufactured fields (300 K offset), double precision.
    if error > 1e-9:
        raise RuntimeError(f"analytic reference failed: {error}")
    verification = verify(case, mesh, reference)
    measurements = Measurements()
    profiler = cProfile.Profile()
    with ExitStack() as stack:
        measurements.install(stack, module)
        profiler.enable()
        start = perf_counter()
        result = solve()
        measured_wall = perf_counter() - start
        profiler.disable()
    parity = float(np.max(np.abs(field(result) - field(reference))))
    if parity != 0:
        raise RuntimeError(f"instrumentation changed solution: {parity}")
    verify(case, mesh, result)
    measurements.wrap("export", export)(case, output, mesh, result, field)
    # Separate memory pass: tracing overhead is absent from wall/phase times.
    tracemalloc.start()
    memory_result = solve()
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    np.testing.assert_array_equal(field(memory_result), field(reference))
    try:
        import resource

        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss *= 1 if sys.platform == "darwin" else 1024
    except ImportError:
        rss = None
    profiler.dump_stats(str(output / "solve.prof"))
    stats = pstats.Stats(profiler)
    top = sorted(stats.stats.items(), key=lambda item: item[1][2], reverse=True)[:20]
    report = {
        "schema": "fem-toolbox-performance-v1",
        "case": case,
        "size": n,
        "steps": steps if case == "transient" else None,
        "nodes": mesh.n_points,
        "cells": mesh.n_tetra4,
        "dofs": mesh.n_points * (4 if case == "electrode" else 1),
        "wall_s": wall,
        "profiled_wall_s": measured_wall,
        "phases_s": dict(measurements.seconds),
        "unclassified_profiled_s": measured_wall
        - sum(t for name, t in measurements.seconds.items() if name != "export"),
        "factorized_systems": measurements.systems,
        "matrix_solves": len(measurements.systems),
        "traced_peak_bytes": peak,
        "process_peak_rss_bytes": rss,
        "max_absolute_error": error,
        "verification": verification,
        "instrumented_parity_error": parity,
        "input_sha256": reference.provenance.input_sha256,
        "environment": environment(),
        "top_self_time": [
            {
                "file": Path(key[0]).name,
                "line": key[1],
                "function": key[2],
                "calls": data[1],
                "self_s": data[2],
            }
            for key, data in top
        ],
    }
    (output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def positive(value):
    integer = int(value)
    if integer < 1:
        raise argparse.ArgumentTypeError("expected positive integer")
    return integer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=CASES)
    parser.add_argument("--size", type=positive, default=2)
    parser.add_argument("--sizes", type=positive, nargs="+", default=[2, 4, 6])
    parser.add_argument("--steps", type=positive, default=5)
    parser.add_argument("--repeats", type=positive, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker, args.size, args.steps, args.output)
        return
    args.output.mkdir(parents=True, exist_ok=False)
    records, summary = [], []
    for case in CASES:
        for n in args.sizes:
            samples = []
            for repeat in range(args.repeats):
                target = args.output / f"{case}-{n}-{repeat}"
                subprocess.run(
                    [
                        sys.executable,
                        str(Path(__file__).resolve()),
                        "--worker",
                        case,
                        "--size",
                        str(n),
                        "--steps",
                        str(args.steps),
                        "--output",
                        str(target),
                    ],
                    check=True,
                    env={**os.environ, **dict.fromkeys(THREAD_VARS, "1")},
                )
                record = json.loads((target / "report.json").read_text())
                records.append(record)
                samples.append(record["wall_s"])
            summary.append(
                {
                    "case": case,
                    "size": n,
                    "median_s": statistics.median(samples),
                    "min_s": min(samples),
                    "max_s": max(samples),
                    "stdev_s": statistics.stdev(samples) if len(samples) > 1 else None,
                }
            )
            print(json.dumps(summary[-1]), flush=True)
    (args.output / "baseline.json").write_text(
        json.dumps({"samples": records, "summary": summary}, indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
