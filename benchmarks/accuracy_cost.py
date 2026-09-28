"""M7 accuracy/cost comparison. Fresh sequential workers; no timing CI gates."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import subprocess
import sys
from importlib.metadata import version
from itertools import product
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy import sparse

# Developer tool only: reuse the baseline mesh, never a test fixture.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmarks.profile_baseline import THREAD_VARS, cube
from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import ElasticMaterial, ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.elasticity import (
    mechanical_dofs,
    solve_elasticity,
)
from fem_toolbox.solvers.heat import solve_heat
from fem_toolbox.solvers.linear_backend import ScipyBackend

EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))


def make_mesh(n, order):
    mesh = cube(n)
    if order == 1:
        return mesh
    if order != 2:
        raise ValueError("order must be 1 or 2")
    points = mesh.points.tolist()
    edges = {}
    cells = []
    for cell in mesh.cells:
        mids = []
        for a, b in EDGES:
            edge = tuple(sorted((int(cell[a]), int(cell[b]))))
            if edge not in edges:
                edges[edge] = len(points)
                points.append(
                    ((mesh.points[edge[0]] + mesh.points[edge[1]]) / 2).tolist()
                )
            mids.append(edges[edge])
        cells.append([*cell, *mids])
    facets = [
        [
            *face,
            *[
                edges[tuple(sorted((int(face[a]), int(face[b]))))]
                for a, b in ((0, 1), (1, 2), (2, 0))
            ],
        ]
        for face in mesh.facets
    ]
    return Mesh(points=points, tetra10=cells, triangle6=facets)


def reference(case):
    if case == "heat":
        coefficients = np.array([1, 0.5, 0.25])
        return (
            lambda p: (1 + p**4 @ coefficients)[:, None],
            lambda p: (4 * p**3 * coefficients)[:, None, :],
        )

    def value(p):
        x, _, z = p.T
        return 0.01 * np.column_stack((-3 * x * x * z, np.zeros(len(p)), x**3))

    def gradient(p):
        x, _, z = p.T
        g = np.zeros((len(p), 3, 3))
        g[:, 0, 0] = -0.06 * x * z
        g[:, 0, 2] = -0.03 * x * x
        g[:, 2, 0] = 0.03 * x * x
        return g

    return value, gradient


def relative_errors(mesh, values, exact, exact_gradient):
    """Independent Duffy/Legendre-6 integration and explicit P1/P2 basis."""
    roots, weights = np.polynomial.legendre.leggauss(6)
    roots = (roots + 1) / 2
    weights = weights / 2
    bary = []
    qw = []
    for i, j, k in product(range(6), repeat=3):
        a, b, c = roots[[i, j, k]]
        x, y, z = a, (1 - a) * b, (1 - a) * (1 - b) * c
        bary.append([1 - x - y - z, x, y, z])
        qw.append(weights[i] * weights[j] * weights[k] * (1 - a) ** 2 * (1 - b))
    bary = np.array(bary)
    qw = np.array(qw)
    shapes = (
        bary
        if mesh.n_tetra4
        else np.column_stack(
            [bary * (2 * bary - 1), *[4 * bary[:, i] * bary[:, j] for i, j in EDGES]]
        )
    )
    total = np.zeros(4)
    for cell in mesh.cells:
        corners = mesh.points[cell[:4]]
        jacobian = (corners[1:] - corners[0]).T
        bg = np.vstack([-np.ones(3), np.eye(3)]) @ np.linalg.inv(jacobian)
        gradients = (
            np.broadcast_to(bg, (len(bary), 4, 3))
            if mesh.n_tetra4
            else np.concatenate(
                [
                    (4 * bary - 1)[:, :, None] * bg[None, :, :],
                    np.stack(
                        [
                            4 * (bary[:, i, None] * bg[j] + bary[:, j, None] * bg[i])
                            for i, j in EDGES
                        ],
                        axis=1,
                    ),
                ],
                axis=1,
            )
        )
        xyz = bary @ corners
        ex = exact(xyz)
        eg = exact_gradient(xyz)
        numerical = shapes @ values[cell]
        ng = np.einsum("qni,nj->qji", gradients, values[cell])
        integrands = np.column_stack(
            [
                np.sum((numerical - ex) ** 2, axis=1),
                np.sum((ng - eg) ** 2, axis=(1, 2)),
                np.sum(ex**2, axis=1),
                np.sum(eg**2, axis=(1, 2)),
            ]
        )
        total += abs(np.linalg.det(jacobian)) * (qw @ integrands)
    return {
        "l2_relative": float(np.sqrt(total[0] / total[2])),
        "h1_relative": float(np.sqrt(total[1] / total[3])),
    }


class MeasuredBackend:
    """Observe the original operator boundary, without patching production code."""

    def __init__(self, backend):
        self.backend = backend
        self.name = backend.name
        self.systems = []

    def factorize(self, operator, **kwargs):
        start = perf_counter()
        factor = self.backend.factorize(operator, **kwargs)
        elapsed = perf_counter() - start
        record = {
            "rows": operator.shape[0],
            "nnz": int(sparse.csc_matrix(operator).nnz),
            "factor_s": elapsed,
            "solve_s": 0.0,
        }
        self.systems.append(record)

        class Factor:
            def solve(self, rhs):
                start = perf_counter()
                result = factor.solve(rhs)
                record["solve_s"] += perf_counter() - start
                record["iterations"] = result.report.iterations
                return result

        return Factor()


def backend_for(name, scalar):
    if name == "scipy":
        return ScipyBackend()
    if scalar:
        paths = list(
            Path("/usr/lib/petscdir").glob(
                f"petsc*/**/*-{scalar}/lib/python3/dist-packages"
            )
        )
        if len(paths) != 1:
            raise RuntimeError(
                f"Expected exactly one {scalar} PETSc build, found {paths}"
            )
        sys.path.insert(0, str(paths[0]))
    from fem_toolbox.solvers.petsc_backend import PetscBackend

    backend = PetscBackend()
    if scalar and backend.scalar_dtype.kind != ("c" if scalar == "complex" else "f"):
        raise RuntimeError("Wrong PETSc scalar build")
    return backend


def rss():
    try:
        import resource

        return int(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            * (1 if sys.platform == "darwin" else 1024)
        )
    except ImportError:
        return None


def worker(args):
    args.output.mkdir(parents=True, exist_ok=False)
    mesh = make_mesh(args.size, args.order)
    exact, gradient = reference(args.worker)
    fixed = np.unique(mesh.facets)
    backend = MeasuredBackend(backend_for(args.backend, args.petsc_scalar))
    if args.worker == "heat":

        def solve():
            return solve_heat(
                mesh,
                ThermalMaterial.isotropic(1),
                [Dirichlet(fixed, exact(mesh.points[fixed]).ravel())],
                source=lambda p: -12 * (p * p @ np.array([1, 0.5, 0.25])),
                backend=backend,
            )
    else:

        def solve():
            return solve_elasticity(
                mesh,
                ElasticMaterial.isotropic(100, 0),
                [
                    Dirichlet(
                        mechanical_dofs(mesh).indices("u", fixed),
                        exact(mesh.points[fixed]).ravel(),
                    )
                ],
                body_force=lambda p: [6 * p[2], 0, 0],
                backend=backend,
            )

    rss_before = rss()
    start = perf_counter()
    result = solve()
    wall = perf_counter() - start
    peak = rss()  # Before reference quadrature and field serialization.
    values = (
        result.temperature[:, None] if args.worker == "heat" else result.displacement
    )
    balance = (
        abs(result.power_balance)
        if args.worker == "heat"
        else max(
            np.max(np.abs(result.force_balance)), np.max(np.abs(result.moment_balance))
        )
    )
    if result.relative_residual > 1e-9 or balance > 1e-9:
        raise RuntimeError("Physical residual/balance gate failed")
    report = {
        "case": args.worker,
        "order": args.order,
        "size": args.size,
        "backend": args.backend,
        "backend_name": backend.name,
        "cells": mesh.n_cells,
        "nodes": mesh.n_points,
        "dofs": int(values.size),
        "wall_s": wall,
        "rss_bytes": peak,
        "rss_before_bytes": rss_before,
        "systems": backend.systems,
        "residual": result.relative_residual,
        "balance": float(balance),
        "input_sha256": result.provenance.input_sha256,
        **relative_errors(mesh, values, exact, gradient),
    }
    np.save(args.output / "field.npy", values)
    (args.output / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )


def comparison(rows, error_targets, time_budgets, memory_budget):
    """Select measured feasible groups; unattained targets stay explicit nulls."""
    keys = sorted({(r["case"], r["backend"], r["order"], r["size"]) for r in rows})
    groups = []
    for case, backend, order, size in keys:
        sample = [
            r
            for r in rows
            if (r["case"], r["backend"], r["order"], r["size"])
            == (case, backend, order, size)
        ]
        times = [r["wall_s"] for r in sample]
        memories = [r["rss_bytes"] for r in sample]
        groups.append(
            {
                "case": case,
                "backend": backend,
                "order": order,
                "size": size,
                "repeats": len(sample),
                "h1_relative": max(r["h1_relative"] for r in sample),
                "rss_max_bytes": None if None in memories else max(memories),
                "wall_s": {
                    "median": statistics.median(times),
                    "min": min(times),
                    "max": max(times),
                    "stdev": statistics.stdev(times) if len(times) > 1 else None,
                },
            }
        )
    accuracy = []
    cost = []
    for case, backend, order in sorted({k[:3] for k in keys}):
        candidates = [
            g
            for g in groups
            if (g["case"], g["backend"], g["order"]) == (case, backend, order)
        ]
        for target in error_targets:
            feasible = [g for g in candidates if g["h1_relative"] <= target]
            accuracy.append(
                {
                    "case": case,
                    "backend": backend,
                    "order": order,
                    "target": target,
                    "selected": min(feasible, key=lambda g: g["wall_s"]["median"])
                    if feasible
                    else None,
                }
            )
        for budget in time_budgets:
            feasible = [
                g
                for g in candidates
                if g["wall_s"]["median"] <= budget
                and g["rss_max_bytes"] is not None
                and g["rss_max_bytes"] <= memory_budget
            ]
            cost.append(
                {
                    "case": case,
                    "backend": backend,
                    "order": order,
                    "seconds": budget,
                    "memory_bytes": memory_budget,
                    "selected": min(feasible, key=lambda g: g["h1_relative"])
                    if feasible
                    else None,
                }
            )
    return {"groups": groups, "accuracy": accuracy, "cost": cost}


def code_hash(root):
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py")):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def coordinate(args):
    import fem_toolbox

    args.output.mkdir(parents=True, exist_ok=False)
    env = {**os.environ, **dict.fromkeys(THREAD_VARS, "1")}
    rows = []
    references = {}
    # Interleave backends and element orders at each size/repeat, run serially.
    for case, n, repeat, order, backend in product(
        ("heat", "bending"), args.sizes, range(args.repeats), (1, 2), ("scipy", "petsc")
    ):
        output = args.output / f"{case}-n{n}-r{repeat}-p{order}-{backend}"
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            case,
            "--size",
            str(n),
            "--order",
            str(order),
            "--backend",
            backend,
            "--output",
            str(output),
        ]
        if args.petsc_scalar:
            command += ["--petsc-scalar", args.petsc_scalar]
        subprocess.run(command, check=True, env=env)
        report = json.loads((output / "report.json").read_text())
        field = np.load(output / "field.npy")
        key = (case, n, order)
        if key not in references:
            references[key] = field
        reference = references[key]
        np.testing.assert_allclose(
            field, reference, rtol=1e-8, atol=1e-10 if case == "heat" else 1e-11
        )
        report["parity_max_absolute"] = float(np.max(np.abs(field - reference)))
        report["repeat"] = repeat
        rows.append(report)
        print(output.name, report["wall_s"], report["h1_relative"], flush=True)
        # Preserve all completed samples if a later worker fails.
        (args.output / "samples.json").write_text(
            json.dumps(rows, indent=2, allow_nan=False) + "\n"
        )
    report = {
        "schema": "fem-toolbox-accuracy-cost-v1",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "logical_cpus": os.cpu_count(),
            "python": platform.python_version(),
            "versions": {n: version(n) for n in ("numpy", "scipy", "meshio")},
            "thread_environment": {n: env[n] for n in THREAD_VARS},
            "numpy_configuration": np.show_config(mode="dicts"),
            "package_path": str(Path(fem_toolbox.__file__).resolve()),
            "production_sha256": code_hash(Path(fem_toolbox.__file__).parent),
            "benchmark_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "mesh_generator_sha256": hashlib.sha256(
                Path(__file__).with_name("profile_baseline.py").read_bytes()
            ).hexdigest(),
            "revision": args.revision,
            "petsc_scalar": args.petsc_scalar,
        },
        "samples": rows,
        **comparison(rows, [0.5, 0.25, 0.125], [1.0, 5.0, 20.0], 512 * 1024**2),
    }
    (args.output / "comparison.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )


def positive(raw):
    number = int(raw)
    if number < 1:
        raise argparse.ArgumentTypeError("expected positive integer")
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", choices=("heat", "bending"))
    parser.add_argument("--size", type=positive, default=2)
    parser.add_argument("--sizes", type=positive, nargs="+", default=[2, 4, 6])
    parser.add_argument("--order", type=int, choices=(1, 2), default=1)
    parser.add_argument("--backend", choices=("scipy", "petsc"), default="scipy")
    parser.add_argument("--petsc-scalar", choices=("real", "complex"))
    parser.add_argument("--repeats", type=positive, default=3)
    parser.add_argument("--revision", default="unspecified")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.worker:
        worker(args)
    else:
        coordinate(args)


if __name__ == "__main__":
    main()
