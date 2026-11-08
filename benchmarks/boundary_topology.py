"""Issue #287 scalar/vectorized timings in separate processes.

Run: PYTHONPATH=src:. .venv/bin/python benchmarks/boundary_topology.py
Mesh construction and result summarization are outside the measured call.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import subprocess
import sys
from itertools import permutations
from pathlib import Path
from time import perf_counter

import numpy as np

from fem_toolbox.mesh import Mesh, inspect_boundary_topology
from tests.reference.legacy_boundary_topology import (
    inspect_boundary_topology as legacy_inspect,
)


def cube(n: int) -> Mesh:
    grid = np.indices((n + 1,) * 3).reshape(3, -1).T
    points = grid / n
    bases = np.indices((n,) * 3).reshape(3, -1).T
    strides = np.array([(n + 1) ** 2, n + 1, 1])
    offsets = []
    for order in permutations(range(3)):
        steps = np.vstack([np.zeros(3, dtype=int), np.eye(3, dtype=int)[list(order)]])
        offsets.append(steps.cumsum(axis=0) @ strides)
    cells = (bases @ strides)[:, None, None] + np.array(offsets)[None, :, :]
    cells = cells.reshape(-1, 4)
    faces = cells[:, [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]].reshape(-1, 3)
    # An exterior face lies on a cube plane. No production topology is used.
    xyz = points[faces]
    outside = np.any(np.all(xyz == 0, axis=1) | np.all(xyz == 1, axis=1), axis=1)
    return Mesh(points=points, tetra4=cells, triangle3=faces[outside])


def worker(n: int, mode: str) -> dict[str, object]:
    mesh = cube(n)
    inspect = legacy_inspect if mode == "scalar" else inspect_boundary_topology
    start = perf_counter()
    result = inspect(mesh, require_complete=True)
    elapsed = perf_counter() - start
    # Hash discrete fields independently of tolerated floating-point differences.
    digest = hashlib.sha256()
    for key, owners in result.facet_incidence.items():
        digest.update(np.asarray((*key, *owners), dtype=np.int64).tobytes())
    for face in result.exterior_faces:
        digest.update(
            np.asarray(
                (
                    *face.nodes,
                    face.cell_index,
                    face.cell_id,
                    face.opposite_node,
                    face.triangle_index,
                    face.triangle_id,
                ),
                dtype=np.int64,
            ).tobytes()
        )
    return dict(
        n=n,
        cells=mesh.n_cells,
        mode=mode,
        seconds=elapsed,
        peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        peak_rss_unit="bytes" if sys.platform == "darwin" else "KiB",
        incidence_sha256=digest.hexdigest(),
        exterior_faces=len(result.exterior_faces),
        total_area=sum(face.area for face in result.exterior_faces),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", nargs="+", type=int, default=[12, 26, 55])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--worker", choices=["scalar", "vectorized"])
    parser.add_argument(
        "--output", type=Path, default=Path("results/boundary_topology.json")
    )
    args = parser.parse_args()
    if min(args.sizes) <= 0 or args.repeats <= 0:
        parser.error("sizes and repeats must be positive")
    if args.worker:
        print(json.dumps(worker(args.sizes[0], args.worker)))
        return
    report = dict(
        python=sys.version,
        numpy=np.__version__,
        platform=platform.platform(),
        processor=platform.processor(),
        repeats=args.repeats,
        base_commit=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        source_sha256=hashlib.sha256(
            Path("src/fem_toolbox/mesh/topology.py").read_bytes()
        ).hexdigest(),
        oracle_sha256=hashlib.sha256(
            Path("tests/reference/legacy_boundary_topology.py").read_bytes()
        ).hexdigest(),
        samples=[],
    )
    env = {
        **os.environ,
        **dict.fromkeys(
            (
                "OMP_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
                "BLIS_NUM_THREADS",
            ),
            "1",
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for n in args.sizes:
        for _ in range(args.repeats):
            pair = []
            for mode in ("scalar", "vectorized"):
                raw = subprocess.check_output(
                    [sys.executable, __file__, "--worker", mode, "--sizes", str(n)],
                    env=env,
                    text=True,
                )
                sample = json.loads(raw)
                pair.append(sample)
                report["samples"].append(sample)
                print(json.dumps(sample), flush=True)
                args.output.write_text(json.dumps(report, indent=2) + "\n")
            assert pair[0]["incidence_sha256"] == pair[1]["incidence_sha256"]
            assert np.isclose(pair[0]["total_area"], pair[1]["total_area"], rtol=2e-14)


if __name__ == "__main__":
    main()
