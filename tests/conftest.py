"""Deterministic structured verification meshes, independent of production code."""

from collections import Counter
from itertools import permutations, product

import numpy as np
import pytest

from fem_toolbox.cfd import FluidMesh2D
from fem_toolbox.mesh import Mesh


@pytest.fixture
def quad_mesh():
    """Authored rectangular/sheared FV meshes; no production mesh generator."""

    def build(nx=4, ny=4, *, shear=0.0):
        points = np.array(
            [
                [i / nx + shear * j / ny, j / ny]
                for j in range(ny + 1)
                for i in range(nx + 1)
            ]
        )

        def node(i, j):
            return j * (nx + 1) + i

        cells = [
            [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)]
            for j in range(ny)
            for i in range(nx)
        ]
        edges = [[node(i, 0), node(i + 1, 0)] for i in range(nx)]
        edges += [[node(nx, j), node(nx, j + 1)] for j in range(ny)]
        edges += [[node(i + 1, ny), node(i, ny)] for i in range(nx)]
        edges += [[node(0, j + 1), node(0, j)] for j in range(ny)]
        return FluidMesh2D(
            points,
            cells,
            edges,
            {
                "bottom": np.arange(nx),
                "right": np.arange(nx, nx + ny),
                "top": np.arange(nx + ny, 2 * nx + ny),
                "left": np.arange(2 * nx + ny, 2 * (nx + ny)),
            },
        )

    return build


@pytest.fixture
def cube_mesh():
    def build(n, *, distort=False):
        points = np.array(list(product(np.linspace(0, 1, n + 1), repeat=3)))
        if distort:
            interior = np.all((points > 0) & (points < 1), axis=1)
            points[interior] += np.random.default_rng(12).uniform(
                -0.12 / n, 0.12 / n, (interior.sum(), 3)
            )

        def index(p):
            return (p[0] * (n + 1) + p[1]) * (n + 1) + p[2]

        cells = []
        for base in product(range(n), repeat=3):
            for order in permutations(range(3)):
                p = np.array(base)
                cell = [index(p)]
                for axis in order:
                    p = p.copy()
                    p[axis] += 1
                    cell.append(index(p))
                cells.append(cell)
        cells = np.array(cells)
        faces = Counter(
            tuple(sorted(np.delete(cell, opposite)))
            for cell in cells
            for opposite in range(4)
        )
        boundary = np.array([face for face, count in faces.items() if count == 1])
        return Mesh(points=points, tetra4=cells, triangle3=boundary)

    return build
