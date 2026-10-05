"""Deterministic structured verification meshes, independent of production code."""

from collections import Counter
from itertools import permutations, product

import numpy as np
import pytest

from femantis.mesh import Mesh


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
