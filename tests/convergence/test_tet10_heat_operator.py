"""Test-only Tet10 assembly and independent quartic manufactured heat solution.

Exercises the element operators; production Model/solver dispatch is #35d.
"""

from itertools import product

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from fem_toolbox.elements import Tet10Geometry
from fem_toolbox.elements.quadratic import conductance_matrix
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.physics.quadratic import source_load


def elevate(mesh):
    points = mesh.points.tolist()
    edges, cells = {}, []
    for cell in mesh.tetra4:
        mids = []
        for a, b in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)):
            edge = tuple(sorted((int(cell[a]), int(cell[b]))))
            if edge not in edges:
                edges[edge] = len(points)
                points.append(
                    ((mesh.points[edge[0]] + mesh.points[edge[1]]) / 2).tolist()
                )
            mids.append(edges[edge])
        cells.append([*cell, *mids])
    return np.array(points), np.array(cells)


def test_four_mesh_quartic_heat_solution_rates_and_balance(cube_mesh):
    roots, weights = np.polynomial.legendre.leggauss(6)
    roots, weights = (roots + 1) / 2, weights / 2
    bary, quadrature = [], []
    for i, j, k in product(range(6), repeat=3):
        a, b, c = roots[[i, j, k]]
        x, y, z = a, (1 - a) * b, (1 - a) * (1 - b) * c
        bary.append([1 - x - y - z, x, y, z])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - a) ** 2 * (1 - b))
    bary, quadrature = np.array(bary), np.array(quadrature)
    coefficients = np.array([1.0, 0.5, 0.25])
    errors = []
    levels = (1, 2, 3, 4)
    for n in levels:
        points, cells = elevate(cube_mesh(n))
        rows, cols, data = [], [], []
        rhs = np.zeros(len(points))
        geometries = []
        for cell in cells:
            g = Tet10Geometry(points[cell])
            geometries.append(g)
            matrix = conductance_matrix(g, ThermalMaterial.isotropic(2))
            rows.extend(np.repeat(cell, 10))
            cols.extend(np.tile(cell, 10))
            data.extend(matrix.ravel())
            rhs[cell] += source_load(g, lambda p: -24 * (p**2 @ coefficients))
        matrix = sparse.coo_matrix(
            (data, (rows, cols)), shape=(len(points), len(points))
        ).tocsr()
        boundary = np.flatnonzero(np.any((points == 0) | (points == 1), axis=1))
        free = np.setdiff1d(np.arange(len(points)), boundary)
        solution = points**4 @ coefficients
        solution[free] = spsolve(
            matrix[free][:, free],
            rhs[free] - matrix[free][:, boundary] @ solution[boundary],
        )
        residual = matrix @ solution - rhs
        assert np.linalg.norm(residual[free]) < 1e-12
        assert abs(residual[boundary].sum() + rhs.sum()) < 2e-12
        assert abs(rhs.sum() + 14) < 2e-12
        l2, h1 = 0.0, 0.0
        for cell, g in zip(cells, geometries, strict=True):
            xyz = bary @ points[cell[:4]]
            values = np.array([g.shape_values(p) @ solution[cell] for p in xyz])
            gradient = np.array([g.shape_gradients(p).T @ solution[cell] for p in xyz])
            l2 += 6 * g.volume * (quadrature @ (values - xyz**4 @ coefficients) ** 2)
            h1 += (
                6
                * g.volume
                * (
                    quadrature
                    @ np.sum((gradient - 4 * xyz**3 * coefficients) ** 2, axis=1)
                )
            )
        errors.append(np.sqrt([l2, h1]))
    rates = (
        np.log(np.array(errors[:-1]) / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("Tet10 heat L2/H1 errors:", errors, "rates:", rates)
    assert np.all((rates[-2:, 0] > 2.7) & (rates[-2:, 0] < 3.3))
    assert np.all((rates[-2:, 1] > 1.8) & (rates[-2:, 1] < 2.2))
