"""Four-mesh interpolation convergence, independently overintegrated (no PDE solve)."""

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet10Geometry


def test_cubic_interpolation_has_l2_order_three_and_h1_order_two(cube_mesh):
    # Independent 5-point product integrates degree-six squared interpolation
    # error exactly, including the collapsed Jacobian (degree <= 8 per axis).
    roots, weights = np.polynomial.legendre.leggauss(5)
    roots, weights = (roots + 1) / 2, weights / 2
    bary, quadrature = [], []
    for i, j, k in product(range(5), repeat=3):
        a, b, c = roots[[i, j, k]]
        x, y, z = a, (1 - a) * b, (1 - a) * (1 - b) * c
        bary.append([1 - x - y - z, x, y, z])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - a) ** 2 * (1 - b))
    bary, quadrature = np.array(bary), np.array(quadrature)
    errors = []
    levels = (1, 2, 3, 4)
    coefficients = np.array([1.0, 0.5, 0.25])
    for n in levels:
        mesh = cube_mesh(n)
        l2, h1 = 0.0, 0.0
        for cell in mesh.tetra4:
            corners = mesh.points[cell]
            nodes = np.vstack(
                [
                    corners,
                    [
                        (corners[a] + corners[b]) / 2
                        for a, b in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
                    ],
                ]
            )
            g = Tet10Geometry(nodes)
            nodal = nodes**3 @ coefficients
            xyz = bary @ corners
            values = np.array([g.shape_values(p) @ nodal for p in xyz])
            gradients = np.array([g.shape_gradients(p).T @ nodal for p in xyz])
            l2 += 6 * g.volume * (quadrature @ (values - xyz**3 @ coefficients) ** 2)
            h1 += (
                6
                * g.volume
                * (
                    quadrature
                    @ np.sum((gradients - 3 * xyz**2 * coefficients) ** 2, axis=1)
                )
            )
        errors.append(np.sqrt([l2, h1]))
    rates = (
        np.log(np.array(errors[:-1]) / errors[1:])
        / np.log(np.array(levels[1:]) / levels[:-1])[:, None]
    )
    print("Tet10 interpolation L2/H1 errors:", errors)
    print("rates:", rates)
    # Cubic data give exact scaling on these self-similar affine meshes;
    # 1e-9 accommodates accumulated floating point cancellation, not PDE error.
    assert_allclose(rates, np.tile([3.0, 2.0], (3, 1)), rtol=0, atol=1e-9)
