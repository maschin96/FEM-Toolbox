"""Positive reference simplex quadrature by a collapsed Gauss product."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.typing import NDArray

from fem_toolbox.elements.tet4 import _frozen_float_array


@dataclass(frozen=True, slots=True)
class QuadratureRule:
    """Points (n, d), weights (n,) and guaranteed total polynomial degree.

    Reference weights integrate d-dimensional reference measure (1/6 for a
    tetrahedron, 1/2 for a triangle). Geometry methods return physical points
    in metres and weights in m³ or m². Arrays returned by the factory functions
    and geometry methods are immutable. No normalization to unit weight.
    """

    points: NDArray[np.float64]
    weights: NDArray[np.float64]
    degree: int


def _rule(
    points: NDArray[np.float64], weights: NDArray[np.float64], degree: int
) -> QuadratureRule:
    return QuadratureRule(
        _frozen_float_array(points, shape=points.shape),
        _frozen_float_array(weights, shape=weights.shape),
        degree,
    )


def tetrahedron_quadrature() -> QuadratureRule:
    """64 positive points on x,y,z >= 0, x+y+z <= 1; exact through degree 5.

    Four Gauss points per collapsed coordinate integrate degrees up to seven;
    the Jacobian adds at most two powers, giving total simplex degree five.
    This includes degree-four products of quadratic shape functions.
    """
    nodes, weights = leggauss(4)
    nodes, weights = (nodes + 1) / 2, weights / 2
    points, values = [], []
    for i, j, k in product(range(4), repeat=3):
        a, b, c = nodes[i], nodes[j], nodes[k]
        points.append((a, (1 - a) * b, (1 - a) * (1 - b) * c))
        values.append(weights[i] * weights[j] * weights[k] * (1 - a) ** 2 * (1 - b))
    return _rule(np.asarray(points), np.asarray(values), 5)


def triangle_quadrature() -> QuadratureRule:
    """Nine positive points on x,y >= 0, x+y <= 1; exact through degree 4."""
    nodes, weights = leggauss(3)
    nodes, weights = (nodes + 1) / 2, weights / 2
    points, values = [], []
    for i, j in product(range(3), repeat=2):
        a, b = nodes[i], nodes[j]
        points.append((a, (1 - a) * b))
        values.append(weights[i] * weights[j] * (1 - a))
    return _rule(np.asarray(points), np.asarray(values), 4)
