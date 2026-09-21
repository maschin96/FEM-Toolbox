"""Consistent mechanical loads; vectors are in global Cartesian axes."""

from collections.abc import Callable
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError

VectorField: TypeAlias = ArrayLike | Callable[[NDArray[np.float64]], ArrayLike]


def body_load(geometry: Tet4Geometry, force: VectorField) -> NDArray[np.float64]:
    """Integrate body force (N/m³) to a (12,) nodal force vector (N).

    Constant forces are exact; callables receive x (3,) in metres. Four-point
    degree-two tetrahedral quadrature is exact for affine force distributions.
    More general spatial loads require a separate quadrature convergence check.
    """
    if not callable(force):
        return np.tile(real_array(force, (3,)) * geometry.volume / 4, 4)
    a = (5 + 3 * np.sqrt(5)) / 20
    b = (5 - np.sqrt(5)) / 20
    shapes = np.full((4, 4), b)
    np.fill_diagonal(shapes, a)
    result = np.zeros((4, 3))
    for shape in shapes:
        result += (
            np.outer(shape, real_array(force(shape @ geometry.coordinates), (3,)))
            * geometry.volume
            / 4
        )
    return result.ravel()


def traction_load(coordinates: ArrayLike, traction: VectorField) -> NDArray[np.float64]:
    """Integrate prescribed global traction (N/m²) on Triangle3 to (9,) N.

    No pressure/normal sign is inferred. Swapping triangle orientation only
    permutes nodal forces. Three-point degree-two quadrature integrates affine
    tractions exactly. Degenerate triangles and invalid vectors are rejected.
    """
    p = real_array(coordinates, (3, 3))
    edges = p[1:] - p[0]
    area = float(np.linalg.norm(np.cross(edges[0], edges[1])) / 2)
    scale = float(
        max(
            np.linalg.norm(edges[0]),
            np.linalg.norm(edges[1]),
            np.linalg.norm(p[2] - p[1]),
        )
    )
    if not np.isfinite(area) or area <= 1e-12 * scale**2:
        raise ValidationError("degenerate traction triangle", field="traction")
    if not callable(traction):
        return np.tile(real_array(traction, (3,)) * area / 3, 3)
    shapes = np.full((3, 3), 1 / 6)
    np.fill_diagonal(shapes, 2 / 3)
    result = np.zeros((3, 3))
    for shape in shapes:
        result += np.outer(shape, real_array(traction(shape @ p), (3,))) * area / 3
    return result.ravel()
