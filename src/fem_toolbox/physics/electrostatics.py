"""Charge integrals for -K_phi phi = surface(D.n) - volume(rho)."""

from collections.abc import Callable
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError

ChargeDensity: TypeAlias = float | Callable[[NDArray[np.float64]], float]


def volume_charge_load(
    geometry: Tet4Geometry, density: ChargeDensity
) -> NDArray[np.float64]:
    """Return -integral N rho (4,) C; rho in C/m³, callable x (3,) in m.

    Four-point degree-two quadrature integrates affine density exactly.
    """
    if not callable(density):
        return np.full(4, -float(real_array(density, ())) * geometry.volume / 4)
    shapes = np.full((4, 4), (5 - np.sqrt(5)) / 20)
    np.fill_diagonal(shapes, (5 + 3 * np.sqrt(5)) / 20)
    result = np.zeros(4)
    for n in shapes:
        result -= (
            n
            * float(real_array(density(n @ geometry.coordinates), ()))
            * geometry.volume
            / 4
        )
    return result


def surface_charge_load(
    coordinates: ArrayLike, outward_displacement: ChargeDensity
) -> NDArray[np.float64]:
    """Return +integral N omega (3,) C; omega = D.n_out in C/m².

    omega is the outward dielectric displacement, not the opposite-signed free
    charge on an adjacent conductor. Three-point degree-two triangle quadrature
    integrates affine density exactly, independent of triangle orientation.
    """
    p = real_array(coordinates, (3, 3))
    area = float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])) / 2)
    length = max(
        float(np.linalg.norm(p[j] - p[i])) for i, j in ((0, 1), (0, 2), (1, 2))
    )
    if not np.isfinite(area) or area <= 1e-12 * length**2:
        raise ValidationError(
            "degenerate electrical boundary triangle", field="surface"
        )
    if not callable(outward_displacement):
        return np.full(3, float(real_array(outward_displacement, ())) * area / 3)
    shapes = np.full((3, 3), 1 / 6)
    np.fill_diagonal(shapes, 2 / 3)
    result = np.zeros(3)
    for n in shapes:
        result += n * float(real_array(outward_displacement(n @ p), ())) * area / 3
    return result
