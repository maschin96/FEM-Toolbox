"""Thermal sources, outward fluxes and linear Robin boundary integration."""

from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError

ScalarField: TypeAlias = float | Callable[[NDArray[np.float64]], float]


def _value(field: ScalarField, point: NDArray[np.float64]) -> float:
    return float(real_array(field(point) if callable(field) else field, ()))


def _triangle(coordinates: ArrayLike) -> tuple[NDArray[np.float64], float]:
    p = real_array(coordinates, (3, 3))
    area = float(np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])) / 2)
    length = max(
        float(np.linalg.norm(p[j] - p[i])) for i, j in ((0, 1), (0, 2), (1, 2))
    )
    if not np.isfinite(area) or area <= 1e-12 * length**2:
        raise ValidationError("degenerate thermal triangle", field="surface")
    return p, area


def source_load(geometry: Tet4Geometry, source: ScalarField) -> NDArray[np.float64]:
    """Integrate Q (W/m³, positive heat generation) to (4,) nodal power in W.

    Callables receive x (3,) in metres. Degree-two four-point quadrature is
    exact for affine sources. Higher variation requires quadrature verification.
    """
    if not callable(source):
        return np.full(4, float(real_array(source, ())) * geometry.volume / 4)
    shapes = np.full((4, 4), (5 - np.sqrt(5)) / 20)
    np.fill_diagonal(shapes, (5 + 3 * np.sqrt(5)) / 20)
    result = np.zeros(4)
    for n in shapes:
        result += n * _value(source, n @ geometry.coordinates) * geometry.volume / 4
    return result


def _surface_load(
    p: NDArray[np.float64], area: float, field: ScalarField
) -> NDArray[np.float64]:
    if not callable(field):
        return np.full(3, float(real_array(field, ())) * area / 3)
    shapes = np.full((3, 3), 1 / 6)
    np.fill_diagonal(shapes, 2 / 3)
    result = np.zeros(3)
    for n in shapes:
        result += n * _value(field, n @ p) * area / 3
    return result


def flux_load(coordinates: ArrayLike, outward_flux: ScalarField) -> NDArray[np.float64]:
    """Return -(integral N q_out) (3,) W; q_out in W/m² is positive OUTWARD.

    Uses degree-two triangle quadrature; exact for affine fluxes. The sign is
    independent of triangle node ordering because q_out is a prescribed scalar.
    """
    p, area = _triangle(coordinates)
    return -_surface_load(p, area, outward_flux)


@dataclass(frozen=True)
class Convection:
    """Robin law q_out=h*(T-T_inf): h >= 0 in W/(m² K), ambient in K.

    h is constant on a facet; ambient is a finite scalar or spatial callable.
    Temperature and ambient use the same reference (absolute K by default).
    """

    coefficient: float
    ambient: ScalarField

    def __post_init__(self) -> None:
        h = float(real_array(self.coefficient, ()))
        if h < 0:
            raise ValidationError(
                "convection coefficient must be nonnegative", field="convection"
            )
        object.__setattr__(self, "coefficient", h)
        if not callable(self.ambient):
            object.__setattr__(self, "ambient", float(real_array(self.ambient, ())))


def convection_terms(
    coordinates: ArrayLike, convection: Convection
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return boundary conductance (3,3) W/K and ambient load (3,) W.

    Constant h gives the exact consistent triangle mass matrix. The load is
    integrated at degree two, exact for affine ambient temperatures.
    """
    p, area = _triangle(coordinates)
    h = convection.coefficient
    matrix = h * area / 12 * (np.ones((3, 3)) + np.eye(3))
    return matrix, h * _surface_load(p, area, convection.ambient)
