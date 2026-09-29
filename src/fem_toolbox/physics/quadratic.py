"""Tet10 and Triangle6 source/boundary integrals with existing SI field signs.

Materials are constant; spatial load callables receive a physical (3,) point in
metres. Volume rules integrate degree-five polynomials, surface rules degree four.
Thus loads are exact through degree three in volume and degree two on surfaces.
More general callables require separate quadrature convergence verification.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet10Geometry, Triangle6Geometry
from fem_toolbox.errors import ValidationError
from fem_toolbox.physics.elasticity import VectorField
from fem_toolbox.physics.heat import Convection, ScalarField

FloatArray = NDArray[np.float64]


def _validate(geometry: object, expected: type) -> None:
    if not isinstance(geometry, expected):
        raise ValidationError(f"expected {expected.__name__}", field="geometry")


def _scalar(
    geometry: Tet10Geometry | Triangle6Geometry, field: ScalarField
) -> FloatArray:
    result = np.zeros(len(geometry.coordinates))
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        value = float(real_array(field(point) if callable(field) else field, ()))
        result += weight * value * geometry.shape_values(point)
    return result


def _vector(
    geometry: Tet10Geometry | Triangle6Geometry, field: VectorField
) -> FloatArray:
    result = np.zeros((len(geometry.coordinates), 3))
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        value = real_array(field(point) if callable(field) else field, (3,))
        result += weight * np.outer(geometry.shape_values(point), value)
    return result.ravel()


def source_load(geometry: Tet10Geometry, source: ScalarField) -> FloatArray:
    """Integral N Q (10,) W; positive generation Q in W/m³."""
    _validate(geometry, Tet10Geometry)
    return _scalar(geometry, source)


def body_load(geometry: Tet10Geometry, force: VectorField) -> FloatArray:
    """Integral N body force (30,) N; Cartesian force density in N/m³."""
    _validate(geometry, Tet10Geometry)
    return _vector(geometry, force)


def volume_charge_load(geometry: Tet10Geometry, density: ScalarField) -> FloatArray:
    """Negative integral N rho (10,) C; charge density in C/m³."""
    _validate(geometry, Tet10Geometry)
    return -_scalar(geometry, density)


def flux_load(geometry: Triangle6Geometry, outward_flux: ScalarField) -> FloatArray:
    """Negative integral N q_out (6,) W; outward heat flux in W/m²."""
    _validate(geometry, Triangle6Geometry)
    return -_scalar(geometry, outward_flux)


def surface_charge_load(
    geometry: Triangle6Geometry, outward_displacement: ScalarField
) -> FloatArray:
    """Positive integral N (D.n_out) (6,) C; input in C/m².

    Adjacent conductor charge has the opposite sign. Orientation never changes
    the sign of a supplied scalar boundary datum.
    """
    _validate(geometry, Triangle6Geometry)
    return _scalar(geometry, outward_displacement)


def traction_load(geometry: Triangle6Geometry, traction: VectorField) -> FloatArray:
    """Integral N prescribed Cartesian traction (18,) N; datum in N/m²."""
    _validate(geometry, Triangle6Geometry)
    return _vector(geometry, traction)


def convection_terms(
    geometry: Triangle6Geometry, convection: Convection
) -> tuple[FloatArray, FloatArray]:
    """Boundary conductance (6,6) W/K and ambient load (6,) W, h>=0.

    Existing Convection carries constant h [W/(m² K)] and ambient K. Matrix is
    exactly integrated; ambient polynomials through degree two are exact.
    """
    _validate(geometry, Triangle6Geometry)
    _validate(convection, Convection)
    matrix = np.zeros((6, 6))
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        n = geometry.shape_values(point)
        matrix += weight * convection.coefficient * np.outer(n, n)
    return matrix, convection.coefficient * _scalar(geometry, convection.ambient)
