"""Internal thermal element dispatch; Tet4 formulas remain the reference path."""

from __future__ import annotations

from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import NDArray

from femantis.elements import (
    Tet4Geometry,
    Tet10Geometry,
    Triangle6Geometry,
    heat,
    quadratic,
    storage,
)
from femantis.errors import ValidationError
from femantis.materials import ThermalMaterial
from femantis.physics import heat as loads
from femantis.physics import quadratic as quadratic_loads

ThermalGeometry: TypeAlias = Tet4Geometry | Tet10Geometry
Lumping: TypeAlias = Literal["consistent", "hrz"]
FloatArray: TypeAlias = NDArray[np.float64]


def conductance_matrix(g: ThermalGeometry, material: ThermalMaterial) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic.conductance_matrix(g, material)
    return heat.conductance_matrix(g, material)


def capacity_matrix(
    g: ThermalGeometry,
    density: float,
    specific_heat: float,
    *,
    lumped: bool = False,
    lumping: Lumping | None = None,
) -> FloatArray:
    if lumping is not None and lumping not in ("consistent", "hrz"):
        raise ValidationError(
            "invalid lumping; choose consistent or hrz", field="lumping"
        )
    if isinstance(g, Tet10Geometry):
        if lumped:
            raise ValidationError(
                "Tet10 lumped=True is invalid; choose lumping='hrz'", field="lumping"
            )
        return quadratic.capacity_matrix(
            g, density, specific_heat, lumping=lumping or "consistent"
        )
    if lumping is not None:
        raise ValidationError(
            "explicit lumping is only supported for Tet10", field="lumping"
        )
    return storage.capacity_matrix(g, density, specific_heat, lumped=lumped)


def source_load(g: ThermalGeometry, source: loads.ScalarField) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic_loads.source_load(g, source)
    return loads.source_load(g, source)


def flux_load(points: FloatArray, flux: loads.ScalarField) -> FloatArray:
    if len(points) == 6:
        return quadratic_loads.flux_load(Triangle6Geometry(points), flux)
    return loads.flux_load(points, flux)


def convection_terms(
    points: FloatArray, condition: loads.Convection
) -> tuple[FloatArray, FloatArray]:
    if len(points) == 6:
        return quadratic_loads.convection_terms(Triangle6Geometry(points), condition)
    return loads.convection_terms(points, condition)


def temperature_gradient(g: ThermalGeometry, temperature: FloatArray) -> FloatArray:
    values = temperature - temperature[0]
    if isinstance(g, Tet10Geometry):
        return np.asarray(g.shape_gradients(g.coordinates[:4].mean(axis=0)).T @ values)
    return np.asarray(g.gradient_matrix @ values)


def contains(g: ThermalGeometry, point: FloatArray) -> bool:
    if isinstance(g, Tet10Geometry):
        corners = g.coordinates[:4]
        coordinates = np.linalg.solve((corners[1:] - corners[0]).T, point - corners[0])
        weights = np.r_[1 - coordinates.sum(), coordinates]
    else:
        weights = g.shape_values(point)
    return bool(np.all(weights >= -1e-12) and np.all(weights <= 1 + 1e-12))
