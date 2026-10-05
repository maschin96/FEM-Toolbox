"""Internal electric dispatch; integrated energy is distinct from field samples."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from femantis.elements import (
    Tet4Geometry,
    Tet10Geometry,
    Triangle6Geometry,
    quadratic,
)
from femantis.elements import electrostatics as linear
from femantis.materials import DielectricMaterial
from femantis.physics import electrostatics as loads
from femantis.physics import quadratic as quadratic_loads

ElectricGeometry = Tet4Geometry | Tet10Geometry
FloatArray = NDArray[np.float64]


def dielectric_matrix(g: ElectricGeometry, material: DielectricMaterial) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic.dielectric_matrix(g, material)
    return linear.dielectric_matrix(g, material)


def volume_charge_load(g: ElectricGeometry, charge: loads.ChargeDensity) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic_loads.volume_charge_load(g, charge)
    return loads.volume_charge_load(g, charge)


def surface_charge_load(points: FloatArray, charge: loads.ChargeDensity) -> FloatArray:
    if len(points) == 6:
        return quadratic_loads.surface_charge_load(Triangle6Geometry(points), charge)
    return loads.surface_charge_load(points, charge)


def electric_field(g: ElectricGeometry, potential: FloatArray) -> FloatArray:
    values = potential - potential[0]
    if isinstance(g, Tet10Geometry):
        return np.asarray(-g.shape_gradients(g.coordinates[:4].mean(axis=0)).T @ values)
    return np.asarray(-g.gradient_matrix @ values)


def stored_energy(
    g: ElectricGeometry, material: DielectricMaterial, potential: FloatArray
) -> float:
    if isinstance(g, Tet10Geometry):
        values = potential - potential[0]
        rule = g.quadrature()
        total = 0.0
        for point, weight in zip(rule.points, rule.weights, strict=True):
            field = -g.shape_gradients(point).T @ values
            total += 0.5 * weight * float(field @ material.permittivity @ field)
        return total
    field = electric_field(g, potential)
    return 0.5 * g.volume * float(field @ material.permittivity @ field)
