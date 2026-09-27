"""Internal static mechanical dispatch and positive Tet10 strain-energy integral."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.elements import (
    Tet4Geometry,
    Tet10Geometry,
    Triangle6Geometry,
    quadratic,
)
from fem_toolbox.elements import elasticity as linear
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion
from fem_toolbox.physics import elasticity as loads
from fem_toolbox.physics import quadratic as quadratic_loads

MechanicalGeometry = Tet4Geometry | Tet10Geometry
FloatArray = NDArray[np.float64]


def stiffness_matrix(g: MechanicalGeometry, material: ElasticMaterial) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic.stiffness_matrix(g, material)
    return linear.stiffness_matrix(g, material)


def strain_matrix(g: MechanicalGeometry) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic.strain_matrix(g, g.coordinates[:4].mean(axis=0))
    return linear.strain_matrix(g)


def body_load(g: MechanicalGeometry, force: loads.VectorField) -> FloatArray:
    if isinstance(g, Tet10Geometry):
        return quadratic_loads.body_load(g, force)
    return loads.body_load(g, force)


def traction_load(points: FloatArray, traction: loads.VectorField) -> FloatArray:
    if len(points) == 6:
        return quadratic_loads.traction_load(Triangle6Geometry(points), traction)
    return loads.traction_load(points, traction)


def quadratic_energy(
    g: Tet10Geometry,
    material: ElasticMaterial,
    displacement: FloatArray,
    *,
    temperature: FloatArray | None = None,
    expansion: ThermalExpansion | None = None,
) -> float:
    # Subtract translation without changing the displacement gradient.
    values = (displacement.reshape(-1, 3) - displacement[:3]).ravel()
    rule = g.quadrature()
    total = 0.0
    for point, weight in zip(rule.points, rule.weights, strict=True):
        strain = quadratic.strain_matrix(g, point) @ values
        if temperature is not None and expansion is not None:
            strain -= expansion.alpha * (
                g.shape_values(point) @ (temperature - expansion.reference_temperature)
            )
        total += 0.5 * weight * float(strain @ material.stiffness @ strain)
    return total


def thermal_terms(
    g: MechanicalGeometry,
    material: ElasticMaterial,
    temperature: FloatArray,
    expansion: ThermalExpansion,
) -> tuple[FloatArray, FloatArray]:
    """Return reported eigenstrain and consistent load; Tet10 report is a sample."""
    if isinstance(g, Tet10Geometry):
        eigenstrain = expansion.alpha * (
            g.shape_values(g.coordinates[:4].mean(axis=0))
            @ (temperature - expansion.reference_temperature)
        )
        return eigenstrain, quadratic.thermal_load(g, material, temperature, expansion)
    eigenstrain = expansion.alpha * (
        float(np.mean(temperature)) - expansion.reference_temperature
    )
    load = g.volume * linear.strain_matrix(g).T @ material.stiffness @ eigenstrain
    return eigenstrain, load
