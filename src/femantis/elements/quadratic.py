"""Constant-material, affine Tet10 operators; no global solver dispatch.

Node-major mechanics uses 30 DOFs, scalar fields 10, piezo [u(30), phi(10)].
Quadrature is exact for these polynomial integrands. All arrays retain input
node ordering. SI units and field signs follow the Tet4 contracts.
"""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.assembly import real_array
from femantis.elements import Tet10Geometry
from femantis.elements.storage import _positive
from femantis.errors import ValidationError
from femantis.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
)

FloatArray = NDArray[np.float64]
Lumping = Literal["consistent", "hrz"]


def _validate(geometry: Tet10Geometry, material: object, expected: type) -> None:
    if not isinstance(geometry, Tet10Geometry):
        raise ValidationError("expected Tet10Geometry", field="geometry")
    if not isinstance(material, expected):
        raise ValidationError(f"expected {expected.__name__}", field="material")


def _strain(gradients: FloatArray) -> FloatArray:
    b = np.zeros((6, 30))
    for node, (x, y, z) in enumerate(gradients):
        b[:, 3 * node : 3 * node + 3] = [
            [x, 0, 0],
            [0, y, 0],
            [0, 0, z],
            [0, z, y],
            [z, 0, x],
            [y, x, 0],
        ]
    return b


def strain_matrix(geometry: Tet10Geometry, point: ArrayLike) -> FloatArray:
    """B (6,30) in 1/m at physical point (3,), engineering Voigt yz,xz,xy."""
    _validate(geometry, geometry, Tet10Geometry)
    return _strain(geometry.shape_gradients(point))


def _scalar_matrix(geometry: Tet10Geometry, tensor: FloatArray) -> FloatArray:
    matrix = np.zeros((10, 10))
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        g = geometry.gradient_matrix(point)
        matrix += weight * (g.T @ tensor @ g)
    return matrix


def conductance_matrix(
    geometry: Tet10Geometry, material: ThermalMaterial
) -> FloatArray:
    """Integral G.T k G (10,10), W/K; constant anisotropic conductivity."""
    _validate(geometry, material, ThermalMaterial)
    return _scalar_matrix(geometry, material.conductivity)


def dielectric_matrix(
    geometry: Tet10Geometry, material: DielectricMaterial
) -> FloatArray:
    """Positive integral G.T kappa G (10,10), F; piezo uses its negative."""
    _validate(geometry, material, DielectricMaterial)
    return _scalar_matrix(geometry, material.permittivity)


def stiffness_matrix(geometry: Tet10Geometry, material: ElasticMaterial) -> FloatArray:
    """Integral B.T C B (30,30), N/m; constant engineering-Voigt stiffness."""
    _validate(geometry, material, ElasticMaterial)
    matrix = np.zeros((30, 30))
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        b = strain_matrix(geometry, point)
        matrix += weight * (b.T @ material.stiffness @ b)
    return matrix


def piezoelectric_matrix(
    geometry: Tet10Geometry, material: PiezoelectricMaterial
) -> FloatArray:
    """Signed symmetric tangent (40,40), [u,phi], with +B.T e.T G coupling.

    Mechanical block is integral B.T C B, dielectric block -G.T kappa G.
    Units: N/m, N/V (C/m), C/V in the corresponding blocks. E=-grad(phi).
    """
    _validate(geometry, material, PiezoelectricMaterial)
    mechanical, coupling, dielectric = (
        np.zeros((30, 30)),
        np.zeros((30, 10)),
        np.zeros((10, 10)),
    )
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        gradients = geometry.shape_gradients(point)
        b, g = _strain(gradients), gradients.T
        mechanical += weight * (b.T @ material.stiffness @ b)
        coupling += weight * (b.T @ material.stress_coupling.T @ g)
        dielectric -= weight * (g.T @ material.permittivity @ g)
    return np.block([[mechanical, coupling], [coupling.T, dielectric]])


def thermal_load(
    geometry: Tet10Geometry,
    material: ElasticMaterial | PiezoelectricMaterial,
    temperature: ArrayLike,
    expansion: ThermalExpansion,
    *,
    electric_thermal_strain: bool = False,
) -> FloatArray:
    """Thermal RHS (30,) N or (40,) [N,C] from nodal temperature (10,) K.

    Integrates B.T C alpha*(N.T-Tref). For piezo, the explicit electric switch
    adds G.T e alpha*(N.T-Tref), matching mechanical-strain D convention.
    Default D uses total strain. No pyroelectricity or thermal feedback.
    Constant expansion in engineering Voigt and reference temperature are SI.
    Invalid geometry, material, temperature, expansion or switch raise ValidationError.
    """
    _validate(geometry, geometry, Tet10Geometry)
    if not isinstance(
        material, (ElasticMaterial, PiezoelectricMaterial)
    ) or not isinstance(expansion, ThermalExpansion):
        raise ValidationError(
            "expected elastic/piezo material and ThermalExpansion", field="thermal"
        )
    if not isinstance(electric_thermal_strain, bool) or (
        electric_thermal_strain and not isinstance(material, PiezoelectricMaterial)
    ):
        raise ValidationError(
            "electric thermal strain requires piezo material and boolean switch",
            field="thermal",
        )
    values = real_array(temperature, (10,)) - expansion.reference_temperature
    result = np.zeros(40 if isinstance(material, PiezoelectricMaterial) else 30)
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        thermal = expansion.alpha * (geometry.shape_values(point) @ values)
        gradients = geometry.shape_gradients(point)
        result[:30] += weight * (_strain(gradients).T @ material.stiffness @ thermal)
        if electric_thermal_strain and isinstance(material, PiezoelectricMaterial):
            result[30:] += weight * (gradients @ material.stress_coupling @ thermal)
    return result


def _storage(
    geometry: Tet10Geometry, coefficient: float, lumping: Lumping
) -> FloatArray:
    _validate(geometry, geometry, Tet10Geometry)
    if lumping not in ("consistent", "hrz"):
        raise ValidationError(
            "lumping must be 'consistent' or 'hrz'; row-sum Tet10 masses are negative",
            field="lumping",
        )
    total = coefficient * geometry.volume
    if not np.isfinite(total) or total <= 0:
        raise ValidationError("storage total is not representable", field="storage")
    matrix = np.zeros((10, 10))
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        n = geometry.shape_values(point)
        matrix += (coefficient * weight) * np.outer(n, n)
    if not np.isfinite(matrix).all() or np.any(np.diag(matrix) <= 0):
        raise ValidationError(
            "storage coefficients underflow or overflow", field="storage"
        )
    if lumping == "hrz":
        diagonal = np.diag(matrix)
        matrix = np.diag((diagonal / diagonal.sum()) * total)
        if np.any(np.diag(matrix) <= 0):
            raise ValidationError(
                "lumped storage coefficients underflow", field="storage"
            )
    return matrix


def mass_matrix(
    geometry: Tet10Geometry, density: float, *, lumping: Lumping = "consistent"
) -> FloatArray:
    """Mechanical (30,30) mass in kg; positive rho [kg/m³], no electric mass.

    Default: exact rho Ni Nj tensor I3. 'hrz' scales the positive consistent
    diagonal to conserve total mass per direction; not P1 row-sum lumping.
    """
    return np.asarray(
        np.kron(_storage(geometry, _positive(density, "density"), lumping), np.eye(3)),
        dtype=np.float64,
    )


def capacity_matrix(
    geometry: Tet10Geometry,
    density: float,
    specific_heat: float,
    *,
    lumping: Lumping = "consistent",
) -> FloatArray:
    """(10,10) capacity in J/K; positive rho [kg/m³], c [J/(kg K)].

    'consistent' integrates rho c Ni Nj; explicit 'hrz' conserves total capacity
    by scaling its positive diagonal. Invalid or unrepresentable values fail.
    """
    coefficient = _positive(density, "density") * _positive(
        specific_heat, "specific_heat"
    )
    return _storage(geometry, coefficient, lumping)
