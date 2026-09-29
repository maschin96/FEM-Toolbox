"""Internal static Tet4/Tet10 piezo dispatch with integrated enthalpy terms."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet10Geometry, quadratic
from fem_toolbox.elements import piezoelectric as linear
from fem_toolbox.elements.piezoelectric import PiezoelectricElementResponse
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial, ThermalExpansion
from fem_toolbox.solvers.mechanical_elements import MechanicalGeometry


def piezoelectric_matrix(
    geometry: MechanicalGeometry, material: PiezoelectricMaterial
) -> NDArray[np.float64]:
    if isinstance(geometry, Tet10Geometry):
        return quadratic.piezoelectric_matrix(geometry, material)
    return linear.piezoelectric_matrix(geometry, material)


def piezoelectric_response(
    geometry: MechanicalGeometry,
    material: PiezoelectricMaterial,
    displacement: ArrayLike,
    potential: ArrayLike,
    *,
    thermal_strain: ArrayLike | None = None,
    temperature: ArrayLike | None = None,
    expansion: ThermalExpansion | None = None,
    electric_thermal_strain: bool = False,
) -> PiezoelectricElementResponse:
    if not isinstance(geometry, Tet10Geometry):
        return linear.piezoelectric_response(
            geometry,
            material,
            displacement,
            potential,
            thermal_strain=thermal_strain,
            electric_thermal_strain=electric_thermal_strain,
        )
    if (temperature is None) != (expansion is None):
        raise ValidationError(
            "temperature and expansion must be supplied together", field="thermal"
        )
    delta = (
        None
        if temperature is None or expansion is None
        else real_array(temperature, (10,)) - expansion.reference_temperature
    )

    def eigenstrain(point: NDArray[np.float64]) -> NDArray[np.float64]:
        if delta is None or expansion is None:
            return np.zeros(6)
        return expansion.alpha * (geometry.shape_values(point) @ delta)

    u = real_array(displacement, (10, 3))
    u = (u - u[0]).ravel()
    phi = real_array(potential, (10,))
    phi = phi - phi[0]
    center = geometry.coordinates[:4].mean(axis=0)
    strain = quadratic.strain_matrix(geometry, center) @ u
    field = -geometry.shape_gradients(center).T @ phi
    mechanical = strain - eigenstrain(center)
    coupling = mechanical if electric_thermal_strain else strain
    stress = material.stiffness @ mechanical - material.stress_coupling.T @ field
    electric_displacement = (
        material.stress_coupling @ coupling + material.permittivity @ field
    )
    energies = np.zeros(3)
    rule = geometry.quadrature()
    for point, weight in zip(rule.points, rule.weights, strict=True):
        eps = quadratic.strain_matrix(geometry, point) @ u
        mechanical = eps - eigenstrain(point)
        coupling = mechanical if electric_thermal_strain else eps
        electric = -geometry.shape_gradients(point).T @ phi
        energies += weight * np.array(
            [
                0.5 * mechanical @ material.stiffness @ mechanical,
                -coupling @ material.stress_coupling.T @ electric,
                -0.5 * electric @ material.permittivity @ electric,
            ]
        )
    return PiezoelectricElementResponse(
        strain,
        field,
        stress,
        electric_displacement,
        float(energies[0]),
        float(energies[1]),
        float(energies[2]),
        float(energies.sum()),
    )
