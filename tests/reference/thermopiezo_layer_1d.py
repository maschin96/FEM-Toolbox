"""Independent continuum reduction for a two-layer thermo-piezoelectric bar.

The 3D cube has unit cross-section. The exact solution is invariant in x/y;
traction, electric displacement and heat flux point in z. This module uses
only constitutive algebra and one-dimensional integration, no FEM operators.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

T_BOTTOM = 310.0  # K
T_TOP = 300.0  # K
T_REFERENCE = 300.0  # K
HEAT_SOURCE = 8.0  # W/m^3
BODY_FORCE = 0.5  # N/m^3, +z
VOLUME_CHARGE = 0.2  # C/m^3
VOLTAGE = 0.7  # V, top minus bottom
TOP_CHARGE = 0.3  # C, adjacent conductor


@dataclass(frozen=True)
class Layer:
    conductivity: float  # W/(m K)
    stiffness: float  # Pa, C^E_zz
    coupling: float  # C/m^2, e_zz
    permittivity: float  # F/m, kappa^S_zz
    expansion: float  # 1/K, alpha_zz


LOWER = Layer(2.0, 12.0, 1.5, 2.5, 0.001)
UPPER = Layer(4.0, 18.0, 2.0, 3.0, 0.0015)


def layer(z: float) -> Layer:
    return LOWER if z < 0.5 else UPPER


def heat_flux(z: float) -> float:
    return _q0() + HEAT_SOURCE * z


def _q0() -> float:
    resistance = 0.5 / LOWER.conductivity + 0.5 / UPPER.conductivity
    first_moment = 0.125 / LOWER.conductivity + 0.375 / UPPER.conductivity
    return (T_BOTTOM - T_TOP - HEAT_SOURCE * first_moment) / resistance


def temperature(z: float | np.ndarray) -> float | np.ndarray:
    z_array = np.asarray(z)
    middle = T_BOTTOM - (_q0() * 0.5 + HEAT_SOURCE * 0.125) / LOWER.conductivity
    below = (
        T_BOTTOM
        - (_q0() * z_array + 0.5 * HEAT_SOURCE * z_array**2) / LOWER.conductivity
    )
    above = (
        middle
        - (_q0() * (z_array - 0.5) + 0.5 * HEAT_SOURCE * (z_array**2 - 0.25))
        / UPPER.conductivity
    )
    return np.where(z_array <= 0.5, below, above)


def temperature_coefficients(material: Layer) -> tuple[float, float, float]:
    slope = -_q0() / material.conductivity
    curvature = -HEAT_SOURCE / (2 * material.conductivity)
    if material is LOWER:
        return T_BOTTOM, slope, curvature
    middle = float(temperature(0.5))
    constant = middle - slope * 0.5 - curvature * 0.25
    return constant, slope, curvature


def inverse_coefficients(material: Layer) -> tuple[float, float, float]:
    denominator = material.stiffness * material.permittivity + material.coupling**2
    return (
        material.permittivity / denominator,
        material.coupling / denominator,
        -material.stiffness / denominator,
    )


def integrate_piecewise(function, end: float = 1.0) -> float:
    """Five-point Gauss integration, exact for all polynomials used here."""
    roots, weights = np.polynomial.legendre.leggauss(5)
    total = 0.0
    for start, stop in ((0.0, min(end, 0.5)), (0.5, end)):
        if stop <= start:
            continue
        z = start + 0.5 * (roots + 1) * (stop - start)
        total += (
            0.5
            * (stop - start)
            * sum(
                weight * function(float(point))
                for point, weight in zip(z, weights, strict=True)
            )
        )
    return float(total)


@dataclass(frozen=True)
class Reference:
    stress_at_bottom: float
    voltage: float
    electric_displacement_at_bottom: float

    def stress(self, z: float) -> float:
        return self.stress_at_bottom - BODY_FORCE * z

    def electric_displacement(self, z: float) -> float:
        return self.electric_displacement_at_bottom + VOLUME_CHARGE * z

    def derivatives(self, z: float) -> tuple[float, float]:
        material = layer(z)
        a, b, g = inverse_coefficients(material)
        sigma = self.stress(z)
        displacement = self.electric_displacement(z)
        u_prime = (
            material.expansion * (float(temperature(z)) - T_REFERENCE)
            + a * sigma
            + b * displacement
        )
        phi_prime = b * sigma + g * displacement
        return u_prime, phi_prime

    def _polynomials(self, material: Layer) -> tuple[np.ndarray, np.ndarray]:
        constant, slope, curvature = temperature_coefficients(material)
        a, b, g = inverse_coefficients(material)
        u_prime = np.array(
            [
                material.expansion * (constant - T_REFERENCE)
                + a * self.stress_at_bottom
                + b * self.electric_displacement_at_bottom,
                material.expansion * slope - a * BODY_FORCE + b * VOLUME_CHARGE,
                material.expansion * curvature,
            ]
        )
        phi_prime = np.array(
            [
                b * self.stress_at_bottom + g * self.electric_displacement_at_bottom,
                -b * BODY_FORCE + g * VOLUME_CHARGE,
                0.0,
            ]
        )
        return u_prime, phi_prime

    def _integrated_field(
        self, z: float | np.ndarray, field: int
    ) -> float | np.ndarray:
        argument = np.asarray(z)
        lower = self._polynomials(LOWER)[field]
        upper = self._polynomials(UPPER)[field]

        def primitive(
            coefficient: np.ndarray, value: np.ndarray | float
        ) -> np.ndarray | float:
            return (
                coefficient[0] * value
                + coefficient[1] * value**2 / 2
                + coefficient[2] * value**3 / 3
            )

        result = np.where(
            argument <= 0.5,
            primitive(lower, argument),
            primitive(lower, 0.5) + primitive(upper, argument) - primitive(upper, 0.5),
        )
        return float(result) if result.ndim == 0 else result

    def displacement(self, z: float) -> float:
        return float(self._integrated_field(z, 0))

    def potential(self, z: float) -> float:
        return float(self._integrated_field(z, 1))

    def displacement_array(self, z: np.ndarray) -> np.ndarray:
        return np.asarray(self._integrated_field(z, 0))

    def potential_array(self, z: np.ndarray) -> np.ndarray:
        return np.asarray(self._integrated_field(z, 1))

    def derivatives_array(self, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lower_u, lower_phi = self._polynomials(LOWER)
        upper_u, upper_phi = self._polynomials(UPPER)

        def evaluate(coefficient: np.ndarray) -> np.ndarray:
            return coefficient[0] + coefficient[1] * z + coefficient[2] * z**2

        return (
            np.where(z < 0.5, evaluate(lower_u), evaluate(upper_u)),
            np.where(z < 0.5, evaluate(lower_phi), evaluate(upper_phi)),
        )

    @property
    def top_charge(self) -> float:
        return -self.electric_displacement(1.0)

    @property
    def heat_energy(self) -> float:
        return integrate_piecewise(lambda z: heat_flux(z) ** 2 / layer(z).conductivity)

    @property
    def enthalpy(self) -> float:
        def density(z: float) -> float:
            material = layer(z)
            u_prime, phi_prime = self.derivatives(z)
            mechanical = u_prime - material.expansion * (
                float(temperature(z)) - T_REFERENCE
            )
            return (
                0.5 * material.stiffness * mechanical**2
                + material.coupling * mechanical * phi_prime
                - 0.5 * material.permittivity * phi_prime**2
            )

        return integrate_piecewise(density)


def reference(
    *, top_voltage: float | None = None, top_charge: float | None = None
) -> Reference:
    if (top_voltage is None) == (top_charge is None):
        raise ValueError("choose top voltage or top charge")

    thermal = integrate_piecewise(
        lambda z: layer(z).expansion * (float(temperature(z)) - T_REFERENCE)
    )
    a = integrate_piecewise(lambda z: inverse_coefficients(layer(z))[0])
    b = integrate_piecewise(lambda z: inverse_coefficients(layer(z))[1])
    g = integrate_piecewise(lambda z: inverse_coefficients(layer(z))[2])
    mechanical_load = integrate_piecewise(
        lambda z: (
            -BODY_FORCE * z * inverse_coefficients(layer(z))[0]
            + VOLUME_CHARGE * z * inverse_coefficients(layer(z))[1]
        )
    )
    electric_load = integrate_piecewise(
        lambda z: (
            -BODY_FORCE * z * inverse_coefficients(layer(z))[1]
            + VOLUME_CHARGE * z * inverse_coefficients(layer(z))[2]
        )
    )
    if top_voltage is not None:
        sigma0, displacement0 = np.linalg.solve(
            np.array([[a, b], [b, g]]),
            np.array([-thermal - mechanical_load, top_voltage - electric_load]),
        )
        voltage = top_voltage
    else:
        assert top_charge is not None
        displacement0 = -top_charge - VOLUME_CHARGE
        sigma0 = (-thermal - mechanical_load - b * displacement0) / a
        voltage = b * sigma0 + g * displacement0 + electric_load
    return Reference(float(sigma0), float(voltage), float(displacement0))
