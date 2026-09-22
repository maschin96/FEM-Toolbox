"""Exactly integrated quasistatic piezoelectric Tet4 operators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.elasticity import strain_matrix
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import PiezoelectricMaterial


def _validate(geometry: Tet4Geometry, material: PiezoelectricMaterial) -> None:
    if not isinstance(geometry, Tet4Geometry):
        raise ValidationError("expected Tet4Geometry", field="geometry")
    if not isinstance(material, PiezoelectricMaterial):
        raise ValidationError("expected PiezoelectricMaterial", field="material")


def _thermal_strain(
    value: ArrayLike | None, electric_thermal_strain: bool
) -> NDArray[np.float64]:
    if not isinstance(electric_thermal_strain, bool):
        raise ValidationError(
            "expected boolean thermal-electric switch", field="thermal"
        )
    if electric_thermal_strain and value is None:
        raise ValidationError(
            "electric thermal strain requires an explicit thermal strain",
            field="thermal",
        )
    return np.zeros(6) if value is None else real_array(value, (6,))


def _readonly(value: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.frombuffer(
        np.ascontiguousarray(value).tobytes(), dtype=np.float64
    ).reshape(value.shape)


def piezoelectric_matrix(
    geometry: Tet4Geometry, material: PiezoelectricMaterial
) -> NDArray[np.float64]:
    """Return the signed symmetric (16,16) local tangent in [u(12), phi(4)] order.

    With B (6,12), G (3,4), and volume V, the blocks are
    ``V B.T C B``, ``V B.T e.T G``, its transpose, and ``-V G.T kappa G``.
    Mechanical rows have N/m for u columns and N/V for phi columns;
    electric rows have C/m for u columns and C/V for phi columns. Invalid
    geometry or material types raise ValidationError. Geometry itself validates
    orientation and degeneracy on construction.
    """
    _validate(geometry, material)
    b = strain_matrix(geometry)
    g = geometry.gradient_matrix
    c = material.stiffness
    e = material.stress_coupling
    k = material.permittivity
    mechanical = b.T @ c @ b
    coupling = b.T @ e.T @ g
    dielectric = -(g.T @ k @ g)
    return np.asarray(
        geometry.volume * np.block([[mechanical, coupling], [coupling.T, dielectric]]),
        dtype=np.float64,
    )


def piezoelectric_thermal_load(
    geometry: Tet4Geometry,
    material: PiezoelectricMaterial,
    thermal_strain: ArrayLike,
    *,
    electric_thermal_strain: bool = False,
) -> NDArray[np.float64]:
    """Return local thermal RHS (16,) for a constant engineering-Voigt strain.

    The mechanical part is ``V B.T C epsilon_th`` in N. If the explicit
    electric switch is enabled, the electric part is ``V G.T e epsilon_th``
    in C; otherwise it is exactly zero. This switch selects whether D uses
    total strain or mechanical strain. The stiffness matrix is unchanged.
    """
    _validate(geometry, material)
    thermal = _thermal_strain(thermal_strain, electric_thermal_strain)
    b = strain_matrix(geometry)
    g = geometry.gradient_matrix
    force = b.T @ material.stiffness @ thermal
    charge = (
        g.T @ material.stress_coupling @ thermal
        if electric_thermal_strain
        else np.zeros(4)
    )
    return np.asarray(
        geometry.volume * np.concatenate([force, charge]), dtype=np.float64
    )


@dataclass(frozen=True, slots=True, eq=False)
class PiezoelectricElementResponse:
    """Constant fields and integrated electric enthalpy terms of one Tet4.

    ``strain`` uses engineering Voigt order (6,), ``electric_field`` and
    ``electric_displacement`` are (3,), and ``stress`` is (6,). Energies are
    in J: mechanical, coupling, negative dielectric, and their sum.
    """

    strain: NDArray[np.float64]
    electric_field: NDArray[np.float64]
    stress: NDArray[np.float64]
    electric_displacement: NDArray[np.float64]
    mechanical_energy: float
    coupling_energy: float
    dielectric_enthalpy: float
    total_enthalpy: float


def piezoelectric_response(
    geometry: Tet4Geometry,
    material: PiezoelectricMaterial,
    displacement: ArrayLike,
    potential: ArrayLike,
    *,
    thermal_strain: ArrayLike | None = None,
    electric_thermal_strain: bool = False,
) -> PiezoelectricElementResponse:
    """Evaluate constant fields and integrated electric enthalpy for one Tet4.

    ``displacement`` has shape (4,3) in m, ``potential`` (4,) in V, and the
    optional constant ``thermal_strain`` (6,) is dimensionless. If enabled,
    D uses strain minus thermal strain; otherwise D uses total strain. In both
    cases stress uses mechanical strain. Invalid shapes or a missing strain for
    the enabled switch raise ValidationError.
    """
    _validate(geometry, material)
    u = real_array(displacement, (4, 3))
    phi = real_array(potential, (4,))
    thermal = _thermal_strain(thermal_strain, electric_thermal_strain)
    strain = strain_matrix(geometry) @ u.ravel()
    field = -(geometry.gradient_matrix @ phi)
    mechanical_strain = strain - thermal
    coupling_strain = mechanical_strain if electric_thermal_strain else strain
    stress = material.stiffness @ mechanical_strain - material.stress_coupling.T @ field
    electric_displacement = (
        material.stress_coupling @ coupling_strain + material.permittivity @ field
    )
    volume = geometry.volume
    mechanical_energy = float(
        0.5 * volume * mechanical_strain @ material.stiffness @ mechanical_strain
    )
    coupling_energy = float(
        -volume * coupling_strain @ material.stress_coupling.T @ field
    )
    dielectric_enthalpy = float(-0.5 * volume * field @ material.permittivity @ field)
    return PiezoelectricElementResponse(
        strain=_readonly(strain),
        electric_field=_readonly(field),
        stress=_readonly(stress),
        electric_displacement=_readonly(electric_displacement),
        mechanical_energy=mechanical_energy,
        coupling_energy=coupling_energy,
        dielectric_enthalpy=dielectric_enthalpy,
        total_enthalpy=mechanical_energy + coupling_energy + dielectric_enthalpy,
    )
