"""Project-wide mathematical, component-ordering, and SI-unit conventions.

These values are executable counterparts of ``Doc/02_MATHEMATICAL_MODEL.md``.
Changing one is an architecture decision and requires corresponding documentation
and verification updates.
"""

from collections.abc import Mapping
from enum import IntEnum
from types import MappingProxyType
from typing import Final

SPATIAL_DIMENSION: Final = 3
SPATIAL_COMPONENTS: Final = ("x", "y", "z")

VOIGT_SIZE: Final = 6
VOIGT_COMPONENTS: Final = ("xx", "yy", "zz", "yz", "xz", "xy")
SHEAR_STRAIN_CONVENTION: Final = "engineering"
ELECTRIC_FIELD_DEFINITION: Final = "E = -grad(phi)"

UNIT_SYSTEM: Final = "SI"
SI_UNITS: Final[Mapping[str, str]] = MappingProxyType(
    {
        "length": "m",
        "displacement": "m",
        "density": "kg/m^3",
        "stress": "Pa",
        "young_modulus": "Pa",
        "temperature_difference": "K",
        "thermal_conductivity": "W/(m K)",
        "specific_heat": "J/(kg K)",
        "electric_potential": "V",
        "electric_field": "V/m",
        "electric_displacement": "C/m^2",
        "permittivity": "F/m",
        "piezoelectric_stress": "C/m^2",
    }
)


class VoigtIndex(IntEnum):
    """Indices for engineering-strain and stress vectors in three dimensions."""

    XX = 0
    YY = 1
    ZZ = 2
    YZ = 3
    XZ = 4
    XY = 5
