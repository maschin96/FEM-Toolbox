"""Contract tests for global mathematical and unit conventions."""

import pytest

from fem_toolbox.conventions import (
    ELECTRIC_FIELD_DEFINITION,
    SHEAR_STRAIN_CONVENTION,
    SI_UNITS,
    SPATIAL_COMPONENTS,
    SPATIAL_DIMENSION,
    UNIT_SYSTEM,
    VOIGT_COMPONENTS,
    VOIGT_SIZE,
    VoigtIndex,
)


def test_spatial_and_voigt_order_match_mathematical_model() -> None:
    """Executable ordering agrees with Doc/02_MATHEMATICAL_MODEL.md."""
    assert SPATIAL_DIMENSION == 3
    assert SPATIAL_COMPONENTS == ("x", "y", "z")
    assert VOIGT_SIZE == 6
    assert VOIGT_COMPONENTS == ("xx", "yy", "zz", "yz", "xz", "xy")
    assert tuple(index.value for index in VoigtIndex) == tuple(range(VOIGT_SIZE))
    assert tuple(index.name.lower() for index in VoigtIndex) == VOIGT_COMPONENTS
    assert SHEAR_STRAIN_CONVENTION == "engineering"
    assert ELECTRIC_FIELD_DEFINITION == "E = -grad(phi)"


def test_si_unit_map_is_explicit_and_immutable() -> None:
    """Core quantities have one internal unit system that cannot drift at runtime."""
    assert UNIT_SYSTEM == "SI"
    assert SI_UNITS["length"] == "m"
    assert SI_UNITS["stress"] == "Pa"
    assert SI_UNITS["electric_potential"] == "V"
    assert SI_UNITS["permittivity"] == "F/m"

    with pytest.raises(TypeError):
        SI_UNITS["length"] = "mm"  # type: ignore[index]
