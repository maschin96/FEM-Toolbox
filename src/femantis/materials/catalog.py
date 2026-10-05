"""Small versioned thermal parameter library from identified primary sources.

Typical room-temperature constants, not temperature-dependent constitutive laws.
See Doc/32_MATERIAL_LIBRARY.md for source tables, conversions and limitations.
"""

from dataclasses import dataclass

from femantis.errors import ValidationError
from femantis.materials.heat import ThermalMaterial, ThermalSource

CATALOG_VERSION = 1


@dataclass(frozen=True)
class ThermalPreset:
    """SI scalar data: density kg/m³, specific_heat J/(kg K), conductivity W/(m K)."""

    key: str
    name: str
    density: float
    specific_heat: float
    conductivity: float
    reference_temperature: float
    source_url: str
    note: str

    def thermal_material(self) -> ThermalMaterial:
        """Build an immutable isotropic material with storage and source metadata."""
        return ThermalMaterial.isotropic(
            self.conductivity,
            density=self.density,
            specific_heat=self.specific_heat,
            source=ThermalSource(
                self.key, self.source_url, self.reference_temperature, self.note
            ),
        )


_PRESETS = (
    ThermalPreset(
        "aluminium_6060_t6_20c",
        "EN AW-6060 T6",
        2700.0,
        901.0,
        200.0,
        293.15,
        "https://www.hydro.com/globalassets/08-about-hydro/hydro-worldwide/austria/"
        "nenzing/alloy-data-sheets/hydro-en-aw-6060.pdf",
        "Hydro 6060 Rev.4, 03/2023, p.2: typical properties at 20°C; "
        "T6 conductivity. Source accessed 2026-09-23.",
    ),
    ThermalPreset(
        "copper_cda_20c",
        "Copper, CDA reference",
        8920.0,
        386.0,
        394.0,
        293.15,
        "https://www.copper.org/publications/pub_list/pdf/a1360.pdf",
        "CDA The Copper Advantage, Table 3, p.9: physical properties of copper "
        "at 20°C. Generic reference copper, not a certified alloy/temper lot.",
    ),
    ThermalPreset(
        "stainless_14301_rt",
        "Stainless steel EN 1.4301 (304)",
        7900.0,
        500.0,
        15.0,
        293.15,
        "https://www.outokumpu.com/-/media/from-sharepoint/documents/product/"
        "outokumpu-steel-grades-properties-global-standards.pdf",
        "Outokumpu Steel Grades, Properties and Global Standards, Table 13, p.13: "
        "room-temperature values; RT represented as 293.15 K, no validity range given.",
    ),
)


def material_catalog() -> tuple[ThermalPreset, ...]:
    """Return the immutable version-1 preset catalog (no network access)."""
    return _PRESETS


def thermal_material(key: str) -> ThermalMaterial:
    """Build a preset by exact key; unknown names raise ValidationError."""
    for preset in _PRESETS:
        if preset.key == key:
            return preset.thermal_material()
    raise ValidationError(f"unknown thermal preset {key!r}", field="material")
