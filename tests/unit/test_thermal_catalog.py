"""Thermal storage and source-backed room-temperature presets."""

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError
from fem_toolbox.materials import ThermalMaterial, assign_region_materials
from fem_toolbox.materials.catalog import material_catalog, thermal_material


def test_presets_create_usable_materials():
    entries = material_catalog()
    assert len(entries) == 3
    for entry in entries:
        material = thermal_material(entry.key)
        assert isinstance(material, ThermalMaterial)
        assert material.density == entry.density
        assert material.specific_heat == entry.specific_heat
        assert entry.source_url.startswith("https://")
        assert entry.reference_temperature == 293.15
    copper = thermal_material("copper_cda_20c")
    assert copper.density == 8920
    assert copper.specific_heat == 386
    np.testing.assert_array_equal(copper.conductivity, np.eye(3) * 394)
    with pytest.raises(ValidationError):
        thermal_material("unknown")


def test_storage_is_optional_and_region_provenance_preserves_it(cube_mesh):
    from fem_toolbox.mesh import Mesh

    base = cube_mesh(1)
    mesh = Mesh(
        points=base.points,
        tetra4=base.tetra4,
        volume_sets={"body": np.arange(base.n_tetra4)},
    )
    old = ThermalMaterial.isotropic(1.0)
    assert old.density is None and old.specific_heat is None
    mat = thermal_material("aluminium_6060_t6_20c")
    assignment = assign_region_materials(mesh, "T", {"body": "al"}, {"al": mat})
    data = assignment.to_provenance_dict()["catalog"]["al"]
    assert data["density"] == 2700
    assert data["specific_heat"] == 901
    assert data["source"]["key"] == "aluminium_6060_t6_20c"


@pytest.mark.parametrize(
    "rho,cp", [(0, 1), (1, -1), (True, 1), (1, np.inf), (1, None), (None, 1)]
)
def test_invalid_thermal_storage(rho, cp):
    with pytest.raises(ValidationError):
        ThermalMaterial.isotropic(1.0, density=rho, specific_heat=cp)
