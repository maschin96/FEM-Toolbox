"""Region maps must resolve every active Tet4 before any solver runs."""

import numpy as np
import pytest

from fem_toolbox.errors import MaterialAssignmentError
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
    assign_region_materials,
)
from fem_toolbox.mesh import Mesh


def two_region_mesh() -> Mesh:
    return Mesh(
        points=np.array(
            [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
            dtype=float,
        ),
        tetra4=np.array([[0, 1, 2, 3], [1, 2, 3, 4]]),
        tetra4_ids=np.array([101, 102]),
        volume_sets={"left": np.array([0]), "right": np.array([1])},
    )


def test_complete_assignment_keeps_cell_order_and_traceable_data() -> None:
    mesh = two_region_mesh()
    left = ThermalMaterial.isotropic(2)
    right = ThermalMaterial.isotropic(4)
    result = assign_region_materials(
        mesh, "T", {"right": "k4", "left": "k2"}, {"k2": left, "k4": right}
    )
    assert result.cell_materials == (left, right)
    assert result.material_names == ("k2", "k4")
    assert result.region_names == (("left",), ("right",))
    assert result.cell_ids == (101, 102)
    assert result.to_provenance_dict()["catalog"]["k2"]["conductivity"][0][0] == 2


def test_identical_overlap_is_allowed_and_records_both_regions() -> None:
    mesh = two_region_mesh()
    mesh = Mesh(
        points=mesh.points,
        tetra4=mesh.tetra4,
        volume_sets={"left": [0], "same": [0], "right": [1]},
    )
    result = assign_region_materials(
        mesh,
        "T",
        {"left": "k2", "same": "k2", "right": "k4"},
        {"k2": ThermalMaterial.isotropic(2), "k4": ThermalMaterial.isotropic(4)},
    )
    assert result.region_names[0] == ("left", "same")


@pytest.mark.parametrize(
    ("regions", "catalog", "kind", "cell_id"),
    [
        ({"left": "k2"}, {"k2": 2}, "missing", 102),
        (
            {"left": "k2", "right": "k4", "ghost": "k2"},
            {"k2": 2, "k4": 4},
            "unknown_region",
            None,
        ),
        (
            {"left": "k2", "right": "k4"},
            {"k2": 2, "k4": 4, "extra": 8},
            "unused_material",
            None,
        ),
    ],
)
def test_invalid_assignment_is_diagnosed(
    regions: dict[str, str], catalog: dict[str, int], kind: str, cell_id: int | None
) -> None:
    with pytest.raises(MaterialAssignmentError) as caught:
        assign_region_materials(
            two_region_mesh(),
            "T",
            regions,
            {name: ThermalMaterial.isotropic(value) for name, value in catalog.items()},
        )
    assert caught.value.kind == kind
    assert caught.value.cell_id == cell_id


def test_conflicting_overlap_reports_both_region_names() -> None:
    mesh = two_region_mesh()
    mesh = Mesh(
        points=mesh.points,
        tetra4=mesh.tetra4,
        tetra4_ids=mesh.tetra4_ids,
        volume_sets={"left": [0], "other": [0], "right": [1]},
    )
    with pytest.raises(MaterialAssignmentError) as caught:
        assign_region_materials(
            mesh,
            "T",
            {"left": "k2", "other": "k4", "right": "k4"},
            {"k2": ThermalMaterial.isotropic(2), "k4": ThermalMaterial.isotropic(4)},
        )
    assert caught.value.kind == "conflict"
    assert caught.value.cell_id == 101
    assert set(caught.value.regions) == {"left", "other"}


@pytest.mark.parametrize(
    ("field", "material"),
    [
        ("u", ElasticMaterial.isotropic(2e6, 0.2)),
        ("phi", DielectricMaterial.isotropic(2e-9)),
        ("expansion", ThermalExpansion([1e-5, 0, 0, 0, 0, 0], 293.15)),
        (
            "piezo",
            PiezoelectricMaterial(np.eye(6) * 2e6, np.zeros((3, 6)), np.eye(3) * 2e-9),
        ),
    ],
)
def test_supported_field_types_and_provenance(field: str, material: object) -> None:
    result = assign_region_materials(
        two_region_mesh(), field, {"left": "one", "right": "one"}, {"one": material}
    )
    assert result.cell_materials == (material, material)
    assert result.to_provenance_dict()["field"] == field


def test_wrong_material_type_is_rejected() -> None:
    with pytest.raises(MaterialAssignmentError) as caught:
        assign_region_materials(
            two_region_mesh(),
            "T",
            {"left": "wrong", "right": "wrong"},
            {"wrong": DielectricMaterial.isotropic(2e-9)},
        )
    assert caught.value.kind == "material_type"


def test_principal_expansion_axes_are_recorded_per_region() -> None:
    axes = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    expansion = ThermalExpansion.from_principal([1e-5, 2e-5, 3e-5], 293.15, axes)
    result = assign_region_materials(
        two_region_mesh(),
        "expansion",
        {"left": "rotated", "right": "rotated"},
        {"rotated": expansion},
    )
    assert (
        result.to_provenance_dict()["catalog"]["rotated"]["axes_local_to_global"]
        == axes.tolist()
    )


def test_empty_mesh_and_nonmapping_specs_are_rejected() -> None:
    with pytest.raises(MaterialAssignmentError) as empty:
        assign_region_materials(Mesh.empty(), "T", {}, {})
    assert empty.value.kind == "empty_mesh"
    with pytest.raises(MaterialAssignmentError) as invalid:
        assign_region_materials(two_region_mesh(), "T", [], {})  # type: ignore[arg-type]
    assert invalid.value.kind == "regions_type"
