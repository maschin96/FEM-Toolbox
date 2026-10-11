"""Source-only signatures must expose added/changed public contracts."""

import json
import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def test_ast_inventory_detects_signature_and_dataclass_field_changes(tmp_path):
    api = runpy.run_path(str(ROOT / "scripts/api_contracts.py"))
    source = tmp_path / "example.py"
    source.write_text("""__all__ = ["Law", "solve"]
class Law:
    mean: float
    def sample(self, size: int = 2) -> list[float]:
        return []
def solve(x: float, *, seed: int = 1) -> float:
    return x
""")
    original = api["module_contract"](source)
    source.write_text(source.read_text().replace("seed: int = 1", "seed: int = 2"))
    assert original != api["module_contract"](source)
    source.write_text(
        source.read_text().replace("mean: float", "mean: float\n    deviation: float")
    )
    assert original != api["module_contract"](source)


def test_uncertainty_signatures_match_documented_contract():
    api = runpy.run_path(str(ROOT / "scripts/api_contracts.py"))
    documented = json.loads((ROOT / "Doc/latex/uncertainty-api.json").read_text())
    actual = api["module_contract"](ROOT / "src/fem_toolbox/uncertainty.py")
    assert documented == actual, "uncertainty API drift: update contract and handbook"


def test_full_inventory_resolves_reexports_and_rejects_missing_module(tmp_path):
    api = runpy.run_path(str(ROOT / "scripts/api_contracts.py"))
    package = tmp_path / "src/fem_toolbox"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text(
        'from fem_toolbox.actual import Law\n__all__ = ["Law"]\n'
    )
    (package / "actual.py").write_text(
        'class Law:\n    """Explicit law."""\n    mean: float\n'
    )
    inventory = api["collect_api"](tmp_path)
    assert (
        inventory["modules"]["fem_toolbox"]["reexports"]["Law"]
        == "fem_toolbox.actual.Law"
    )
    assert (
        inventory["modules"]["fem_toolbox.actual"]["definitions"]["Law"]["fields"][
            "mean"
        ]["annotation"]
        == "float"
    )
    (package / "actual.py").unlink()
    with pytest.raises(ValueError, match="reexport"):
        api["collect_api"](tmp_path)


def test_repository_api_inventory_has_no_signature_or_missing_reference_drift():
    api = runpy.run_path(str(ROOT / "scripts/api_contracts.py"))
    actual = api["collect_api"](ROOT)
    documented = json.loads((ROOT / "Doc/latex/api-inventory.json").read_text())
    assert actual == documented, "API inventory drift: regenerate and review contracts"
    mapping = json.loads((ROOT / "Doc/latex/api-map.json").read_text())
    api["validate_api_map"](ROOT, actual, mapping)
    broken = json.loads(json.dumps(mapping))
    broken["modules"].pop(next(iter(broken["modules"])))
    with pytest.raises(ValueError, match="coverage"):
        api["validate_api_map"](ROOT, actual, broken)


def test_rendered_api_chapter_is_generated_from_exact_inventory():
    api = runpy.run_path(str(ROOT / "scripts/api_contracts.py"))
    inventory = json.loads((ROOT / "Doc/latex/api-inventory.json").read_text())
    mapping = json.loads((ROOT / "Doc/latex/api-map.json").read_text())
    assert (ROOT / "Doc/latex/chapters/api/inventory.tex").read_text() == api[
        "render_inventory"
    ](inventory, mapping)


def test_canonical_stable_entry_cannot_disappear_from_inventory_map():
    api = runpy.run_path(str(ROOT / "scripts/api_contracts.py"))
    inventory = api["collect_api"](ROOT)
    mapping = json.loads((ROOT / "Doc/latex/api-map.json").read_text())
    mapping["modules"]["fem_toolbox.model"]["stable_symbols"].remove("Model")
    with pytest.raises(ValueError, match="canonical"):
        api["validate_api_map"](ROOT, inventory, mapping)
