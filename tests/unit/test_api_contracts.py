"""Source-only signatures must expose added/changed public contracts."""

import json
import runpy
from pathlib import Path

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
