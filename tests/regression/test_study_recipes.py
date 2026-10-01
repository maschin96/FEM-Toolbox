"""Execute the documented recipes, including direct and analytic comparisons."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]


@pytest.mark.parametrize("recipe,count", [("heat", 2), ("piezo", 4)])
def test_recipe_acceptance(tmp_path, recipe, count):
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "examples/study_recipes/run.py"),
            recipe,
            "--output",
            str(tmp_path / recipe),
        ],
        check=True,
    )
    report = json.loads((tmp_path / recipe / "verification.json").read_text())
    assert len(report["cases"]) == count
    assert report["resume_attempts_unchanged"]
    assert report["input_change_invalidated"]
    assert report["failed_cases"] == 1
    assert report["direct_and_analytic_checks_passed"]
