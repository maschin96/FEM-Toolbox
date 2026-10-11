"""Keep the independently runnable handbook inputs tied to the candidate code."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BOOK = ROOT / "Doc/latex"


def test_handbook_inputs_are_complete_and_identical_to_versioned_examples():
    manifest = json.loads((BOOK / "recipes.json").read_text())
    originals = {
        str(path.relative_to(ROOT))
        for path in (ROOT / "examples").rglob("*")
        if path.is_file() and path.suffix in {".py", ".toml", ".msh", ".csv"}
    }
    assert set(manifest["inputs"]) == originals
    for source, digest in manifest["inputs"].items():
        original = (ROOT / source).read_bytes()
        embedded = (BOOK / source).read_bytes()
        assert embedded == original, f"handbook source drift: {source}"
        assert hashlib.sha256(embedded).hexdigest() == digest


def test_every_executable_example_has_a_documented_recipe_and_pdf_listing():
    manifest = json.loads((BOOK / "recipes.json").read_text())
    recipes = manifest["recipes"]
    assert len({recipe["id"] for recipe in recipes}) == len(recipes)
    scripts = {source for source in manifest["inputs"] if source.endswith(".py")}
    assert {recipe["script"] for recipe in recipes} == scripts
    text = "\n".join(
        path.read_text() for path in (BOOK / "chapters/application").glob("*.tex")
    )
    for recipe in recipes:
        assert recipe["command"][0:2] == ["python", "-I"]
        assert recipe["command"][2] == recipe["script"]
        assert recipe["expectation"] and recipe["limits"] and recipe["troubleshooting"]
        assert recipe["artifacts"]
        assert recipe["script"] in recipe["inputs"]
        assert set(recipe["inputs"]) <= set(manifest["inputs"])
        assert recipe["runtime_budget_seconds"] == 90
        for command in recipe.get("continuation_commands", []):
            assert command[:3] == recipe["command"][:3]
        assert f"\\label{{recipe:{recipe['id']}}}" in text
        assert f"\\lstinputlisting{{{recipe['script']}}}" in text
