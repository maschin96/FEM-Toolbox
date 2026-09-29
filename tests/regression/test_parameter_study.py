"""Documented installed-style study example and deterministic continuation."""

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parents[2]


def test_study_example_and_resume(tmp_path: Path) -> None:
    command = [
        sys.executable,
        str(ROOT / "examples/parameter_study/run.py"),
        "--output",
        str(tmp_path),
    ]
    env = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
    subprocess.run([*command, "--max-cases", "1"], cwd=tmp_path, env=env, check=True)
    catalog = tmp_path / "catalog.json"
    assert [c["status"] for c in json.loads(catalog.read_text())["cases"]] == [
        "succeeded",
        "pending",
    ]
    subprocess.run(command, cwd=tmp_path, env=env, check=True)
    before = catalog.read_bytes()
    subprocess.run(command, cwd=tmp_path, env=env, check=True)
    assert catalog.read_bytes() == before
    assert [c["attempts"] for c in json.loads(before)["cases"]] == [1, 1]
