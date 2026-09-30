"""Exercise the public Model-to-JPEG example without a display."""

import runpy
from pathlib import Path

import pytest

pytest.importorskip("matplotlib", reason="requires visualization extra")


def test_example(tmp_path):
    from PIL import Image

    example = Path(__file__).parents[2] / "examples/field_visualization/run.py"
    runpy.run_path(str(example))["run"](tmp_path)
    images = list(tmp_path.glob("*.jpg"))
    assert len(images) == 6
    for path in images:
        with Image.open(path) as image:
            assert image.format == "JPEG"
            assert image.size == (1200, 900)
