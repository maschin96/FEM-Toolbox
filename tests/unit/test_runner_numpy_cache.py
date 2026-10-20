"""A runner dependency cache must never silently accept changed wheel bytes."""

import hashlib
import json
import runpy
from pathlib import Path

import pytest

API = runpy.run_path(str(Path(__file__).parents[2] / "ci/build_runner_numpy.py"))


def test_cache_accepts_exact_receipt_and_rejects_wheel_corruption(tmp_path):
    wheel = tmp_path / "numpy-test.whl"
    wheel.write_bytes(b"wheel bytes")
    receipt = {
        "source": {"sha256": "source"},
        "wheel": {
            "name": wheel.name,
            "sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        },
    }
    (tmp_path / "BUILD_INFO.json").write_text(json.dumps(receipt))
    assert API["cached_wheel"](tmp_path, "source") == wheel
    wheel.write_bytes(b"changed wheel")
    with pytest.raises(ValueError, match="hash"):
        API["cached_wheel"](tmp_path, "source")


def test_cache_rejects_wrong_source_or_escaping_wheel_path(tmp_path):
    receipt = {
        "source": {"sha256": "wrong"},
        "wheel": {"name": "../escape.whl", "sha256": "hash"},
    }
    (tmp_path / "BUILD_INFO.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="source"):
        API["cached_wheel"](tmp_path, "source")
    receipt["source"]["sha256"] = "source"
    (tmp_path / "BUILD_INFO.json").write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="path"):
        API["cached_wheel"](tmp_path, "source")
