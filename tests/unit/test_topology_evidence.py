"""Candidate measurements remain tied to the actual numerical implementation."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).parents[2]


def test_m11_evidence_has_exact_sources_and_preserves_no_go():
    archive = ROOT / "Doc/evidence/m11"
    receipts = json.loads((archive / "environment.json").read_text())
    for name, digest in receipts["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    for name, digest in receipts["evidence_sha256"].items():
        assert hashlib.sha256((archive / name).read_bytes()).hexdigest() == digest, name
    final = json.loads((archive / "verification.json").read_text())
    initial = json.loads((archive / "unscaled_no_go.json").read_text())
    assert final["complete_matrix"] and final["gates_passed"]
    assert len(final["rows"]) == 6
    assert all(row["geometry_gate_passed"] for row in final["rows"])
    assert not initial["gates_passed"]
    assert receipts["independent_human_review"] == "not_available"
