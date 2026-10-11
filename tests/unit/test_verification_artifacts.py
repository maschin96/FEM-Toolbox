"""Audit documentation evidence without rerunning historical timing measurements."""

import copy
import json
import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]


def api():
    return runpy.run_path(str(ROOT / "scripts/verification_artifacts.py"))


def baseline():
    return json.loads(
        (ROOT / "benchmarks/baselines/m7-37c-linux-arm64.json").read_text()
    )


def test_complete_baseline_has_independently_recomputed_statistics_and_missing_goals():
    report = api()["audit_baseline"](baseline())
    assert len(report["groups"]) == 24
    assert report["sample_count"] == 72
    assert sum(r["selected"] is None for r in report["accuracy"]) == 4
    rows = [r for r in report["groups"] if r["case"] == "heat" and r["order"] == 1]
    assert min(r["h1_relative"] for r in rows) > 0.125


def test_missing_duplicate_failed_or_altered_historical_results_are_rejected():
    audit = api()["audit_baseline"]
    original = baseline()
    for mutation in ("missing", "duplicate", "failed", "altered"):
        value = copy.deepcopy(original)
        if mutation == "missing":
            value["samples"].pop()
        elif mutation == "duplicate":
            value["samples"][1] = value["samples"][0]
        elif mutation == "failed":
            value["samples"][0]["balance"] = 0.001
        else:
            value["groups"][0]["wall_s"]["median"] *= 2
        with pytest.raises(ValueError):
            audit(value)


def test_handbook_evidence_cannot_omit_supported_physics_or_reference():
    validate = api()["validate_map"]
    mapping = json.loads((ROOT / "Doc/latex/verification-map.json").read_text())
    validate(ROOT, mapping)
    broken = copy.deepcopy(mapping)
    broken["entries"] = [r for r in broken["entries"] if r["id"] != "harmonic"]
    with pytest.raises(ValueError, match="coverage"):
        validate(ROOT, broken)
    broken = copy.deepcopy(mapping)
    broken["entries"][0]["tests"][0] = "tests/absent.py"
    with pytest.raises(ValueError, match="reference"):
        validate(ROOT, broken)


def test_generated_evidence_matches_versioned_sources_and_figures():
    api()["check_artifacts"](ROOT)
