"""Claims require matching evidence, not a closed issue or a branch name."""

import copy
import json
import runpy
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).parents[2]


def api():
    return runpy.run_path(str(ROOT / "scripts/capability_catalog.py"))


def catalog():
    return json.loads((ROOT / "Doc/evidence/catalog/catalog.json").read_text())


def test_published_candidate_and_missing_measurements_stay_distinct():
    rows = api()["validate_catalog"](ROOT, catalog())
    by_id = {row["id"]: row for row in rows}
    assert by_id["heat"]["released_version"] == "1.4.0"
    assert by_id["adaptive_heat"]["released_version"] is None
    assert by_id["quad_fv"]["status"] == "branch_candidate"
    assert by_id["experimental_piezo"]["verification"] == "missing_evidence"
    assert all(row["experimental_validation"] == "not_performed" for row in rows)


@pytest.mark.parametrize(
    "key,value", [("released_version", "1.5.0"), ("commit", "0" * 40)]
)
def test_false_release_provenance_is_rejected(key, value):
    data = catalog()
    data["entries"][0][key] = value
    with pytest.raises(ValueError, match="release provenance"):
        api()["validate_catalog"](ROOT, data)


def test_branch_claim_cannot_be_promoted_by_an_issue_state():
    data = catalog()
    row = next(r for r in data["entries"] if r["id"] == "quad_fv")
    row.update(status="released", released_version="1.4.0", release="v1.4.0")
    with pytest.raises(ValueError, match="release"):
        api()["validate_catalog"](ROOT, data)


@pytest.mark.parametrize("kind", ["hash", "missing", "escape"])
def test_missing_mismatched_or_external_evidence_is_rejected(kind):
    data = catalog()
    evidence = data["entries"][0]["evidence"][0]
    evidence.update(
        {"sha256": "0" * 64}
        if kind == "hash"
        else {"path": "absent.json" if kind == "missing" else "../outside.json"}
    )
    with pytest.raises(ValueError, match="evidence"):
        api()["validate_catalog"](ROOT, data)


def test_released_scope_cannot_exceed_archived_handbook_scope():
    data = catalog()
    data["entries"][0]["scope"] = "Arbitrary turbulent fluid and solid models"
    with pytest.raises(ValueError, match="release scope"):
        api()["validate_catalog"](ROOT, data)


def test_unsupported_schema_and_duplicate_identity_fail_closed():
    data = catalog()
    data["schema"] = 2
    with pytest.raises(ValueError, match="schema"):
        api()["validate_catalog"](ROOT, data)
    data = catalog()
    data["entries"].append(copy.deepcopy(data["entries"][0]))
    with pytest.raises(ValueError, match="duplicate"):
        api()["validate_catalog"](ROOT, data)


def test_reader_exposes_limits_and_missing_experimental_evidence():
    data = catalog()
    text = api()["render_catalog"](api()["validate_catalog"](ROOT, data))
    assert "branch_candidate" in text
    assert "missing_evidence" in text
    assert "not_performed" in text
    assert "Tet4" in text
    assert "Review" in text


def test_omitted_release_capability_and_unbacked_verification_are_rejected():
    data = catalog()
    data["entries"].pop(0)
    with pytest.raises(ValueError, match="coverage"):
        api()["validate_catalog"](ROOT, data)
    data = catalog()
    data["entries"][0]["evidence"] = []
    with pytest.raises(ValueError, match="missing evidence"):
        api()["validate_catalog"](ROOT, data)


def test_package_hash_must_match_published_verification_receipt():
    data = catalog()
    data["releases"]["v1.4.0"]["source_distribution"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="evidence hash"):
        api()["validate_catalog"](ROOT, data)


def test_cli_reads_catalog_outside_checkout_and_reports_missing_file(tmp_path):
    command = [sys.executable, str(ROOT / "scripts/capability_catalog.py")]
    result = subprocess.run(
        [*command, "--format", "json"], cwd=tmp_path, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    assert len(json.loads(result.stdout)) == 20
    result = subprocess.run(
        [*command, "--catalog", "absent.json"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "absent.json" in result.stderr
