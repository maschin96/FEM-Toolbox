"""Audit versioned capability claims and render their evidence and limits."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from email.parser import BytesParser
from pathlib import Path
from typing import Any

CATALOG = "Doc/evidence/catalog/catalog.json"


def read_evidence(root: Path, receipt: dict[str, Any]) -> bytes:
    """Bind each local reference to exact bytes inside the source tree."""
    name = receipt.get("path")
    digest = receipt.get("sha256")
    if not isinstance(name, str) or not isinstance(digest, str):
        raise ValueError("evidence requires path and sha256")
    path = root / name
    if (
        Path(name).is_absolute()
        or not path.resolve().is_relative_to(root.resolve())
        or not path.is_file()
    ):
        raise ValueError(f"evidence missing/invalid path: {name}")
    content = path.read_bytes()
    if hashlib.sha256(content).hexdigest() != digest:
        raise ValueError(f"evidence hash mismatch: {name}")
    return content


def validate_catalog(root: Path, data: dict[str, Any]) -> list[dict[str, Any]]:
    """Reject unsupported provenance; retain explicit missing-evidence records."""
    if (
        not isinstance(data, dict)
        or type(data.get("schema")) is not int
        or data["schema"] != 1
    ):
        raise ValueError("unsupported capability catalog schema")
    entries = data.get("entries")
    releases = data.get("releases")
    if not isinstance(entries, list) or not entries or not isinstance(releases, dict):
        raise ValueError("catalog requires entries and releases")
    release_scopes = {}
    for name, release in releases.items():
        manifest = json.loads(read_evidence(root, release["manifest"]))
        verification = json.loads(read_evidence(root, release["verification"]))
        mapping_bytes = read_evidence(root, release["capability_map"])
        mapping = json.loads(mapping_bytes)
        package = release["source_distribution"]
        read_evidence(root, package)
        with tarfile.open(root / package["path"]) as archive:
            metadata_file = archive.extractfile(
                f"fem_toolbox-{manifest['version']}/PKG-INFO"
            )
            if metadata_file is None:
                raise ValueError(f"release package metadata missing: {name}")
            metadata = BytesParser().parsebytes(metadata_file.read())
        if (
            manifest.get("schema") != 1
            or manifest.get("status") != "release"
            or manifest.get("source_state") != "clean"
            or name != f"v{manifest.get('version')}"
            or verification.get("version") != manifest["version"]
            or verification.get("commit") != manifest["commit"]
            or verification.get("source_state") != "clean"
            or verification.get("ci", {}).get("conclusion") != "success"
            or verification.get("ci", {}).get("head_sha") != manifest["commit"]
            or metadata["Name"] != "fem-toolbox"
            or metadata["Version"] != manifest["version"]
            or verification.get("package_artifacts", {}).get(Path(package["path"]).name)
            != package["sha256"]
            or verification.get("handbook_manifest_sha256")
            != release["manifest"]["sha256"]
            or hashlib.sha256(mapping_bytes).hexdigest()
            != manifest.get("inputs", {}).get("Doc/latex/verification-map.json")
        ):
            raise ValueError(f"release provenance mismatch: {name}")
        release_scopes[name] = (
            manifest,
            {row["id"]: row for row in mapping["entries"]},
        )
    ids = set()
    for row in entries:
        identity = row.get("id")
        if not isinstance(identity, str) or not re.fullmatch(
            r"[a-z][a-z0-9_]*", identity
        ):
            raise ValueError("invalid capability identity")
        if identity in ids:
            raise ValueError(f"duplicate capability: {identity}")
        ids.add(identity)
        for field in ("title", "scope", "criterion", "limits", "review"):
            if not isinstance(row.get(field), str) or not row[field].strip():
                raise ValueError(f"{identity}: missing {field}")
        if row.get("status") not in {"released", "branch_candidate", "planned"}:
            raise ValueError(f"{identity}: invalid publication status")
        if row.get("verification") not in {"verified", "missing_evidence"}:
            raise ValueError(f"{identity}: invalid verification status")
        if row.get("experimental_validation") != "not_performed":
            raise ValueError(f"{identity}: schema 1 has no measurement-validation gate")
        evidence = row.get("evidence")
        if not isinstance(evidence, list):
            raise ValueError(f"{identity}: invalid evidence")
        if row["verification"] == "verified" and not evidence:
            raise ValueError(f"{identity}: missing evidence for verification claim")
        for receipt in evidence:
            read_evidence(root, receipt)
        commit = row.get("commit")
        if row["status"] != "planned" and (
            not isinstance(commit, str) or not re.fullmatch(r"[0-9a-f]{40}", commit)
        ):
            raise ValueError(f"{identity}: missing source commit")
        if row["status"] == "released":
            release = row.get("release")
            if release not in release_scopes:
                raise ValueError(f"{identity}: missing release provenance")
            manifest, scopes = release_scopes[release]
            if (
                row.get("released_version") != manifest["version"]
                or commit != manifest["commit"]
            ):
                raise ValueError(f"{identity}: release provenance mismatch")
            if (
                identity not in scopes
                or row["scope"] != scopes[identity]["scope"]
                or row["criterion"] != scopes[identity]["criterion"]
                or row["verification"] != "verified"
            ):
                raise ValueError(f"{identity}: release scope mismatch")
            paths = {receipt["path"] for receipt in evidence}
            if not {scopes[identity]["contract"], *scopes[identity]["tests"]} <= paths:
                raise ValueError(f"{identity}: missing release evidence")
        elif row.get("released_version") is not None or row.get("release") is not None:
            raise ValueError(f"{identity}: candidate/planned has release provenance")
    for name, (_, scopes) in release_scopes.items():
        captured = {row["id"] for row in entries if row.get("release") == name}
        if captured != set(scopes):
            raise ValueError(f"release coverage missing/extra capabilities: {name}")
    return entries


def render_catalog(rows: list[dict[str, Any]]) -> str:
    """Connect acceptance quantities, references, scope and review in a read view."""

    def cell(value: str) -> str:
        return value.replace("|", "\\|").replace("\n", " ")

    lines = [
        "# Nachweiskatalog (Schema 1)",
        "",
        "Verifikation belegt den angegebenen Umfang; experimentelle Validierung wird separat ausgewiesen.",
        "",
        "| Fähigkeit / Version | Status | Kennwert / Kriterium | Nachweis | Umfang / Grenze | Review / experimentelle Validierung |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        references = (
            "; ".join(
                f"[{r.get('source_path', r['path'])}](https://github.com/maschin96/FEM-Toolbox/blob/{row['commit']}/{r.get('source_path', r['path'])})"
                for r in row["evidence"]
            )
            or "missing_evidence"
        )
        values = [
            f"{row['title']} / {row['released_version'] or 'keine freigegebene Version'}",
            f"{row['status']} / {row['verification']}",
            row["criterion"],
            references,
            f"{row['scope']}; {row['limits']}",
            f"{row['review']}; {row['experimental_validation']}",
        ]
        lines.append("| " + " | ".join(cell(value) for value in values) + " |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=Path(CATALOG))
    parser.add_argument("--format", choices=("markdown", "json"), default="markdown")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    try:
        rows = validate_catalog(root, json.loads((root / args.catalog).read_text()))
    except (OSError, ValueError, KeyError, TypeError, tarfile.TarError) as error:
        parser.exit(2, f"capability catalog: {error}\n")
    print(
        json.dumps(rows, indent=2, ensure_ascii=False)
        if args.format == "json"
        else render_catalog(rows),
        end="\n",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
