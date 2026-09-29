"""Inspect a FEM Toolbox result manifest using only the Python standard library."""

from __future__ import annotations

import argparse
import json
from hashlib import sha256
from pathlib import Path


def read_bundle(path: Path) -> dict:
    """Read versions 1/2 and verify the companion VTU checksum."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("bundle_schema_version") not in (1, 2):
        raise ValueError("unsupported result bundle schema")
    files = data["files"]
    name = files["fields"]
    if (
        not isinstance(name, str)
        or Path(name).name != name
        or not name.endswith(".vtu")
    ):
        raise ValueError("invalid VTU filename")
    actual = sha256((path.parent / name).read_bytes()).hexdigest()
    if actual != files["fields_sha256"]:
        raise ValueError("VTU checksum mismatch")
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    bundle = read_bundle(args.manifest)
    print(
        json.dumps(
            {
                "physics": bundle["physics"],
                "input_sha256": bundle["provenance"]["input_sha256"],
                "metrics": bundle["metrics"],
            },
            sort_keys=True,
            allow_nan=False,
        )
    )
