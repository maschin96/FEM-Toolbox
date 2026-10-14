"""Build current NumPy for the runner CPU and record source/wheel provenance."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
import urllib.request
import zipfile
from email.parser import BytesParser
from pathlib import Path


def run(*arguments: str) -> None:
    subprocess.run([sys.executable, "-m", "pip", *arguments], check=True)


def main() -> None:
    output = Path(os.environ["RUNNER_TEMP"]) / "numpy-wheels"
    source = output / "source"
    source.mkdir(parents=True)
    # Let pip resolve the same current version as the ordinary binary install,
    # without importing it or compiling metadata with NumPy's default ISA.
    selection = output / "selection"
    run(
        "download",
        "--no-deps",
        "--only-binary=:all:",
        "--dest",
        str(selection),
        "numpy",
    )
    selected = list(selection.glob("numpy-*.whl"))
    if len(selected) != 1:
        raise ValueError("exactly one resolved NumPy wheel required")
    with zipfile.ZipFile(selected[0]) as archive_metadata:
        names = [
            name
            for name in archive_metadata.namelist()
            if name.endswith(".dist-info/METADATA")
        ]
        if len(names) != 1:
            raise ValueError("exactly one NumPy metadata record required")
        version = BytesParser().parsebytes(archive_metadata.read(names[0]))["Version"]
    with urllib.request.urlopen(
        f"https://pypi.org/pypi/numpy/{version}/json", timeout=60
    ) as response:
        metadata = json.load(response)
    sdists = [record for record in metadata["urls"] if record["packagetype"] == "sdist"]
    if len(sdists) != 1:
        raise ValueError("exactly one official NumPy source archive required")
    record = sdists[0]
    archive = source / Path(record["filename"]).name
    with urllib.request.urlopen(record["url"], timeout=60) as response:
        archive.write_bytes(response.read())
    if hashlib.sha256(archive.read_bytes()).hexdigest() != record["digests"]["sha256"]:
        raise ValueError("NumPy source does not match official PyPI hash")
    run(
        "wheel",
        "--no-deps",
        "--no-cache-dir",
        "--wheel-dir",
        str(output),
        "--config-settings=setup-args=-Dcpu-baseline=none",
        "--config-settings=setup-args=-Dcpu-dispatch=none",
        "--config-settings=setup-args=-Dallow-noblas=true",
        "--config-settings=compile-args=-j2",
        str(archive),
    )
    wheels = list(output.glob("numpy-*.whl"))
    if len(wheels) != 1:
        raise ValueError("exactly one runner-compatible NumPy wheel required")
    wheel = wheels[0]
    run("install", str(wheel))
    subprocess.run(
        [
            sys.executable,
            "-I",
            "-c",
            "import numpy as n; "
            "n.testing.assert_allclose(n.linalg.eigvalsh(n.diag([1., 2.])), [1., 2.]); "
            "n.testing.assert_allclose(n.eye(3) @ n.arange(3.), n.arange(3.)); "
            "print(n.__version__); n.show_config()",
        ],
        check=True,
    )
    receipt = {
        "schema": 1,
        "python": platform.python_version(),
        "source": {
            "name": archive.name,
            "sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        },
        "wheel": {
            "name": wheel.name,
            "sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        },
        "cpu_baseline": "none",
        "cpu_dispatch": "none",
        "blas": "numpy bundled fallback",
        "scope": "CI dependency only; no NumPy wheel is part of the FEM release",
    }
    info = output / "BUILD_INFO.json"
    info.write_text(json.dumps(receipt, indent=2) + "\n")
    (output / "SHA256SUMS").write_text(
        "".join(
            f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n"
            for path in (wheel, info)
        )
    )


if __name__ == "__main__":
    main()
