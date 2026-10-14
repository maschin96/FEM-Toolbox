"""Build the handbook with a pinned, isolated TeX toolchain (stdlib only)."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import runpy
import shutil
import subprocess
import tarfile
import tempfile
import tomllib
from pathlib import Path
from typing import Any

BOOKS = ("main", "theory", "application", "api", "verification")


def publication_metadata(
    release: bool, state: str, version: str, commit: str
) -> tuple[str, str, str]:
    """Label releases explicitly and reject incomplete release provenance."""
    if not release:
        return "candidate", "KANDIDAT", "Gesamtabnahme"
    if (
        state != "clean"
        or not re.fullmatch(r"\d+\.\d+\.\d+", version)
        or not re.fullmatch(r"[0-9a-f]{40}", commit)
    ):
        raise ValueError("release provenance requires clean source, version and commit")
    return "release", "RELEASESTAND", "Verifizierte Dokumentation"


def check_log(text: str) -> list[str]:
    """Reject unresolved content; preserve layout warnings for human review."""
    fatal = re.search(
        r"^!|(?:Reference|Citation).*undefined|undefined references|"
        r"I didn't find a database|There were undefined citations|"
        r"Rerun to get cross-references right|multiply[ -]defined|"
        r"Fatal error|Emergency stop",
        text,
        re.MULTILINE,
    )
    if fatal:
        raise ValueError(f"documentation build: {fatal.group(0)}")
    return [line for line in text.splitlines() if "full \\" in line]


def required_sources(root: Path) -> set[str]:
    """Inventory Markdown, supporting assets and every public API group."""
    sources = {"README.md", "INSTALL.md"}
    sources.update(str(p.relative_to(root)) for p in (root / "Doc").rglob("*.md"))
    sources.update(
        str(p.relative_to(root))
        for p in (root / "Doc").rglob("*")
        if p.is_file() and "latex" not in p.parts
    )
    package = root / "src/fem_toolbox"
    sources.update(str(p.relative_to(root)) for p in package.glob("*.py"))
    sources.update(str(p.relative_to(root)) for p in package.rglob("__init__.py"))
    return sources


def validate_matrix(root: Path, matrix: dict[str, Any]) -> None:
    """Fail closed on newly uncovered sources and dangling ownership entries."""
    entries = matrix["entries"]
    sources = [entry["source"] for entry in entries]
    if len(sources) != len(set(sources)):
        raise ValueError("duplicate coverage source")
    missing = required_sources(root) - set(sources)
    if missing:
        raise ValueError(f"uncovered sources: {sorted(missing)}")
    for entry in entries:
        for field in ("source", "canonical", "target"):
            path = root / entry[field]
            if not path.resolve().is_relative_to(root.resolve()) or not path.is_file():
                raise ValueError(f"missing/invalid {field}: {entry[field]}")
        if not entry["owner"] or entry["status"] not in {
            "open",
            "sample",
            "historical",
        }:
            raise ValueError(f"missing ownership/status: {entry['source']}")

    for topic in matrix.get("planned", []):
        target = root / topic["target"]
        if (
            not target.resolve().is_relative_to(root.resolve())
            or not target.is_file()
            or not topic["owner"]
        ):
            raise ValueError(f"invalid planned topic: {topic}")


def tex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(character, character) for character in text)


def run(command: list[str], cwd: Path, log: Path | None = None) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    output = result.stdout + result.stderr
    if log is not None:
        log.write_text(output, encoding="utf-8")
    if result.returncode:
        raise ValueError(
            f"documentation build: {command[0]} failed; {log}\n{output[-1500:]}"
        )
    return output


def compile_book(work: Path, book: str) -> list[str]:
    command = [
        "pdflatex",
        "-no-shell-escape",
        "-halt-on-error",
        "-interaction=nonstopmode",
        f"{book}.tex",
    ]
    run(command, work, work / f"{book}-pass1.txt")
    bibliography = run(["bibtex", book], work, work / f"{book}-bibtex.txt")
    check_log(bibliography)
    run(["makeindex", f"{book}.idx"], work, work / f"{book}-index.txt")
    run(command, work, work / f"{book}-pass2.txt")
    run(command, work, work / f"{book}-pass3.txt")
    return check_log(
        (work / f"{book}.log").read_text(encoding="utf-8", errors="replace")
    )


def negative_checks(work: Path, output: Path) -> dict[str, str]:
    results = {}
    for name, invalid in {
        "missing-file": r"\input{deliberately-absent}",
        "missing-reference": r"\ref{deliberately-absent}",
        "missing-citation": r"\cite{deliberately-absent}",
    }.items():
        with tempfile.TemporaryDirectory(dir=work.parent) as directory:
            fixture = Path(directory) / "source"
            shutil.copytree(work, fixture)
            chapter = fixture / "chapters/theory.tex"
            chapter.write_text(chapter.read_text() + "\n" + invalid + "\n")
            try:
                compile_book(fixture, "theory")
            except ValueError as error:
                results[name] = str(error).splitlines()[0]
            else:
                raise ValueError(f"documentation build: negative fixture {name} passed")
            for logfile in fixture.glob("theory*.txt"):
                shutil.copy2(logfile, output / f"{name}-{logfile.name}")
    return results


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def native_build(root: Path, output: Path, args: argparse.Namespace) -> None:
    status, label, subtitle = publication_metadata(
        args.release, args.state, args.version, args.commit
    )
    output.mkdir(parents=True, exist_ok=True)
    # Never ship a prior success manifest after a failed rebuild.
    (output / "manifest.json").unlink(missing_ok=True)
    source = root / "Doc/latex"
    matrix = json.loads((source / "coverage.json").read_text())
    validate_matrix(root, matrix)
    os.environ["SOURCE_DATE_EPOCH"] = args.epoch
    os.environ["FORCE_SOURCE_DATE"] = "1"
    with tempfile.TemporaryDirectory(dir=output) as directory:
        work = Path(directory) / "source"
        shutil.copytree(source, work)
        metadata = (
            rf"\newcommand{{\ToolboxVersion}}{{{tex_escape(args.version)}}}"
            + "\n"
            + rf"\newcommand{{\SourceCommit}}{{{tex_escape(args.commit)}}}"
            + "\n"
            + rf"\newcommand{{\SourceState}}{{{tex_escape(args.state)}}}"
            + "\n"
            + rf"\newcommand{{\PublicationLabel}}{{{label}}}"
            + "\n"
            + rf"\newcommand{{\PublicationSubtitle}}{{{subtitle}}}"
            + "\n"
            + r"\newif\ifReleasePublication"
            + "\n"
            + (
                r"\ReleasePublicationtrue"
                if args.release
                else r"\ReleasePublicationfalse"
            )
            + "\n"
        )
        (work / "metadata.tex").write_text(metadata)
        shutil.copy2(root / "LICENSE", work / "LICENSE")
        # Package the exact inputs before any TeX intermediates are generated.
        with tarfile.open(output / "handbook-source.tar.gz", "w:gz") as archive:
            for path in sorted(work.rglob("*")):
                if path.is_file():
                    archive.add(path, arcname="handbook/" + str(path.relative_to(work)))
        artifacts = {}
        for book in BOOKS:
            try:
                warnings = compile_book(work, book)
            except ValueError:
                for logfile in work.glob(f"{book}*"):
                    if logfile.suffix in {".txt", ".log", ".blg"}:
                        shutil.copy2(logfile, output / logfile.name)
                raise
            shutil.copy2(work / f"{book}.pdf", output / f"{book}.pdf")
            text = run(["pdftotext", f"{book}.pdf", "-"], work)
            for marker in (
                args.version,
                args.commit,
                label,
                "Wärme",
                "Literatur",
                "Sachregister",
            ):
                if marker not in text:
                    raise ValueError(
                        f"documentation build: {book}: missing PDF text {marker}"
                    )
            if book in {"main", "theory"} and not all(
                symbol in text for symbol in ("∇", "∫")
            ):
                raise ValueError(
                    f"documentation build: {book}: mathematical Unicode missing"
                )
            if b"/Outlines" not in (work / f"{book}.pdf").read_bytes():
                raise ValueError(f"documentation build: {book}: missing PDF bookmarks")
            fonts = run(["pdffonts", f"{book}.pdf"], work)
            font_rows = fonts.splitlines()[2:]
            if not font_rows or any(
                not re.search(r"\byes\s+yes\s+yes\b", row) for row in font_rows
            ):
                raise ValueError(
                    f"documentation build: {book}: fonts not embedded/searchable\n{fonts}"
                )
            for suffix in ("pdf", "log"):
                shutil.copy2(work / f"{book}.{suffix}", output / f"{book}.{suffix}")
            (output / f"{book}.txt").write_text(text)
            artifacts[book] = {
                "sha256": sha256(output / f"{book}.pdf"),
                "warnings": warnings,
                "pdfinfo": run(["pdfinfo", f"{book}.pdf"], work),
                "fonts": fonts,
            }
        negatives = negative_checks(work, output) if args.self_test else {}
        with tempfile.TemporaryDirectory(dir=output) as extracted_directory:
            extracted = Path(extracted_directory)
            run(
                [
                    "tar",
                    "-xzf",
                    str(output / "handbook-source.tar.gz"),
                    "-C",
                    str(extracted),
                ],
                root,
            )
            run(
                ["make", "-C", str(extracted / "handbook")],
                root,
                output / "source-rebuild.log",
            )
            for book in BOOKS:
                rebuilt_text = run(
                    ["pdftotext", f"{book}.pdf", "-"], extracted / "handbook"
                )
                if rebuilt_text != (output / f"{book}.txt").read_text():
                    raise ValueError(
                        f"documentation build: source archive differs: {book}"
                    )
    inputs = {
        str(p.relative_to(root)): sha256(p)
        for p in sorted(source.rglob("*"))
        if p.is_file()
    }
    manifest = {
        "schema": 1,
        "status": status,
        "version": args.version,
        "commit": args.commit,
        "source_state": args.state,
        "source_date_epoch": args.epoch,
        "reproducibility": "content; byte identity not claimed across architectures/toolchains",
        "inputs": inputs,
        "artifacts": artifacts,
        "negative_checks": negatives,
        "source_archive_verified": True,
        "source_archive_sha256": sha256(output / "handbook-source.tar.gz"),
        "toolchain": run(["dpkg-query", "-W"], root),
        "engines": {
            tool: run(
                [tool, "-v" if tool == "pdftotext" else "--version"], root
            ).splitlines()[0]
            for tool in ("pdflatex", "bibtex", "pdftotext")
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    )
    checksums = "".join(
        f"{sha256(p)}  {p.name}\n"
        for p in sorted(output.iterdir())
        if p.is_file()
        and p.name not in {"SHA256SUMS", "build-console.log", "toolchain-build.log"}
    )
    (output / "SHA256SUMS").write_text(checksums)
    print(f"Handbook: {len(BOOKS)} PDFs, {len(negatives)} negative checks; {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--native", action="store_true", help="inside pinned toolchain only"
    )
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--release", action="store_true", help="label clean release inputs"
    )
    parser.add_argument("--check-matrix", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("build/handbook"))
    parser.add_argument("--commit", default="")
    parser.add_argument("--version", default="")
    parser.add_argument("--state", default="")
    parser.add_argument("--epoch", default="")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runpy.run_path(str(root / "scripts/api_contracts.py"))["validate_reference"](root)
    runpy.run_path(str(root / "scripts/verification_artifacts.py"))["check_artifacts"](
        root
    )
    catalog_api = runpy.run_path(str(root / "scripts/capability_catalog.py"))
    catalog_api["validate_catalog"](
        root, json.loads((root / catalog_api["CATALOG"]).read_text())
    )
    matrix = json.loads((root / "Doc/latex/coverage.json").read_text())
    validate_matrix(root, matrix)
    if args.check_matrix:
        print(f"Coverage: {len(matrix['entries'])} owned sources")
        return
    output = args.output.resolve()
    if args.native:
        if not all((args.commit, args.version, args.state, args.epoch)):
            parser.error("native build requires explicit provenance")
        native_build(root, output, args)
        return
    args.commit = run(["git", "rev-parse", "HEAD"], root).strip()
    args.epoch = run(["git", "show", "-s", "--format=%ct", "HEAD"], root).strip()
    args.state = (
        "dirty" if run(["git", "status", "--porcelain"], root).strip() else "clean"
    )
    args.version = tomllib.loads((root / "pyproject.toml").read_text())["project"][
        "version"
    ]
    publication_metadata(args.release, args.state, args.version, args.commit)
    output.mkdir(parents=True, exist_ok=True)
    run(
        [
            "docker",
            "build",
            "-f",
            "Doc/latex/Dockerfile",
            "-t",
            "fem-handbook:206",
            ".",
        ],
        root,
        output / "toolchain-build.log",
    )
    command = [
        "docker",
        "run",
        "--rm",
        "--network",
        "none",
        "--user",
        f"{os.getuid()}:{os.getgid()}",
        "-v",
        f"{root}:/work:ro",
        "-v",
        f"{output}:/output",
        "fem-handbook:206",
        "python3",
        "scripts/build_handbook.py",
        "--native",
        "--output",
        "/output",
        "--commit",
        args.commit,
        "--version",
        args.version,
        "--state",
        args.state,
        "--epoch",
        args.epoch,
    ]
    if args.self_test:
        command.append("--self-test")
    if args.release:
        command.append("--release")
    print(run(command, root, output / "build-console.log"))


if __name__ == "__main__":
    main()
