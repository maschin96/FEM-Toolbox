"""Behavioral checks for documentation gates without a TeX dependency."""

import json
import runpy
from pathlib import Path

import pytest

API = runpy.run_path(str(Path(__file__).parents[2] / "scripts/build_handbook.py"))


@pytest.mark.parametrize(
    "message",
    [
        "LaTeX Warning: Reference `missing' on page 1 undefined on input line 2.",
        "LaTeX Warning: Citation `unknown' on page 1 undefined on input line 4.",
        'Warning--I didn\'t find a database entry for "unknown"',
        "! LaTeX Error: File `absent.tex' not found.",
        "LaTeX Warning: There were undefined references.",
        "! Emergency stop.",
    ],
)
def test_fatal_diagnostics_are_rejected(message):
    with pytest.raises(ValueError, match="documentation build"):
        API["check_log"](message)


def test_layout_warnings_remain_visible_but_are_not_fatal():
    warning = "Overfull \\hbox (3.0pt too wide) in paragraph at lines 12--14"
    assert API["check_log"](warning) == [warning]
    assert API["check_log"]("Output written on main.pdf (10 pages).") == []


def test_matrix_detects_uncovered_sources_and_missing_targets(tmp_path):
    (tmp_path / "Doc/latex/chapters").mkdir(parents=True)
    (tmp_path / "src/fem_toolbox/mesh").mkdir(parents=True)
    for name in ["README.md", "INSTALL.md", "Doc/README.md"]:
        (tmp_path / name).write_text("source")
    (tmp_path / "src/fem_toolbox/mesh/__init__.py").write_text("__all__ = ()")
    target = "Doc/latex/chapters/api.tex"
    (tmp_path / target).write_text("chapter")
    matrix = {
        "entries": [
            {
                "source": name,
                "target": target,
                "canonical": name,
                "owner": "API maintainers",
                "status": "open",
            }
            for name in [
                "README.md",
                "INSTALL.md",
                "Doc/README.md",
                "src/fem_toolbox/mesh/__init__.py",
            ]
        ]
    }
    API["validate_matrix"](tmp_path, matrix)
    (tmp_path / "Doc/new.md").write_text("new content")
    with pytest.raises(ValueError, match="uncovered"):
        API["validate_matrix"](tmp_path, matrix)
    matrix["entries"].append(
        {
            "source": "Doc/new.md",
            "target": "absent.tex",
            "canonical": "Doc/new.md",
            "owner": "docs",
            "status": "open",
        }
    )
    with pytest.raises(ValueError, match="target"):
        API["validate_matrix"](tmp_path, matrix)


def test_repository_matrix_is_complete():
    root = Path(__file__).parents[2]
    matrix = json.loads((root / "Doc/latex/coverage.json").read_text())
    API["validate_matrix"](root, matrix)


@pytest.mark.parametrize(
    "message",
    [
        "LaTeX Warning: Label `eq:a' multiply defined.",
        "LaTeX Warning: Rerun to get cross-references right.",
    ],
)
def test_unstable_and_duplicate_references_are_rejected(message):
    with pytest.raises(ValueError):
        API["check_log"](message)


def test_planned_topics_have_existing_targets_and_owners():
    root = Path(__file__).parents[2]
    matrix = json.loads((root / "Doc/latex/coverage.json").read_text())
    matrix["planned"][0]["target"] = "missing.tex"
    with pytest.raises(ValueError, match="planned"):
        API["validate_matrix"](root, matrix)
