"""Source-bound V&V tables and historical benchmark figures; no solver imports."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from itertools import product
from pathlib import Path

BASELINE = "benchmarks/baselines/m7-37c-linux-arm64.json"
MAP = "Doc/latex/verification-map.json"
DATA = "Doc/latex/verification-data.json"
TABLE = "Doc/latex/chapters/verification/measured.tex"
MATRIX = "Doc/latex/chapters/verification/matrix.tex"
FIGURE = "Doc/latex/assets/verification-accuracy-cost.pdf"
UQ = "Doc/latex/evidence/uq-linear.json"
UQ_FIGURE = "Doc/latex/assets/verification-uncertainty.pdf"
UQ_TABLE = "Doc/latex/chapters/verification/uq-measured.tex"
RAW_COPY = "Doc/latex/evidence/m7-37c-linux-arm64.json"
REQUIRED = {
    "heat",
    "transient_heat",
    "electrostatics",
    "elasticity",
    "thermomechanics",
    "piezoelectric",
    "thermopiezoelectric",
    "modal",
    "newmark",
    "harmonic",
    "geometry",
    "backends",
    "indicator",
    "studies",
    "surrogates",
    "uncertainty",
    "cost",
}


def read(root, name):
    return json.loads((root / name).read_text())


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_baseline(data):
    """Recalculate archived statistics and selections using the fixed M7 contract."""
    if data["schema"] != "fem-toolbox-accuracy-cost-v1":
        raise ValueError("benchmark schema")
    rows = data["samples"]
    expected = set(
        product(("heat", "bending"), ("scipy", "petsc"), (1, 2), (2, 4, 6), range(3))
    )
    identities = [
        (r["case"], r["backend"], r["order"], r["size"], r["repeat"]) for r in rows
    ]
    if len(identities) != 72 or set(identities) != expected:
        raise ValueError("benchmark missing/duplicate measurement")
    for row in rows:
        for field in (
            "wall_s",
            "h1_relative",
            "l2_relative",
            "residual",
            "balance",
            "parity_max_absolute",
        ):
            if isinstance(row[field], bool) or not math.isfinite(row[field]):
                raise ValueError("benchmark nonfinite measurement")
        if row["wall_s"] <= 0 or row["rss_bytes"] is None or row["rss_bytes"] <= 0:
            raise ValueError("benchmark missing time/RSS")
        if abs(row["residual"]) >= 1e-9 or abs(row["balance"]) >= 1e-9:
            raise ValueError("benchmark failed physical check")
    groups = []
    for key in sorted({r[:4] for r in identities}):
        samples = [
            row
            for identity, row in zip(identities, rows, strict=True)
            if identity[:4] == key
        ]
        times = [r["wall_s"] for r in samples]
        groups.append(
            dict(
                zip(("case", "backend", "order", "size"), key, strict=True),
                repeats=3,
                h1_relative=max(r["h1_relative"] for r in samples),
                rss_max_bytes=max(r["rss_bytes"] for r in samples),
                wall_s={
                    "median": statistics.median(times),
                    "min": min(times),
                    "max": max(times),
                    "stdev": statistics.stdev(times),
                },
            )
        )
    accuracy, cost = [], []
    for key in sorted({r[:3] for r in identities}):
        candidates = [g for g in groups if (g["case"], g["backend"], g["order"]) == key]
        base = dict(zip(("case", "backend", "order"), key, strict=True))
        for target in (0.5, 0.25, 0.125):
            feasible = [g for g in candidates if g["h1_relative"] <= target]
            accuracy.append(
                {
                    **base,
                    "target": target,
                    "selected": min(feasible, key=lambda g: g["wall_s"]["median"])
                    if feasible
                    else None,
                }
            )
        for seconds in (1.0, 5.0, 20.0):
            feasible = [
                g
                for g in candidates
                if g["wall_s"]["median"] <= seconds
                and g["rss_max_bytes"] <= 512 * 1024**2
            ]
            cost.append(
                {
                    **base,
                    "seconds": seconds,
                    "memory_bytes": 512 * 1024**2,
                    "selected": min(feasible, key=lambda g: g["h1_relative"])
                    if feasible
                    else None,
                }
            )
    for name, value in (("groups", groups), ("accuracy", accuracy), ("cost", cost)):
        if data[name] != value:
            raise ValueError(f"benchmark altered {name}")
    return {
        "sample_count": len(rows),
        "groups": groups,
        "accuracy": accuracy,
        "cost": cost,
    }


def validate_map(root, mapping):
    identities = [r["id"] for r in mapping["entries"]]
    if (
        mapping["schema"] != 1
        or set(identities) != REQUIRED
        or len(identities) != len(REQUIRED)
    ):
        raise ValueError("verification coverage")
    for entry in mapping["entries"]:
        if not entry["tests"] or not entry["criterion"] or not entry["scope"]:
            raise ValueError("verification empty contract")
        for name in [entry["contract"], *entry["tests"]]:
            path = root / name
            if not path.resolve().is_relative_to(root.resolve()) or not path.is_file():
                raise ValueError(f"verification reference: {name}")


def escape(text):
    return "".join(
        {"&": r"\&", "%": r"\%", "_": r"\_", "#": r"\#"}.get(c, c) for c in text
    )


def render_matrix(mapping):
    lines = [
        "% Generated by scripts/verification_artifacts.py; review sources before regeneration."
    ]
    for row in mapping["entries"]:
        lines += [
            r"\subsection{" + escape(row["title"]) + "}",
            r"\label{sec:vv-" + row["id"] + "}",
            escape(row["scope"]) + ". " + escape(row["criterion"]),
            r"Fachvertrag: \path{" + row["contract"] + "}.",
            r"\begin{itemize}",
        ]
        lines += [r"\item \path{" + path + "}" for path in row["tests"]]
        lines += [r"\end{itemize}"]
    return "\n".join(lines) + "\n"


def render_table(report):
    lines = [
        "% Generated from unchanged M7 raw samples, not fresh timings.",
        r"\begin{longtable}{llrrrrr}",
        r"\caption{Archivierte M7-Messung: Median $\pm$ Stichprobenstandardabweichung, RSS-Maximum.}\label{tab:vv-cost}\\",
        r"Fall & Element & $N$ & $H^1$ rel. & SciPy / s & PETSc / s & RSS / MiB\\\hline",
        r"\endfirsthead",
        r"Fall & Element & $N$ & $H^1$ rel. & SciPy / s & PETSc / s & RSS / MiB\\\hline",
        r"\endhead",
    ]
    for case, order, size in product(("heat", "bending"), (1, 2), (2, 4, 6)):
        pairs = [
            next(
                g
                for g in report["groups"]
                if (g["case"], g["order"], g["size"], g["backend"])
                == (case, order, size, backend)
            )
            for backend in ("scipy", "petsc")
        ]
        a, b = pairs
        times = [
            f"${g['wall_s']['median']:.4f}\\pm{g['wall_s']['stdev']:.4f}$"
            for g in pairs
        ]
        lines.append(
            f"{'Wärme' if case == 'heat' else 'Biegung'} & Tet{4 if order == 1 else 10} & {size} & {max(a['h1_relative'], b['h1_relative']):.5f} & {times[0]} & {times[1]} & {a['rss_max_bytes'] / 1024**2:.1f}/{b['rss_max_bytes'] / 1024**2:.1f}"
            + r"\\"
        )
    return "\n".join([*lines, r"\end{longtable}"]) + "\n"


def render_uq(reference):
    summary = reference["summary"]
    lines = [
        r"% Generated; review sources before regenerating.",
        r"\begin{tabular}{lrr}\toprule Größe & Analytisch & Stichprobe\\\midrule",
        f"Mittelwert / K & 8 & {summary['mean']:.8f}" + r"\\",
        f"Varianz / K$^2$ & 13 & {summary['variance']:.8f}" + r"\\",
    ]
    for probability, value in zip(
        ("0.025", "0.5", "0.975"), summary["quantiles"], strict=True
    ):
        lines.append(f"Quantil {escape(probability)} / K & --- & {value:.8f}" + r"\\")
    return "\n".join([*lines, r"\bottomrule\end{tabular}"]) + "\n"


def source_report(root):
    mapping = read(root, MAP)
    validate_map(root, mapping)
    report = audit_baseline(read(root, BASELINE))
    names = {
        BASELINE,
        MAP,
        RAW_COPY,
        "scripts/verification_artifacts.py",
        "benchmarks/accuracy_cost.py",
        "Doc/55_TET10_ACCEPTANCE.md",
        "Doc/65_VALIDATION_DATA_CONTRACT.md",
        "Doc/66_M9_DATA_SEARCH_AND_RELEASE_OPTIONS.md",
    }
    names.update(
        {UQ, "src/fem_toolbox/uncertainty.py", "tests/unit/test_uncertainty.py"}
    )
    for entry in mapping["entries"]:
        names.update([entry["contract"], *entry["tests"]])
    return {
        "schema": 1,
        "measurement_date": "2026-09-28",
        "experimental_validation": "not_performed",
        "benchmark_environment": read(root, BASELINE)["environment"],
        "sources_sha256": {name: digest(root / name) for name in sorted(names)},
        "figure_sha256": digest(root / FIGURE),
        "uq_figure_sha256": digest(root / UQ_FIGURE),
        "uq_reference": read(root, UQ),
        **report,
    }


def draw_figure(root, report):
    # Optional authoring dependency only. Check/build/archive need stdlib + frozen PDF.
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4), constrained_layout=True)
    for ax, case in zip(axes, ("heat", "bending"), strict=True):
        for order, color in ((1, "#215d8a"), (2, "#af5c16")):
            for backend, style, marker in (("scipy", "-", "o"), ("petsc", "--", "s")):
                rows = sorted(
                    [
                        g
                        for g in report["groups"]
                        if (g["case"], g["order"], g["backend"])
                        == (case, order, backend)
                    ],
                    key=lambda r: r["size"],
                )
                ax.plot(
                    [r["wall_s"]["median"] for r in rows],
                    [r["h1_relative"] for r in rows],
                    style,
                    marker=marker,
                    color=color,
                    label=f"Tet{4 if order == 1 else 10} / {backend}",
                )
        ax.set(
            xscale="log",
            yscale="log",
            xlabel="Erster Solverlauf / s (Median)",
            ylabel="Relativer H¹-Halbnormfehler",
            title="Wärmegradient" if case == "heat" else "Synthetische Biegung",
        )
        ax.grid(alpha=0.2, which="both")
        ax.legend(fontsize=7)
    fig.savefig(
        root / FIGURE,
        metadata={
            "Creator": "FEM Toolbox / verification_artifacts.py",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)

    reference = read(root, UQ)
    prefixes = reference["summary"]["convergence"]
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), constrained_layout=True)
    for ax, field, exact, unit in zip(
        axes, ("mean", "variance"), (8, 13), ("K", "K²"), strict=True
    ):
        ax.plot(
            [r["size"] for r in prefixes],
            [r[field] for r in prefixes],
            "o-",
            color="#215d8a",
            markersize=3,
            label="PCG64 / Seed 641",
        )
        ax.axhline(exact, color="#af5c16", linestyle="--", label="Analytische Referenz")
        ax.set(
            xscale="log",
            xlabel="Gemeinsame Stichprobenzahl n",
            ylabel=f"{'Mittel' if field == 'mean' else 'Varianz'} / {unit}",
        )
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)
    fig.savefig(
        root / UQ_FIGURE,
        metadata={
            "Creator": "FEM Toolbox / verification_artifacts.py",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    plt.close(fig)


def collect_uq(root):
    """Reproduce the existing independently defined correlated linear test case."""
    import sys

    sys.path.insert(0, str(root / "src"))
    import numpy as np

    from fem_toolbox.uncertainty import (
        Normal,
        SamplingPlan,
        UncertainParameter,
        local_sensitivity,
        summarize_samples,
    )

    parameters = tuple(
        UncertainParameter(
            name,
            f"loads.{name}.value",
            "1",
            Normal(mean, sd),
            "Synthetic analytic verification; no measured material",
        )
        for name, mean, sd in (("x", 1, 2), ("y", -2, 3))
    )
    plan = SamplingPlan(
        parameters, size=100_000, seed=641, correlation=((1, 0.5), (0.5, 1))
    )
    samples = plan.sample()
    values = 4 + 2 * samples[:, 0] - samples[:, 1]
    summary = summarize_samples(values, unit="K", numerical_error_bound=0.01)
    if abs(summary["mean"] - 8) >= 6 * np.sqrt(13 / len(values)) or abs(
        summary["variance"] - 13
    ) >= 6 * 13 * np.sqrt(2 / (len(values) - 1)):
        raise ValueError("UQ existing analytical verification gate")
    sensitivity = local_sensitivity(
        lambda x: 4 + 2 * x[0] - x[1], center=[1, -2], steps=[0.01, 0.02]
    )
    np.testing.assert_allclose(sensitivity["gradient"], [2, -1], atol=1e-12, rtol=0)
    result = {
        "scope": "synthetic linear probability reference, no experiment",
        "formula": "4+2X-Y",
        "expected_mean_K": 8,
        "expected_variance_K2": 13,
        "sampling": plan.to_document(),
        "sample_sha256": hashlib.sha256(
            np.asarray(samples, dtype="<f8").tobytes()
        ).hexdigest(),
        "summary": summary,
        "sensitivity": sensitivity,
    }
    (root / UQ).parent.mkdir(parents=True, exist_ok=True)
    (root / UQ).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")


def check_artifacts(root):
    if (root / RAW_COPY).read_bytes() != (root / BASELINE).read_bytes():
        raise ValueError("verification raw archive drift")
    actual = source_report(root)
    if actual != read(root, DATA):
        raise ValueError("verification evidence drift: review sources and regenerate")
    if (root / TABLE).read_text() != render_table(actual):
        raise ValueError("verification table drift")
    if (root / MATRIX).read_text() != render_matrix(read(root, MAP)):
        raise ValueError("verification matrix drift")
    if (root / UQ_TABLE).read_text() != render_uq(read(root, UQ)):
        raise ValueError("verification UQ table drift")


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    if args.write:
        mapping = read(root, MAP)
        validate_map(root, mapping)
        report = audit_baseline(read(root, BASELINE))
        (root / TABLE).parent.mkdir(parents=True, exist_ok=True)
        collect_uq(root)
        (root / RAW_COPY).write_bytes((root / BASELINE).read_bytes())
        draw_figure(root, report)
        (root / DATA).write_text(
            json.dumps(source_report(root), indent=2, ensure_ascii=False) + "\n"
        )
        (root / TABLE).write_text(render_table(report))
        (root / MATRIX).write_text(render_matrix(mapping))
        (root / UQ_TABLE).write_text(render_uq(read(root, UQ)))
    check_artifacts(root)
    print(
        "V&V: 17 evidence groups / 72 historical measurements / sources and figure verified"
    )


if __name__ == "__main__":
    main()
