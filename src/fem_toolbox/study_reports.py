"""Deterministic SI quantity tables and reports from stored stationary studies."""

from __future__ import annotations

import csv
import io
import json
import re
from collections.abc import Callable, Sequence
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Literal

import numpy as np

from fem_toolbox.studies import _inside, _json, _load_catalog, _valid_result


@dataclass(frozen=True)
class Quantity:
    """One scalar output definition; indices are global nodes or volume-cell rows.

    Metric fields use dotted bundle paths (e.g. heat.power_balance_W). Point/cell
    fields use exported array names. Select a vector component by its metadata
    name or 'norm'. No selection means all active nodes/volume cells; no padding
    or smoothing. scalar requires exactly one selected value. Tolerances, if
    supplied, evaluate abs(value) <= tolerance in this quantity's SI unit.
    """

    name: str
    source: Literal["metric", "point", "cell"]
    field: str
    unit: str
    definition: str
    reduction: Literal["scalar", "min", "max", "mean", "sum", "max_abs"] = "scalar"
    component: str = "scalar"
    indices: tuple[int, ...] | None = None
    absolute_tolerance: float | None = None

    def __post_init__(self) -> None:
        if not re.fullmatch("[A-Za-z][A-Za-z0-9_]*", self.name) or self.name in {
            "case_id",
            "status",
        }:
            raise ValueError("quantity name must be a non-reserved identifier")
        if (
            self.source not in {"metric", "point", "cell"}
            or not self.field
            or not self.unit
            or not self.definition
        ):
            raise ValueError(
                "quantity source, field, SI unit and definition are required"
            )
        if self.reduction not in {"scalar", "min", "max", "mean", "sum", "max_abs"}:
            raise ValueError("unknown reduction")
        if self.indices is not None:
            if (
                not self.indices
                or any(type(i) is not int or i < 0 for i in self.indices)
                or len(set(self.indices)) != len(self.indices)
            ):
                raise ValueError("indices must be unique nonnegative integers")
            object.__setattr__(self, "indices", tuple(self.indices))
        if self.absolute_tolerance is not None and (
            not np.isfinite(self.absolute_tolerance) or self.absolute_tolerance < 0
        ):
            raise ValueError("balance tolerance must be finite and nonnegative")


def _metric_unit(path: str) -> str:
    name = path.split(".")[-1]
    if name.startswith("relative_residual"):
        return "1"
    for suffix, unit in [
        ("_N_m", "N m"),
        ("_C_V", "C V"),
        ("_W", "W"),
        ("_J", "J"),
        ("_N", "N"),
        ("_C", "C"),
        ("_V", "V"),
    ]:
        if name.endswith(suffix):
            return unit
    raise ValueError("metric has no supported SI unit contract")


def _reduce(values: Any, q: Quantity, components: list[str]) -> float:
    a = np.asarray(values, dtype=float)
    if not np.isfinite(a).all() or not a.size:
        raise ValueError("empty or nonfinite quantity")
    if a.ndim == 2:
        if q.component == "norm":
            a = np.linalg.norm(a, axis=1)
        elif q.component in components:
            a = a[:, components.index(q.component)]
        else:
            raise ValueError("select an explicit vector/tensor component or norm")
    elif a.ndim == 1 and q.component != "scalar":
        raise ValueError("scalar quantity has no requested component")
    if q.reduction == "scalar":
        if a.size != 1:
            raise ValueError("scalar reduction requires one value")
        return float(a.reshape(-1)[0])
    operations: dict[str, Callable[[Any], Any]] = {
        "min": np.min,
        "max": np.max,
        "mean": np.mean,
        "sum": np.sum,
        "max_abs": lambda v: np.max(np.abs(v)),
    }
    result = float(operations[q.reduction](a))
    if not np.isfinite(result):
        raise ValueError("nonfinite reduction")
    return result


def _quantity(
    bundle: dict[str, Any], mesh: Any, q: Quantity
) -> tuple[float, dict[str, Any]]:
    if q.source == "metric":
        value: Any = bundle["metrics"]
        for key in q.field.split("."):
            value = value[int(key)] if isinstance(value, list) else value[key]
        if q.unit != _metric_unit(q.field):
            raise ValueError("metric unit mismatch")
        a = np.atleast_1d(value)
        if q.component != "scalar":
            if a.shape != (3,) or q.component not in ("x", "y", "z", "norm"):
                raise ValueError("invalid metric component")
            a = np.array(
                [
                    np.linalg.norm(a)
                    if q.component == "norm"
                    else a[("x", "y", "z").index(q.component)]
                ]
            )
        if q.indices is not None:
            a = a[list(q.indices)]
        # Component selection has already occurred for vector metrics.
        scalar = Quantity(q.name, q.source, q.field, q.unit, q.definition, q.reduction)
        return _reduce(a, scalar, ["scalar"]), {
            "kind": "solver_metric",
            "source": q.field,
        }
    location = (
        "point"
        if q.source == "point"
        else ("tetra10" if "tetra10" in bundle["fields"] else "tetra4")
    )
    meta = bundle["fields"][location][q.field]
    if meta["unit"] != q.unit:
        raise ValueError("field unit mismatch")
    if q.source == "point":
        active = meta["active_node_ids"]
        selected = active if q.indices is None else list(q.indices)
        if not set(selected).issubset(active):
            raise ValueError("selection includes inactive nodes")
        values = mesh.point_data[q.field][selected]
    else:
        cell_type = "tetra10" if location == "tetra10" else "tetra"
        values = mesh.cell_data_dict[q.field][cell_type]
        if q.indices is not None:
            values = values[list(q.indices)]
    return _reduce(values, q, meta["components"]), {
        "kind": meta["kind"],
        "location": location,
        "source": q.field,
    }


def export_study_table(
    catalog_path: str | Path, output: str | Path, quantities: Sequence[Quantity]
) -> Path:
    """Write values.csv, metadata.json and report.md into a new output directory.

    All active catalog cases are retained; unavailable quantities are null with
    reasons, encoded as empty CSV cells. Existing directories are never replaced.
    Invalid bundles remain visible; models are never solved. Returned path points
    to metadata.json. A failed write can leave an incomplete directory, which
    must not be consumed until metadata.json and the CSV hash validate.
    """
    import meshio  # type: ignore[import-untyped]

    path = Path(catalog_path).resolve()
    catalog = _load_catalog(path)
    if len({q.name for q in quantities}) != len(quantities):
        raise ValueError("duplicate quantity names")
    columns: dict[str, Any] = {}
    for case in catalog["cases"]:
        for name, p in case["parameters"].items():
            key = "parameter:" + name
            meta = {"role": "parameter", "unit": p["unit"], "definition": p["target"]}
            if key in columns and columns[key] != meta:
                raise ValueError("inconsistent parameter metadata")
            columns[key] = meta
    columns = dict(sorted(columns.items()))
    for q in quantities:
        columns[q.name] = {"role": "quantity", **asdict(q)}
    rows = []
    for case in catalog["cases"]:
        row: dict[str, Any] = {
            "case_id": case["case_id"],
            "status": case["status"],
            "values": {},
            "missing": {},
            "origins": {},
            "fingerprints": case["fingerprints"],
            "versions": case["versions"],
            "result": case["result"],
            "solver_status": "not_available",
            "balance_status": "not_checked",
            "balance_checks": {},
            "mesh_convergence": "not_checked",
            "error": case["error"],
        }
        for name in columns:
            if name.startswith("parameter:"):
                p = case["parameters"].get(name.removeprefix("parameter:"))
                row["values"][name] = None if p is None else p["value"]
                if p is None:
                    row["missing"][name] = "parameter not available"
        bundle = mesh = None
        if case["status"] == "succeeded":
            if not _valid_result(path.parent, case):
                row.update(status="invalid", error="missing or corrupt bundle")
            else:
                try:
                    bundle = json.loads(
                        _inside(path.parent, case["result"]["manifest"]).read_text()
                    )
                    mesh = meshio.read(_inside(path.parent, case["result"]["fields"]))
                    row["solver_status"] = bundle["solver_status"]
                    row["provenance"] = bundle["provenance"]
                except (OSError, ValueError, KeyError) as error:
                    bundle = None
                    row.update(status="invalid", error=str(error))
        for q in quantities:
            try:
                if bundle is None:
                    raise ValueError(row["error"] or f"case {row['status']}")
                value, origin = _quantity(bundle, mesh, q)
                row["values"][q.name] = value
                row["origins"][q.name] = origin
                if q.absolute_tolerance is not None:
                    row["balance_checks"][q.name] = (
                        "passed" if abs(value) <= q.absolute_tolerance else "failed"
                    )
            except (ValueError, KeyError, IndexError, TypeError) as error:
                row["values"][q.name] = None
                row["missing"][q.name] = str(error)
                if q.absolute_tolerance is not None:
                    row["balance_checks"][q.name] = "not_available"
        checks = list(row["balance_checks"].values())
        if checks:
            row["balance_status"] = (
                "failed"
                if "failed" in checks
                else "not_available"
                if "not_available" in checks
                else "passed"
            )
        rows.append(row)
    table: dict[str, Any] = {
        "table_schema_version": 1,
        "columns": columns,
        "rows": rows,
        "catalog_sha256": sha256(path.read_bytes()).hexdigest(),
        "limitations": [
            "No mesh convergence or experimental validation implied",
            "Reductions are unweighted; no smoothing or integration of field samples",
        ],
    }
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(["case_id", "status", *columns])
    for row in rows:
        writer.writerow(
            [
                row["case_id"],
                row["status"],
                *[_cell(row["values"][name]) for name in columns],
            ]
        )
    csv_bytes = stream.getvalue().encode("utf-8")
    table["csv_sha256"] = sha256(csv_bytes).hexdigest()
    target = Path(output)
    target.mkdir(parents=True, exist_ok=False)
    (target / "values.csv").write_bytes(csv_bytes)
    report = ["# Stationary study comparison", "", "Units and definitions:", ""]
    report += [
        f"- `{name}` [{meta['unit']}]: {meta['definition']}"
        for name, meta in columns.items()
    ]
    report += [
        "",
        "| Case | Status | Solver | Balance | Mesh convergence | Values / missing reasons |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        values = "; ".join(
            f"{name}={row['values'][name] if row['values'][name] is not None else 'unavailable: ' + row['missing'][name]}"
            for name in columns
        )
        report.append(
            "| "
            + " | ".join(
                str(v).replace("|", "\\|").replace("\n", " ")
                for v in [
                    row["case_id"],
                    row["status"],
                    row["solver_status"],
                    row["balance_status"],
                    row["mesh_convergence"],
                    values,
                ]
            )
            + " |"
        )
    report += [
        "",
        "See metadata.json for case fingerprints, result references, reductions and provenance.",
        *table["limitations"],
    ]
    (target / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (target / "metadata.json").write_text(_json(table) + "\n", encoding="utf-8")
    return target / "metadata.json"


def _cell(value: Any) -> str:
    return "" if value is None else _json(value)


def read_study_table(directory: str | Path) -> dict[str, Any]:
    """Read schema-1 metadata and verify the matching CSV, without FEM execution.

    Raises ValueError for unsupported/corrupt data, OSError for missing files.
    Returned JSON data includes column contracts, per-case values and provenance.
    """
    path = Path(directory)
    data = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    if (
        type(data.get("table_schema_version")) is not int
        or data["table_schema_version"] != 1
    ):
        raise ValueError("unsupported table schema")
    raw = (path / "values.csv").read_bytes()
    if sha256(raw).hexdigest() != data["csv_sha256"]:
        raise ValueError("table CSV hash mismatch")
    # JSON sorting does not preserve insertion order: recover CSV column order.
    csv_rows = list(csv.reader(io.StringIO(raw.decode("utf-8"))))
    if not csv_rows or csv_rows[0][:2] != ["case_id", "status"]:
        raise ValueError("invalid CSV header")
    names = csv_rows[0][2:]
    if (
        len(set(names)) != len(names)
        or set(names) != set(data["columns"])
        or len(csv_rows) != len(data["rows"]) + 1
    ):
        raise ValueError("CSV/metadata columns or rows disagree")
    ids = set()
    for actual, row in zip(csv_rows[1:], data["rows"], strict=True):
        if row["case_id"] in ids:
            raise ValueError("duplicate table case IDs")
        ids.add(row["case_id"])
        expected = [
            row["case_id"],
            row["status"],
            *[_cell(row["values"][n]) for n in names],
        ]
        if actual != expected:
            raise ValueError("CSV/metadata values disagree")
    return dict(data)
