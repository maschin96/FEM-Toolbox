"""Inspectable polynomial surrogates for real scalar FEM quantities.

No FEM solve, ML dependency, hyperparameter search or uncertainty claim is made.
All arrays use sample rows and explicitly ordered parameter/output columns.
"""

from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis import __version__


def _canonical(value: object) -> str:
    try:
        return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ValueError("provenance/artifact must contain finite JSON data") from error


def _digest(value: object) -> str:
    return sha256(_canonical(value).encode("utf-8")).hexdigest()


def _array(value: ArrayLike, *, ndim: int, name: str) -> NDArray[np.float64]:
    raw = np.asarray(value)
    if raw.dtype.kind not in "iuf" or raw.ndim != ndim or not raw.size:
        raise ValueError(f"{name} must be a nonempty real numeric {ndim}D array")
    a = np.asarray(raw, dtype=np.float64)
    if not np.isfinite(a).all():
        raise ValueError(f"{name} must be finite")
    return np.frombuffer(a.tobytes(), dtype=np.float64).reshape(a.shape)


@dataclass(frozen=True)
class Variable:
    """Explicit column identity, SI unit and physical definition (no conversion)."""

    name: str
    unit: str
    definition: str

    def __post_init__(self) -> None:
        for value in (self.name, self.unit, self.definition):
            if (
                not isinstance(value, str)
                or not value.strip()
                or value.strip() != value
            ):
                raise ValueError(
                    "variable name, SI unit and definition must be nonempty strings"
                )


def _variables(values: Sequence[Variable], count: int) -> tuple[Variable, ...]:
    result = tuple(values)
    if len(result) != count or any(not isinstance(v, Variable) for v in result):
        raise ValueError("variable metadata must match column count")
    if len({v.name for v in result}) != count:
        raise ValueError("duplicate variable names")
    return result


class SurrogateData:
    """Validated owned data X[n,p], Y[n,q], IDs and mandatory source provenance.

    X/Y must be finite real matrices, never complex/bool/object or empty. Data
    and metadata are detached from caller mutation; array storage is read-only.
    Every row needs a unique nonempty case ID. Units are explicit SI declarations,
    not a unit-conversion engine. Raises ValueError for incompatible data.
    """

    def __init__(
        self,
        X: ArrayLike,
        Y: ArrayLike,
        case_ids: Sequence[str],
        parameters: Sequence[Variable],
        targets: Sequence[Variable],
        *,
        provenance: Mapping[str, object],
    ) -> None:
        self._X = _array(X, ndim=2, name="X")
        self._Y = _array(Y, ndim=2, name="Y")
        self.case_ids = tuple(case_ids)
        if len(self._X) != len(self._Y) or len(self.case_ids) != len(self._X):
            raise ValueError("X, Y and case IDs must have the same sample count")
        if any(not isinstance(i, str) or not i.strip() for i in self.case_ids) or len(
            set(self.case_ids)
        ) != len(self.case_ids):
            raise ValueError("case IDs must be unique nonempty strings")
        self.parameters = _variables(parameters, self._X.shape[1])
        self.targets = _variables(targets, self._Y.shape[1])
        if set(v.name for v in self.parameters) & set(v.name for v in self.targets):
            raise ValueError("parameter and target names must be distinct")
        if not isinstance(provenance, Mapping) or not provenance:
            raise ValueError("source provenance is required")
        self._provenance = _canonical(dict(provenance))

    @classmethod
    def from_table(
        cls,
        directory: str | Path,
        *,
        parameter_columns: Sequence[str],
        target_columns: Sequence[str],
        exclude_invalid: bool = False,
    ) -> SurrogateData:
        """Load selected scalar columns from a verified #60 CSV/JSON table.

        Failed/unavailable cases, explicitly failed balance checks and missing,
        nonfinite or nonscalar selected values are rejected by default. Explicit
        exclusion records every excluded case and all reasons in provenance.
        Unchecked balances remain visibly unchecked, never promoted to passed.
        """
        from femantis.study_reports import read_study_table

        if type(exclude_invalid) is not bool:
            raise ValueError("exclude_invalid must be boolean")
        table = read_study_table(directory)
        names = [*parameter_columns, *target_columns]
        if not parameter_columns or not target_columns or len(set(names)) != len(names):
            raise ValueError("select distinct nonempty parameter/target columns")
        metadata = []
        for selected, role in (
            (parameter_columns, "parameter"),
            (target_columns, "quantity"),
        ):
            variables = []
            for name in selected:
                if (
                    name not in table["columns"]
                    or table["columns"][name]["role"] != role
                ):
                    raise ValueError(f"column {name!r} missing or incompatible role")
                column = table["columns"][name]
                variables.append(Variable(name, column["unit"], column["definition"]))
            metadata.append(variables)
        X, Y, ids, exclusions, sources = [], [], [], [], {}
        for row in table["rows"]:
            reasons = []
            if row["status"] != "succeeded" or row["solver_status"] != "converged":
                reasons.append(
                    f"case status {row['status']}, solver {row['solver_status']}"
                )
            if row["balance_status"] == "failed":
                reasons.append("explicit balance checks failed")
            for name in names:
                value = row["values"].get(name)
                if type(value) not in (int, float) or not math.isfinite(value):
                    reasons.append(
                        f"{name}: {row['missing'].get(name, 'not a finite scalar')}"
                    )
            if reasons:
                exclusions.append({"case_id": row["case_id"], "reasons": reasons})
                continue
            ids.append(row["case_id"])
            X.append([row["values"][name] for name in parameter_columns])
            Y.append([row["values"][name] for name in target_columns])
            sources[row["case_id"]] = {
                key: row[key]
                for key in (
                    "fingerprints",
                    "versions",
                    "result",
                    "origins",
                    "balance_status",
                    "mesh_convergence",
                )
            }
        if exclusions and not exclude_invalid:
            raise ValueError(
                "invalid table cases require explicit exclusion: "
                + _canonical(exclusions)
            )
        if not ids:
            raise ValueError("no valid rows after table selection")
        return cls(
            X,
            Y,
            ids,
            metadata[0],
            metadata[1],
            provenance={
                "table_schema_version": table["table_schema_version"],
                "table_sha256": _digest(table),
                "csv_sha256": table["csv_sha256"],
                "catalog_sha256": table["catalog_sha256"],
                "selected_columns": names,
                "rows": sources,
                "exclusions": exclusions,
            },
        )

    @property
    def X(self) -> NDArray[np.float64]:
        """Read-only parameters in declared SI units, shape [n,p]."""
        return self._X

    @property
    def Y(self) -> NDArray[np.float64]:
        """Read-only targets in declared SI units, shape [n,q]."""
        return self._Y

    @property
    def provenance(self) -> dict[str, Any]:
        """Detached source metadata, including any table exclusion log."""
        return dict(json.loads(self._provenance))

    @property
    def fingerprint(self) -> str:
        """SHA-256 of values, row/column identity, definitions and provenance."""
        return _digest(
            {
                "X": self.X.tolist(),
                "Y": self.Y.tolist(),
                "ids": self.case_ids,
                "parameters": [asdict(v) for v in self.parameters],
                "targets": [asdict(v) for v in self.targets],
                "provenance": self.provenance,
            }
        )


@dataclass(frozen=True)
class DataSplit:
    """Explicit train/test ID partition; seed records optional random splitting."""

    train_ids: tuple[str, ...]
    test_ids: tuple[str, ...]
    seed: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "train_ids", tuple(self.train_ids))
        object.__setattr__(self, "test_ids", tuple(self.test_ids))
        if self.seed is not None and (type(self.seed) is not int or self.seed < 0):
            raise ValueError("seed must be a nonnegative integer or None")


def split_cases(data: SurrogateData, *, test_fraction: float, seed: int) -> DataSplit:
    """Seeded split of distinct X groups; duplicate inputs always stay together.

    Fraction applies to groups, rounded up, with at least one group on each side.
    Order within each side follows the original dataset. No fitting is performed.
    """
    if (
        type(seed) is not int
        or seed < 0
        or not math.isfinite(test_fraction)
        or not 0 < test_fraction < 1
    ):
        raise ValueError("need nonnegative integer seed and 0 < test_fraction < 1")
    groups = list(dict.fromkeys(tuple(row) for row in data.X))
    if len(groups) < 2:
        raise ValueError("split needs at least two distinct parameter points")
    ntest = min(len(groups) - 1, math.ceil(len(groups) * test_fraction))
    test = {
        groups[i] for i in np.random.default_rng(seed).permutation(len(groups))[:ntest]
    }
    return DataSplit(
        tuple(
            i
            for i, x in zip(data.case_ids, data.X, strict=True)
            if tuple(x) not in test
        ),
        tuple(
            i for i, x in zip(data.case_ids, data.X, strict=True) if tuple(x) in test
        ),
        seed,
    )


def _split(data: SurrogateData, split: DataSplit) -> tuple[list[int], list[int]]:
    train, test = split.train_ids, split.test_ids
    if (
        not train
        or not test
        or len(set(train)) != len(train)
        or len(set(test)) != len(test)
    ):
        raise ValueError("train/test must both be nonempty without duplicate IDs")
    if set(train) & set(test) or set(train) | set(test) != set(data.case_ids):
        raise ValueError("train/test IDs must be a disjoint complete dataset partition")
    index = {name: i for i, name in enumerate(data.case_ids)}
    tr, te = [index[i] for i in train], [index[i] for i in test]
    if {tuple(x) for x in data.X[tr]} & {tuple(x) for x in data.X[te]}:
        raise ValueError(
            "identical parameter points cannot cross train/test boundaries"
        )
    return tr, te


def _basis(p: int, degree: int) -> list[list[int]]:
    exponents = [[0] * p]
    for i in range(p):
        row = [0] * p
        row[i] = 1
        exponents.append(row)
    if degree == 2:
        for i in range(p):
            for j in range(i, p):
                row = [0] * p
                row[i] += 1
                row[j] += 1
                exponents.append(row)
    return exponents


def _design(
    X: NDArray[np.float64],
    center: NDArray[np.float64],
    scale: NDArray[np.float64],
    basis: list[list[int]],
) -> NDArray[np.float64]:
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        z = (X - center) / scale
        a = np.column_stack([np.prod(z ** np.asarray(e), axis=1) for e in basis])
    if not np.isfinite(a).all():
        raise ValueError("nonfinite polynomial design; check input scales")
    return a


def _metrics(
    expected: NDArray[np.float64], actual: NDArray[np.float64]
) -> dict[str, list[float]]:
    with np.errstate(over="ignore", invalid="ignore"):
        error = actual - expected
        values = {
            "mae": np.mean(np.abs(error), axis=0),
            "rmse": np.sqrt(np.mean(error**2, axis=0)),
            "max_abs": np.max(np.abs(error), axis=0),
        }
    if any(not np.isfinite(v).all() for v in values.values()):
        raise ValueError("nonfinite error metric; check output scales")
    return {key: value.tolist() for key, value in values.items()}


@dataclass(frozen=True)
class Prediction:
    """Values [n,q] in target units and an outside-training-box flag per row."""

    values: NDArray[np.float64]
    extrapolated: NDArray[np.bool_]
    parameters: tuple[Variable, ...]
    targets: tuple[Variable, ...]


class PolynomialSurrogate:
    """Polynomial fitted by fit_surrogate; prediction needs only its artifact.

    Min/max checks are box checks, not guarantees of accuracy inside the box.
    Test errors are empirical errors, never uncertainty intervals/error bounds.
    """

    def __init__(self, document: Mapping[str, Any]) -> None:
        self._document: dict[str, Any] = json.loads(_canonical(dict(document)))
        _validate_artifact(self._document)

    @property
    def report(self) -> dict[str, Any]:
        """Detached training/test errors, split, source and version provenance."""
        return dict(json.loads(_canonical(self._document["report"])))

    def to_document(self) -> dict[str, Any]:
        """Return an inspectable detached JSON document."""
        return dict(json.loads(_canonical(self._document)))

    def save(self, path: str | Path, *, overwrite: bool = False) -> None:
        """Atomically save schema-1 JSON/checksum; default refuses overwriting.

        Parent must exist. No pickle, executable data or training files required.
        Overwrite uses atomic replacement; exclusive publication uses a hard link.
        """
        if type(overwrite) is not bool:
            raise ValueError("overwrite must be boolean")
        target = Path(path)
        if target.suffix.lower() != ".json":
            raise ValueError("surrogate artifact needs .json suffix")
        document = self.to_document()
        content = _canonical({"model": document, "artifact_sha256": _digest(document)})
        temporary = None
        try:
            with NamedTemporaryFile(
                mode="w", encoding="utf-8", dir=target.parent, delete=False
            ) as stream:
                temporary = Path(stream.name)
                stream.write(content + "\n")
                stream.flush()
                os.fsync(stream.fileno())
            if overwrite:
                os.replace(temporary, target)
            else:
                os.link(temporary, target)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)

    @classmethod
    def load(cls, path: str | Path) -> PolynomialSurrogate:
        """Validate checksum, version and complete model structure from JSON.

        Raises ValueError for damaged/incompatible artifacts; filesystem errors
        propagate. A checksum detects accidental damage, not hostile tampering.
        """
        try:
            wrapper = json.loads(Path(path).read_text(encoding="utf-8"))
            if not isinstance(wrapper, dict) or set(wrapper) != {
                "model",
                "artifact_sha256",
            }:
                raise ValueError("invalid artifact wrapper")
            if _digest(wrapper["model"]) != wrapper["artifact_sha256"]:
                raise ValueError("surrogate artifact checksum mismatch")
            return cls(wrapper["model"])
        except (KeyError, TypeError) as error:
            raise ValueError("invalid surrogate artifact") from error

    def predict(self, X: ArrayLike, *, allow_extrapolation: bool = False) -> Prediction:
        """Predict one [p] or several [n,p] points; always return [n,q] values.

        Outside-box input raises ValueError by default. Explicit permission
        returns per-row flags. Columns must follow the exposed parameter order.
        Complex, nonfinite, empty and dimensionally incompatible arrays fail.
        """
        raw = np.asarray(X)
        if raw.ndim == 1:
            raw = raw.reshape(1, -1)
        x = _array(raw, ndim=2, name="prediction X")
        doc = self._document
        if (
            x.shape[1] != len(doc["parameters"])
            or type(allow_extrapolation) is not bool
        ):
            raise ValueError("prediction columns or extrapolation flag invalid")
        low, high = np.asarray(doc["minimum"]), np.asarray(doc["maximum"])
        flags = np.any((x < low) | (x > high), axis=1)
        if flags.any() and not allow_extrapolation:
            raise ValueError(
                "prediction outside training parameter box; explicit permission required"
            )
        a = _design(
            x, np.asarray(doc["center"]), np.asarray(doc["scale"]), doc["basis"]
        )
        with np.errstate(over="ignore", invalid="ignore"):
            values = a @ np.asarray(doc["coefficients"])
        result = _array(values, ndim=2, name="prediction")
        return Prediction(
            result,
            np.frombuffer(flags.tobytes(), dtype=np.bool_),
            tuple(Variable(**v) for v in doc["parameters"]),
            tuple(Variable(**v) for v in doc["targets"]),
        )


def fit_surrogate(
    data: SurrogateData, split: DataSplit, *, degree: int = 1, ridge: float = 0.0
) -> PolynomialSurrogate:
    """Fit total-degree 1/2 polynomial with train-only box scaling and intercept.

    Ridge solves ||A B - Y||² + ridge*||B[1:]||² by augmented least squares;
    intercept is unpenalized. No normal equations or automatic tuning. At least
    as many distinct training points as basis terms are required, even for ridge.
    Constant parameters are rejected (remove them explicitly); constant outputs
    are allowed. Unregularized rank deficiency raises ValueError. Test data is
    used only for the returned evaluation report, never scaling or fitting.
    """
    if type(degree) is not int or degree not in (1, 2):
        raise ValueError("degree must be integer 1 or 2")
    if (
        isinstance(ridge, bool)
        or not isinstance(ridge, (int, float))
        or not math.isfinite(ridge)
        or ridge < 0
    ):
        raise ValueError("ridge must be finite and nonnegative")
    tr, te = _split(data, split)
    train = data.X[tr]
    low, high = train.min(axis=0), train.max(axis=0)
    center, scale = low / 2 + high / 2, high / 2 - low / 2
    if np.any(scale <= 0) or not np.isfinite(scale).all():
        raise ValueError(
            "constant or unusably scaled training parameters; remove explicitly"
        )
    basis = _basis(train.shape[1], degree)
    if len({tuple(x) for x in train}) < len(basis):
        raise ValueError(
            "underdetermined fit: need at least as many distinct training points as basis terms"
        )
    a = _design(train, center, scale, basis)
    b = data.Y[tr]
    rank = int(np.linalg.matrix_rank(a))
    if rank < len(basis) and ridge == 0:
        raise ValueError("rank-deficient unregularized polynomial design")
    if ridge:
        penalty = np.sqrt(ridge) * np.eye(len(basis))
        penalty[0, 0] = 0
        augmented = np.vstack([a, penalty])
        targets = np.vstack([b, np.zeros((len(basis), b.shape[1]))])
    else:
        augmented, targets = a, b
    coefficients, _, solved_rank, _ = np.linalg.lstsq(augmented, targets, rcond=None)
    if solved_rank < len(basis):
        raise ValueError(
            "numerically rank-deficient least-squares system; reduce ridge or rescale"
        )
    if not np.isfinite(coefficients).all():
        raise ValueError("nonfinite fit coefficients")
    test_design = _design(data.X[te], center, scale, basis)
    flags = np.any((data.X[te] < low) | (data.X[te] > high), axis=1)
    report = {
        "dataset_sha256": data.fingerprint,
        "case_sha256": {
            i: _digest({"X": x.tolist(), "Y": y.tolist()})
            for i, x, y in zip(data.case_ids, data.X, data.Y, strict=True)
        },
        "split": asdict(split),
        "degree": degree,
        "ridge": float(ridge),
        "design_rank": rank,
        "training": _metrics(b, a @ coefficients),
        "test": _metrics(data.Y[te], test_design @ coefficients),
        "test_extrapolated": flags.tolist(),
        "targets": [asdict(v) for v in data.targets],
        "parameters": [asdict(v) for v in data.parameters],
        "training_box": {"minimum": low.tolist(), "maximum": high.tolist()},
        "source": data.provenance,
        "versions": {
            "toolbox": __version__,
            "numpy": version("numpy"),
            "surrogate_schema": 1,
        },
        "limitations": [
            "Box membership does not guarantee accuracy",
            "Test errors are not uncertainty bounds or experimental validation",
        ],
    }
    return PolynomialSurrogate(
        {
            "surrogate_schema_version": 1,
            "method": "scaled_polynomial_lstsq",
            "degree": degree,
            "ridge": float(ridge),
            "basis": basis,
            "center": center.tolist(),
            "scale": scale.tolist(),
            "minimum": low.tolist(),
            "maximum": high.tolist(),
            "coefficients": coefficients.tolist(),
            "parameters": [asdict(v) for v in data.parameters],
            "targets": [asdict(v) for v in data.targets],
            "report": report,
        }
    )


def _validate_artifact(doc: dict[str, Any]) -> None:
    try:
        expected = {
            "surrogate_schema_version",
            "method",
            "degree",
            "ridge",
            "basis",
            "center",
            "scale",
            "minimum",
            "maximum",
            "coefficients",
            "parameters",
            "targets",
            "report",
        }
        if (
            set(doc) != expected
            or type(doc["surrogate_schema_version"]) is not int
            or doc["surrogate_schema_version"] != 1
        ):
            raise ValueError("unsupported surrogate schema or keys")
        degree = doc["degree"]
        if (
            type(degree) is not int
            or degree not in (1, 2)
            or doc["method"] != "scaled_polynomial_lstsq"
        ):
            raise ValueError("invalid polynomial method/degree")
        if (
            type(doc["ridge"]) not in (int, float)
            or not math.isfinite(doc["ridge"])
            or doc["ridge"] < 0
        ):
            raise ValueError("invalid ridge")
        parameters = _variables(
            [Variable(**v) for v in doc["parameters"]], len(doc["parameters"])
        )
        targets = _variables(
            [Variable(**v) for v in doc["targets"]], len(doc["targets"])
        )
        if (
            not parameters
            or not targets
            or {v.name for v in parameters} & {v.name for v in targets}
        ):
            raise ValueError("invalid parameter/target identities")
        p, q = len(parameters), len(targets)
        basis = _basis(p, degree)
        if doc["basis"] != basis or any(
            type(e) is not int for row in doc["basis"] for e in row
        ):
            raise ValueError("invalid basis ordering")
        coefficients = _array(doc["coefficients"], ndim=2, name="coefficients")
        if coefficients.shape != (len(basis), q):
            raise ValueError("coefficient shape disagrees with basis/targets")
        vectors = {
            key: _array(doc[key], ndim=1, name=key)
            for key in ("center", "scale", "minimum", "maximum")
        }
        if any(len(v) != p for v in vectors.values()):
            raise ValueError("invalid parameter shape")
        low, high = vectors["minimum"], vectors["maximum"]
        if (
            np.any(low >= high)
            or not np.array_equal(vectors["center"], low / 2 + high / 2)
            or not np.array_equal(vectors["scale"], high / 2 - low / 2)
            or np.any(vectors["scale"] <= 0)
        ):
            raise ValueError("inconsistent training box/scaling")
        report = doc["report"]
        split = DataSplit(**report["split"])
        ids = (*split.train_ids, *split.test_ids)
        if (
            not split.train_ids
            or not split.test_ids
            or any(not isinstance(i, str) or not i for i in ids)
            or len(set(ids)) != len(ids)
        ):
            raise ValueError("invalid stored split")
        if set(report["case_sha256"]) != set(ids):
            raise ValueError("case provenance disagrees with split")
        for digest in [report["dataset_sha256"], *report["case_sha256"].values()]:
            if (
                not isinstance(digest, str)
                or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)
            ):
                raise ValueError("invalid provenance fingerprint")
        if (
            report["parameters"] != doc["parameters"]
            or report["targets"] != doc["targets"]
            or report["degree"] != degree
            or report["ridge"] != doc["ridge"]
            or report["training_box"]
            != {"minimum": doc["minimum"], "maximum": doc["maximum"]}
        ):
            raise ValueError("report/model metadata disagree")
        if type(report["design_rank"]) is not int or not 1 <= report[
            "design_rank"
        ] <= len(basis):
            raise ValueError("invalid design rank")
        if len(split.train_ids) < len(basis) or (
            doc["ridge"] == 0 and report["design_rank"] != len(basis)
        ):
            raise ValueError("invalid training rank/count")
        for part in ("training", "test"):
            if set(report[part]) != {"mae", "rmse", "max_abs"}:
                raise ValueError("missing error metrics")
            for values in report[part].values():
                a = _array(values, ndim=1, name="error metric")
                if len(a) != q or np.any(a < 0):
                    raise ValueError("invalid error metric")
        flags = report["test_extrapolated"]
        if len(flags) != len(split.test_ids) or any(type(v) is not bool for v in flags):
            raise ValueError("invalid test extrapolation flags")
        versions = report["versions"]
        if not isinstance(versions, dict) or set(versions) != {
            "toolbox",
            "numpy",
            "surrogate_schema",
        }:
            raise ValueError("missing or unknown version provenance")
        if (
            any(
                not isinstance(versions[k], str) or not versions[k]
                for k in ("toolbox", "numpy")
            )
            or type(versions["surrogate_schema"]) is not int
        ):
            raise ValueError("invalid version provenance")
        if (
            not isinstance(report["source"], dict)
            or not report["source"]
            or report["versions"]["surrogate_schema"] != 1
        ):
            raise ValueError("invalid source/version provenance")
    except (KeyError, TypeError, AttributeError, OverflowError) as error:
        raise ValueError(f"invalid surrogate artifact structure: {error}") from error
