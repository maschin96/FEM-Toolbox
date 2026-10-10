"""Explicit Monte Carlo input uncertainty, separate from numerical model errors.

Laws describe user assumptions, never inferred measurement distributions. Rows
are joint samples, columns follow parameter declaration order. No optimizer,
inverse identification or experimental validation is performed.
"""

from __future__ import annotations

import json
import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import asdict, dataclass
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import ndtr  # type: ignore[import-untyped]

from fem_toolbox.studies import Study
from fem_toolbox.study_reports import Quantity, export_study_table, read_study_table

__all__ = [
    "Normal",
    "SamplingPlan",
    "UncertainParameter",
    "UncertaintyStudy",
    "Uniform",
    "local_sensitivity",
    "summarize_samples",
]


def _finite(value: object, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, float, np.integer, np.floating)
    ):
        raise ValueError(f"{name} must be a finite real scalar")
    number = float(value)
    if not np.isfinite(number):
        raise ValueError(f"{name} must be finite")
    return number


def _vector(values: ArrayLike, name: str) -> NDArray[np.float64]:
    raw = np.asarray(values)
    if raw.ndim != 1 or not raw.size or raw.dtype.kind not in "iuf":
        raise ValueError(f"{name} must be a nonempty real vector")
    result = np.asarray(raw, dtype=np.float64)
    if not np.isfinite(result).all():
        raise ValueError(f"{name} must be finite")
    return result.copy()


@dataclass(frozen=True)
class Normal:
    """Unbounded normal law with finite mean and positive standard deviation.

    Values are in the parameter's declared SI unit. Negative/otherwise
    unphysical samples are not clipped; a model may reject them explicitly.
    """

    mean: float
    standard_deviation: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "mean", _finite(self.mean, "mean"))
        sd = _finite(self.standard_deviation, "standard_deviation")
        if sd <= 0:
            raise ValueError("standard_deviation must be positive")
        object.__setattr__(self, "standard_deviation", sd)


@dataclass(frozen=True)
class Uniform:
    """Uniform probability law on a finite increasing SI interval.

    An interval alone does not justify a uniform law; the caller must state
    why this probability assumption is appropriate for the intended study.
    """

    lower: float
    upper: float

    def __post_init__(self) -> None:
        lower, upper = _finite(self.lower, "lower"), _finite(self.upper, "upper")
        if lower >= upper or not np.isfinite(upper - lower):
            raise ValueError("uniform bounds must be increasing with finite width")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True)
class UncertainParameter:
    """Scalar M8 target, explicit SI unit, law and mandatory justification."""

    name: str
    target: str
    unit: str
    distribution: Normal | Uniform
    justification: str

    def __post_init__(self) -> None:
        import re

        if not isinstance(self.name, str) or not re.fullmatch(
            r"[A-Za-z][A-Za-z0-9_]*", self.name
        ):
            raise ValueError("parameter name must be an ASCII identifier")
        for value in (self.target, self.unit, self.justification):
            if (
                not isinstance(value, str)
                or not value.strip()
                or value != value.strip()
            ):
                raise ValueError("target, SI unit and justification are required")
        if not isinstance(self.distribution, (Normal, Uniform)):
            raise ValueError("unsupported distribution")


class SamplingPlan:
    """Detached joint Monte Carlo specification with PCG64 and explicit seed.

    ``size`` is an integer >=2, ``seed`` a nonnegative integer. Default inputs
    are independent. An optional real correlation[p,p] must be symmetric, PSD
    and have unit diagonal (absolute roundoff tolerance 1e-12). Nonidentity
    correlation is supported only for all-normal inputs. Singular normal
    correlation is valid; eigenvalues in [-1e-12,0] are rounded to zero.
    ``sample`` returns immutable float64 rows [size,p], detached from callers.
    No distribution fitting, truncation or automatic unit conversion.
    """

    def __init__(
        self,
        parameters: Sequence[UncertainParameter],
        *,
        size: int,
        seed: int,
        correlation: ArrayLike | None = None,
    ) -> None:
        parameters = tuple(parameters)
        if not parameters or any(
            not isinstance(p, UncertainParameter) for p in parameters
        ):
            raise ValueError("parameters must be nonempty UncertainParameter records")
        if len({p.name for p in parameters}) != len(parameters) or len(
            {p.target for p in parameters}
        ) != len(parameters):
            raise ValueError("duplicate parameter names or targets")
        if type(size) is not int or size < 2 or type(seed) is not int or seed < 0:
            raise ValueError(
                "size must be an integer >=2 and seed a nonnegative integer"
            )
        c = np.eye(len(parameters)) if correlation is None else np.asarray(correlation)
        if (
            c.shape != (len(parameters), len(parameters))
            or c.dtype.kind not in "iuf"
            or not np.isfinite(c).all()
        ):
            raise ValueError("correlation must be a finite real square matrix")
        c = np.asarray(c, dtype=float)
        if not np.allclose(c, c.T, atol=1e-12, rtol=0) or not np.allclose(
            np.diag(c), 1, atol=1e-12, rtol=0
        ):
            raise ValueError("correlation must be symmetric with unit diagonal")
        c = (c + c.T) / 2
        if np.linalg.eigvalsh(c).min() < -1e-12:
            raise ValueError("correlation must be positive semidefinite")
        if not np.array_equal(c, np.eye(len(parameters))) and not all(
            isinstance(p.distribution, Normal) for p in parameters
        ):
            raise ValueError("dependent sampling requires all-normal distributions")
        self._parameters = parameters
        self._size, self._seed = size, seed
        self._correlation = tuple(tuple(float(v) for v in row) for row in c)

    @property
    def parameters(self) -> tuple[UncertainParameter, ...]:
        return self._parameters

    @property
    def size(self) -> int:
        return self._size

    @property
    def seed(self) -> int:
        return self._seed

    def to_document(self) -> dict[str, Any]:
        """Return detached schema-1 provenance, including laws and RNG method."""
        return {
            "sampling_schema_version": 1,
            "method": "Monte_Carlo",
            "generator": "PCG64/standard_normal/normal_CDF_uniform/eigh_PSD",
            "size": self.size,
            "seed": self.seed,
            "numpy_version": np.__version__,
            "scipy_version": version("scipy"),
            "correlation": [list(row) for row in self._correlation],
            "parameters": [
                {
                    **asdict(p),
                    "distribution": {
                        "kind": type(p.distribution).__name__.lower(),
                        **asdict(p.distribution),
                    },
                }
                for p in self.parameters
            ],
        }

    def sample(self) -> NDArray[np.float64]:
        """Generate joint samples; same specification/runtime gives same bytes."""
        z = np.random.Generator(np.random.PCG64(self.seed)).standard_normal(
            (self.size, len(self.parameters))
        )
        c = np.asarray(self._correlation)
        if not np.array_equal(c, np.eye(len(self.parameters))):
            eigenvalues, vectors = np.linalg.eigh(c)
            z = z @ (vectors * np.sqrt(np.maximum(eigenvalues, 0))).T
        values = np.empty_like(z)
        for j, p in enumerate(self.parameters):
            law = p.distribution
            if isinstance(law, Normal):
                values[:, j] = law.mean + law.standard_deviation * z[:, j]
            else:
                values[:, j] = law.lower + (law.upper - law.lower) * ndtr(z[:, j])
        if not np.isfinite(values).all():
            raise ValueError("sample overflow; revise input scales")
        return np.frombuffer(values.tobytes(), dtype=np.float64).reshape(values.shape)


def summarize_samples(
    values: ArrayLike,
    *,
    unit: str,
    numerical_error_bound: float | None,
    quantile_probabilities: Sequence[float] = (0.025, 0.5, 0.975),
) -> dict[str, Any]:
    """Summarize finite real output rows; never omit failed/nonfinite samples.

    Variance is unbiased (ddof=1); quantiles use NumPy linear interpolation.
    Quantile bands are output probability estimates, not confidence intervals
    or guaranteed bounds. Prefix moments document sample convergence, not a
    proof that sampling error is negligible. ``numerical_error_bound`` is an
    externally justified deterministic bound in ``unit``, or None if unknown;
    it is reported separately and is never added to statistical variance.
    """
    y = _vector(values, "output samples")
    if len(y) < 2 or not isinstance(unit, str) or not unit.strip():
        raise ValueError("at least two samples and an explicit SI unit are required")
    q = _vector(quantile_probabilities, "quantile probabilities")
    if np.any(q < 0) or np.any(q > 1) or np.any(np.diff(q) <= 0):
        raise ValueError("quantile probabilities must increase in [0,1]")
    bound = (
        None
        if numerical_error_bound is None
        else _finite(numerical_error_bound, "numerical_error_bound")
    )
    if bound is not None and bound < 0:
        raise ValueError("numerical_error_bound must be nonnegative")
    sizes = sorted(
        {2, len(y), *(2**i for i in range(1, len(y).bit_length()) if 2**i <= len(y))}
    )
    result = {
        "size": len(y),
        "unit": unit,
        "mean": float(np.mean(y)),
        "variance": float(np.var(y, ddof=1)),
        "quantile_probabilities": q.tolist(),
        "quantiles": np.quantile(y, q, method="linear").tolist(),
        "numerical_error_bound": bound,
        "convergence": [
            {
                "size": n,
                "mean": float(np.mean(y[:n])),
                "variance": float(np.var(y[:n], ddof=1)),
            }
            for n in sizes
        ],
    }
    if not np.isfinite([result["mean"], result["variance"]]).all():
        raise ValueError("summary overflow; revise output scales")
    return result


def local_sensitivity(
    evaluate: Callable[[NDArray[np.float64]], float],
    *,
    center: ArrayLike,
    steps: ArrayLike,
) -> dict[str, Any]:
    """Central finite-difference gradient at explicit SI center and positive steps.

    Each derivative has output-unit/input-unit; it is local, not a global
    variance decomposition. The caller must establish step convergence and
    admissibility of center +/- step. Callback exceptions propagate. Exactly
    2*p scalar evaluations; callback inputs are detached copies.
    """
    x, h = _vector(center, "center"), _vector(steps, "steps")
    if x.shape != h.shape or np.any(h <= 0):
        raise ValueError("positive steps must match center")
    gradient = []
    for j in range(len(x)):
        plus, minus = x.copy(), x.copy()
        plus[j] += h[j]
        minus[j] -= h[j]
        if (
            not np.isfinite([plus[j], minus[j]]).all()
            or plus[j] == x[j]
            or minus[j] == x[j]
        ):
            raise ValueError("steps must produce distinct finite points")
        a, b = _finite(evaluate(plus), "response"), _finite(evaluate(minus), "response")
        gradient.append(_finite((a - b) / (2 * h[j]), "gradient"))
    return {
        "method": "central_difference",
        "center": x.tolist(),
        "steps": h.tolist(),
        "gradient": gradient,
    }


class UncertaintyStudy:
    """Joint scalar sampling on the public stationary M8 Study/Quantity APIs.

    The relative ``model`` TOML and its mesh must remain inside ``root``. Each
    joint row runs as a one-case Study, so correlations are retained without a
    Cartesian product. M8 owns validation, retry, provenance and hashed result
    reuse. Reusing an output directory requires identical sampling provenance;
    changed model/mesh bytes invalidate cached M8 results. Normal-law physical
    failures stay visible and suppress all uncertainty statistics (no silent
    clipping or conditional subset). Results are not an experimental claim.
    """

    def __init__(self, plan: SamplingPlan, *, model: str, root: str | Path) -> None:
        if not isinstance(plan, SamplingPlan):
            raise ValueError("plan must be a SamplingPlan")
        self.plan = plan
        self.model = model
        self.root = Path(root).resolve()

    def _study(self, values: NDArray[np.float64]) -> Study:
        return Study(
            {
                "study_schema_version": 1,
                "model": self.model,
                "parameters": [
                    {
                        "name": p.name,
                        "target": p.target,
                        "unit": p.unit,
                        "values": [float(v)],
                    }
                    for p, v in zip(self.plan.parameters, values, strict=True)
                ],
            },
            root=self.root,
        )

    def run(
        self,
        output: str | Path,
        *,
        quantities: Sequence[Quantity],
        numerical_error_bounds: Mapping[str, float | None],
        sensitivity_center: ArrayLike | None = None,
        sensitivity_steps: ArrayLike | None = None,
    ) -> dict[str, Any]:
        """Write sampling.json and atomic uncertainty.json; return detached data.

        Every quantity requires an explicit deterministic numerical-error bound
        in its SI unit or None (unknown), never a default zero. Optional local
        sensitivities require center[p] and positive steps[p] together; each
        quantity derivative uses its output unit divided by the input unit.
        Caller establishes admissibility and step convergence. Sensitivities
        are suppressed if a perturbation fails; sampling statistics may still
        be complete. Tables/case receipts and fingerprints remain inspectable.

        An exclusive lock rejects concurrent writers. Interrupted case runs
        retain their M8 recovery state; rerunning resumes successes. After a
        hard kill, verify no writer remains before manually removing the lock.
        Schema errors propagate, physical failures appear in the report.
        """
        quantities = tuple(quantities)
        if (
            not quantities
            or any(not isinstance(q, Quantity) for q in quantities)
            or len({q.name for q in quantities}) != len(quantities)
        ):
            raise ValueError("quantities must be nonempty with unique names")
        if set(numerical_error_bounds) != {q.name for q in quantities}:
            raise ValueError("numerical_error_bounds must cover exactly all quantities")
        bounds: dict[str, float | None] = {}
        for name, value in numerical_error_bounds.items():
            bounds[name] = (
                None if value is None else _finite(value, "numerical_error_bounds")
            )
            bound = bounds[name]
            if bound is not None and bound < 0:
                raise ValueError("numerical_error_bounds must be nonnegative")
        sensitivity_values: list[NDArray[np.float64]] = []
        center = steps = None
        if sensitivity_center is not None or sensitivity_steps is not None:
            if sensitivity_center is None or sensitivity_steps is None:
                raise ValueError("sensitivity center and steps are required together")
            center = _vector(sensitivity_center, "sensitivity center")
            steps = _vector(sensitivity_steps, "sensitivity steps")
            if (
                center.shape != (len(self.plan.parameters),)
                or steps.shape != center.shape
                or np.any(steps <= 0)
            ):
                raise ValueError(
                    "sensitivity steps must be positive and match parameter count"
                )
            for j in range(len(center)):
                for sign in (1, -1):
                    row = center.copy()
                    row[j] += sign * steps[j]
                    if not np.isfinite(row).all() or row[j] == center[j]:
                        raise ValueError(
                            "sensitivity steps must produce distinct finite points"
                        )
                    sensitivity_values.append(row)
        samples = self.plan.sample()
        # Resolve all schema/unit contracts before creating output or solving.
        studies = [self._study(row) for row in [*samples, *sensitivity_values]]
        for study in studies:
            study.resolve()
        path = Path(output).resolve()
        path.mkdir(parents=True, exist_ok=True)
        lock = path / ".uncertainty.lock"
        with lock.open("x") as stream:
            stream.write(str(os.getpid()))
        try:
            document = self.plan.to_document()
            specification = path / "sampling.json"
            if specification.exists():
                if json.loads(specification.read_text(encoding="utf-8")) != document:
                    raise ValueError(
                        "output belongs to a different sampling specification"
                    )
            else:
                _write_json(specification, document)
            # Never leave a prior successful summary after an interrupted rebuild.
            (path / "uncertainty.json").unlink(missing_ok=True)
            records = []
            for index, study in enumerate(studies):
                group = "sample" if index < len(samples) else "sensitivity"
                ordinal = index if group == "sample" else index - len(samples)
                directory = path / f"{group}-{ordinal:06d}"
                catalog = study.run(directory)
                catalog_path = directory / "catalog.json"
                digest = sha256(catalog_path.read_bytes()).hexdigest()
                table_key = sha256(
                    _canonical(
                        {
                            "catalog": digest,
                            "quantities": [asdict(q) for q in quantities],
                        }
                    ).encode()
                ).hexdigest()
                table = directory / f"table-{table_key}"
                if not table.exists():
                    export_study_table(catalog_path, table, quantities)
                data = read_study_table(table)
                if data["catalog_sha256"] != digest:
                    raise ValueError("table receipt does not match current catalog")
                for q in quantities:
                    if _canonical(data["columns"].get(q.name)) != _canonical(
                        {"role": "quantity", **asdict(q)}
                    ):
                        raise ValueError("table quantity contract mismatch")
                if len(data["rows"]) != 1 or len(catalog["cases"]) != 1:
                    raise ValueError("joint row must have exactly one M8 case")
                row = data["rows"][0]
                records.append(
                    {
                        "index": ordinal,
                        "case_id": row["case_id"],
                        "catalog": str(catalog_path.relative_to(path)),
                        "catalog_sha256": digest,
                        "table": str(table.relative_to(path)),
                        "table_metadata_sha256": sha256(
                            (table / "metadata.json").read_bytes()
                        ).hexdigest(),
                        "values": {q.name: row["values"][q.name] for q in quantities},
                        "row": row,
                    }
                )
            sampled = records[: len(samples)]

            def available(record: dict[str, Any]) -> bool:
                return record["row"]["status"] == "succeeded" and all(
                    record["values"][q.name] is not None for q in quantities
                )

            complete = all(available(record) for record in sampled)
            summaries = (
                {
                    q.name: summarize_samples(
                        [r["values"][q.name] for r in sampled],
                        unit=q.unit,
                        numerical_error_bound=bounds[q.name],
                    )
                    for q in quantities
                }
                if complete
                else {}
            )
            sensitivity = {}
            perturbed = records[len(samples) :]
            if perturbed and all(available(record) for record in perturbed):
                assert center is not None and steps is not None
                for q in quantities:
                    gradient = [
                        (
                            perturbed[2 * j]["values"][q.name]
                            - perturbed[2 * j + 1]["values"][q.name]
                        )
                        / (2 * steps[j])
                        for j in range(len(center))
                    ]
                    if not np.isfinite(gradient).all():
                        raise ValueError("sensitivity overflow")
                    sensitivity[q.name] = {
                        "method": "central_difference",
                        "center": center.tolist(),
                        "steps": steps.tolist(),
                        "gradient": gradient,
                        "gradient_units": [
                            f"({q.unit})/({p.unit})" for p in self.plan.parameters
                        ],
                    }
            report = {
                "uncertainty_schema_version": 1,
                "complete": complete,
                "sampling": document,
                "sample_sha256": sha256(samples.astype("<f8").tobytes()).hexdigest(),
                "quantities": [asdict(q) for q in quantities],
                "numerical_error_bounds": bounds,
                "samples": sampled,
                "summaries": summaries,
                "sensitivity": sensitivity,
                "sensitivity_samples": perturbed,
                "sensitivity_status": "not_requested"
                if not perturbed
                else "complete"
                if sensitivity
                else "incomplete",
                "scope": "input uncertainty under declared laws; no experimental validation",
            }
            _write_json(path / "uncertainty.json", report)
            return dict(json.loads(_canonical(report)))  # detached JSON scalar types
        finally:
            lock.unlink()


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _write_json(path: Path, value: object) -> None:
    temporary: Path | None = None
    try:
        with NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, delete=False
        ) as stream:
            temporary = Path(stream.name)
            stream.write(_canonical(value) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
