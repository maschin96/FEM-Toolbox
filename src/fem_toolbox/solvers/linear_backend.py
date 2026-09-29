"""Internal linear backend boundary; physical reduction stays with adapters.

Operators are free, electrode-condensed systems in their original SI units.
Fields partition their local rows; scale is positive real, x=S*y and SAS*y=S*b.
Multiple RHS columns and real/complex data are supported without scalar loss.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, Literal, Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.linalg import cho_factor, cho_solve  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.errors import ValidationError

Numeric = NDArray[Any]
Kind = Literal["general", "positive_definite"]


def _snapshot(a: Numeric) -> Numeric:
    return np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape)


@dataclass(frozen=True)
class SolveReport:
    """Unscaled residuals; each norm ratio uses || |A||x|+|b| || on that field.

    Matrix RHS uses the Frobenius norm. Residual has RHS shape and equation
    units (N, W or signed C as supplied). Direct methods have no iterations.
    Failure raises ValidationError rather than returning a false success.
    """

    backend: str
    method: str
    converged: bool
    reason: str
    iterations: int | None
    relative_residual: float
    field_relative_residuals: Mapping[str, float]
    residual: Numeric


@dataclass(frozen=True)
class LinearSolution:
    """Owned immutable solution and residual report in original unknown units."""

    values: Numeric
    report: SolveReport


class FactorizedOperator(Protocol):
    def solve(self, rhs: ArrayLike) -> LinearSolution:
        """Solve finite (n,) or (n,nrhs) RHS; reject invalid or failed solves."""
        ...


class LinearBackend(Protocol):
    @property
    def name(self) -> str: ...

    def factorize(
        self,
        operator: Any,
        *,
        fields: Mapping[str, ArrayLike],
        scale: ArrayLike | None = None,
        kind: Kind = "general",
        pivot_tolerance: float = 0.0,
        relative_pivot: bool = False,
    ) -> FactorizedOperator:
        """Snapshot a square finite operator and its complete field partition."""
        ...


class ScipyBackend:
    """Sparse LU or dense positive-definite Cholesky; no optional dependencies."""

    name = "scipy"

    def factorize(
        self,
        operator: Any,
        *,
        fields: Mapping[str, ArrayLike],
        scale: ArrayLike | None = None,
        kind: Kind = "general",
        pivot_tolerance: float = 0.0,
        relative_pivot: bool = False,
    ) -> FactorizedOperator:
        return _ScipyFactor(
            operator, fields, scale, kind, pivot_tolerance, relative_pivot
        )


class _LinearData:
    """Owned validated data shared by backend implementations; no factorization."""

    def __init__(
        self,
        operator: Any,
        fields: Mapping[str, ArrayLike],
        scale: ArrayLike | None,
        kind: Kind,
        pivot_tolerance: float,
        relative_pivot: bool,
    ) -> None:
        try:
            a = sparse.csc_matrix(operator, copy=True)
        except (TypeError, ValueError) as error:
            raise ValidationError(
                "invalid backend operator", field="operator"
            ) from error
        if (
            a.shape[0] != a.shape[1]
            or a.dtype.kind not in "iufc"
            or not np.isfinite(a.data).all()
        ):
            raise ValidationError(
                "expected square finite numeric operator", field="operator"
            )
        self.a = a.astype(np.complex128 if a.dtype.kind == "c" else np.float64)
        n = a.shape[0]
        self.fields = {}
        for name, ids in fields.items():
            raw = np.asarray(ids)
            if (
                not isinstance(name, str)
                or not name.strip()
                or raw.ndim != 1
                or (raw.size and raw.dtype.kind not in "iu")
            ):
                raise ValidationError("invalid backend field indices", field="fields")
            if np.any(raw < 0) or np.any(raw >= n):
                raise ValidationError(
                    "backend field index outside operator", field="fields"
                )
            self.fields[name] = np.array(raw, dtype=np.int64)
        combined = (
            np.concatenate(list(self.fields.values()))
            if self.fields
            else np.empty(0, dtype=np.int64)
        )
        if not np.array_equal(np.sort(combined), np.arange(n)):
            raise ValidationError(
                "backend fields must partition operator rows exactly", field="fields"
            )
        s = np.ones(n) if scale is None else np.asarray(scale)
        if (
            s.shape != (n,)
            or s.dtype.kind not in "iuf"
            or not np.isfinite(s).all()
            or np.any(s <= 0)
        ):
            raise ValidationError(
                "expected positive finite real scaling vector", field="scale"
            )
        self.scale = np.array(s, dtype=np.float64)
        if not np.isfinite(pivot_tolerance) or pivot_tolerance < 0:
            raise ValidationError("invalid pivot tolerance", field="solver")
        self.kind = kind
        if kind not in ("general", "positive_definite"):
            raise ValidationError("unsupported factorization kind", field="solver")
        scaling = sparse.diags(self.scale)
        self.scaled = (scaling @ self.a @ scaling).tocsc()
        if not np.isfinite(self.scaled.data).all():
            raise ValidationError("nonfinite scaled operator", field="scale")

    def prepare_rhs(self, rhs: ArrayLike) -> tuple[Numeric, Numeric, Numeric]:
        b = np.asarray(rhs)
        if (
            b.ndim not in (1, 2)
            or b.shape[0] != self.a.shape[0]
            or b.dtype.kind not in "iufc"
            or not np.isfinite(b).all()
        ):
            raise ValidationError(
                "expected finite backend RHS with operator row count", field="rhs"
            )
        b = np.array(b, dtype=np.result_type(self.a.dtype, b.dtype, np.float64))
        scale = self.scale if b.ndim == 1 else self.scale[:, None]
        scaled_rhs = scale * b
        if not np.isfinite(scaled_rhs).all():
            raise ValidationError("nonfinite scaled RHS", field="scale")
        return b, scale, scaled_rhs

    def solution(
        self,
        b: Numeric,
        x: Numeric,
        *,
        backend: str,
        method: str,
        reason: str = "direct_solve",
        iterations: int | None = None,
    ) -> LinearSolution:
        residual = np.asarray(self.a @ x - b)
        denominator = np.asarray(abs(self.a) @ np.abs(x) + np.abs(b))
        if (
            not np.isfinite(x).all()
            or not np.isfinite(residual).all()
            or not np.isfinite(denominator).all()
        ):
            raise ValidationError(
                "nonfinite backend solution or residual", field="solver"
            )

        def relative(ids: NDArray[np.int64]) -> float:
            # Scale both norms to avoid overflow/underflow of squared SI magnitudes.
            d = denominator[ids]
            r = residual[ids]
            maximum = float(np.max(d)) if d.size else 0.0
            if maximum == 0:
                return 0.0 if not np.any(r) else float("inf")
            return float(np.linalg.norm(r / maximum) / np.linalg.norm(d / maximum))

        ratios = {name: relative(ids) for name, ids in self.fields.items()}
        total = relative(np.arange(len(b), dtype=np.int64))
        if total > 1e-9 or any(value > 1e-9 for value in ratios.values()):
            raise ValidationError(
                "unscaled backend field residual exceeds 1e-9", field="solver"
            )
        report = SolveReport(
            backend,
            method,
            True,
            reason,
            iterations,
            total,
            MappingProxyType(ratios),
            _snapshot(residual),
        )
        return LinearSolution(_snapshot(x), report)


class _ScipyFactor(_LinearData):
    def __init__(
        self,
        operator: Any,
        fields: Mapping[str, ArrayLike],
        scale: ArrayLike | None,
        kind: Kind,
        pivot_tolerance: float,
        relative_pivot: bool,
    ) -> None:
        super().__init__(operator, fields, scale, kind, pivot_tolerance, relative_pivot)
        self.factor: Any = None
        n = self.a.shape[0]
        scaled = self.scaled
        if not n:
            return
        try:
            if kind == "positive_definite":
                dense = scaled.toarray()
                if not np.allclose(dense, dense.conj().T, rtol=1e-12, atol=0):
                    raise ValidationError(
                        "positive definite operator must be Hermitian", field="operator"
                    )
                self.factor = cho_factor(dense)
            else:
                self.factor = splu(scaled)
                pivots = np.abs(self.factor.U.diagonal())
                bound = pivot_tolerance * (
                    float(pivots.max()) if relative_pivot else 1.0
                )
                if not np.isfinite(pivots).all() or np.min(pivots) <= bound:
                    raise ValidationError(
                        "singular or ill-conditioned backend operator", field="solver"
                    )
        except (RuntimeError, np.linalg.LinAlgError) as error:
            message = (
                "non-positive definite operator"
                if kind == "positive_definite"
                else "singular backend operator"
            )
            raise ValidationError(message, field="solver") from error

    def solve(self, rhs: ArrayLike) -> LinearSolution:
        b, scale, scaled_rhs = self.prepare_rhs(rhs)
        if not len(b):
            x = b.copy()
        else:
            try:
                if self.kind == "positive_definite":
                    y = cho_solve(self.factor, scaled_rhs)
                elif self.a.dtype.kind != "c" and b.dtype.kind == "c":
                    y = self.factor.solve(scaled_rhs.real) + 1j * self.factor.solve(
                        scaled_rhs.imag
                    )
                else:
                    y = self.factor.solve(scaled_rhs)
                x = np.asarray(scale * y)
            except (RuntimeError, np.linalg.LinAlgError) as error:
                raise ValidationError("backend solve failed", field="solver") from error
        return self.solution(
            b,
            x,
            backend="scipy",
            method="cholesky" if self.kind == "positive_definite" else "sparse_lu",
        )


def solve_linear(
    operator: Any,
    rhs: ArrayLike,
    *,
    fields: Mapping[str, ArrayLike],
    scale: ArrayLike | None = None,
    kind: Kind = "general",
    pivot_tolerance: float = 0.0,
    relative_pivot: bool = False,
    backend: LinearBackend | None = None,
) -> LinearSolution:
    """Factor and solve once; use backend.factorize explicitly for reuse."""
    selected = ScipyBackend() if backend is None else backend
    return selected.factorize(
        operator,
        fields=fields,
        scale=scale,
        kind=kind,
        pivot_tolerance=pivot_tolerance,
        relative_pivot=relative_pivot,
    ).solve(rhs)
