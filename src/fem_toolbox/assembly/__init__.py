"""Physics-independent COO assembly and preserved original linear systems.

Local vectors already include integration and the weak-form sign. Matrix and
vector entries use the SI units of their respective equations.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]

from fem_toolbox.errors import ValidationError


def indices(value: ArrayLike, size: int | None = None) -> NDArray[np.int64]:
    """Validate a one-dimensional sequence of nonnegative integer indices."""
    a = np.asarray(value)
    if a.ndim != 1 or (a.size and a.dtype.kind not in "iu"):
        raise ValidationError("expected integer index vector", field="indices")
    if a.size and (np.any(a < 0) or np.any(a > np.iinfo(np.int64).max)):
        raise ValidationError("index outside int64 range", field="indices")
    if size is not None and np.any(a >= size):
        raise ValidationError("index outside system", field="indices")
    return np.array(a, dtype=np.int64)


def real_array(value: ArrayLike, shape: tuple[int, ...]) -> NDArray[np.float64]:
    """Copy finite real data with the required shape; reject complex values."""
    a = np.asarray(value)
    if a.shape != shape or a.dtype.kind not in "iuf" or not np.isfinite(a).all():
        raise ValidationError(
            f"expected finite real array of shape {shape}", field="values"
        )
    return np.array(a, dtype=np.float64)


@dataclass(frozen=True)
class MatrixBlock:
    """Local matrix (nrow, ncol) with separate global row/column DOFs."""

    rows: ArrayLike
    columns: ArrayLike
    values: ArrayLike


@dataclass(frozen=True)
class VectorBlock:
    """Integrated local load (ndof,); repeated DOFs accumulate."""

    dofs: ArrayLike
    values: ArrayLike


class LinearSystem:
    """Owned sparse matrix and RHS; accessors return independent copies."""

    def __init__(self, matrix: sparse.spmatrix, rhs: ArrayLike) -> None:
        if (
            not sparse.issparse(matrix)
            or matrix.ndim != 2
            or matrix.shape[0] != matrix.shape[1]
        ):
            raise ValidationError("expected square sparse matrix", field="matrix")
        if matrix.dtype.kind not in "iuf" or not np.isfinite(matrix.tocoo().data).all():
            raise ValidationError("expected finite real matrix", field="matrix")
        self._matrix = matrix.astype(np.float64).copy()
        self._rhs = real_array(rhs, (matrix.shape[0],))

    @property
    def matrix(self) -> sparse.spmatrix:
        """Independent sparse matrix copy."""
        return self._matrix.copy()

    @property
    def rhs(self) -> NDArray[np.float64]:
        """Independent load vector copy, shape (size,)."""
        return self._rhs.copy()

    @property
    def size(self) -> int:
        """Number of global DOFs."""
        return int(self._rhs.size)

    def residual(self, solution: ArrayLike) -> NDArray[np.float64]:
        """Return K x - f in equation units, including constrained reactions."""
        x = real_array(solution, (self.size,))
        return np.asarray(self._matrix @ x - self._rhs, dtype=np.float64)


def assemble(
    size: int,
    matrices: Iterable[MatrixBlock],
    vectors: Iterable[VectorBlock] = (),
    *,
    format: Literal["csr", "csc"] = "csr",
) -> LinearSystem:
    """Sum local blocks through COO into CSR/CSC; no inactive DOFs are added.

    Invalid shapes, indices or nonfinite values raise ValidationError. Empty
    systems and repeated indices are supported. No symmetry is assumed.
    """
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise ValidationError("size must be a nonnegative integer", field="size")
    if format not in ("csr", "csc"):
        raise ValidationError("format must be csr or csc", field="format")
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for block in matrices:
        r, c = indices(block.rows, size), indices(block.columns, size)
        a = real_array(block.values, (len(r), len(c)))
        rows.extend(np.repeat(r, len(c)).tolist())
        columns.extend(np.tile(c, len(r)).tolist())
        values.extend(a.ravel().tolist())
    matrix = sparse.coo_matrix((values, (rows, columns)), shape=(size, size)).asformat(
        format
    )
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    rhs = np.zeros(size)
    for vector in vectors:
        dofs = indices(vector.dofs, size)
        np.add.at(rhs, dofs, real_array(vector.values, (len(dofs),)))
    return LinearSystem(matrix, rhs)
