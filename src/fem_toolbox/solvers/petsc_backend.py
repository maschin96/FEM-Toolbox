"""Optional one-process PETSc direct and FieldSplit backends, imported only on selection.

A real PETSc build handles complex operators by exact real block embedding and
complex RHS for real operators by two solves. Field residuals use the original
complex system. MPI-distributed solves are outside the supported scope.
"""

from __future__ import annotations

import importlib
import weakref
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike
from scipy import sparse  # type: ignore[import-untyped]

from fem_toolbox.errors import ValidationError
from fem_toolbox.solvers.linear_backend import (
    Kind,
    LinearSolution,
    Numeric,
    _LinearData,
)


def _release(petsc: Any, resources: list[Any]) -> None:
    if not petsc.Sys.isFinalized():
        for resource in reversed(resources):
            resource.destroy()
    resources.clear()


class PetscBackend:
    """Experimental native sequential LU; requires petsc4py and double PETSc.

    Instantiation diagnoses missing/incompatible installations. Importing this
    module does not load petsc4py. Factors may be reused and explicitly closed;
    finalization also releases their native resources. No PETSc options database
    is read by the solver and no alternative factor package is selected.
    """

    def __init__(self) -> None:
        try:
            self._petsc = importlib.import_module("petsc4py.PETSc")
        except (ImportError, OSError, ValueError) as error:
            raise ValidationError(
                "optional petsc4py runtime unavailable or incompatible; see Doc/38_PETSC_BACKEND.md",
                field="backend",
            ) from error
        petsc = self._petsc
        self.scalar_dtype = np.dtype(petsc.ScalarType)
        if self.scalar_dtype not in (np.dtype(np.float64), np.dtype(np.complex128)):
            raise ValidationError(
                "PETSc requires float64 or complex128 scalars", field="backend"
            )
        if petsc.COMM_WORLD.getSize() != 1:
            raise ValidationError(
                "PETSc backend supports one process only", field="backend"
            )
        version = ".".join(str(v) for v in petsc.Sys.getVersion())
        self.name = f"petsc:{version}:{self.scalar_dtype.name}:preonly-lu"
        self._iteration: _IterationOptions | None = None

    def factorize(
        self,
        operator: Any,
        *,
        fields: Mapping[str, ArrayLike],
        scale: ArrayLike | None = None,
        kind: Kind = "general",
        pivot_tolerance: float = 0.0,
        relative_pivot: bool = False,
    ) -> _PetscFactor:
        """Snapshot/validate original data and factor the scaled sequential system."""
        return _PetscFactor(
            self._petsc,
            self.name,
            operator,
            fields,
            scale,
            kind,
            pivot_tolerance,
            relative_pivot,
            self._iteration,
        )


@dataclass(frozen=True)
class _IterationOptions:
    rtol: float
    atol: float
    max_iterations: int
    restart: int


class PetscFieldSplitBackend(PetscBackend):
    """FGMRES with full Schur block factorization and A11 preconditioning.

    Two nonempty fields use native FieldSplit with explicit local index sets;
    a single remaining field uses LU preconditioning. Both block inverses are
    native sequential LU. Physical adapters must remove rigid/gauge nullspaces;
    no loads are projected and no missing boundary conditions are invented.
    Tolerances apply to the scaled KSP residual; original field residuals must
    additionally satisfy the shared contract. Limits apply per native RHS solve.
    """

    def __init__(
        self,
        *,
        rtol: float = 1e-13,
        atol: float = 0.0,
        max_iterations: int = 500,
        restart: int = 50,
    ) -> None:
        for label, value in (("rtol", rtol), ("atol", atol)):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not np.isfinite(value)
                or value < 0
            ):
                raise ValidationError(
                    "expected finite nonnegative tolerance", field=label
                )
        if not 0 < rtol < 1:
            raise ValidationError("expected 0 < rtol < 1", field="rtol")
        for label, value in (("max_iterations", max_iterations), ("restart", restart)):
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValidationError("expected positive integer", field=label)
        super().__init__()
        self._iteration = _IterationOptions(rtol, atol, max_iterations, restart)
        self.name = (
            self.name.rsplit(":", 1)[0]
            + f":fgmres-schur:rtol={rtol}:atol={atol}"
            + f":max_iterations={max_iterations}:restart={restart}"
        )


class _PetscFactor(_LinearData):
    def __init__(
        self,
        petsc: Any,
        name: str,
        operator: Any,
        fields: Mapping[str, ArrayLike],
        scale: ArrayLike | None,
        kind: Kind,
        pivot_tolerance: float,
        relative_pivot: bool,
        iteration: _IterationOptions | None,
    ) -> None:
        super().__init__(operator, fields, scale, kind, pivot_tolerance, relative_pivot)
        self.petsc, self.name = petsc, name
        self.closed = False
        self.iteration = iteration
        self.iterations = 0
        self.pivot_tolerance = pivot_tolerance
        self.relative_pivot = relative_pivot
        self.method = "preonly_lu" if iteration is None else "fgmres_lu"
        self.resources: list[Any] = []
        self._finalizer = weakref.finalize(self, _release, petsc, self.resources)
        self.native_complex = np.dtype(petsc.ScalarType).kind == "c"
        self.real_block = not self.native_complex and self.a.dtype.kind == "c"
        if not self.a.shape[0]:
            return
        if kind == "positive_definite":
            # This hint has a real validation contract, even though PETSc uses
            # LU for both kinds. Only existing bounded dense M6 blocks use it.
            dense = self.scaled.toarray()
            if (
                not np.allclose(dense, dense.conj().T, rtol=1e-12, atol=0)
                or np.min(np.linalg.eigvalsh(dense)) <= 0
            ):
                raise ValidationError(
                    "PETSc positive definite operator required", field="operator"
                )
        a = self.scaled
        if self.real_block:
            a = sparse.bmat([[a.real, -a.imag], [a.imag, a.real]], format="csr")
        a = a.tocsr()
        # SeqAIJ symbolic LU requires structural diagonal entries, including zeros.
        a.setdiag(a.diagonal())
        a.sum_duplicates()
        a.sort_indices()
        limit = np.iinfo(petsc.IntType).max
        if max(a.shape[0], a.nnz) > limit:
            raise ValidationError(
                "operator exceeds PETSc index range", field="operator"
            )
        self.csr = (
            np.asarray(a.indptr, dtype=petsc.IntType),
            np.asarray(a.indices, dtype=petsc.IntType),
            np.asarray(a.data, dtype=petsc.ScalarType),
        )
        try:
            self.matrix = petsc.Mat().createAIJ(
                size=a.shape, csr=self.csr, comm=petsc.COMM_SELF
            )
            self.resources.append(self.matrix)
            self.matrix.assemble()
            self.ksp = petsc.KSP().create(comm=petsc.COMM_SELF)
            self.resources.append(self.ksp)
            # FieldSplit initializes sub-KSPs from options internally. Isolate
            # their prefixes from other solvers without changing global options.
            self.ksp.setOptionsPrefix(f"fem_toolbox_{id(self):x}_")
            self.ksp.setOperators(self.matrix)
            if iteration is None:
                self.ksp.setType("preonly")
                self.ksp.setErrorIfNotConverged(True)
                self._configure_lu(self.ksp)
            else:
                self._configure_iteration()
        except petsc.Error as error:
            self.close()
            raise ValidationError(
                f"PETSc factorization failed: {error}", field="solver"
            ) from error
        except ValidationError:
            self.close()
            raise

    def _configure_lu(self, ksp: Any) -> None:
        petsc = self.petsc
        pc = ksp.getPC()
        pc.setType("lu")
        pc.setFactorSolverType("petsc")
        pc.setFactorOrdering("rcm", nzdiag=np.finfo(np.float64).tiny)
        pc.setFactorShift(petsc.Mat.FactorShiftType.NONE)
        pc.setFactorPivot(zeropivot=0.0)
        ksp.setUp()
        diagonal = pc.getFactorMatrix().getDiagonal()
        try:
            pivots = np.abs(diagonal.getArray(readonly=True)).copy()
        finally:
            diagonal.destroy()
        bound = (
            self.pivot_tolerance * (float(pivots.max()) if self.relative_pivot else 1.0)
            if self.kind == "general"
            else 0.0
        )
        if not np.isfinite(pivots).all() or np.min(pivots) <= bound:
            raise ValidationError(
                "singular or ill-conditioned PETSc pivots; check constraints and nullspaces",
                field="solver",
            )

    def _configure_iteration(self) -> None:
        petsc, ksp, options = self.petsc, self.ksp, self.iteration
        assert options is not None
        fields = [(name, ids) for name, ids in self.fields.items() if len(ids)]
        if len(fields) > 2:
            raise ValidationError(
                "Schur backend supports at most two nonempty fields", field="fields"
            )
        # u first for the supported piezo operator, independent of mapping order.
        fields.sort(key=lambda item: (item[0] != "u", item[0]))
        ksp.setType("fgmres")
        ksp.setPCSide(petsc.PC.Side.RIGHT)
        ksp.setNormType(petsc.KSP.NormType.UNPRECONDITIONED)
        ksp.setTolerances(
            rtol=options.rtol,
            atol=options.atol,
            divtol=1e8,
            max_it=options.max_iterations,
        )
        ksp.setGMRESRestart(options.restart)
        ksp.setErrorIfNotConverged(False)
        if len(fields) == 1:
            self._configure_lu(ksp)
            return
        pc = ksp.getPC()
        pc.setType("fieldsplit")
        pc.setFieldSplitType(petsc.PC.CompositeType.SCHUR)
        pc.setFieldSplitSchurFactType(petsc.PC.SchurFactType.FULL)
        pc.setFieldSplitSchurPreType(petsc.PC.SchurPreType.A11)
        splits = []
        for name, ids in fields:
            if self.real_block:
                ids = np.concatenate((ids, ids + self.a.shape[0]))
            index = petsc.IS().createGeneral(
                np.asarray(ids, dtype=petsc.IntType), comm=petsc.COMM_SELF
            )
            self.resources.append(index)
            splits.append((name, index))
        pc.setFieldSplitIS(*splits)
        ksp.setUp()
        for sub in pc.getFieldSplitSubKSP():
            sub.setType("preonly")
            sub.setErrorIfNotConverged(True)
            self._configure_lu(sub)
        self.method = "fgmres_fieldsplit_schur"

    def close(self) -> None:
        """Release native objects once; subsequent solves are explicitly invalid."""
        self.closed = True
        self._finalizer()

    def _reason_name(self, reason: int) -> str:
        return next(
            (
                name
                for name in dir(self.petsc.KSP.ConvergedReason)
                if name.startswith(("CONVERGED_", "DIVERGED_"))
                and getattr(self.petsc.KSP.ConvergedReason, name) == reason
            ),
            f"petsc_reason_{reason}",
        )

    def _native_solve(self, rhs: Numeric) -> tuple[Numeric, int]:
        petsc = self.petsc
        vector = solution = None
        try:
            vector = petsc.Vec().createWithArray(
                np.asarray(rhs, dtype=petsc.ScalarType), comm=petsc.COMM_SELF
            )
            solution = self.matrix.createVecRight()
            self.ksp.solve(vector, solution)
            reason = int(self.ksp.getConvergedReason())
            count = int(self.ksp.getIterationNumber())
            self.iterations += count
            if reason <= 0:
                raise ValidationError(
                    f"PETSc {self._reason_name(reason)} (reason={reason}), "
                    f"iterations={count}, max_iterations="
                    f"{self.iteration.max_iterations if self.iteration else 0}, "
                    f"scaled_residual={self.ksp.getResidualNorm():.6g}; "
                    "check constraints, nullspaces and preconditioner",
                    field="solver",
                )
            return solution.getArray(readonly=True).copy(), reason
        except petsc.Error as error:
            raise ValidationError(
                f"PETSc solve failed: {error}", field="solver"
            ) from error
        finally:
            if solution is not None:
                solution.destroy()
            if vector is not None:
                vector.destroy()

    def solve(self, rhs: ArrayLike) -> LinearSolution:
        """Solve original RHS columns with native or explicitly embedded scalars."""
        if self.closed:
            raise ValidationError("PETSc factor is closed", field="solver")
        self.iterations = 0
        b, scale, scaled_rhs = self.prepare_rhs(rhs)
        x = np.empty_like(b)
        reason = 0
        if len(b):
            columns = scaled_rhs[:, None] if b.ndim == 1 else scaled_rhs
            values = []
            y: Numeric
            for column in columns.T:
                if self.real_block:
                    raw, reason = self._native_solve(
                        np.concatenate((column.real, column.imag))
                    )
                    y = raw[: len(b)] + 1j * raw[len(b) :]
                elif not self.native_complex and column.dtype.kind == "c":
                    real, reason = self._native_solve(column.real)
                    imag, reason = self._native_solve(column.imag)
                    y = real + 1j * imag
                else:
                    y, reason = self._native_solve(column)
                if b.dtype.kind != "c" and np.iscomplexobj(y):
                    if np.any(y.imag != 0):
                        raise ValidationError(
                            "PETSc returned complex values for real system",
                            field="solver",
                        )
                    y = y.real
                values.append(y)
            # A valid (n,0) RHS needs no KSP solves.
            if values:
                solved = np.column_stack(values)
                x = scale * (solved[:, 0] if b.ndim == 1 else solved)
        method = self.method + ("_real_block" if self.real_block else "")
        return self.solution(
            b,
            x,
            backend=self.name,
            method=method,
            reason=self._reason_name(reason)
            if self.iteration and reason
            else (f"petsc_reason_{reason}" if reason else "empty_system"),
            iterations=self.iterations if self.iteration is not None else None,
        )
