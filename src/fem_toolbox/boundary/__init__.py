"""Symmetry-preserving Dirichlet elimination on a copy of a linear system."""

from collections.abc import Iterable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy import sparse  # type: ignore[import-untyped]

from fem_toolbox.assembly import LinearSystem, indices, real_array
from fem_toolbox.errors import ValidationError


@dataclass(frozen=True)
class Dirichlet:
    """Global DOFs and prescribed SI field values, each shape (n,).

    Resolve field selections with DofMap.indices; evaluate spatial/time callables
    before constructing this algebraic condition. Duplicate values must agree
    exactly; no implicit unit-dependent conflict tolerance is applied.
    """

    dofs: ArrayLike
    values: ArrayLike


def eliminate_dirichlet(
    system: LinearSystem, conditions: Iterable[Dirichlet]
) -> LinearSystem:
    """Return P K P + I_fixed and P(f-K g)+g; preserve original reactions.

    P selects free DOFs and g holds prescribed values (zero elsewhere). Both
    rows and columns are eliminated, so symmetry is retained. No solve occurs.
    """
    assigned: dict[int, float] = {}
    for condition in conditions:
        dofs = indices(condition.dofs, system.size)
        values = real_array(condition.values, (len(dofs),))
        for dof, value in zip(dofs, values, strict=True):
            if int(dof) in assigned and assigned[int(dof)] != value:
                raise ValidationError(
                    f"conflicting Dirichlet values at DOF {dof}", field="Dirichlet"
                )
            assigned[int(dof)] = float(value)
    g = np.zeros(system.size)
    free = np.ones(system.size)
    for dof, value in assigned.items():
        g[dof], free[dof] = value, 0
    matrix = system.matrix
    rhs = free * (system.rhs - matrix @ g) + g
    projection = sparse.diags(free, format="csr")
    constrained = projection @ matrix @ projection + sparse.diags(1 - free)
    return LinearSystem(constrained.asformat(matrix.format), rhs)
