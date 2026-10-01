"""Affine interior patch: analytic transport and leading donor-cell error."""

import numpy as np
import pytest

from fem_toolbox.cfd import CartesianFluidGrid, ClosedWalls
from fem_toolbox.cfd.mac import assemble, convection


@pytest.mark.parametrize("scheme", ["donor-cell", "linear-upwind", "central"])
def test_affine_divergence_free_convective_patch(scheme):
    for n in (8, 16, 32, 64):
        h = 1 / n
        system = assemble(CartesianFluidGrid((n, n), (1, 1)), ClosedWalls(), 1)
        q = np.zeros(len(system.weights))
        # u=x, v=-y. The interior patch is divergence-free and (u.grad)u=(x,y).
        for j, i in zip(*np.nonzero(system.uid >= 0), strict=True):
            q[system.uid[j, i]] = i * h
        for j, i in zip(*np.nonzero(system.vid >= 0), strict=True):
            q[system.vid[j, i]] = -j * h
        matrix, constant = convection(system, q, 1, scheme)
        transport = np.asarray(matrix @ q + constant)
        j = i = n // 2
        donor_error = h / 2 if scheme == "donor-cell" else 0
        assert transport[system.uid[j, i]] == pytest.approx(
            i * h - donor_error, rel=0, abs=1e-12
        )
        assert transport[system.vid[j, i]] == pytest.approx(
            j * h + donor_error, rel=0, abs=1e-12
        )
