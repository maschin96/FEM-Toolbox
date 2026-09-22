"""Stand-alone dense 1D piezo bar FEM, independent of production modules."""

import numpy as np


def solve_bar(layers, force, voltage):
    """Solve H=C u'²/2+e u'φ'-κ φ'²/2 with interleaved nodal unknowns.

    Each layer tuple is (C33, e33, kappa33). Length and cross section are
    unity. u=phi=0 at z=0, phi=voltage and force=force at z=1. The dense
    solver has no dependency on production geometry, assembly, or BC code.
    Returns nodal u, nodal phi, high-electrode charge, and electric enthalpy.
    """
    n = len(layers)
    h = 1 / n
    matrix = np.zeros((2 * (n + 1), 2 * (n + 1)))
    derivative = np.array([[-1 / h, 0, 1 / h, 0], [0, -1 / h, 0, 1 / h]])
    for i, (c, e, k) in enumerate(layers):
        tangent = np.array([[c, e], [e, -k]])
        local = h * derivative.T @ tangent @ derivative
        dofs = np.arange(2 * i, 2 * i + 4)
        matrix[np.ix_(dofs, dofs)] += local
    rhs = np.zeros(2 * (n + 1))
    rhs[-2] = force
    fixed = np.array([0, 1, 2 * n + 1])
    solution = np.zeros(len(rhs))
    solution[fixed[-1]] = voltage
    free = np.setdiff1d(np.arange(len(rhs)), fixed)
    solution[free] = np.linalg.solve(
        matrix[np.ix_(free, free)],
        rhs[free] - matrix[np.ix_(free, fixed)] @ solution[fixed],
    )
    reaction = matrix @ solution - rhs
    return (
        solution[::2],
        solution[1::2],
        -reaction[-1],
        0.5 * solution @ matrix @ solution,
    )
