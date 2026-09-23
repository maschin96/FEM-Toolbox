"""Longitudinal conducting bar: continuum reference independent of Tet4 FEM."""

import numpy as np
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs


def conducting_bar(cube_mesh, n, short):
    unit = cube_mesh(n)
    mesh = Mesh(
        points=unit.points * [0.01, 0.01, 0.02],
        tetra4=unit.tetra4,
        triangle3=unit.triangle3,
    )
    e = np.zeros((3, 6))
    e[2, 2] = 5
    mat = PiezoelectricMaterial(
        np.diag([2e9, 2e9, 1e9, 5e8, 5e8, 5e8]), e, np.eye(3) * 1e-8
    )
    dofs = piezoelectric_dofs(mesh)
    bottom = np.flatnonzero(mesh.points[:, 2] == 0)
    fixed = np.concatenate(
        (
            dofs.indices("u", np.arange(mesh.n_points), ["x", "y"]),
            dofs.indices("u", bottom, ["z"]),
        )
    )
    z = mesh.points[mesh.triangle3, 2]
    electrodes = [
        Electrode("ground", np.flatnonzero(np.all(z == 0, axis=1)), potential=0),
        Electrode(
            "top",
            np.flatnonzero(np.all(z == 0.02, axis=1)),
            **({"potential": 0} if short else {"charge": 0}),
        ),
    ]
    return mesh, mat, [Dirichlet(fixed, np.zeros(len(fixed)))], electrodes


def test_open_and_short_continuum_frequency_and_mesh_convergence(cube_mesh):
    # u(0)=0, T(L)=0, rho u_tt=T_z, D_z=0, CD=C+e²/kappa.
    # Open: D=0 -> lambda=pi/2. Short: int(E dz)=0 ->
    # D=e*u(L)/L and lambda*cos(lambda)-alpha*sin(lambda)=0.
    cd = 1e9 + 25 / 1e-8
    alpha = (25 / 1e-8) / cd
    short_root = brentq(lambda x: x / np.tan(x) - alpha, 0.01, np.pi / 2)
    frequencies = []
    for short, root in [(False, np.pi / 2), (True, short_root)]:
        reference = root * np.sqrt(cd / 7500) / (2 * np.pi * 0.02)
        errors = []
        for n in (1, 2, 4):
            mesh, mat, bc, electrodes = conducting_bar(cube_mesh, n, short)
            result = solve_piezoelectric_modes(
                mesh, mat, 7500, bc, electrodes=electrodes, count=1
            )
            errors.append(abs(result.frequencies_hz[0] / reference - 1))
            assert max(result.relative_residual_u) < 1e-9
            assert max(result.relative_residual_phi) < 1e-9
            if short:
                assert_allclose(result.electrode_potentials, 0, atol=1e-13)
            else:
                # Charge cancellation roundoff relative to the modal voltage capacitance.
                scale = 1e-8 * 0.01**2 / 0.02 * abs(result.electrode_potentials[0, 1])
                assert abs(result.electrode_charges[0, 1]) < 1e-10 * scale
        assert errors[2] < errors[1] < errors[0]
        assert errors[2] < 0.02
        frequencies.append(result.frequencies_hz[0])
    assert frequencies[0] > frequencies[1]
