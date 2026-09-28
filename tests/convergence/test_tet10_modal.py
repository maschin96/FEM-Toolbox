"""Independent longitudinal bar continuum spectra on four axial Tet10 meshes."""

from collections import Counter
from itertools import combinations

import numpy as np
from numpy.testing import assert_allclose
from scipy.optimize import brentq

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs
from tests.reference.quadratic_mesh import elevate


def bar(cube_mesh, n, short):
    unit = cube_mesh(1)
    points = np.array(
        [
            [x * 0.01, y * 0.01, z * 0.02 / n]
            for x in range(2)
            for y in range(2)
            for z in range(n + 1)
        ]
    )
    cells = []
    for j in range(n):
        mapping = [
            int(x) * 2 * (n + 1) + int(y) * (n + 1) + j + int(z)
            for x, y, z in unit.points
        ]
        cells.extend([[mapping[i] for i in cell] for cell in unit.cells])
    cells = np.array(cells)
    faces = Counter(
        tuple(sorted(face)) for cell in cells for face in combinations(cell, 3)
    )
    mesh = elevate(
        Mesh(
            points=points,
            tetra4=cells,
            triangle3=[face for face, count in faces.items() if count == 1],
        )
    )
    e = np.zeros((3, 6))
    e[2, 2] = 5
    mat = PiezoelectricMaterial(
        np.diag([2e9, 2e9, 1e9, 5e8, 5e8, 5e8]), e, np.eye(3) * 1e-8
    )
    dofs = piezoelectric_dofs(mesh)
    bottom = np.flatnonzero(mesh.points[:, 2] == 0)
    fixed = np.concatenate(
        [
            dofs.indices("u", np.arange(mesh.n_points), ["x", "y"]),
            dofs.indices("u", bottom, ["z"]),
        ]
    )
    z = mesh.points[mesh.facets, 2]
    electrodes = [
        Electrode("ground", np.flatnonzero(np.all(z == 0, axis=1)), potential=0),
        Electrode(
            "top",
            np.flatnonzero(np.all(z == 0.02, axis=1)),
            **({"potential": 0} if short else {"charge": 0}),
        ),
    ]
    return mesh, mat, [Dirichlet(fixed, np.zeros(len(fixed)))], electrodes


def test_four_mesh_open_short_spectra_and_electrode_constraints(cube_mesh):
    cd = 1e9 + 25 / 1e-8
    alpha = (25 / 1e-8) / cd
    root = brentq(lambda x: x / np.tan(x) - alpha, 0.01, np.pi / 2)
    for short, lam in [(False, np.pi / 2), (True, root)]:
        exact = lam / 0.02 * np.sqrt(cd / 7500)
        errors = []
        for n in (1, 2, 4, 8):
            mesh, mat, bc, electrodes = bar(cube_mesh, n, short)
            result = solve_piezoelectric_modes(
                mesh, mat, 7500, bc, electrodes=electrodes, count=1
            )
            errors.append(abs(result.angular_frequencies[0] / exact - 1))
            assert max(result.relative_residual_u) < 1e-9
            assert max(result.relative_residual_phi) < 1e-9
            top = np.flatnonzero(mesh.points[:, 2] == 0.02)
            assert_allclose(
                result.potential[0, top],
                result.electrode_potentials[0, 1],
                rtol=0,
                atol=1e-8,
            )
            if short:
                assert_allclose(result.electrode_potentials, 0, rtol=0, atol=1e-13)
            else:
                scale = 1e-8 * 0.01**2 / 0.02 * abs(result.electrode_potentials[0, 1])
                assert abs(result.electrode_charges[0, 1]) < 1e-10 * scale
        rates = np.log(np.array(errors[:-1]) / errors[1:]) / np.log(2)
        print("Tet10 modal short=", short, "errors=", errors, "rates=", rates)
        assert np.all(np.diff(errors) < 0)
        assert np.all((rates[-2:] > 3.5) & (rates[-2:] < 4.5))
        assert errors[-1] < 2e-6
