"""Four axial meshes against the independent driven piezo bar continuum."""

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
from tests.convergence.test_tet10_modal import bar


def driven_bar(cube_mesh, n, short):
    mesh, mat, bc, electrodes = bar(cube_mesh, n, short)
    weights = np.zeros(mesh.n_points)
    for face in mesh.facets[np.all(mesh.points[mesh.facets, 2] == 0.02, axis=1)]:
        p = mesh.points[face[:3]]
        area = np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0])) / 2
        weights[face[3:]] += area / 3
    loads = np.zeros(4 * mesh.n_points, dtype=complex)
    traction = 1000 + 2000j
    loads[np.arange(mesh.n_points) * 3 + 2] = weights * traction
    return mesh, mat, bc, electrodes, loads, weights


def test_four_mesh_harmonic_bar_open_short(cube_mesh):
    length = 0.02
    area = 0.01**2
    cd = 1e9 + 25 / 1e-8
    lam = 0.6
    wave = lam / length
    omega = wave * np.sqrt(cd / 7500)
    for short in (False, True):
        denominator = cd * wave * np.cos(lam) - (
            25 / 1e-8 / length * np.sin(lam) if short else 0
        )
        tip = (1000 + 2000j) / denominator * np.sin(lam)
        errors = []
        for n in (1, 2, 4, 8):
            mesh, mat, bc, electrodes, loads, weights = driven_bar(cube_mesh, n, short)
            result = solve_harmonic_piezoelectric(
                mesh,
                mat,
                7500,
                bc,
                omega / (2 * np.pi),
                electrodes=electrodes,
                loads=loads,
            )
            average = weights @ result.displacement[:, 2] / area
            errors.append(abs(average / tip - 1))
            top = np.flatnonzero(mesh.points[:, 2] == length)
            assert_allclose(
                result.potential[top], result.potential[top[0]], rtol=0, atol=1e-10
            )
            if short:
                assert_allclose(result.potential[top], 0, rtol=0, atol=1e-12)
            else:
                scale = 1e-8 * area / length * abs(result.potential[top[0]])
                assert abs(result.residual_phi[top].sum()) < 1e-10 * scale
            assert result.relative_residual_u < 1e-9
            assert result.relative_residual_phi < 1e-9
        rates = np.log2(np.array(errors[:-1]) / errors[1:])
        print("short=", short, "errors=", errors, "rates=", rates)
        assert np.all(np.diff(errors) < 0)
        assert np.all((rates[-2:] > 3.5) & (rates[-2:] < 4.5))
        assert errors[-1] < 2e-6
