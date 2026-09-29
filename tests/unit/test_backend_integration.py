"""Backend substitution at the physical reduction boundary."""

import numpy as np
import pytest

from fem_toolbox.boundary import Dirichlet, Electrode
from fem_toolbox.materials import PiezoelectricMaterial, ThermalMaterial
from fem_toolbox.solvers.heat import solve_heat
from fem_toolbox.solvers.linear_backend import ScipyBackend
from fem_toolbox.solvers.piezo_harmonic import solve_harmonic_piezoelectric
from fem_toolbox.solvers.piezo_transient import solve_transient_piezoelectric
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric
from fem_toolbox.solvers.transient_heat import solve_transient_heat


class RecordingBackend:
    name = "recording-scipy"

    def __init__(self):
        self.calls = []
        self.reports = []

    def factorize(self, operator, **kwargs):
        self.calls.append((operator.copy(), kwargs))
        factor = ScipyBackend().factorize(operator, **kwargs)
        reports = self.reports

        class Factor:
            def solve(self, rhs):
                result = factor.solve(rhs)
                reports.append(result.report)
                return result

        return Factor()


def test_heat_backend_selection_and_time_steps(cube_mesh):
    mesh = cube_mesh(2)
    x = mesh.points[:, 0]
    fixed = np.flatnonzero((x == 0) | (x == 1))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    backend = RecordingBackend()
    result = solve_heat(
        mesh, material, [Dirichlet(fixed, 300 + x[fixed])], backend=backend
    )
    np.testing.assert_allclose(result.temperature, 300 + x, rtol=0, atol=1e-10)
    assert backend.calls[0][0].shape == (9, 9)
    np.testing.assert_array_equal(backend.calls[0][1]["fields"]["T"], np.arange(9))
    assert result.provenance.solver_parameters["backend"] == backend.name
    solve_transient_heat(mesh, material, [0, 0.5, 1], 300, source=12, backend=backend)
    assert len(backend.calls) == 3


@pytest.mark.parametrize("floating", [False, True])
def test_piezo_fields_are_local_after_electrode_condensation(cube_mesh, floating):
    mesh = cube_mesh(2)
    x = mesh.points[:, 0]
    left = np.flatnonzero(np.all(x[mesh.triangle3] == 0, axis=1))
    right = np.flatnonzero(np.all(x[mesh.triangle3] == 1, axis=1))
    dofs = piezoelectric_dofs(mesh)
    ufix = dofs.indices("u", np.flatnonzero(x == 0))
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.5
    material = PiezoelectricMaterial(np.eye(6) * 10, coupling, np.eye(3) * 2)
    backend = RecordingBackend()
    result = solve_piezoelectric(
        mesh,
        material,
        [Dirichlet(ufix, np.zeros(len(ufix)))],
        electrodes=[
            Electrode("left", left, potential=0),
            Electrode("right", right, charge=0)
            if floating
            else Electrode("right", right, potential=1),
        ],
        backend=backend,
    )
    matrix, kwargs = backend.calls[0]
    assert len(kwargs["fields"]["u"]) == 54
    assert len(kwargs["fields"]["phi"]) == (10 if floating else 9)
    np.testing.assert_array_equal(
        kwargs["fields"]["phi"], np.arange(54, matrix.shape[0])
    )
    assert max(result.relative_residual_u, result.relative_residual_phi) < 1e-9


def test_m6_linear_subproblems_use_selected_backend(cube_mesh):
    mesh = cube_mesh(1)
    n = mesh.n_points
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 0.5
    material = PiezoelectricMaterial(np.eye(6) * 10, coupling, np.eye(3) * 2)
    # One supported mechanical unknown and one electric unknown, no electrodes.
    free = [3, 3 * n + 1]
    fixed = np.setdiff1d(np.arange(4 * n), free)
    bc = [Dirichlet(fixed, np.zeros(len(fixed)))]
    backend = RecordingBackend()
    f = np.zeros(4 * n, dtype=complex)
    f[3] = 1 + 2j
    result = solve_harmonic_piezoelectric(
        mesh, material, 6.0, bc, 0.1, loads=f, backend=backend
    )
    assert result.displacement[1, 0].imag != 0
    assert len(backend.calls) == 3  # electric map, mechanical phasor, electric RHS
    assert all(set(k["fields"]) in ({"u"}, {"phi"}) for _, k in backend.calls)
    backend = RecordingBackend()
    solve_transient_piezoelectric(
        mesh, material, 6.0, bc, [0, 0.1, 0.2], backend=backend
    )
    assert len(backend.calls) == 5  # electric map/inverse, mass, two Newmark operators
    assert len(backend.reports) >= 6


def test_elastic_and_electric_backends_preserve_affine_fields(cube_mesh):
    from fem_toolbox.materials import DielectricMaterial, ElasticMaterial
    from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity
    from fem_toolbox.solvers.electrostatics import solve_electrostatics

    mesh = cube_mesh(2)
    boundary = np.flatnonzero(np.any((mesh.points == 0) | (mesh.points == 1), axis=1))
    backend = RecordingBackend()
    u = 0.01 * mesh.points
    bc = Dirichlet(mechanical_dofs(mesh).indices("u", boundary), u[boundary].ravel())
    elastic = solve_elasticity(
        mesh, ElasticMaterial.isotropic(100, 0.25), [bc], backend=backend
    )
    np.testing.assert_allclose(elastic.displacement, u, rtol=0, atol=1e-14)
    assert set(backend.calls[0][1]["fields"]) == {"u"}
    phi = mesh.points[:, 0]
    electric = solve_electrostatics(
        mesh,
        DielectricMaterial.isotropic(2),
        [Dirichlet(boundary, phi[boundary])],
        backend=backend,
    )
    np.testing.assert_allclose(electric.potential, phi, rtol=0, atol=1e-14)
    assert set(backend.calls[1][1]["fields"]) == {"phi"}
