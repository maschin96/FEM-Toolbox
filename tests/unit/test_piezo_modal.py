"""Independent scalar oscillator references for mechanical/electrical reduction."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.errors import ElectricGaugeError, ValidationError
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezo_dynamics import solve_piezoelectric_modes


def oscillator(short=False):
    mesh = Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], tetra4=[[0, 1, 2, 3]]
    )
    coupling = np.zeros((3, 6))
    coupling[0, 0] = 2
    mat = PiezoelectricMaterial(
        np.diag([10, 20, 30, 4, 5, 6]), coupling, np.diag([3, 4, 5])
    )
    fixed = np.setdiff1d(np.arange(16), [3] if short else [3, 13])
    return mesh, mat, [Dirichlet(fixed, np.zeros(len(fixed)))]


@pytest.mark.parametrize("short", [False, True])
def test_one_mechanical_dof_frequency_and_mass_normalization(short):
    mesh, mat, bc = oscillator(short)
    result = solve_piezoelectric_modes(mesh, mat, 6.0, bc, count=1)
    # Unit right tetrahedron: V=1/6, int(N1²)=V/10, Bxx=1, Gx=1.
    mass = 6 / 60
    stiffness = (10 if short else 10 + 4 / 3) / 6
    omega = np.sqrt(stiffness / mass)
    assert_allclose(result.angular_frequencies, [omega], rtol=2e-14)
    assert_allclose(result.frequencies_hz, [omega / (2 * np.pi)], rtol=2e-14)
    assert_allclose(result.displacement[0, 1, 0], 1 / np.sqrt(mass), rtol=2e-14)
    assert_allclose(
        result.potential[0, 1], 0 if short else 2 / 3 / np.sqrt(mass), atol=1e-14
    )
    assert np.max(result.relative_residual_u) < 1e-13
    assert np.max(result.relative_residual_phi) < 1e-13
    assert result.displacement.shape == (1, 4, 3)
    assert result.provenance.solver_parameters["electrical_mass"] == 0


@pytest.mark.parametrize("density", [0, -1, np.nan, True, 1 + 1j, [1, 2]])
def test_invalid_density(density):
    mesh, mat, bc = oscillator()
    with pytest.raises(ValidationError):
        solve_piezoelectric_modes(mesh, mat, density, bc, count=1)


@pytest.mark.parametrize("count", [0, -1, 2, True, 1.5])
def test_invalid_mode_count(count):
    mesh, mat, bc = oscillator()
    with pytest.raises(ValidationError):
        solve_piezoelectric_modes(mesh, mat, 6.0, bc, count=count)


def test_missing_gauge_and_nonzero_modal_conditions_rejected():
    mesh, mat, bc = oscillator()
    with pytest.raises(ElectricGaugeError):
        solve_piezoelectric_modes(
            mesh, mat, 6.0, [Dirichlet(np.arange(12), np.zeros(12))]
        )
    with pytest.raises(ValidationError, match="homogeneous"):
        solve_piezoelectric_modes(mesh, mat, 6.0, [*bc, Dirichlet([13], [1.0])])


def test_result_owns_readonly_arrays_and_stable_provenance():
    mesh, mat, bc = oscillator()
    result = solve_piezoelectric_modes(mesh, mat, [6.0], bc, count=1)
    again = solve_piezoelectric_modes(mesh, mat, 6.0, bc, count=1)
    assert result.provenance == again.provenance
    with pytest.raises(ValueError):
        result.displacement[0, 0, 0] = 3


def test_two_modes_orthogonal_with_independent_diagonal_reference():
    mesh, mat, _ = oscillator()
    fixed = np.setdiff1d(np.arange(16), [3, 7])
    result = solve_piezoelectric_modes(
        mesh, mat, 6.0, [Dirichlet(fixed, np.zeros(len(fixed)))], count=2
    )
    assert_allclose(result.angular_frequencies**2, [10 / 0.6, 20 / 0.6], rtol=3e-14)
    u = result.displacement.reshape(2, 12)
    # Different Cartesian components have zero consistent-mass cross term.
    assert_allclose(0.1 * u @ u.T, np.eye(2), atol=3e-14)
    assert_allclose(result.potential, 0, atol=1e-15)


def test_active_node_order_and_geometry_density_scaling():
    mesh, mat, bc = oscillator()
    reference = solve_piezoelectric_modes(mesh, mat, 6.0, bc, count=1)
    shifted = Mesh(
        points=np.vstack(([9, 9, 9], 2 * mesh.points)), tetra4=mesh.tetra4 + 1
    )
    result = solve_piezoelectric_modes(shifted, mat, 24.0, bc, count=1)
    assert_allclose(result.node_ids, [1, 2, 3, 4])
    assert_allclose(
        result.angular_frequencies, reference.angular_frequencies / 4, rtol=2e-14
    )
    assert result.provenance != reference.provenance


def test_fully_fixed_and_rigid_mechanical_systems_rejected():
    mesh, mat, _ = oscillator()
    with pytest.raises(ValidationError, match="no free mechanical"):
        solve_piezoelectric_modes(
            mesh, mat, 6.0, [Dirichlet(np.arange(16), np.zeros(16))], count=1
        )
    with pytest.raises(ValidationError, match="rigid"):
        solve_piezoelectric_modes(mesh, mat, 6.0, [Dirichlet([12], [0.0])], count=1)


def test_dense_limit_fails_before_assembly(cube_mesh):
    _, mat, _ = oscillator()
    with pytest.raises(ValidationError, match="2048"):
        solve_piezoelectric_modes(cube_mesh(8), mat, 6.0, [], count=1)


def test_nonzero_electrode_perturbation_rejected(cube_mesh):
    from fem_toolbox.boundary import Electrode

    _, mat, _ = oscillator()
    with pytest.raises(ValidationError, match="homogeneous electrodes"):
        solve_piezoelectric_modes(
            cube_mesh(1),
            mat,
            6.0,
            [],
            electrodes=[Electrode("bias", [0], potential=1)],
            count=1,
        )
