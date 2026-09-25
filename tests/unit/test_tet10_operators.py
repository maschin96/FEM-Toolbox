"""Quadratic element operators against analytic energy and simplex moments."""

from math import factorial

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.elements import Tet10Geometry, Triangle6Geometry
from fem_toolbox.materials import (
    DielectricMaterial,
    ElasticMaterial,
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
)

CORNERS = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
POINTS = np.vstack([CORNERS, [(CORNERS[a] + CORNERS[b]) / 2 for a, b in EDGES]])


def material():
    e = np.zeros((3, 6))
    e[0, 0] = 2
    return PiezoelectricMaterial(
        np.diag([10.0, 20, 30, 4, 5, 6]), e, np.diag([3.0, 4, 5])
    )


def test_scalar_quadratic_energy_and_mechanical_rigid_modes():
    from fem_toolbox.elements.quadratic import (
        conductance_matrix,
        dielectric_matrix,
        stiffness_matrix,
    )

    g = Tet10Geometry(POINTS)
    t = POINTS[:, 0] ** 2
    heat = conductance_matrix(g, ThermalMaterial.isotropic(2))
    electric = dielectric_matrix(g, DielectricMaterial.isotropic(2))
    # grad(x²) k grad(x²) = 8 x², integral x² = 1/60.
    assert_allclose(t @ heat @ t, 8 / 60, atol=2e-15)
    assert_allclose(heat, electric, atol=0, rtol=0)
    assert_allclose(heat @ np.ones(10), 0, atol=3e-15)
    k = stiffness_matrix(g, ElasticMaterial(np.eye(6)))
    assert_allclose(k, k.T, atol=2e-15)
    eigenvalues = np.linalg.eigvalsh(k)
    assert np.count_nonzero(eigenvalues > 1e-12) == 24
    for axis in np.eye(3):
        assert_allclose(k @ np.tile(axis, 10), 0, atol=3e-15)
        assert_allclose(k @ np.cross(axis, POINTS).ravel(), 0, atol=3e-15)


def test_storage_against_independent_barycentric_polynomial_integrals():
    from fem_toolbox.elements.quadratic import capacity_matrix, mass_matrix

    # Represent each basis as a polynomial in all four barycentric coordinates.
    polynomials = []
    for i in range(4):
        linear = [0] * 4
        linear[i] = 1
        square = [0] * 4
        square[i] = 2
        polynomials.append({tuple(linear): -1, tuple(square): 2})
    for a, b in EDGES:
        powers = [0] * 4
        powers[a] = powers[b] = 1
        polynomials.append({tuple(powers): 4})
    expected = np.zeros((10, 10))
    for i, p in enumerate(polynomials):
        for j, q in enumerate(polynomials):
            for a, ca in p.items():
                for b, cb in q.items():
                    powers = np.array(a) + b
                    expected[i, j] += (
                        ca
                        * cb
                        * np.prod([factorial(int(k)) for k in powers])
                        / factorial(int(powers.sum()) + 3)
                    )
    g = Tet10Geometry(POINTS)
    c = capacity_matrix(g, 2, 3)
    assert_allclose(c, 6 * expected, atol=2e-16, rtol=3e-14)
    assert np.linalg.eigvalsh(c).min() > 0
    assert_allclose(mass_matrix(g, 2), np.kron(2 * expected, np.eye(3)), atol=2e-16)
    # Row sums at corners are negative: never silently use P1 row-sum lumping.
    assert np.all(c.sum(axis=1)[:4] < 0)
    diagonal = capacity_matrix(g, 2, 3, lumping="hrz")
    assert np.all(np.diag(diagonal) > 0)
    assert_allclose(diagonal, np.diag(np.diag(c) / np.trace(c)), atol=3e-16)
    assert_allclose(diagonal.sum(), 1, atol=3e-16)


def test_piezo_affine_patch_matches_surface_force_and_displacement_flux():
    from fem_toolbox.elements.quadratic import piezoelectric_matrix
    from fem_toolbox.physics.quadratic import surface_charge_load, traction_load

    g = Tet10Geometry(POINTS)
    u = np.zeros((10, 3))
    u[:, 0] = 0.2 * POINTS[:, 0]
    phi = 0.3 * POINTS[:, 0]
    matrix = piezoelectric_matrix(g, material())
    assert_allclose(matrix, matrix.T, atol=4e-15)
    assert_allclose(matrix[30:, :30], matrix[:30, 30:].T, atol=0)
    rhs = np.zeros(40)
    stress = np.diag([2.6, 0, 0])
    displacement = np.array([-0.5, 0, 0])
    for face in g.faces:
        nodes = np.array(face.local_nodes)
        surface = Triangle6Geometry(POINTS[nodes])
        forces = traction_load(surface, stress @ face.unit_normal).reshape(6, 3)
        rhs[:30].reshape(10, 3)[nodes] += forces
        rhs[30:][nodes] += surface_charge_load(surface, displacement @ face.unit_normal)
    state = np.r_[u.ravel(), phi]
    assert_allclose(matrix @ state, rhs, atol=2e-15)
    assert_allclose(
        0.5 * state @ matrix @ state,
        (0.5 * 10 * 0.2**2 + 2 * 0.2 * 0.3 - 0.5 * 3 * 0.3**2) / 6,
        atol=2e-16,
    )


def test_thermal_coupling_free_expansion_and_electric_switch():
    from fem_toolbox.elements.quadratic import piezoelectric_matrix, thermal_load

    g = Tet10Geometry(POINTS)
    expansion = ThermalExpansion([0.1, 0, 0, 0, 0, 0], 300)
    # T=300+x, u_x=.05*x² gives epsilon_xx=.1*x exactly.
    temperature = 300 + POINTS[:, 0]
    u = np.zeros((10, 3))
    u[:, 0] = 0.05 * POINTS[:, 0] ** 2
    tangent = piezoelectric_matrix(g, material())
    force = thermal_load(
        g, material(), temperature, expansion, electric_thermal_strain=True
    )
    assert_allclose(tangent @ np.r_[u.ravel(), np.zeros(10)], force, atol=3e-15)
    mechanical_only = thermal_load(g, material(), temperature, expansion)
    assert_allclose(mechanical_only[:30], force[:30], atol=0)
    assert_allclose(mechanical_only[30:], 0, atol=0)


def test_signed_loads_and_robin_constant_field_balance():
    from fem_toolbox.physics.heat import Convection
    from fem_toolbox.physics.quadratic import (
        body_load,
        convection_terms,
        flux_load,
        source_load,
        surface_charge_load,
        traction_load,
        volume_charge_load,
    )

    g = Tet10Geometry(POINTS)
    # Integral Ni: -V/20 at corners, V/5 at mid-edge nodes.
    integral = np.r_[np.full(4, -1 / 120), np.full(6, 1 / 30)]
    assert_allclose(source_load(g, 2), 2 * integral, atol=2e-16)
    assert_allclose(volume_charge_load(g, 2), -2 * integral, atol=2e-16)
    assert_allclose(
        body_load(g, [1, 2, 3]).reshape(10, 3),
        np.outer(integral, [1, 2, 3]),
        atol=3e-16,
    )
    face = Triangle6Geometry(POINTS[[0, 1, 2, 4, 5, 6]])
    surface = np.r_[np.zeros(3), np.full(3, 1 / 6)]
    assert_allclose(flux_load(face, 2), -2 * surface, atol=3e-16)
    assert_allclose(surface_charge_load(face, 2), 2 * surface, atol=3e-16)
    assert_allclose(
        traction_load(face, [1, 2, 3]).reshape(6, 3),
        np.outer(surface, [1, 2, 3]),
        atol=4e-16,
    )
    h, f = convection_terms(face, Convection(2, 300))
    assert_allclose(h @ np.full(6, 300), f, atol=6e-14)
    assert np.linalg.eigvalsh(h).min() > 0


def test_invalid_storage_and_nonfinite_fields_rejected():
    from fem_toolbox.elements.quadratic import capacity_matrix
    from fem_toolbox.errors import ValidationError
    from fem_toolbox.physics.quadratic import source_load

    g = Tet10Geometry(POINTS)
    for density in (0, -1, np.nan, np.inf):
        with pytest.raises(ValidationError):
            capacity_matrix(g, density, 2)
    with pytest.raises(ValidationError, match="lumping"):
        capacity_matrix(g, 1, 1, lumping="row_sum")
    with pytest.raises(ValidationError):
        source_load(g, lambda p: np.nan)


@pytest.mark.parametrize("floating", [False, True])
def test_electrode_reduction_includes_all_six_nodes(floating):
    from scipy import sparse

    from fem_toolbox.assembly import LinearSystem
    from fem_toolbox.boundary import Electrode
    from fem_toolbox.fields import DofMap, Field
    from fem_toolbox.mesh import Mesh
    from fem_toolbox.solvers.piezoelectric import _electrode_reduction

    face = [0, 2, 1, 6, 5, 4]
    mesh = Mesh(
        points=POINTS,
        tetra10=[np.arange(10)],
        triangle6=[face],
        facet_sets={"port": [0]},
    )
    dofs = DofMap([Field("phi", ("phi",), tuple(range(10)), "V")])
    system = LinearSystem(sparse.eye(10, format="csr"), np.zeros(10))
    electrode = (
        Electrode("port", "port", charge=2)
        if floating
        else Electrode("port", [0], potential=3)
    )
    reduction = _electrode_reduction(mesh, dofs, system, [], [electrode])
    assert set(reduction.data[0][1]) == set(face)
    if floating:
        assert len(set(reduction.index[face])) == 1
        assert reduction.system.matrix.shape == (5, 5)
        assert_allclose(reduction.system.rhs.sum(), -2, atol=0)
    else:
        assert set(reduction.conditions[0].dofs) == set(face)


def test_quadratic_coupled_patch_with_anisotropy_and_body_charge():
    from fem_toolbox.elements.quadratic import piezoelectric_matrix
    from fem_toolbox.physics.quadratic import (
        body_load,
        surface_charge_load,
        traction_load,
        volume_charge_load,
    )

    # Independent polynomial strain, stress and D and their analytic divergence.
    c = np.diag([10.0, 20, 30, 4, 5, 6]) + np.ones((6, 6))
    e = np.arange(18).reshape(3, 6) / 10
    kappa = np.diag([3.0, 4, 5]) + np.ones((3, 3))
    mat = PiezoelectricMaterial(c, e, kappa)
    g = Tet10Geometry(POINTS)
    x, y, z = POINTS.T
    u = np.column_stack([x * x + y * z, y * y + x * z, z * z + x * y])
    phi = x * x + y * y + z * z
    eps_derivative = np.vstack([2 * np.eye(3), 2 * np.eye(3)])
    stress_derivative = c @ eps_derivative + 2 * e.T
    d_derivative = e @ eps_derivative - 2 * kappa
    force = -np.array(
        [
            stress_derivative[0, 0] + stress_derivative[5, 1] + stress_derivative[4, 2],
            stress_derivative[5, 0] + stress_derivative[1, 1] + stress_derivative[3, 2],
            stress_derivative[4, 0] + stress_derivative[3, 1] + stress_derivative[2, 2],
        ]
    )
    rhs = np.r_[body_load(g, force), volume_charge_load(g, np.trace(d_derivative))]

    def stress_tensor(p):
        s = stress_derivative @ p
        return np.array([[s[0], s[5], s[4]], [s[5], s[1], s[3]], [s[4], s[3], s[2]]])

    for face in g.faces:
        nodes = np.array(face.local_nodes)
        surface = Triangle6Geometry(POINTS[nodes])
        normal = face.unit_normal
        rhs[:30].reshape(10, 3)[nodes] += traction_load(
            surface, lambda p, n=normal: stress_tensor(p) @ n
        ).reshape(6, 3)
        rhs[30:][nodes] += surface_charge_load(
            surface, lambda p, n=normal: n @ d_derivative @ p
        )
    matrix = piezoelectric_matrix(g, mat)
    state = np.r_[u.ravel(), phi]
    assert_allclose(matrix @ state, rhs, atol=2e-14, rtol=2e-14)
    # Negative input orientation changes only the local permutation of DOFs.
    order = [0, 2, 1, 3, 6, 5, 4, 7, 9, 8]
    changed = piezoelectric_matrix(Tet10Geometry(POINTS[order]), mat)
    dofs = np.r_[
        (3 * np.array(order)[:, None] + np.arange(3)).ravel(), 30 + np.array(order)
    ]
    assert_allclose(changed, matrix[np.ix_(dofs, dofs)], atol=3e-14, rtol=3e-14)


def test_consistent_storage_reproduces_affine_kinetic_energy():
    from fem_toolbox.elements.quadratic import (
        capacity_matrix,
        mass_matrix,
        stiffness_matrix,
        thermal_load,
    )

    g = Tet10Geometry(POINTS)
    velocity = POINTS.copy()
    # rho/2 integral(x²+y²+z²) = 2/2 * 3/60.
    assert_allclose(
        0.5 * velocity.ravel() @ mass_matrix(g, 2) @ velocity.ravel(),
        1 / 20,
        atol=2e-16,
    )
    assert_allclose(np.ones(10) @ capacity_matrix(g, 2, 3) @ np.ones(10), 1, atol=4e-16)
    elastic = ElasticMaterial(np.eye(6))
    expansion = ThermalExpansion([0.1, 0, 0, 0, 0, 0], 300)
    u = np.zeros((10, 3))
    u[:, 0] = 0.05 * POINTS[:, 0] ** 2
    assert_allclose(
        stiffness_matrix(g, elastic) @ u.ravel(),
        thermal_load(g, elastic, 300 + POINTS[:, 0], expansion),
        atol=2e-16,
    )


def test_thermal_and_load_validation_is_explicit():
    from fem_toolbox.elements.quadratic import (
        capacity_matrix,
        conductance_matrix,
        thermal_load,
    )
    from fem_toolbox.errors import ValidationError
    from fem_toolbox.physics.heat import Convection
    from fem_toolbox.physics.quadratic import body_load, convection_terms

    g = Tet10Geometry(POINTS)
    expansion = ThermalExpansion(np.zeros(6), 300)
    for temperature in (np.zeros(4), np.full(10, np.nan), np.zeros(10, dtype=complex)):
        with pytest.raises(ValidationError):
            thermal_load(g, material(), temperature, expansion)
    for switch in (1, "yes", None):
        with pytest.raises(ValidationError):
            thermal_load(
                g,
                material(),
                np.full(10, 300),
                expansion,
                electric_thermal_strain=switch,
            )
    with pytest.raises(ValidationError):
        thermal_load(
            g,
            ElasticMaterial(np.eye(6)),
            np.full(10, 300),
            expansion,
            electric_thermal_strain=True,
        )
    with pytest.raises(ValidationError):
        conductance_matrix(g, material())
    with pytest.raises(ValidationError):
        body_load(g, lambda p: [1, 2])
    for heat in (0, -1, np.inf):
        with pytest.raises(ValidationError):
            capacity_matrix(g, 1, heat)
    with pytest.raises(ValidationError, match="representable"):
        capacity_matrix(g, 1e308, 1e308)
    face = Triangle6Geometry(POINTS[[0, 1, 2, 4, 5, 6]])
    h, f = convection_terms(face, Convection(0, 300))
    assert_allclose(h, 0, atol=0)
    assert_allclose(f, 0, atol=0)
