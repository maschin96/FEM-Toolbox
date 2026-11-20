"""Independent affine, tensor and conservation expectations for Quad FV."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain, FluidMesh2D, QuadFVGeometry
from fem_toolbox.errors import ValidationError


@pytest.mark.parametrize("shear", [0.0, np.tan(np.pi / 9)])
def test_affine_scalar_vector_and_face_gradient(quad_mesh, shear):
    mesh = quad_mesh(shear=shear)
    fv = QuadFVGeometry(mesh, width_m=2.5)
    x = mesh.cell_centres_m
    b = mesh.face_centres_m[mesh.boundary_face_indices]

    def scalar(x):
        return 2 + x @ np.array([3.0, -4.0])

    gradient = fv.gradient(scalar(x), scalar(b))
    np.testing.assert_allclose(
        gradient, np.tile([3, -4], (len(x), 1)), atol=1e-11, rtol=0
    )
    np.testing.assert_allclose(
        fv.face_values(scalar(x), scalar(b)),
        scalar(mesh.face_centres_m),
        atol=1e-11,
        rtol=0,
    )
    fg = fv.face_gradient(scalar(x), scalar(b))
    np.testing.assert_allclose(
        fg.gradient, np.tile([3, -4], (len(mesh.faces), 1)), atol=1e-11, rtol=0
    )
    np.testing.assert_allclose(
        fg.compact_normal + fg.nonorthogonal_normal,
        np.einsum("ij,ij->i", fg.gradient, mesh.face_normals),
        atol=1e-11,
        rtol=0,
    )
    a = np.array([[2.0, 3.0], [-1.0, 4.0]])

    def u(x):
        return x @ a.T + [1, 2]

    np.testing.assert_allclose(
        fv.gradient(u(x), u(b)), np.tile(a, (len(x), 1, 1)), atol=1e-11, rtol=0
    )


def test_full_variable_stress_and_global_force(quad_mesh):
    mesh = quad_mesh(shear=np.tan(np.pi / 9))
    fv = QuadFVGeometry(mesh, width_m=3)
    x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]
    a = np.array([[2.0, 3.0], [-1.0, 4.0]])

    def u(x):
        return x @ a.T

    def mu(x):
        return 2 + x[:, 0] + 2 * x[:, 1]

    def p(x):
        return 5 + 3 * x[:, 0] - x[:, 1]

    result = fv.stress_flux(
        u(x),
        u(b),
        viscosity_Pa_s=mu(x),
        boundary_viscosity_Pa_s=mu(b),
        pressure_Pa=p(x),
        boundary_pressure_Pa=p(b),
    )
    sigma = mu(mesh.face_centres_m)[:, None, None] * (a + a.T) - p(mesh.face_centres_m)[
        :, None, None
    ] * np.eye(2)
    np.testing.assert_allclose(result.stress_Pa, sigma, atol=1e-11, rtol=0)
    # div(mu*(A+A^T)-pI) = (A+A^T)*(1,2) - (3,-1) = (5,19).
    expected = mesh.cell_areas_m2[:, None] * 3 * np.array([5, 19])
    np.testing.assert_allclose(result.cell_force_N, expected, atol=1e-11, rtol=0)
    np.testing.assert_allclose(
        result.cell_force_N.sum(axis=0), [15, 57], atol=1e-11, rtol=0
    )


@pytest.mark.parametrize("rotation", [False, True])
def test_rigid_motion_has_zero_stress(quad_mesh, rotation):
    mesh = quad_mesh(shear=0.3)
    fv = QuadFVGeometry(mesh, width_m=0.1)
    x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]

    def u(x):
        return (
            np.column_stack((-x[:, 1], x[:, 0]))
            if rotation
            else np.tile([2, 3], (len(x), 1))
        )

    r = fv.stress_flux(u(x), u(b), viscosity_Pa_s=1, boundary_viscosity_Pa_s=1)
    np.testing.assert_allclose(r.stress_Pa, 0, atol=1e-11, rtol=0)


def test_common_flux_conservation(quad_mesh):
    mesh = quad_mesh(2, 1)
    fv = QuadFVGeometry(mesh, width_m=2)
    f = int(np.flatnonzero(mesh.face_neighbours >= 0)[0])
    flux = np.zeros(len(mesh.faces))
    flux[f] = 7
    np.testing.assert_array_equal(fv.cell_balance(flux), [7, -7])
    flux[mesh.boundary_face_indices] = np.arange(6)
    assert fv.cell_balance(flux).sum() == pytest.approx(15, abs=1e-14)


def test_compact_pressure_nullspace_checkerboard_and_outlet(quad_mesh):
    mesh = quad_mesh()
    fv = QuadFVGeometry(mesh, width_m=0.2)
    matrix = fv.pressure_correction_matrix(momentum_diagonal_kg_s=2, density_kg_m3=3)
    eigen = np.linalg.eigvalsh(matrix.toarray())
    assert abs(eigen[0]) < 1e-13
    assert eigen[1] > 1e-4
    np.testing.assert_allclose(matrix @ np.ones(16), 0, atol=1e-13, rtol=0)
    checker = np.array([(-1.0) ** (i + j) for j in range(4) for i in range(4)])
    # rho*A*(V/a)/d = 3*.05*.00625/.25 = .00375 per edge.
    # Checkerboard jumps twice its local sign across each adjacent edge.
    degree = np.array(
        [2 + int(0 < i < 3) + int(0 < j < 3) for j in range(4) for i in range(4)]
    )
    np.testing.assert_allclose(
        matrix @ checker, 0.0075 * degree * checker, atol=1e-13, rtol=0
    )
    outlet = fv.pressure_correction_matrix(
        momentum_diagonal_kg_s=2,
        density_kg_m3=3,
        pressure_boundary_indices=mesh.boundary_sets["right"],
    )
    assert np.linalg.eigvalsh(outlet.toarray())[0] > 1e-4


def test_affine_rhie_chow_preserves_velocity_and_pressure(quad_mesh):
    mesh = quad_mesh(shear=np.tan(np.pi / 9))
    fv = QuadFVGeometry(mesh, width_m=0.4)
    x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]

    def u(x):
        return np.column_stack((x[:, 1] + 1, -x[:, 0]))

    def p(x):
        return 2 + x[:, 0] - 3 * x[:, 1]

    flow = fv.rhie_chow_flux(
        u(x),
        u(b),
        p(x),
        p(b),
        momentum_diagonal_kg_s=1,
        density_kg_m3=2,
        pressure_boundary_indices=mesh.boundary_sets["right"],
    )
    expected = (
        2
        * 0.4
        * mesh.face_lengths_m
        * np.einsum("ij,ij->i", u(mesh.face_centres_m), mesh.face_normals)
    )
    np.testing.assert_allclose(flow, expected, atol=1e-11, rtol=0)
    np.testing.assert_allclose(fv.cell_balance(flow), 0, atol=1e-11, rtol=0)


def test_invalid_diagonal_and_resources(quad_mesh):
    mesh = quad_mesh()
    with pytest.raises(ValidationError, match="limit"):
        QuadFVGeometry(mesh, width_m=1, max_cells=2)
    fv = QuadFVGeometry(mesh, width_m=1)
    with pytest.raises(ValidationError, match="positive"):
        fv.pressure_correction_matrix(momentum_diagonal_kg_s=0, density_kg_m3=1)


def periodic_domain(mesh):
    return FluidDomain.from_document(
        {
            "fluid_domain_schema_version": 1,
            "mesh": "unused.json",
            "width_m": 0.3,
            "fluid": {"density_kg_m3": 2, "viscosity_Pa_s": 1},
            "pressure_reference_cell": 0,
            "boundaries": {
                "bottom": {"type": "wall"},
                "top": {"type": "wall"},
                "left": {
                    "type": "periodic",
                    "partner": "right",
                    "translation_m": [1, 0],
                    "pressure_drop_Pa": 1,
                },
                "right": {
                    "type": "periodic",
                    "partner": "left",
                    "translation_m": [-1, 0],
                    "pressure_drop_Pa": -1,
                },
            },
        },
        mesh=mesh,
    )


def test_periodic_shared_flux_stress_and_pressure_nullspace(quad_mesh):
    mesh = quad_mesh(shear=0.3)
    domain = periodic_domain(mesh)
    fv = QuadFVGeometry.from_domain(domain)
    x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]
    u = np.column_stack((1 + x[:, 1], np.zeros(len(x))))
    ub = np.column_stack((1 + b[:, 1], np.zeros(len(b))))
    p, pb = x[:, 1], b[:, 1]
    flow = fv.rhie_chow_flux(u, ub, p, pb, momentum_diagonal_kg_s=1, density_kg_m3=2)
    for f, q in fv.periodic_partners.items():
        assert flow[f] == -flow[q]
    np.testing.assert_allclose(fv.cell_balance(flow), 0, atol=1e-11, rtol=0)
    np.testing.assert_allclose(
        fv.gradient(p, pb), np.tile([0, 1], (16, 1)), atol=1e-11, rtol=0
    )
    force = fv.stress_flux(u, ub, viscosity_Pa_s=1, boundary_viscosity_Pa_s=1)
    np.testing.assert_allclose(force.cell_force_N, 0, atol=1e-11, rtol=0)
    matrix = fv.pressure_correction_matrix(momentum_diagonal_kg_s=1, density_kg_m3=2)
    e = np.linalg.eigvalsh(matrix.toarray())
    assert abs(e[0]) < 1e-13 and e[1] > 1e-4
    assert len(fv.active_faces) == len(mesh.faces) - 4
    with pytest.raises(ValueError):
        flow[next(iter(fv.periodic_partners))] = 0


@pytest.mark.parametrize(
    "pairs",
    [
        {0: 1},
        {0: 0},
        {True: 3, 3: True},
        {0: 99, 99: 0},
        {0: 1, 1: 0},
        {0.0: 3, 3: 0.0},
    ],
)
def test_invalid_periodic_pairing(quad_mesh, pairs):
    with pytest.raises(ValidationError, match="periodic"):
        QuadFVGeometry(quad_mesh(), width_m=1, periodic_pairs=pairs)


@pytest.mark.parametrize("width", [0, -1, True, float("inf"), 1e-323])
def test_invalid_width(quad_mesh, width):
    with pytest.raises(ValidationError):
        QuadFVGeometry(quad_mesh(), width_m=width)


@pytest.mark.parametrize(
    "values,boundary",
    [
        ([1], [1] * 16),
        ([1] * 16, [1]),
        (np.zeros((16, 0)), np.zeros((16, 0))),
        ([True] * 16, [1] * 16),
        ([float("nan")] * 16, [1] * 16),
        (np.ones((16, 2)), np.ones((16, 1))),
    ],
)
def test_invalid_reconstruction_fields(quad_mesh, values, boundary):
    fv = QuadFVGeometry(quad_mesh(), width_m=1)
    with pytest.raises(ValidationError):
        fv.gradient(values, boundary)


def test_operator_arrays_and_inputs_are_owned(quad_mesh):
    mesh = quad_mesh()
    fv = QuadFVGeometry(mesh, width_m=2)
    values, boundary = np.ones(16), np.ones(16)
    gradient = fv.gradient(values, boundary)
    values[:] = 0
    np.testing.assert_allclose(gradient, 0, atol=1e-14, rtol=0)
    for a in (
        fv.cell_volumes_m3,
        fv.face_areas_m2,
        fv.neighbours,
        fv.active_faces,
        gradient,
    ):
        with pytest.raises(ValueError):
            a.setflags(write=True)
    matrix = fv.pressure_correction_matrix(momentum_diagonal_kg_s=1, density_kg_m3=1)
    matrix.data[:] = 0
    other = fv.pressure_correction_matrix(momentum_diagonal_kg_s=1, density_kg_m3=1)
    assert np.any(other.data != 0)


@pytest.mark.parametrize(
    "indices", [np.empty((0, 2), dtype=int), [1, 1], [True], [1.0], [-1], [100]]
)
def test_invalid_pressure_boundary_rows(quad_mesh, indices):
    fv = QuadFVGeometry(quad_mesh(), width_m=1)
    with pytest.raises(ValidationError):
        fv.pressure_correction_matrix(
            momentum_diagonal_kg_s=1, density_kg_m3=1, pressure_boundary_indices=indices
        )


def test_periodic_face_cannot_be_pressure_boundary_or_unpaired_flux(quad_mesh):
    mesh = quad_mesh()
    fv = QuadFVGeometry.from_domain(periodic_domain(mesh))
    with pytest.raises(ValidationError, match="periodic"):
        fv.pressure_correction_matrix(
            momentum_diagonal_kg_s=1,
            density_kg_m3=1,
            pressure_boundary_indices=mesh.boundary_sets["left"],
        )
    flux = np.zeros(len(mesh.faces))
    flux[next(iter(fv.periodic_partners))] = 1
    with pytest.raises(ValidationError, match="opposite"):
        fv.cell_balance(flux)


def test_checkerboard_is_visible_in_actual_rhie_chow_flux(quad_mesh):
    mesh = quad_mesh(6, 6)
    fv = QuadFVGeometry(mesh, width_m=1)
    checker = np.array([(-1.0) ** (i + j) for j in range(6) for i in range(6)])
    boundary = checker[mesh.face_owners[mesh.boundary_face_indices]]
    flow = fv.rhie_chow_flux(
        np.zeros((36, 2)),
        np.zeros((24, 2)),
        checker,
        boundary,
        momentum_diagonal_kg_s=1,
        density_kg_m3=1,
    )
    # Interior checkerboard cell gradients vanish, so compact jump is exposed.
    inside = (mesh.face_neighbours >= 0) & np.all(
        (mesh.face_centres_m > 1 / 6) & (mesh.face_centres_m < 5 / 6), axis=1
    )
    assert np.max(np.abs(flow[inside])) == pytest.approx(1 / 18, abs=1e-13)
    assert np.linalg.norm(fv.cell_balance(flow)) > 0.1


def test_independent_nonorthogonal_split(quad_mesh):
    mesh = quad_mesh(3, 3, shear=np.tan(np.pi / 9))
    fv = QuadFVGeometry(mesh, width_m=1)
    x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]
    fg = fv.face_gradient(x[:, 0], b[:, 0])
    # On a horizontal face, centres differ by (shear*h,h), normal=(0,1).
    # Compact normal slope is shear, cancelled by the deferred tangential part.
    faces = (mesh.face_neighbours >= 0) & (mesh.face_normals[:, 1] > 0.99)
    np.testing.assert_allclose(
        fg.compact_normal[faces], np.tan(np.pi / 9), atol=1e-13, rtol=0
    )
    np.testing.assert_allclose(
        fg.nonorthogonal_normal[faces], -np.tan(np.pi / 9), atol=1e-13, rtol=0
    )


@pytest.mark.parametrize(
    "diagonal,density",
    [
        (np.ones(2), 1),
        (-1, 1),
        (True, 1),
        (1, 0),
        (1, True),
        (1, float("inf")),
        (1e-323, 1),
    ],
)
def test_invalid_pressure_mobility(quad_mesh, diagonal, density):
    fv = QuadFVGeometry(quad_mesh(), width_m=1)
    with np.errstate(over="ignore"):
        with pytest.raises(ValidationError):
            fv.pressure_correction_matrix(
                momentum_diagonal_kg_s=diagonal, density_kg_m3=density
            )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"viscosity_Pa_s": 0},
        {"boundary_viscosity_Pa_s": -1},
        {"pressure_Pa": np.ones(16)},
        {"pressure_Pa": np.ones((16, 2)), "boundary_pressure_Pa": np.ones((16, 2))},
    ],
)
def test_invalid_stress_data(quad_mesh, kwargs):
    fv = QuadFVGeometry(quad_mesh(), width_m=1)
    args = {"viscosity_Pa_s": 1, "boundary_viscosity_Pa_s": 1}
    args.update(kwargs)
    with pytest.raises(ValidationError):
        fv.stress_flux(np.zeros((16, 2)), np.zeros((16, 2)), **args)


def test_skewed_faces_preserve_affine_values(quad_mesh):
    base = quad_mesh(5, 5)
    points = base.points.copy()
    interior = np.all((points > 0) & (points < 1), axis=1)
    points[interior] += np.random.default_rng(229).uniform(
        -0.008, 0.008, (interior.sum(), 2)
    )
    mesh = FluidMesh2D(points, base.quad4, base.line2, base.boundary_sets)
    assert mesh.max_skewness > 0.001
    fv = QuadFVGeometry(mesh, width_m=1)
    x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]
    value = 1 + 2 * x[:, 0] - 3 * x[:, 1]
    boundary = 1 + 2 * b[:, 0] - 3 * b[:, 1]
    np.testing.assert_allclose(
        fv.gradient(value, boundary), np.tile([2, -3], (25, 1)), atol=1e-11, rtol=0
    )
    expected = 1 + 2 * mesh.face_centres_m[:, 0] - 3 * mesh.face_centres_m[:, 1]
    np.testing.assert_allclose(
        fv.face_values(value, boundary), expected, atol=1e-11, rtol=0
    )


def test_pressure_mobility_uses_same_unrelaxed_diagonal(quad_mesh):
    mesh = quad_mesh(2, 1)
    fv = QuadFVGeometry(mesh, width_m=2)
    # V=(1,1), a=(2,4), d=(.5,.25), df=.375, rho=3, A=2, dn=.5.
    matrix = fv.pressure_correction_matrix(
        momentum_diagonal_kg_s=[2, 4], density_kg_m3=3
    )
    np.testing.assert_allclose(
        matrix.toarray(), [[4.5, -4.5], [-4.5, 4.5]], atol=1e-13, rtol=0
    )
    p = np.array([1.0, -1.0])
    pb = p[mesh.face_owners[mesh.boundary_face_indices]]
    flux = fv.rhie_chow_flux(
        np.zeros((2, 2)),
        np.zeros((6, 2)),
        p,
        pb,
        momentum_diagonal_kg_s=[2, 4],
        density_kg_m3=3,
    )
    face = int(np.flatnonzero(mesh.face_neighbours >= 0)[0])
    assert flux[face] == pytest.approx(4.5, abs=1e-13)


@pytest.mark.parametrize(
    "which", ["mesh", "domain", "pairs", "stress_velocity", "rc_pressure"]
)
def test_invalid_object_and_components(quad_mesh, which):
    mesh = quad_mesh()
    with pytest.raises(ValidationError):
        if which == "mesh":
            QuadFVGeometry(None, width_m=1)
        elif which == "domain":
            QuadFVGeometry.from_domain(None)
        elif which == "pairs":
            QuadFVGeometry(mesh, width_m=1, periodic_pairs=[])
        elif which == "stress_velocity":
            QuadFVGeometry(mesh, width_m=1).stress_flux(
                np.ones(16), np.ones(16), viscosity_Pa_s=1, boundary_viscosity_Pa_s=1
            )
        else:
            QuadFVGeometry(mesh, width_m=1).rhie_chow_flux(
                np.ones((16, 2)),
                np.ones((16, 2)),
                np.ones((16, 2)),
                np.ones((16, 2)),
                momentum_diagonal_kg_s=1,
                density_kg_m3=1,
            )
