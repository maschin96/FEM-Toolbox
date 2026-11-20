"""Independent invariants of the conservative flow reconstruction."""

import numpy as np

from fem_toolbox.cfd import FluidDomain, solve_quad_flow
from fem_toolbox.cfd.quad_system import QuadSystem


def closed_domain(mesh):
    return FluidDomain.from_document(
        {
            "fluid_domain_schema_version": 1,
            "mesh": "unused.json",
            "width_m": 0.2,
            "fluid": {"density_kg_m3": 1, "viscosity_Pa_s": 1},
            "boundaries": {name: {"type": "wall"} for name in mesh.boundary_sets},
            "pressure_reference_cell": 0,
        },
        mesh=mesh,
    )


def test_zero_flow_and_all_result_arrays_are_immutable(quad_mesh):
    result = solve_quad_flow(
        closed_domain(quad_mesh(4, 4)), velocity_scale_m_s=1, length_scale_m=1
    )
    np.testing.assert_allclose(result.velocity_m_s, 0, atol=1e-14, rtol=0)
    np.testing.assert_allclose(result.pressure_Pa, 0, atol=1e-14, rtol=0)
    assert result.power_balance_W == 0
    for name in (
        "velocity_m_s",
        "pressure_Pa",
        "mass_flux_kg_s",
        "face_velocity_m_s",
        "mass_residual_kg_m3_s",
        "momentum_residual_N_m3",
        "nonorthogonal_correction_residual_N_m3",
        "wall_traction_Pa",
        "global_momentum_balance_N",
    ):
        assert not getattr(result, name).flags.writeable


def test_linear_upwind_face_values_are_bounded_by_neighbour_cells(quad_mesh):
    mesh = quad_mesh(6, 6)
    system = QuadSystem(closed_domain(mesh), np.ones(36), np.ones(24))
    random = np.random.default_rng(230)
    velocity = random.normal(size=(36, 2))
    flux = random.normal(size=len(mesh.faces))
    _, transported = system.transport(flux, velocity.ravel(), "linear-upwind")
    values = transported.reshape(-1, 2)
    faces = system.fv.active_faces[system.fv.neighbours[system.fv.active_faces] >= 0]
    owner = velocity[mesh.face_owners[faces]]
    neighbour = velocity[system.fv.neighbours[faces]]
    assert np.all(values[faces] >= np.minimum(owner, neighbour) - 1e-14)
    assert np.all(values[faces] <= np.maximum(owner, neighbour) + 1e-14)


def test_gauged_coupled_operator_has_no_extra_pressure_nullmode(quad_mesh):
    from scipy.sparse import bmat

    system = QuadSystem(closed_domain(quad_mesh(4, 4)), np.ones(16), np.ones(16))
    diagonal = system.diffusion.matrix.diagonal().reshape(16, 2)
    pressure = system.pressure_flux(diagonal)
    np.testing.assert_allclose(
        system.pressure_momentum.matrix @ np.ones(16), 0, atol=1e-12, rtol=0
    )
    alternate = (-1.0) ** np.indices((4, 4)).sum(axis=0).ravel()
    assert np.linalg.norm(system.balance @ pressure.apply(alternate)) > 1e-6
    operator = bmat(
        [
            [system.diffusion.matrix, system.pressure_momentum.matrix],
            [system.mass_velocity.matrix, system.balance @ pressure.matrix],
        ],
        format="lil",
    )
    operator[32, :] = 0
    operator[32, 32] = 1
    # Independent dense singular-value check, including all velocity modes.
    values = np.linalg.svd(operator.toarray(), compute_uv=False)
    assert values[-1] / values[0] > 1e-10


def test_variable_viscosity_full_stress_matches_analytic_source(quad_mesh):
    for shear in (0, np.tan(np.deg2rad(20))):
        mesh = quad_mesh(4, 8, shear=shear)
        domain = FluidDomain.from_document(
            {
                "fluid_domain_schema_version": 1,
                "mesh": "unused.json",
                "width_m": 0.2,
                "fluid": {"density_kg_m3": 1, "viscosity_Pa_s": 1},
                "boundaries": {
                    "bottom": {"type": "wall"},
                    "top": {"type": "wall"},
                    "left": {
                        "type": "periodic",
                        "partner": "right",
                        "translation_m": [1, 0],
                        "pressure_drop_Pa": 0,
                    },
                    "right": {
                        "type": "periodic",
                        "partner": "left",
                        "translation_m": [-1, 0],
                        "pressure_drop_Pa": 0,
                    },
                },
                "pressure_reference_cell": 0,
            },
            mesh=mesh,
        )
        y = mesh.cell_centres_m[:, 1]
        boundary_y = mesh.face_centres_m[mesh.boundary_face_indices, 1]
        system = QuadSystem(domain, 1 + y, 1 + boundary_y)
        velocity = np.column_stack((y * (1 - y), np.zeros(len(y))))
        # -d/dy[(1+y)*(1-2y)] = 1+4y, not -mu*laplace(U).
        source = np.column_stack((1 + 4 * y, np.zeros(len(y))))
        np.testing.assert_allclose(
            system.diffusion.apply(velocity.ravel()).reshape(-1, 2)
            / system.fv.cell_volumes_m3[:, None],
            source,
            atol=1e-10,
            rtol=0,
        )
