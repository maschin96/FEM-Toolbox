"""Unchanged Doc/70/71 accuracy gates; independent analytic point sources."""

import numpy as np
import pytest

from fem_toolbox.cfd import FluidDomain, solve_quad_flow

pytestmark = pytest.mark.cfd_convergence


def problem(mesh, *, rho=1, periodic=False, closed=False):
    boundaries = {"bottom": {"type": "wall"}, "top": {"type": "wall"}}
    if periodic:
        boundaries.update(
            left={
                "type": "periodic",
                "partner": "right",
                "translation_m": [1, 0],
                "pressure_drop_Pa": 2,
            },
            right={
                "type": "periodic",
                "partner": "left",
                "translation_m": [-1, 0],
                "pressure_drop_Pa": -2,
            },
        )
    elif closed:
        boundaries.update(left={"type": "wall"}, right={"type": "wall"})
    else:
        x = mesh.face_centres_m[mesh.boundary_face_indices[mesh.boundary_sets["left"]]]
        boundaries.update(
            left={
                "type": "inlet",
                "velocity_m_s": np.column_stack(
                    (x[:, 1] * (1 - x[:, 1]), np.zeros(len(x)))
                ).tolist(),
                "nu_tilde_m2_s": 0,
            },
            right={"type": "pressure_outlet", "pressure_Pa": 0},
        )
    doc = dict(
        fluid_domain_schema_version=1,
        mesh="unused.json",
        width_m=1,
        fluid={"density_kg_m3": rho, "viscosity_Pa_s": 1},
        boundaries=boundaries,
    )
    if periodic or closed:
        doc["pressure_reference_cell"] = 0
    return FluidDomain.from_document(doc, mesh=mesh)


def exact_mms(points, density):
    # psi=sin²(pi*x)sin²(pi*y), U=(psi_y,-psi_x), div U identically zero.
    x, y = points.T
    pi = np.pi
    cx, cy, sx, sy = (
        np.cos(2 * pi * x),
        np.cos(2 * pi * y),
        np.sin(2 * pi * x),
        np.sin(2 * pi * y),
    )
    u, v = pi / 2 * (1 - cx) * sy, -pi / 2 * sx * (1 - cy)
    ux, uy = pi**2 * sx * sy, pi**2 * (1 - cx) * cy
    vx, vy = -(pi**2) * cx * (1 - cy), -(pi**2) * sx * sy
    fx = -2 * pi**3 * (2 * cx - 1) * sy + 2 * pi * cx * sy + density * (u * ux + v * uy)
    fy = -2 * pi**3 * sx * (1 - 2 * cy) + 2 * pi * sx * cy + density * (u * vx + v * vy)
    return np.column_stack((u, v)), sx * sy, np.column_stack((fx, fy))


def test_quad_poiseuille_four_meshes(quad_mesh, record_property):
    errors, pressure_errors = [], []
    for n in (8, 16, 32, 64):
        print(f"Poiseuille: starting {2 * n}x{n} cells", flush=True)
        mesh = quad_mesh(2 * n, n, size_m=(2, 1))
        r = solve_quad_flow(problem(mesh), velocity_scale_m_s=1, length_scale_m=2)
        y = mesh.cell_centres_m[:, 1]
        errors.append(float(np.max(np.abs(r.velocity_m_s[:, 0] - y * (1 - y)))))
        p = r.pressure_Pa.reshape(n, 2 * n).mean(axis=0)
        pressure_errors.append(float(np.sqrt(np.mean((np.diff(p) * n + 2) ** 2))))
        print(
            f"Poiseuille {2 * n}x{n}: U error={errors[-1]:.6g}, "
            f"pressure-gradient error={pressure_errors[-1]:.6g}, "
            f"outer={r.iterations}, SIMPLE={r.simple_steps}",
            flush=True,
        )
        assert abs(r.global_mass_balance_kg_s) < 1e-9 * r.inlet_mass_flow_kg_s
    record_property("velocity_errors", errors)
    record_property("pressure_gradient_errors", pressure_errors)
    assert errors[-1] < 0.001, errors
    assert pressure_errors[-1] < 0.02, pressure_errors
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 0.8), errors


@pytest.mark.parametrize(
    "scheme,rho,limit,order",
    [("central", 0.1, 0.01, 1.8), ("linear-upwind", 1.0, 0.1, 0.7)],
)
def test_quad_nonlinear_mms_four_meshes(
    quad_mesh, record_property, scheme, rho, limit, order
):
    errors, p_errors = [], []
    for n in (16, 32, 64, 128):
        print(f"MMS {scheme}: starting {n}x{n} cells", flush=True)
        mesh = quad_mesh(n, n)
        exact, p, force = exact_mms(mesh.cell_centres_m, rho)
        r = solve_quad_flow(
            problem(mesh, rho=rho, closed=True),
            velocity_scale_m_s=1,
            length_scale_m=1,
            force_N_m3=force,
            convection_scheme=scheme,
        )
        errors.append(
            float(np.sqrt(np.mean(np.sum((r.velocity_m_s - exact) ** 2, axis=1))))
        )
        p_errors.append(float(np.sqrt(np.mean((r.pressure_Pa - p + p[0]) ** 2))))
        print(
            f"MMS {scheme} {n}x{n}: U error={errors[-1]:.6g}, "
            f"p error={p_errors[-1]:.6g}, outer={r.iterations}, "
            f"SIMPLE={r.simple_steps}",
            flush=True,
        )
        assert np.max(np.abs(r.mass_residual_kg_m3_s)) < 1e-8 * rho
    record_property("velocity_errors", errors)
    record_property("pressure_errors", p_errors)
    assert errors[-1] < limit, errors
    assert p_errors[-1] < limit, p_errors
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > order), errors
    assert np.all(np.log2(np.array(p_errors[:-1]) / p_errors[1:]) > order), p_errors


@pytest.mark.parametrize("warp", [0.0, 0.08])
def test_quad_variable_mu_four_meshes(quad_mesh, record_property, warp):
    errors = []
    for n in (16, 32, 64, 128):
        print(f"Variable mu (warp={warp}): starting {n}x{n} cells", flush=True)
        mesh = quad_mesh(n, n, warp=warp)
        x, b = mesh.cell_centres_m, mesh.face_centres_m[mesh.boundary_face_indices]
        force = np.column_stack(
            (-1 + 4 * x[:, 1], np.zeros(len(x)))
        )  # + periodic driving (2,0).
        r = solve_quad_flow(
            problem(mesh, periodic=True),
            velocity_scale_m_s=1,
            length_scale_m=1,
            force_N_m3=force,
            viscosity_Pa_s=1 + x[:, 1],
            boundary_viscosity_Pa_s=1 + b[:, 1],
            convection_scheme="central",
            velocity_relaxation=1,
            pressure_relaxation=1,
        )
        errors.append(
            float(
                np.sqrt(np.mean((r.velocity_m_s[:, 0] - x[:, 1] * (1 - x[:, 1])) ** 2))
            )
        )
        print(
            f"Variable mu {n}x{n}: U error={errors[-1]:.6g}, "
            f"outer={r.iterations}, SIMPLE={r.simple_steps}",
            flush=True,
        )
    record_property("velocity_errors", errors)
    assert errors[-1] < 0.001, errors
    if warp == 0:
        # Independently checked polynomial-exact full-stress patch: refine
        # below iteration/roundoff error, where a log-ratio is undefined.
        assert max(errors) < 1e-9, errors
    else:
        # Same physical square, source and four mesh sizes. Smoothly distorted
        # quads expose measurable face-quadrature truncation, not roundoff.
        assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 1.8), errors
