"""Elastically guided lid: bidirectional stationary FSI verification (#92).

Parameters, references and tolerances come from Doc/71 (#92 row) and were
fixed before any coupled result was computed:

* laminar: Q=1e-3 m^3/s, H0=0.1 m, L=b=1 m, mu=rho=1, k_s=1e4 N/m,
  p_external=p_out=0. Root of d(H0+d)^3=6 mu L^2 Q/k_s; finest lid displacement
  and pressure drop within 1 %, last mesh change <0.2 %.
* SA-noft2: nu=(0.05 m)(0.01 m/s)/587.19, k_s=1 N/m, Q from the independent
  1D SA BVP at H0 with f_x=0.002 N/m^3; the BVP gives Delta p at the final
  height. Within 2 %, last mesh change <0.4 %, y+<=1.

Scales are chosen before the run and documented in Doc/93: Lref=H0 so the
Doc/71 deformation and coupling gates refer to the channel height. Laminar
Uref=1 m/s makes rho Uref^2 = 1 Pa commensurate with the 12 Pa Stokes
pressure load (rho Ubulk^2=1e-4 Pa would place the load gate at 1e-12 of the
load). SA Uref=0.2 m/s=20 u_tau as in the established small SI RANS check.
"""

from itertools import pairwise

import numpy as np
import pytest

from fem_toolbox.coupling import StationaryFsiProblem, solve_stationary_fsi
from tests.reference.fsi_lid import (
    LaminarLid,
    SaLid,
    graded_wall_nodes,
    lid_channel,
    lid_structure,
)

pytestmark = pytest.mark.cfd_convergence

LAMINAR_LEVELS = ((10, 4, 5), (20, 8, 10), (40, 16, 20), (80, 32, 40))
SA_LEVELS = ((32, 4), (48, 6), (64, 8), (96, 12))


def result_flow(result):
    return getattr(result.fluid, "flow", result.fluid)


def pressure_drop(result, length_m=1.0):
    """Least-squares linear p(x) over all cells, Delta p = p(0) - p(L)."""
    x = result.mesh_motion.domain.mesh.cell_centres_m[:, 0]
    pressure = np.asarray(result_flow(result).pressure_Pa)
    return -np.polyfit(x, pressure, 1)[0] * length_m


def audit_and_measure(result, reference, height_m):
    """Shared balance gates; returns lid displacement and pressure-drop data."""
    wet = np.unique(result.interface._triangles)
    lid = np.asarray(result.structure_displacement_m)[wet]
    geometry = np.asarray(result.geometry_displacement_m)[wet, 1].mean()
    flow = result_flow(result)
    domain = result.mesh_motion.domain
    # Doc/71 transfer gates on every coupling step.
    for row in result.history:
        assert row.transfer_relative_error <= 1e-10
        assert row.maximum_displacement_ratio <= 0.02
        assert row.maximum_displacement_gradient <= 0.02
        assert row.minimum_area_ratio >= 0.5
    # Interface virtual work with the final structural displacement.
    lookup = {int(r): i for i, r in enumerate(domain.wall_quality.boundary_indices)}
    rows = [lookup[int(r)] for r in domain.mesh.boundary_sets["top"]]
    traction = np.asarray(flow.wall_traction_Pa)[rows]
    mapping = result.interface
    forces = np.zeros_like(mapping.quadrature_points_m)
    forces[:, :2] = (
        -traction[mapping.segment_indices] * mapping.quadrature_areas_m2[:, None]
    )
    work = mapping.audit(forces, result.structure_displacement_m)
    assert work.maximum_relative_error <= 1e-10
    # Structure equilibrium and guided reactions balance the lid load.
    final = result.history[-1]
    # The 1e9 plate/spring modulus ratio needed for a rigid lid limits the
    # reaction reconstruction K u - f to ~1e-7 of the load; the solver's own
    # backward-error gate (<=1e-10) is unchanged and enforced in the loop.
    load = abs(final.fluid_force_N[1])
    assert np.max(np.abs(result.structure.force_balance)) <= 1e-6 * load
    assert result.structure.relative_residual <= 1e-10
    np.testing.assert_allclose(
        final.structure_reaction_N[1], -final.fluid_force_N[1], rtol=1e-6
    )
    # Fluid balances: fixed inlet flow leaves through the outlet.
    q, settings = reference.flow_m3_s, result.settings
    assert abs(flow.inlet_mass_flow_kg_s) == pytest.approx(q, rel=1e-9)
    mass_gate = 1e-10 * settings["velocity_scale_m_s"] * settings["reference_area_m2"]
    assert abs(flow.global_mass_balance_kg_s) <= mass_gate
    assert abs(abs(flow.outlet_mass_flow_kg_s) - q) <= mass_gate
    # Lid: nearly rigid (span-invariant gate already enforced by the map).
    assert np.ptp(lid[:, 1]) <= 1e-6 * abs(lid[:, 1].mean())
    assert np.max(np.abs(lid[:, [0, 2]])) == 0
    return {
        "lid_displacement_m": float(lid[:, 1].mean()),
        "geometry_height_m": float(height_m + geometry),
        "pressure_drop_Pa": float(pressure_drop(result)),
        "one_way_displacement_m": float(result.history[0].displacement_residual_m),
        "iterations": result.iterations,
        "transfer_relative_error": max(
            r.transfer_relative_error for r in result.history
        ),
        "virtual_work_relative_error": work.maximum_relative_error,
        "minimum_area_ratio": result.mesh_motion.minimum_area_ratio,
        "final_displacement_residual_m": final.displacement_residual_m,
        "final_load_change_N": final.load_change_N,
    }


def check_sequence(levels, metrics, reference, *, tolerance, change, record_property):
    d_ref = reference.coupled_displacement_m
    d1_ref = reference.one_way_displacement_m
    for level, m in zip(levels, metrics, strict=True):
        m["displacement_error"] = m["lid_displacement_m"] / d_ref - 1
        m["pressure_drop_reference_Pa"] = reference.pressure_drop_Pa(
            m["geometry_height_m"]
        )
        m["pressure_drop_error"] = (
            m["pressure_drop_Pa"] / m["pressure_drop_reference_Pa"] - 1
        )
        m["one_way_error"] = m["one_way_displacement_m"] / d1_ref - 1
        for name, value in m.items():
            record_property(f"{'x'.join(map(str, level))}_{name}", value)
        print(f"FSI lid {level}: {m}", flush=True)
    finest, previous = metrics[-1], metrics[-2]
    assert abs(finest["displacement_error"]) < tolerance
    assert abs(finest["pressure_drop_error"]) < tolerance
    for key in ("lid_displacement_m", "pressure_drop_Pa"):
        assert abs(finest[key] / previous[key] - 1) < change
    # Discretization error decreases monotonically with refinement.
    errors = [abs(m["displacement_error"]) for m in metrics]
    assert all(b < a for a, b in pairwise(errors))
    # One-way and bidirectional effect: the coupled lid opens less than one-way.
    analytic_effect = 1 - d_ref / d1_ref
    numeric_effect = 1 - finest["lid_displacement_m"] / finest["one_way_displacement_m"]
    record_property("analytic_bidirectional_effect", analytic_effect)
    record_property("numeric_bidirectional_effect", numeric_effect)
    print(f"bidirectional effect analytic={analytic_effect} numeric={numeric_effect}")
    assert abs(finest["one_way_error"]) < tolerance
    assert analytic_effect > tolerance
    assert numeric_effect == pytest.approx(analytic_effect, rel=0.05)


def test_laminar_guided_lid_four_mesh_pairs(record_property):
    reference = LaminarLid()
    metrics = []
    for nx, ny, sx in LAMINAR_LEVELS:
        print(f"laminar lid: starting fluid {nx}x{ny}, structure {sx}", flush=True)
        domain = lid_channel(
            nx=nx,
            y_nodes_m=np.linspace(0, reference.height_m, ny + 1),
            inlet_profile=reference.inlet_profile,
        )
        model, wet = lid_structure(
            nx=sx,
            height_m=reference.height_m,
            spring_stiffness_N_m=reference.spring_stiffness_N_m,
        )
        result = solve_stationary_fsi(
            StationaryFsiProblem(
                fluid_domain=domain,
                wall_name="top",
                structure=model,
                facet_set=wet,
                length_scale_m=reference.height_m,
                velocity_scale_m_s=1.0,
                relaxation=1.0,
            )
        )
        assert result.converged
        metrics.append(audit_and_measure(result, reference, reference.height_m))
    check_sequence(
        LAMINAR_LEVELS,
        metrics,
        reference,
        tolerance=0.01,
        change=0.002,
        record_property=record_property,
    )


def test_sa_guided_lid_four_mesh_pairs(record_property):
    reference = SaLid()
    base = reference.base
    nu = reference.kinematic_viscosity_m2_s
    record_property("reference_flow_m3_s", reference.flow_m3_s)
    record_property("reference_re_tau", base.re_tau)
    print(f"SA lid reference Q={reference.flow_m3_s}, Re_tau={base.re_tau}")
    metrics = []
    for n, nx in SA_LEVELS:
        print(
            f"SA lid: starting {n} cells/half, {nx} along, structure {nx}", flush=True
        )
        domain = lid_channel(
            nx=nx,
            y_nodes_m=graded_wall_nodes(n, reference.height_m),
            inlet_profile=reference.inlet_profile,
            viscosity_Pa_s=reference.density_kg_m3 * nu,
        )
        model, wet = lid_structure(
            nx=nx,
            height_m=reference.height_m,
            spring_stiffness_N_m=reference.spring_stiffness_N_m,
        )
        # Initial guess only: developed BVP nu_tilde, positive in every cell.
        start = np.maximum(base.evaluate(domain.mesh.cell_centres_m[:, 1])[1], 3 * nu)
        result = solve_stationary_fsi(
            StationaryFsiProblem(
                fluid_domain=domain,
                wall_name="top",
                structure=model,
                facet_set=wet,
                length_scale_m=reference.height_m,
                velocity_scale_m_s=0.2,
                fluid_model="sa-noft2",
                initial_nu_tilde_m2_s=start,
                fluid_options={"sa_relaxation": 0.1},
                relaxation=1.0,
            )
        )
        assert result.converged
        assert np.max(result.fluid.wall_y_plus) <= 1
        assert np.all(np.asarray(result.fluid.nu_tilde_m2_s) >= 0)
        m = audit_and_measure(result, reference, reference.height_m)
        m["maximum_y_plus"] = float(np.max(result.fluid.wall_y_plus))
        metrics.append(m)
    # Independent BVP at the final height confirms Delta p ~ H^-3 at fixed Q.
    height = metrics[-1]["geometry_height_m"]
    similarity = reference.pressure_drop_Pa(height) / (
        reference.pressure_drive_N_m3 * (reference.height_m / height) ** 3
    )
    record_property("bvp_similarity_ratio", similarity)
    assert similarity == pytest.approx(1, abs=1e-6)
    check_sequence(
        SA_LEVELS,
        metrics,
        reference,
        tolerance=0.02,
        change=0.004,
        record_property=record_property,
    )
