"""Independent MKM model comparison; canonical external data stay outside Git."""

import os
from pathlib import Path

import numpy as np
import pytest

from tests.convergence.test_rans_channel import run_case
from tests.reference.dns_channel import (
    compare_dns_profile,
    load_dns_profile,
    numerical_channel_profile,
)

pytestmark = pytest.mark.cfd_convergence


def test_rans_mkm_model_budget_and_last_mesh_change(record_property):
    value = os.environ.get("M12_DNS_PROFILE")
    if value is None:
        pytest.skip(
            "External DNS profile not supplied; CI sa-dns job downloads and verifies canonical input"
        )
    profile = load_dns_profile(Path(value))
    record_property("dns_sha256", profile.sha256)
    record_property("dns_re_tau", profile.re_tau)
    predictions = []
    bulks = []
    for n in (128, 256):
        print(f"DNS model comparison: starting {n} cells/half", flush=True)
        domain, reference, result, _ = run_case(n)
        mesh = domain.mesh
        rows = np.concatenate(
            [
                mesh.boundary_sets[name]
                for name, bc in domain.boundaries.items()
                if bc["type"] == "wall"
            ]
        )
        bottom = np.isin(rows, mesh.boundary_sets["bottom"])
        normals = mesh.boundary_normals[rows]
        traction = result.flow.wall_traction_Pa
        tangential = (
            traction - np.einsum("ij,ij->i", traction, normals)[:, None] * normals
        )
        areas = mesh.face_lengths_m[mesh.boundary_face_indices[rows]] * domain.width_m
        shear = float(
            np.linalg.norm(tangential[bottom], axis=1)
            @ areas[bottom]
            / areas[bottom].sum()
        )
        inputs = dict(
            cell_y_m=mesh.cell_centres_m[:, 1],
            velocity_m_s=result.flow.velocity_m_s[:, 0],
            cell_volumes_m3=mesh.cell_areas_m2 * domain.width_m,
            half_height_m=1,
            wall_shear_Pa=shear,
            viscosity_Pa_s=domain.fluid.viscosity_Pa_s,
        )
        metrics = compare_dns_profile(
            profile, friction_velocity_m_s=reference.friction_velocity_m_s, **inputs
        )
        assert metrics["relative_rms_velocity_plus"] <= 0.05
        assert metrics["relative_bulk_velocity"] <= 0.05
        interpolant = numerical_channel_profile(**inputs)
        predictions.append(
            interpolant(profile.y_over_delta) / reference.friction_velocity_m_s
        )
        bulks.append(metrics["bulk_velocity_m_s"])
        for key, metric in metrics.items():
            record_property(f"n{n}_{key}", metric)
        print(f"DNS model comparison {n}: {metrics}", flush=True)
    change = float(
        np.linalg.norm(predictions[1] - predictions[0]) / np.linalg.norm(predictions[1])
    )
    assert change < 0.01 and change <= 0.05 / 5
    assert abs(bulks[1] / bulks[0] - 1) < 0.01
    record_property("last_mesh_relative_profile_change", change)
    print(f"DNS last mesh change: {change}", flush=True)
