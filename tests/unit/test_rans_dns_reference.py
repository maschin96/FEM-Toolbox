"""Offline source-integrity and interpolation tests; no test downloads data."""

import hashlib

import numpy as np
import pytest

from tests.reference.dns_channel import compare_dns_profile, load_dns_profile


def profile_file(tmp_path, *, re_header="587.19", normalization="U_tau, h"):
    text = f"# Re_tau = {re_header}\n# Normalization: {normalization}\n"
    text += "0 0 0 0 0 0 0\n0.5 293.595 0.375 0 0 0 0\n1 587.19 0.5 0 0 0 0\n"
    path = tmp_path / "synthetic.means"
    path.write_text(text)
    source = {
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "rows": 3,
        "re_tau": 587.19,
        "normalization": "U_tau, h",
    }
    return path, source


def test_source_checksum_units_and_immutable_loaded_data(tmp_path):
    path, source = profile_file(tmp_path)
    data = load_dns_profile(path, source=source)
    assert data.re_tau == 587.19
    assert data.sha256 == source["sha256"]
    np.testing.assert_array_equal(data.y_over_delta, [0, 0.5, 1])
    assert data.bulk_velocity_plus == pytest.approx(1 / 3)
    with pytest.raises(ValueError):
        data.velocity_plus.setflags(write=True)
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError, match="SHA"):
        load_dns_profile(path, source=source)


@pytest.mark.parametrize("change", [{"re_header": "590"}, {"normalization": "m/s, m"}])
def test_wrong_reynolds_or_normalization_is_not_silently_reinterpreted(
    tmp_path, change
):
    path, source = profile_file(tmp_path, **change)
    with pytest.raises(ValueError):
        load_dns_profile(path, source=source)


def test_numerical_interpolation_preserves_parabolic_wall_and_center_conditions(
    tmp_path,
):
    path, source = profile_file(tmp_path)
    data = load_dns_profile(path, source=source)
    y = (np.arange(6) + 0.5) / 3
    velocity = 0.5 * y * (2 - y)
    metrics = compare_dns_profile(
        data,
        cell_y_m=y,
        velocity_m_s=velocity,
        cell_volumes_m3=np.full(6, 1 / 3),
        half_height_m=1,
        friction_velocity_m_s=1,
        wall_shear_Pa=1,
        viscosity_Pa_s=1,
    )
    assert metrics["relative_rms_velocity_plus"] < 1e-14
    assert metrics["bulk_velocity_m_s"] == pytest.approx(float(velocity.mean()))
    assert metrics["relative_bulk_velocity"] > 0


def test_coordinate_roundoff_between_streamwise_cells_does_not_create_extra_rows(
    tmp_path,
):
    path, source = profile_file(tmp_path)
    data = load_dns_profile(path, source=source)
    base = (np.arange(6) + 0.5) / 3
    columns = np.tile([-1, 0, 1], len(base))
    y = np.repeat(base, 3) + columns * np.finfo(float).eps
    velocity = np.repeat(0.5 * base * (2 - base), 3) + columns * 1e-12
    metrics = compare_dns_profile(
        data,
        cell_y_m=y,
        velocity_m_s=velocity,
        cell_volumes_m3=np.full(len(y), 1 / 9),
        half_height_m=1,
        friction_velocity_m_s=1,
        wall_shear_Pa=1,
        viscosity_Pa_s=1,
    )
    assert metrics["relative_rms_velocity_plus"] < 1e-10
