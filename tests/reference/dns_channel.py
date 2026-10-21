"""Offline, checksum-bound MKM source reader and numerical-profile comparison."""

import hashlib
import io
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline


def _readonly(value):
    a = np.asarray(value, dtype=np.float64)
    return np.frombuffer(a.tobytes(), dtype=np.float64).reshape(a.shape)


@dataclass(frozen=True)
class DNSChannelProfile:
    y_over_delta: np.ndarray
    velocity_plus: np.ndarray
    re_tau: float
    sha256: str

    @property
    def bulk_velocity_plus(self):
        return float(simpson(self.velocity_plus, x=self.y_over_delta))


def load_dns_profile(path, *, source=None):
    """No network access or automatic source replacement; reject hash/unit drift."""
    if source is None:
        manifest = (
            Path(__file__).resolve().parents[2]
            / "Doc/evidence/m12/reference_sources.json"
        )
        source = json.loads(manifest.read_text())["mkm_channel_590"]
    path = Path(path)
    if path.stat().st_size > 1_000_000:
        raise ValueError("DNS profile exceeds input resource limit")
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != source["sha256"]:
        raise ValueError(f"DNS SHA-256 mismatch: {digest}")
    text = raw.decode("ascii")
    match = re.search(r"^#\s*Re_tau\s*=\s*([0-9.eE+-]+)\s*$", text, re.M)
    if match is None or not np.isclose(
        float(match[1]), source["re_tau"], rtol=0, atol=1e-10
    ):
        raise ValueError("DNS Reynolds number does not match declared reference")
    if source["normalization"] != "U_tau, h" or "Normalization: U_tau, h" not in text:
        raise ValueError("DNS normalization must be U_tau, h")
    values = np.loadtxt(io.BytesIO(raw), comments="#")
    if (
        values.shape != (source["rows"], 7)
        or source["rows"] < 3
        or not np.all(np.isfinite(values))
    ):
        raise ValueError("DNS profile needs the declared complete seven-column table")
    y, yp, u = values[:, :3].T
    if (
        np.any(np.diff(y) <= 0)
        or abs(y[0]) > 1e-10
        or abs(y[-1] - 1) > 1e-10
        or np.any(u < 0)
    ):
        raise ValueError("invalid normalized DNS coordinates/mean velocity")
    # The original text columns are rounded to a few decimal digits.
    if not np.allclose(yp, y * source["re_tau"], rtol=2e-4, atol=1e-4):
        raise ValueError("DNS y+ is inconsistent with normalized wall distance")
    return DNSChannelProfile(
        _readonly(y), _readonly(u), float(source["re_tau"]), digest
    )


def numerical_channel_profile(
    *,
    cell_y_m,
    velocity_m_s,
    cell_volumes_m3,
    half_height_m,
    wall_shear_Pa,
    viscosity_Pa_s,
):
    """Cubic numerical interpolant with actual wall shear and centre symmetry.

    A symmetric quartic through the three nearest centre-cell values supplies
    the centre value. No independent DNS/BVP values enter this reconstruction.
    """
    y, u, vol = (
        np.asarray(a, dtype=float) for a in (cell_y_m, velocity_m_s, cell_volumes_m3)
    )
    if y.ndim != 1 or u.shape != y.shape or vol.shape != y.shape or not y.size:
        raise ValueError("matching nonempty numerical cell vectors required")
    if not all(
        np.isfinite(x) and x > 0 for x in (half_height_m, wall_shear_Pa, viscosity_Pa_s)
    ):
        raise ValueError("positive physical interpolation parameters required")
    if (
        not np.all(np.isfinite([y, u, vol]))
        or np.any(vol <= 0)
        or np.any(y <= 0)
        or np.any(y >= 2 * half_height_m)
    ):
        raise ValueError("finite interior cell data and positive volumes required")
    # Polygon centroids can differ by a few ulps between streamwise cells of
    # the same physical row. Exact np.unique would create near-coincident
    # interpolation knots. Group only at a machine-precision geometry scale.
    row_tolerance = (
        128 * np.finfo(float).eps * max(half_height_m, float(np.max(np.abs(y))))
    )
    order = np.argsort(y, kind="stable")
    ordered = y[order]
    inverse = np.cumsum(np.r_[True, np.diff(ordered) > row_tolerance]) - 1
    weights = np.bincount(inverse, weights=vol[order])
    rows = np.bincount(inverse, weights=vol[order] * ordered) / weights
    means = np.bincount(inverse, weights=vol[order] * u[order]) / weights
    lower = np.full(len(rows), np.inf)
    upper = np.full(len(rows), -np.inf)
    np.minimum.at(lower, inverse, ordered)
    np.maximum.at(upper, inverse, ordered)
    if np.any(upper - lower > 4 * row_tolerance):
        raise ValueError("unresolved physical row separation; interpolation rejected")
    half = rows < half_height_m - row_tolerance
    ys, us = rows[half], means[half]
    if len(ys) < 3:
        raise ValueError("at least three independent half-channel rows required")
    at_center = np.flatnonzero(np.abs(rows - half_height_m) <= row_tolerance)
    if len(at_center):
        center = float(means[at_center[0]])
    else:
        distance = half_height_m - ys[-3:]
        scaled = distance / distance.max()
        center = float(
            np.linalg.solve(
                np.column_stack((np.ones(3), scaled**2, scaled**4)), us[-3:]
            )[0]
        )
    x = np.r_[0, ys, half_height_m]
    values = np.r_[0, us, center]
    return CubicSpline(x, values, bc_type=((1, wall_shear_Pa / viscosity_Pa_s), (1, 0)))


def compare_dns_profile(
    profile,
    *,
    cell_y_m,
    velocity_m_s,
    cell_volumes_m3,
    half_height_m,
    friction_velocity_m_s,
    wall_shear_Pa,
    viscosity_Pa_s,
):
    if not np.isfinite(friction_velocity_m_s) or friction_velocity_m_s <= 0:
        raise ValueError("positive reference friction velocity required")
    interpolant = numerical_channel_profile(
        cell_y_m=cell_y_m,
        velocity_m_s=velocity_m_s,
        cell_volumes_m3=cell_volumes_m3,
        half_height_m=half_height_m,
        wall_shear_Pa=wall_shear_Pa,
        viscosity_Pa_s=viscosity_Pa_s,
    )
    predicted = (
        interpolant(profile.y_over_delta * half_height_m) / friction_velocity_m_s
    )
    relative_rms = float(
        np.linalg.norm(predicted - profile.velocity_plus)
        / np.linalg.norm(profile.velocity_plus)
    )
    volume = np.asarray(cell_volumes_m3)
    bulk = float(np.asarray(velocity_m_s) @ volume / volume.sum())
    dns_bulk = profile.bulk_velocity_plus * friction_velocity_m_s
    return {
        "relative_rms_velocity_plus": relative_rms,
        "relative_bulk_velocity": abs(bulk / dns_bulk - 1),
        "bulk_velocity_m_s": bulk,
        "dns_bulk_velocity_plus": profile.bulk_velocity_plus,
        "reference_friction_velocity_m_s": friction_velocity_m_s,
    }
