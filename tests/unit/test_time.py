from decimal import Decimal, localcontext

import numpy as np
import pytest

from fem_toolbox.errors import ValidationError
from fem_toolbox.time import Ramp, Step, TimeSeries, integration_times


def test_step_ramp_and_event_grid():
    step = Step(0.4, before=2.0, after=6.0)
    assert step(0.4) == 6.0
    assert step.value(0.4, side="left") == 2.0
    ramp = Ramp(0.2, 0.8, start_value=1.0, end_value=4.0)
    assert ramp(0.5) == pytest.approx(2.5)
    assert ramp(0.0) == 1.0 and ramp(1.0) == 4.0
    np.testing.assert_array_equal(
        integration_times([0.0, 1.0], [step, ramp]), [0.0, 0.2, 0.4, 0.8, 1.0]
    )


def test_measured_curve_roundtrip(tmp_path):
    path = tmp_path / "data.csv"
    path.write_text("time_s,power_W\n0,10\n1,30\n2,20\n")
    curve = TimeSeries.from_csv(path, value_column="power_W", unit="W")
    assert curve(0.25) == 15.0
    assert curve.unit == "W"
    with pytest.raises(ValidationError):
        curve(3.0)
    hold = TimeSeries([0.0, 1.0, 2.0], [10.0, 30.0, 20.0], interpolation="previous")
    assert hold(1.0) == 30.0
    assert hold.value(1.0, side="left") == 10.0
    with pytest.raises(ValueError):
        curve.times[0] = -1


@pytest.mark.parametrize("times", [[0, 0], [1, 0], [0, np.nan], [0], [True, False]])
def test_invalid_times(times):
    with pytest.raises(ValidationError):
        integration_times(times)


def test_invalid_series_and_ramp():
    with pytest.raises(ValidationError):
        TimeSeries([0, 1], [0, np.inf])
    with pytest.raises(ValidationError):
        Ramp(1.0, 1.0, 0.0, 1.0)
    with pytest.raises(ValidationError):
        TimeSeries([0, 1], [0, 1], interpolation="spline")


@pytest.mark.parametrize(
    ("time", "expected"),
    [(-5e307, 0.25), (0.0, 0.5), (5e307, 0.75)],
)
def test_ramp_interpolates_when_finite_time_span_overflows(time, expected):
    ramp = Ramp(-1e308, 1e308, 0.0, 1.0)
    # Fractions 1/4, 1/2 and 3/4 follow from the symmetric interval,
    # independently of the floating-point implementation. Eight epsilons
    # cover the few normalized operations for a dimensionless O(1) result.
    assert ramp(time) == pytest.approx(
        expected, rel=0, abs=8 * np.finfo(np.float64).eps
    )


@pytest.mark.parametrize("fraction", [0.1, 0.3, 0.5, 0.7, 0.9])
def test_ramp_constant_at_largest_finite_amplitude(fraction):
    maximum = np.finfo(np.float64).max
    ramp = Ramp(0.0, 1.0, maximum, maximum)
    # A constant endpoint value must be returned exactly: no arithmetic
    # approximation is needed and infinity is never a valid result.
    assert ramp(fraction) == maximum


@pytest.mark.parametrize(
    ("before", "after", "fraction", "expected_fraction"),
    [
        (1.0, -1.0, 0.25, 0.5),
        (-1.0, 1.0, 0.75, 0.5),
        (1.0, 0.5, 0.5, 0.75),
        (-1.0, -0.5, 0.5, -0.75),
        (0.5, 1.0, 0.75, 0.875),
        (-0.5, -1.0, 0.25, -0.625),
    ],
)
def test_ramp_large_finite_amplitudes(before, after, fraction, expected_fraction):
    maximum = np.finfo(np.float64).max
    ramp = Ramp(0.0, 1.0, before * maximum, after * maximum)
    result = ramp(fraction)
    assert np.isfinite(result)
    # Normalize before comparison: absolute tolerance 8 eps corresponds to
    # a relative O(eps) amplitude error without overflowing the test itself.
    assert result / maximum == pytest.approx(
        expected_fraction, rel=0, abs=8 * np.finfo(np.float64).eps
    )


@pytest.mark.parametrize("side", ["left", "right"])
def test_ramp_extreme_endpoints_and_holds_are_exact(side):
    maximum = np.finfo(np.float64).max
    ramp = Ramp(-1e308, 1e308, maximum, -maximum)
    for time, expected in [
        (-maximum, maximum),
        (-1e308, maximum),
        (1e308, -maximum),
        (maximum, -maximum),
    ]:
        assert ramp.value(time, side=side) == expected
    assert ramp.value(0.0, side=side) == 0.0


@pytest.mark.parametrize(
    ("start", "end", "time"),
    [
        (1e308, np.nextafter(1e308, np.inf), 1e308),
        (0.0, 2 * np.nextafter(0.0, 1.0), np.nextafter(0.0, 1.0)),
    ],
)
def test_ramp_resolvable_small_time_intervals(start, end, time):
    ramp = Ramp(start, end, 2.0, 6.0)
    expected = 2.0 if time == start else 4.0
    assert ramp(time) == expected


@pytest.mark.parametrize("invalid", [np.inf, -np.inf, np.nan])
def test_ramp_rejects_nonfinite_parameters_and_query_times(invalid):
    for name in ("start", "end", "start_value", "end_value"):
        parameters = dict(start=0.0, end=1.0, start_value=2.0, end_value=3.0)
        parameters[name] = invalid
        with pytest.raises(ValidationError):
            Ramp(**parameters)
    with pytest.raises(ValidationError):
        Ramp(0.0, 1.0, 2.0, 3.0).value(invalid)


@pytest.mark.parametrize(
    ("start", "end", "time", "before", "after"),
    [
        (-1e308, 1e308, 3e307, 8e307, 1.6e308),
        (-1.6e308, 8e307, 0.0, -1.7e308, 1.2e308),
        (-4e307, 8e307, 1e307, -1e308, -5e307),
        (-1.7e308, 1.7e308, -6e307, 1.7e308, -1.7e308),
        (0.0, 1.0, 0.125, 1.7e308, 1.6e308),
        (0.0, 1.0, 0.875, -1.6e308, -1.7e308),
    ],
)
def test_ramp_matches_independent_high_precision_reference(
    start, end, time, before, after
):
    with localcontext() as context:
        context.prec = 100
        a, b, t, u, v = map(Decimal.from_float, (start, end, time, before, after))
        reference = u + (t - a) * (v - u) / (b - a)
        scale = max(abs(before), abs(after))
        expected = float(reference / Decimal.from_float(scale))
    result = Ramp(start, end, before, after)(time)
    assert np.isfinite(result)
    # Decimal evaluates the analytical straight line at the exact binary64
    # input values. Eight epsilons of the endpoint-amplitude scale cover the
    # small fixed number of float64 operations, including cancellation.
    assert result / scale == pytest.approx(
        expected, rel=0, abs=8 * np.finfo(np.float64).eps
    )
