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
