"""Deterministic SI time grids and scalar prescribed histories (no expressions)."""

from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, TypeAlias

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError

Side = Literal["left", "right"]


def _scalar(value: float) -> float:
    return float(real_array(value, ()))


def _times(value: ArrayLike) -> NDArray[np.float64]:
    raw = np.asarray(value)
    if raw.ndim != 1 or len(raw) < 2:
        raise ValidationError("need at least two times", field="times")
    times = real_array(raw, raw.shape)
    if np.any(np.diff(times) <= 0) or not np.isfinite(np.diff(times)).all():
        raise ValidationError("times must be strictly increasing", field="times")
    return times


@dataclass(frozen=True)
class Step:
    """Right-continuous scalar step; time in s, before/after in caller's SI unit."""

    time: float
    before: float
    after: float

    def __post_init__(self) -> None:
        for name in ("time", "before", "after"):
            object.__setattr__(self, name, _scalar(getattr(self, name)))

    @property
    def knots(self) -> tuple[float, ...]:
        return (self.time,)

    def value(self, time: float, *, side: Side = "right") -> float:
        t = _scalar(time)
        return (
            self.before
            if t < self.time or (t == self.time and side == "left")
            else self.after
        )

    def __call__(self, time: float) -> float:
        return self.value(time)


@dataclass(frozen=True)
class Ramp:
    """Linear ramp on [start,end] s, held constant before and after."""

    start: float
    end: float
    start_value: float
    end_value: float

    def __post_init__(self) -> None:
        for name in ("start", "end", "start_value", "end_value"):
            object.__setattr__(self, name, _scalar(getattr(self, name)))
        if self.end <= self.start:
            raise ValidationError("ramp end must exceed start", field="time")

    @property
    def knots(self) -> tuple[float, ...]:
        return (self.start, self.end)

    def value(self, time: float, *, side: Side = "right") -> float:
        fraction = float(
            np.clip((_scalar(time) - self.start) / (self.end - self.start), 0, 1)
        )
        return (1 - fraction) * self.start_value + fraction * self.end_value

    def __call__(self, time: float) -> float:
        return self.value(time)


@dataclass(frozen=True, init=False, eq=False)
class TimeSeries:
    """Immutable measured scalar samples. Times in seconds, values in stated SI unit.

    Linear interpolation is default; previous is right-continuous zero-order
    hold. Extrapolation defaults to error; hold must be selected explicitly.
    Duplicate/nonfinite times, nonfinite values and unknown modes are rejected.
    """

    times: NDArray[np.float64]
    values: NDArray[np.float64]
    interpolation: Literal["linear", "previous"]
    extrapolation: Literal["error", "hold"]
    unit: str

    def __init__(
        self,
        times: ArrayLike,
        values: ArrayLike,
        *,
        interpolation: Literal["linear", "previous"] = "linear",
        extrapolation: Literal["error", "hold"] = "error",
        unit: str = "1",
    ) -> None:
        t = _times(times)
        v = real_array(values, t.shape)
        if interpolation not in ("linear", "previous") or extrapolation not in (
            "error",
            "hold",
        ):
            raise ValidationError("unknown interpolation/extrapolation", field="series")
        if not isinstance(unit, str) or not unit.strip():
            raise ValidationError("unit must be specified", field="unit")
        object.__setattr__(self, "times", np.frombuffer(t.tobytes(), dtype=np.float64))
        object.__setattr__(self, "values", np.frombuffer(v.tobytes(), dtype=np.float64))
        object.__setattr__(self, "interpolation", interpolation)
        object.__setattr__(self, "extrapolation", extrapolation)
        object.__setattr__(self, "unit", unit)

    @property
    def knots(self) -> tuple[float, ...]:
        return tuple(float(t) for t in self.times)

    def value(self, time: float, *, side: Side = "right") -> float:
        t = _scalar(time)
        if self.extrapolation == "error" and not self.times[0] <= t <= self.times[-1]:
            raise ValidationError("time outside measured range", field="series")
        if self.interpolation == "linear":
            return float(np.interp(t, self.times, self.values))
        index = int(np.searchsorted(self.times, t, side=side)) - 1
        return float(self.values[max(0, min(index, len(self.values) - 1))])

    def __call__(self, time: float) -> float:
        return self.value(time)

    @classmethod
    def from_csv(
        cls,
        path: str | Path,
        *,
        value_column: str,
        unit: str,
        time_column: str = "time_s",
        delimiter: str = ",",
        interpolation: Literal["linear", "previous"] = "linear",
        extrapolation: Literal["error", "hold"] = "error",
    ) -> TimeSeries:
        """Read numeric UTF-8 CSV columns, explicitly in s and the given SI unit."""
        try:
            with Path(path).open(encoding="utf-8-sig", newline="") as stream:
                reader = csv.DictReader(stream, delimiter=delimiter)
                if not reader.fieldnames or len(set(reader.fieldnames)) != len(
                    reader.fieldnames
                ):
                    raise ValueError("missing or duplicate CSV headers")
                rows = list(reader)
                return cls(
                    [float(r[time_column]) for r in rows],
                    [float(r[value_column]) for r in rows],
                    unit=unit,
                    interpolation=interpolation,
                    extrapolation=extrapolation,
                )
        except (OSError, KeyError, TypeError, ValueError, csv.Error) as error:
            raise ValidationError(
                f"invalid measured CSV: {error}", field="series"
            ) from error


TimeSignal: TypeAlias = float | Step | Ramp | TimeSeries


def signal_value(signal: TimeSignal, time: float, *, side: Side = "right") -> float:
    """Evaluate a finite scalar history at a specified side of an event."""
    if isinstance(signal, (Step, Ramp, TimeSeries)):
        return signal.value(time, side=side)
    return _scalar(signal)


def integration_times(
    times: ArrayLike,
    signals: Iterable[TimeSignal] = (),
) -> NDArray[np.float64]:
    """Validate times and insert all interior signal knots, without rounding."""
    grid = _times(times)
    events: list[float] = []
    for signal in signals:
        signal_value(signal, float(grid[0]))
        signal_value(signal, float(grid[-1]))
        if isinstance(signal, (Step, Ramp, TimeSeries)):
            events.extend(t for t in signal.knots if grid[0] < t < grid[-1])
    return np.unique(np.concatenate((grid, events)))
