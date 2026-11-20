"""Validated uniform fluid mesh, SI material and explicit boundary contracts."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.stokes import _owned
from fem_toolbox.errors import ValidationError


@dataclass(frozen=True, slots=True)
class CartesianFluidGrid:
    """Uniform rectangle: cells=(nx,ny)>=2, size_m=(Lx,Ly)>0 in metres."""

    cells: tuple[int, int]
    size_m: tuple[float, float]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.cells, (list, tuple))
            or len(self.cells) != 2
            or any(type(n) is not int or n < 2 for n in self.cells)
        ):
            raise ValidationError("two integer cell counts >=2 required", field="cells")
        if not isinstance(self.size_m, (list, tuple)) or len(self.size_m) != 2:
            raise ValidationError("two SI lengths required", field="size_m")
        lengths = tuple(_scalar(v, "size_m") for v in self.size_m)
        if min(lengths) <= 0:
            raise ValidationError("domain lengths must be positive", field="size_m")
        hx, hy = lengths[0] / self.cells[0], lengths[1] / self.cells[1]
        scales = [hx * hx, hy * hy, hx * hy]
        if not np.all(np.isfinite(scales)) or min(scales) <= 0:
            raise ValidationError("unresolved cell spacing/area", field="size_m")
        object.__setattr__(self, "cells", tuple(self.cells))
        object.__setattr__(self, "size_m", lengths)


@dataclass(frozen=True, slots=True)
class IncompressibleFluid:
    """Constant density kg/m³ and dynamic viscosity Pa s, both positive."""

    density_kg_m3: float
    viscosity_Pa_s: float

    def __post_init__(self) -> None:
        for name in ("density_kg_m3", "viscosity_Pa_s"):
            value = _scalar(getattr(self, name), name)
            if value <= 0:
                raise ValidationError(f"{name} must be positive", field=name)
            object.__setattr__(self, name, value)


@dataclass(frozen=True, slots=True)
class ClosedWalls:
    """Zero normal speed; constant bottom/top/left/right tangential speeds m/s."""

    tangential_velocities_m_s: tuple[float, float, float, float] = (0, 0, 0, 0)

    def __post_init__(self) -> None:
        if (
            not isinstance(self.tangential_velocities_m_s, (list, tuple))
            or len(self.tangential_velocities_m_s) != 4
        ):
            raise ValidationError("four tangential wall speeds required", field="walls")
        object.__setattr__(
            self,
            "tangential_velocities_m_s",
            tuple(_scalar(v, "walls") for v in self.tangential_velocities_m_s),
        )


@dataclass(frozen=True, slots=True, init=False)
class ChannelBoundary:
    """Left u(y)>=0 inlet, v=0; right prescribed pressure and zero du/dn.

    Profile shape is (ny,), owned and read-only; top/bottom are no-slip.
    Reverse flow at the right outlet is unsupported and explicitly rejected.
    """

    inlet_velocity_m_s: NDArray[np.float64]
    outlet_pressure_Pa: float

    def __init__(
        self, inlet_velocity_m_s: ArrayLike, outlet_pressure_Pa: float = 0
    ) -> None:
        try:
            profile = np.asarray(inlet_velocity_m_s, dtype=np.float64)
        except (ValueError, TypeError, OverflowError) as error:
            raise ValidationError("invalid inlet profile", field="inlet") from error
        if (
            profile.ndim != 1
            or len(profile) < 2
            or not np.all(np.isfinite(profile))
            or np.any(profile < 0)
        ):
            raise ValidationError(
                "finite nonnegative 1D inlet profile required", field="inlet"
            )
        object.__setattr__(self, "inlet_velocity_m_s", _owned(profile))
        object.__setattr__(
            self,
            "outlet_pressure_Pa",
            _scalar(outlet_pressure_Pa, "outlet_pressure_Pa"),
        )
