"""Constant anisotropic thermal conductivity in global Cartesian axes."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError


@dataclass(frozen=True)
class ThermalSource:
    """Traceable preset identity and source; no claim of a validity interval."""

    key: str
    source_url: str
    reference_temperature: float
    note: str


@dataclass(frozen=True, init=False, eq=False)
class ThermalMaterial:
    """Immutable symmetric positive-definite (3,3) conductivity in W/(m K).

    Optional density [kg/m³] and specific_heat [J/(kg K)] must be supplied
    together, finite and positive, for transient analysis.

    Symmetry tolerance is 1e-12 times the largest coefficient. Accepted roundoff
    asymmetry is averaged; eigenvalue ratios <= 1e-12 are rejected. Input arrays
    are copied. Material-axis rotation and temperature dependence are not implicit.
    """

    conductivity: NDArray[np.float64]

    density: float | None
    specific_heat: float | None
    source: ThermalSource | None

    def __init__(
        self,
        conductivity: ArrayLike,
        *,
        density: float | None = None,
        specific_heat: float | None = None,
        source: ThermalSource | None = None,
    ) -> None:
        if (density is None) != (specific_heat is None):
            raise ValidationError(
                "provide both density and specific_heat", field="storage"
            )
        if density is not None and specific_heat is not None:
            density = float(real_array(density, ()))
            specific_heat = float(real_array(specific_heat, ()))
            if density <= 0 or specific_heat <= 0:
                raise ValidationError(
                    "density and specific_heat must be positive", field="storage"
                )
        if source is not None and not isinstance(source, ThermalSource):
            raise ValidationError("expected ThermalSource", field="source")
        object.__setattr__(self, "density", density)
        object.__setattr__(self, "specific_heat", specific_heat)
        object.__setattr__(self, "source", source)
        k = real_array(conductivity, (3, 3))
        scale = float(np.max(np.abs(k)))
        if scale == 0 or np.max(np.abs(k - k.T)) > 1e-12 * scale:
            raise ValidationError(
                "conductivity must be symmetric and nonzero", field="conductivity"
            )
        k = k / 2 + k.T / 2
        eigenvalues = np.linalg.eigvalsh(k / scale)
        if eigenvalues[0] <= 1e-12 * eigenvalues[-1]:
            raise ValidationError(
                "conductivity must be positive definite and well resolved",
                field="conductivity",
            )
        object.__setattr__(
            self,
            "conductivity",
            np.frombuffer(k.tobytes(), dtype=np.float64).reshape(3, 3),
        )

    @classmethod
    def isotropic(
        cls,
        conductivity: float,
        *,
        density: float | None = None,
        specific_heat: float | None = None,
        source: ThermalSource | None = None,
    ) -> "ThermalMaterial":
        """Construct from a positive finite scalar conductivity in W/(m K)."""
        k = float(real_array(conductivity, ()))
        if k <= 0:
            raise ValidationError("conductivity must be positive", field="conductivity")
        return cls(
            np.eye(3) * k, density=density, specific_heat=specific_heat, source=source
        )
