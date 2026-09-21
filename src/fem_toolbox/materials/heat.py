"""Constant anisotropic thermal conductivity in global Cartesian axes."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError


@dataclass(frozen=True, init=False, eq=False)
class ThermalMaterial:
    """Immutable symmetric positive-definite (3,3) conductivity in W/(m K).

    Symmetry tolerance is 1e-12 times the largest coefficient. Accepted roundoff
    asymmetry is averaged; eigenvalue ratios <= 1e-12 are rejected. Input arrays
    are copied. Material-axis rotation and temperature dependence are not implicit.
    """

    conductivity: NDArray[np.float64]

    def __init__(self, conductivity: ArrayLike) -> None:
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
    def isotropic(cls, conductivity: float) -> "ThermalMaterial":
        """Construct from a positive finite scalar conductivity in W/(m K)."""
        k = float(real_array(conductivity, ()))
        if k <= 0:
            raise ValidationError("conductivity must be positive", field="conductivity")
        return cls(np.eye(3) * k)
