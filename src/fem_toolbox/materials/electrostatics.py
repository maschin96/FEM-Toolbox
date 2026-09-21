"""Constant anisotropic absolute dielectric permittivity in global Cartesian axes."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError


@dataclass(frozen=True, init=False, eq=False)
class DielectricMaterial:
    """Immutable symmetric positive-definite (3,3) permittivity in F/m.

    Symmetry tolerance is 1e-12 times the largest coefficient. Accepted roundoff
    asymmetry is averaged; eigenvalue ratios <= 1e-12 are rejected. Input arrays
    are copied. Material-axis rotation and field dependence are not implicit.
    """

    permittivity: NDArray[np.float64]

    def __init__(self, permittivity: ArrayLike) -> None:
        k = real_array(permittivity, (3, 3))
        scale = float(np.max(np.abs(k)))
        if scale == 0 or np.max(np.abs(k - k.T)) > 1e-12 * scale:
            raise ValidationError(
                "permittivity must be symmetric and nonzero", field="permittivity"
            )
        k = k / 2 + k.T / 2
        eigenvalues = np.linalg.eigvalsh(k / scale)
        if eigenvalues[0] <= 1e-12 * eigenvalues[-1]:
            raise ValidationError(
                "permittivity must be positive definite and well resolved",
                field="permittivity",
            )
        object.__setattr__(
            self,
            "permittivity",
            np.frombuffer(k.tobytes(), dtype=np.float64).reshape(3, 3),
        )

    @classmethod
    def isotropic(cls, permittivity: float) -> "DielectricMaterial":
        """Construct from a positive finite scalar permittivity in F/m."""
        k = float(real_array(permittivity, ()))
        if k <= 0:
            raise ValidationError("permittivity must be positive", field="permittivity")
        return cls(np.eye(3) * k)
