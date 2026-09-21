"""Validated, immutable material parameters in global axes."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials.heat import ThermalMaterial as ThermalMaterial


@dataclass(frozen=True, init=False, eq=False)
class ElasticMaterial:
    """Symmetric positive-definite 6x6 stiffness in Pa.

    Voigt order is xx, yy, zz, yz, xz, xy with engineering shear strains.
    Anisotropic data must already be expressed in the global coordinate frame.
    Symmetry tolerance is 1e-12 times the largest coefficient; accepted rounding
    asymmetry is averaged. Eigenvalues <= 1e-12 times the largest eigenvalue are
    rejected as numerically unstable. Storage is copied onto immutable bytes.
    """

    stiffness: NDArray[np.float64]

    def __init__(self, stiffness: ArrayLike) -> None:
        c = real_array(stiffness, (6, 6))
        scale = float(np.max(np.abs(c)))
        if scale == 0 or np.max(np.abs(c - c.T)) > 1e-12 * scale:
            raise ValidationError(
                "stiffness must be symmetric and nonzero", field="stiffness"
            )
        c = c / 2 + c.T / 2
        eigenvalues = np.linalg.eigvalsh(c / scale)
        if eigenvalues[0] <= 1e-12 * eigenvalues[-1]:
            raise ValidationError(
                "stiffness must be positive definite and well resolved",
                field="stiffness",
            )
        object.__setattr__(
            self,
            "stiffness",
            np.frombuffer(c.tobytes(), dtype=np.float64).reshape(6, 6),
        )

    @classmethod
    def isotropic(cls, young: float, poisson: float) -> "ElasticMaterial":
        """Construct from Young's modulus (Pa) and Poisson ratio (dimensionless).

        Requires E > 0 and -1 < nu < 0.49; nearly incompressible Tet4 inputs are
        rejected to avoid presenting a locking-prone formulation as suitable.
        """
        e, nu = real_array([young, poisson], (2,))
        if e <= 0 or not -1 < nu < 0.49:
            raise ValidationError(
                "require E > 0 and -1 < nu < 0.49 for Tet4", field="elasticity"
            )
        mu = e / (2 * (1 + nu))
        lam = e * nu / ((1 + nu) * (1 - 2 * nu))
        c = np.diag([2 * mu] * 3 + [mu] * 3)
        c[:3, :3] += lam
        return cls(c)
