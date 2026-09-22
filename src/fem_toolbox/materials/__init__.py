"""Validated, immutable material parameters in global axes."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.errors import ValidationError
from fem_toolbox.materials.electrostatics import (
    DielectricMaterial as DielectricMaterial,
)
from fem_toolbox.materials.heat import ThermalMaterial as ThermalMaterial
from fem_toolbox.materials.piezoelectric import (
    PiezoelectricMaterial as PiezoelectricMaterial,
)


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


@dataclass(frozen=True, init=False, eq=False)
class ThermalExpansion:
    """Global engineering-Voigt expansion (6,) in 1/K and reference T in K.

    Material-axis tensors must be rotated into global axes by the caller.
    ``source_axes`` records the local-to-global rotation only for values
    constructed through ``from_principal``; direct inputs are already global.
    """

    alpha: NDArray[np.float64]
    reference_temperature: float
    source_axes: NDArray[np.float64] | None

    def __init__(self, alpha: ArrayLike, reference_temperature: float) -> None:
        a = real_array(alpha, (6,))
        reference = float(real_array(reference_temperature, ()))
        if reference < 0:
            raise ValidationError(
                "reference temperature must be nonnegative K", field="temperature"
            )
        object.__setattr__(self, "alpha", np.frombuffer(a.tobytes(), dtype=np.float64))
        object.__setattr__(self, "reference_temperature", reference)
        object.__setattr__(self, "source_axes", None)

    @classmethod
    def from_principal(
        cls, principal_alpha: ArrayLike, reference_temperature: float, axes: ArrayLike
    ) -> "ThermalExpansion":
        """Rotate principal expansion (3,) in 1/K; axes columns are local axes."""
        values = real_array(principal_alpha, (3,))
        rotation = real_array(axes, (3, 3))
        if not np.allclose(
            rotation.T @ rotation, np.eye(3), rtol=0, atol=1e-12
        ) or not np.isclose(np.linalg.det(rotation), 1, rtol=0, atol=1e-12):
            raise ValidationError(
                "material axes must form a proper orthonormal frame", field="axes"
            )
        tensor = rotation @ np.diag(values) @ rotation.T
        material = cls(
            [
                tensor[0, 0],
                tensor[1, 1],
                tensor[2, 2],
                2 * tensor[1, 2],
                2 * tensor[0, 2],
                2 * tensor[0, 1],
            ],
            reference_temperature,
        )
        object.__setattr__(
            material,
            "source_axes",
            np.frombuffer(rotation.tobytes(), dtype=np.float64).reshape(3, 3),
        )
        return material


from fem_toolbox.materials import regions as _regions  # noqa: E402

MaterialAssignment = _regions.MaterialAssignment
assign_region_materials = _regions.assign_region_materials
