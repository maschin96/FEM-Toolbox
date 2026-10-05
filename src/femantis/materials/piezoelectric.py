"""Linear piezoelectric material in the project's stress-charge convention."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from femantis.assembly import real_array
from femantis.errors import ValidationError

_PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))


def _owned(value: NDArray[np.float64]) -> NDArray[np.float64]:
    """Store a truly immutable copy, including when callers retain their input."""
    return np.frombuffer(
        np.ascontiguousarray(value).tobytes(), dtype=np.float64
    ).reshape(value.shape)


def _positive_symmetric(
    value: ArrayLike, shape: tuple[int, int], field: str
) -> NDArray[np.float64]:
    matrix = real_array(value, shape)
    scale = float(np.max(np.abs(matrix)))
    if scale == 0 or np.max(np.abs(matrix - matrix.T)) > 1e-12 * scale:
        raise ValidationError(f"{field} must be symmetric and nonzero", field=field)
    matrix = (matrix + matrix.T) / 2
    eigenvalues = np.linalg.eigvalsh(matrix / scale)
    if eigenvalues[0] <= 1e-12 * eigenvalues[-1]:
        raise ValidationError(
            f"{field} must be positive definite and well resolved", field=field
        )
    return matrix


def _axes(value: ArrayLike | None) -> NDArray[np.float64]:
    axes = np.eye(3) if value is None else real_array(value, (3, 3))
    if not np.allclose(axes.T @ axes, np.eye(3), rtol=0, atol=1e-12) or not np.isclose(
        np.linalg.det(axes), 1, rtol=0, atol=1e-12
    ):
        raise ValidationError(
            "axes must be a proper orthonormal frame with local axes as columns",
            field="axes",
        )
    return axes


def _voigt_rotation(axes: NDArray[np.float64], *, strain: bool) -> NDArray[np.float64]:
    """Map a local Voigt vector to global, with engineering strain shears."""
    transform = np.zeros((6, 6))
    for column, (i, j) in enumerate(_PAIRS):
        tensor = np.zeros((3, 3))
        tensor[i, j] = 0.5 if strain and i != j else 1.0
        tensor[j, i] = tensor[i, j]
        rotated = axes @ tensor @ axes.T
        transform[:, column] = [
            rotated[a, b] * (2 if strain and a != b else 1) for a, b in _PAIRS
        ]
    return transform


@dataclass(frozen=True, init=False, eq=False)
class PiezoelectricMaterial:
    """Immutable stress-charge tensors in global Cartesian axes.

    ``stiffness`` is C^E (6,6) Pa, ``stress_coupling`` is e (3,6) C/m²,
    and ``permittivity`` is kappa^S (3,3) F/m. Strain uses engineering shear
    and Voigt order xx, yy, zz, yz, xz, xy. Electric components are x, y, z.
    The constitutive laws are sigma=C^E epsilon-e.T E and
    D=e epsilon+kappa^S E. ``axes`` columns are local basis vectors in global
    coordinates; omitted axes mean input tensors already use global axes.
    Invalid shape, symmetry, stability, axes, or nonfinite data raise
    ValidationError. All arrays are owned read-only copies.
    """

    stiffness: NDArray[np.float64]
    stress_coupling: NDArray[np.float64]
    permittivity: NDArray[np.float64]
    axes: NDArray[np.float64]
    source_form: Literal["stress_charge", "strain_charge"]
    _source: Mapping[str, NDArray[np.float64]]

    def __init__(
        self,
        stiffness: ArrayLike,
        stress_coupling: ArrayLike,
        permittivity: ArrayLike,
        *,
        axes: ArrayLike | None = None,
    ) -> None:
        source_c = real_array(stiffness, (6, 6))
        source_e = real_array(stress_coupling, (3, 6))
        source_k = real_array(permittivity, (3, 3))
        c = _positive_symmetric(stiffness, (6, 6), "stiffness")
        e = source_e
        k = _positive_symmetric(permittivity, (3, 3), "permittivity")
        rotation = _axes(axes)
        stress_to_global = _voigt_rotation(rotation, strain=False)
        strain_to_local = _voigt_rotation(rotation.T, strain=True)
        c_global = stress_to_global @ c @ strain_to_local
        e_global = rotation @ e @ strain_to_local
        k_global = rotation @ k @ rotation.T
        object.__setattr__(
            self,
            "stiffness",
            _owned(_positive_symmetric(c_global, (6, 6), "stiffness")),
        )
        object.__setattr__(
            self, "stress_coupling", _owned(real_array(e_global, (3, 6)))
        )
        object.__setattr__(
            self,
            "permittivity",
            _owned(_positive_symmetric(k_global, (3, 3), "permittivity")),
        )
        object.__setattr__(self, "axes", _owned(rotation))
        object.__setattr__(self, "source_form", "stress_charge")
        object.__setattr__(
            self,
            "_source",
            MappingProxyType(
                {
                    "stiffness": _owned(source_c),
                    "stress_coupling": _owned(source_e),
                    "permittivity": _owned(source_k),
                }
            ),
        )

    @classmethod
    def from_strain_charge(
        cls,
        compliance: ArrayLike,
        strain_coupling: ArrayLike,
        permittivity_at_stress: ArrayLike,
        *,
        axes: ArrayLike | None = None,
    ) -> PiezoelectricMaterial:
        """Convert s^E (1/Pa), d (C/N), kappa^T (F/m) to stress-charge.

        Requires both s^E and the resulting kappa^S = kappa^T-d C^E d.T
        to be positive definite. Conversion uses linear solves, not an explicit
        inverse, and retains the supplied local tensors in provenance.
        """
        source_s = real_array(compliance, (6, 6))
        source_k = real_array(permittivity_at_stress, (3, 3))
        s = _positive_symmetric(compliance, (6, 6), "compliance")
        d = real_array(strain_coupling, (3, 6))
        k_t = _positive_symmetric(permittivity_at_stress, (3, 3), "permittivity")
        c = np.linalg.solve(s, np.eye(6))
        e = d @ c
        k_s = k_t - d @ c @ d.T
        material = cls(c, e, k_s, axes=axes)
        object.__setattr__(material, "source_form", "strain_charge")
        object.__setattr__(
            material,
            "_source",
            MappingProxyType(
                {
                    "compliance": _owned(source_s),
                    "strain_coupling": _owned(d),
                    "permittivity_at_stress": _owned(source_k),
                }
            ),
        )
        return material

    def as_strain_charge(
        self,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Return detached global s^E, d, kappa^T tensors in SI units."""
        s = np.asarray(np.linalg.solve(self.stiffness, np.eye(6)), dtype=np.float64)
        d = np.asarray(self.stress_coupling @ s, dtype=np.float64)
        k_t = np.asarray(
            self.permittivity + d @ self.stress_coupling.T, dtype=np.float64
        )
        return s, d, k_t

    def to_provenance_dict(self) -> dict[str, object]:
        """Return detached JSON-compatible local inputs and global tensors."""
        return {
            "source_form": self.source_form,
            "axes_local_to_global": self.axes.tolist(),
            "source": {name: value.tolist() for name, value in self._source.items()},
            "global": {
                "stiffness": self.stiffness.tolist(),
                "stress_coupling": self.stress_coupling.tolist(),
                "permittivity": self.permittivity.tolist(),
            },
            "voigt_order": ["xx", "yy", "zz", "yz", "xz", "xy"],
            "shear_strain": "engineering",
            "electric_field": "E = -grad(phi)",
            "units": {
                "stiffness": "Pa",
                "stress_coupling": "C/m^2",
                "permittivity": "F/m",
            },
            "source_units": (
                {
                    "compliance": "1/Pa",
                    "strain_coupling": "C/N",
                    "permittivity_at_stress": "F/m",
                }
                if self.source_form == "strain_charge"
                else {
                    "stiffness": "Pa",
                    "stress_coupling": "C/m^2",
                    "permittivity": "F/m",
                }
            ),
        }
