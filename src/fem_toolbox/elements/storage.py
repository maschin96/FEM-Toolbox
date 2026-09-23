"""Exactly integrated, constant-coefficient Tet4 storage operators."""

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.assembly import real_array
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError


def _positive(value: float, name: str) -> float:
    result = float(real_array(value, ()))
    if result <= 0:
        raise ValidationError(f"{name} must be positive", field=name)
    return result


def _storage(total: float, lumped: bool) -> NDArray[np.float64]:
    if not isinstance(lumped, bool):
        raise ValidationError("lumped must be a boolean", field="lumped")
    if not np.isfinite(total) or total <= 0:
        raise ValidationError("storage total is not representable", field="storage")
    coefficient = total / (4 if lumped else 20)
    if coefficient == 0:
        raise ValidationError("storage coefficients underflow", field="storage")
    return coefficient * (np.eye(4) if lumped else np.ones((4, 4)) + np.eye(4))


def mass_matrix(
    geometry: Tet4Geometry, density: float, *, lumped: bool = False
) -> NDArray[np.float64]:
    """Return mechanical mass (12,12) in kg, node-major xyz local order.

    Density is constant on the element, finite and positive in kg/m³.
    The default integrates rho Ni Nj exactly; lumped=True selects row-sum
    lumping. Invalid inputs or unrepresentable storage raise ValidationError.
    Electrical DOFs carry no mass and are not part of this operator.
    """
    scalar = _storage(_positive(density, "density") * geometry.volume, lumped)
    return np.asarray(np.kron(scalar, np.eye(3)), dtype=np.float64)


def capacity_matrix(
    geometry: Tet4Geometry,
    density: float,
    specific_heat: float,
    *,
    lumped: bool = False,
) -> NDArray[np.float64]:
    """Return thermal capacity (4,4) in J/K in original local node order.

    Density [kg/m³] and specific_heat [J/(kg K)] must be finite positive
    constants on this element. Default: exact integral rho c Ni Nj;
    lumped=True selects row-sum lumping. Invalid inputs and unrepresentable
    storage raise ValidationError. No temperature dependence is implied.
    """
    rho = _positive(density, "density")
    c = _positive(specific_heat, "specific_heat")
    return _storage(rho * c * geometry.volume, lumped)
