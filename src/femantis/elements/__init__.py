"""Reference finite-element geometry and operators."""

from femantis.elements.quadrature import (
    QuadratureRule,
    tetrahedron_quadrature,
    triangle_quadrature,
)
from femantis.elements.tet4 import (
    DEFAULT_DEGENERACY_TOLERANCE,
    Tet4Face,
    Tet4Geometry,
    Tet4Quality,
)
from femantis.elements.tet10 import Tet10Face, Tet10Geometry, Triangle6Geometry

__all__ = (
    "DEFAULT_DEGENERACY_TOLERANCE",
    "QuadratureRule",
    "Tet4Face",
    "Tet4Geometry",
    "Tet4Quality",
    "Tet10Face",
    "Tet10Geometry",
    "Triangle6Geometry",
    "tetrahedron_quadrature",
    "triangle_quadrature",
)
