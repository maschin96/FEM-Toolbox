"""Exactly integrated Tet4 positive dielectric operator."""

import numpy as np
from numpy.typing import NDArray

from femantis.elements import Tet4Geometry
from femantis.materials import DielectricMaterial


def dielectric_matrix(
    geometry: Tet4Geometry, material: DielectricMaterial
) -> NDArray[np.float64]:
    """Return V G.T k G (4,4), in F, in original local node order."""
    g = geometry.gradient_matrix
    return np.asarray(
        geometry.volume * (g.T @ material.permittivity @ g), dtype=np.float64
    )
