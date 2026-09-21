"""Exactly integrated Tet4 thermal conductance operator."""

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ThermalMaterial


def conductance_matrix(
    geometry: Tet4Geometry, material: ThermalMaterial
) -> NDArray[np.float64]:
    """Return V G.T k G (4,4), in W/K, in original local node order."""
    g = geometry.gradient_matrix
    return np.asarray(
        geometry.volume * (g.T @ material.conductivity @ g), dtype=np.float64
    )
