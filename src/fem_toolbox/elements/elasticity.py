"""Linear Tet4 mechanical operators in engineering-shear Voigt convention."""

import numpy as np
from numpy.typing import NDArray

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.materials import ElasticMaterial


def strain_matrix(geometry: Tet4Geometry) -> NDArray[np.float64]:
    """Return B (6,12), in 1/m, using the geometry's original node order."""
    b = np.zeros((6, 12))
    for node, (x, y, z) in enumerate(geometry.shape_gradients):
        b[:, 3 * node : 3 * node + 3] = [
            [x, 0, 0],
            [0, y, 0],
            [0, 0, z],
            [0, z, y],
            [z, 0, x],
            [y, x, 0],
        ]
    return b


def stiffness_matrix(
    geometry: Tet4Geometry, material: ElasticMaterial
) -> NDArray[np.float64]:
    """Return exactly integrated V B.T C B (12,12), in N/m."""
    b = strain_matrix(geometry)
    return np.asarray(
        geometry.volume * (b.T @ material.stiffness @ b), dtype=np.float64
    )
