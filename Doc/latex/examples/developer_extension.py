"""Local teaching extension, using a research load hook; no production API change."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class LinearHeatSource:
    """Q(x)=offset+gradient.x; SI offset W/m^3, gradient W/m^4, point m."""

    offset: float
    gradient: tuple[float, float, float]

    def __post_init__(self) -> None:
        values = np.asarray(self.gradient, dtype=float)
        if (
            values.shape != (3,)
            or not np.isfinite(values).all()
            or not np.isfinite(self.offset)
        ):
            raise ValueError("finite offset and gradient[3] are required")
        object.__setattr__(self, "gradient", tuple(float(v) for v in values))

    def __call__(self, point: NDArray[np.float64]) -> float:
        return float(self.offset + np.dot(self.gradient, point))


def verify() -> None:
    from numpy.testing import assert_allclose

    from fem_toolbox.elements import Tet4Geometry
    from fem_toolbox.physics.heat import source_load

    geometry = Tet4Geometry(np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]))
    values = source_load(geometry, LinearHeatSource(2, (3, 0, 0)))
    # Exact barycentric integrals: V/20 off-diagonal, V/10 diagonal.
    assert_allclose(values, [13 / 120, 2 / 15, 13 / 120, 13 / 120], atol=2e-15, rtol=0)
    assert_allclose(values.sum(), 11 / 24, atol=2e-15, rtol=0)
    print("Linear source extension: independent integral passed")


if __name__ == "__main__":
    verify()
