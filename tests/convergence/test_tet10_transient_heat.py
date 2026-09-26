"""Independent manufactured time solution on quadratic elements."""

import numpy as np

from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.solvers.transient_heat import solve_transient_heat
from tests.reference.quadratic_mesh import elevate


def test_tet10_backward_euler_order_and_storage(cube_mesh):
    mesh = elevate(cube_mesh(1))
    material = ThermalMaterial.isotropic(2, density=3, specific_heat=4)
    errors = []
    for steps in (2, 4, 8, 16):
        # T=300+t², C=12, Q=24t, insulated domain. BE endpoint error=dt.
        result = solve_transient_heat(
            mesh,
            material,
            np.linspace(0, 1, steps + 1),
            300.0,
            source=lambda x, t: 24 * t,
        )
        error = np.max(abs(result.frames[-1].temperature - 301.0))
        errors.append(error)
        np.testing.assert_allclose(error, 1 / steps, rtol=0, atol=5e-11)
        assert max(abs(b.balance) for b in result.balances) < 2e-10
    rates = np.log2(np.array(errors[:-1]) / errors[1:])
    assert np.all((rates > 0.99) & (rates < 1.01))
