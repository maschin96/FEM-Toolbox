"""2x3x4 m block, k=2 W/(m K), hot face 400 K, Robin h=4, ambient=300 K.

Expected: T=400-40*x K; outward flux=80 W/m²; transferred power=960 W.
Run after installing fem-toolbox; output is steady_heat.vtu.
"""

from itertools import permutations, product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.io import write_vtu
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.heat import Convection
from fem_toolbox.solvers.heat import solve_heat, thermal_dofs

points = np.array(list(product((0.0, 1.0), repeat=3))) * [2, 3, 4]
cells = []
for order in permutations(range(3)):
    node = np.zeros(3, dtype=int)
    cell = [0]
    for axis in order:
        node[axis] = 1
        cell.append(int(node @ [4, 2, 1]))
    cells.append(cell)
mesh = Mesh(points=points, tetra4=cells, triangle3=[[4, 6, 7], [4, 5, 7]])
left = np.flatnonzero(points[:, 0] == 0)
dofs = thermal_dofs(mesh)
result = solve_heat(
    mesh,
    ThermalMaterial.isotropic(2),
    [Dirichlet(dofs.indices("T", left), np.full(len(left), 400.0))],
    convection={0: Convection(4, 300), 1: Convection(4, 300)},
)
np.testing.assert_allclose(
    result.temperature, 400 - 40 * points[:, 0], atol=3e-12, rtol=0
)
np.testing.assert_allclose(result.convection_outflow, 960, atol=1e-10, rtol=0)
np.testing.assert_allclose(result.power_balance, 0, atol=1e-10, rtol=0)
print(
    f"Temperature: {result.temperature.min():.6f} to {result.temperature.max():.6f} K"
)
print(f"Convective outflow: {result.convection_outflow:.6f} W")
print(f"Power balance: {result.power_balance:.3e} W")
write_vtu(
    "steady_heat.vtu",
    mesh,
    point_data={"temperature_K": result.temperature, "reaction_W": result.reactions},
    cell_data={
        "heat_flux_W_m2": {"tetra4": result.heat_flux, "triangle3": np.zeros((2, 3))},
        "temperature_gradient_K_m": {
            "tetra4": result.gradient,
            "triangle3": np.zeros((2, 3)),
        },
    },
)
