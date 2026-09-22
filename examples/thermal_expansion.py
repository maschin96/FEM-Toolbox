"""One-metre cube at 313.15 K: free anisotropic thermal expansion.

Reference T=293.15 K, alpha=(10,20,0) microstrain/K. Expected displacement
u=(0.0002*x, 0.0004*y, 0) m, zero stress and thermal power balance.
Run after installing fem-toolbox; output is thermal_expansion.vtu.
"""

from itertools import permutations, product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.io import write_vtu
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.heat import thermal_dofs
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics

points = np.array(list(product((0.0, 1.0), repeat=3)))
cells = []
for order in permutations(range(3)):
    node = np.zeros(3, dtype=int)
    cell = [0]
    for axis in order:
        node[axis] = 1
        cell.append(int(node @ [4, 2, 1]))
    cells.append(cell)
mesh = Mesh(points=points, tetra4=cells)
thermal_bc = [Dirichlet(thermal_dofs(mesh).indices("T", [0]), [313.15])]

# Anchor six rigid motions without preventing the expected expansion.
x_node = int(np.flatnonzero(np.all(points == [1, 0, 0], axis=1))[0])
y_node = int(np.flatnonzero(np.all(points == [0, 1, 0], axis=1))[0])
fixed = [0, 1, 2, 3 * x_node + 1, 3 * x_node + 2, 3 * y_node + 2]
mechanical_bc = [Dirichlet(fixed, np.zeros(len(fixed)))]

result = solve_thermomechanics(
    mesh,
    ThermalMaterial.isotropic(2),
    ElasticMaterial.isotropic(100e9, 0.25),
    ThermalExpansion([10e-6, 20e-6, 0, 0, 0, 0], 293.15),
    thermal_bc,
    mechanical_bc,
)
expected = points * [0.0002, 0.0004, 0]
np.testing.assert_allclose(result.heat.temperature, 313.15, atol=1e-11, rtol=0)
np.testing.assert_allclose(result.elasticity.displacement, expected, atol=1e-17, rtol=0)
np.testing.assert_allclose(result.elasticity.stress, 0, atol=1e-5, rtol=0)
np.testing.assert_allclose(result.heat.power_balance, 0, atol=1e-9, rtol=0)
np.testing.assert_allclose(result.elasticity.force_balance, 0, atol=1e-5, rtol=0)
print(f"Maximum displacement: {np.max(np.linalg.norm(expected, axis=1)):.9f} m")
print(
    f"Maximum stress norm: {np.max(np.linalg.norm(result.elasticity.stress, axis=1)):.3e} Pa"
)
print(f"Heat power balance: {result.heat.power_balance:.3e} W")
write_vtu(
    "thermal_expansion.vtu",
    mesh,
    point_data={
        "temperature_K": result.heat.temperature,
        "displacement_m": result.elasticity.displacement,
    },
    cell_data={"stress_Pa": {"tetra4": result.elasticity.stress}},
)
