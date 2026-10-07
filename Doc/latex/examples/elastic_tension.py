"""Unit cube in uniaxial tension; run after installing fem-toolbox.

E=1000 Pa, nu=0.25, traction=10 Pa: ux at x=1 is 0.01 m.
Output: elastic_tension.vtu with displacement, raw element strain/stress.
"""

from itertools import permutations, product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.io import write_vtu
from fem_toolbox.materials import ElasticMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.elasticity import mechanical_dofs, solve_elasticity

points = np.array(list(product((0.0, 1.0), repeat=3)))
cells = []
for order in permutations(range(3)):
    node = np.zeros(3, dtype=int)
    cell = [0]
    for axis in order:
        node[axis] = 1
        cell.append(int(node @ [4, 2, 1]))
    cells.append(cell)
mesh = Mesh(points=points, tetra4=cells, triangle3=[[4, 6, 7], [4, 5, 7]])
dofs = mechanical_dofs(mesh)
conditions = []
for axis, component in enumerate(("x", "y", "z")):
    nodes = np.flatnonzero(points[:, axis] == 0)
    conditions.append(
        Dirichlet(dofs.indices("u", nodes, [component]), np.zeros(len(nodes)))
    )
result = solve_elasticity(
    mesh,
    ElasticMaterial.isotropic(1000.0, 0.25),
    conditions,
    tractions={0: [10.0, 0.0, 0.0], 1: [10.0, 0.0, 0.0]},
)
np.testing.assert_allclose(
    result.displacement, points * [0.01, -0.0025, -0.0025], atol=1e-15
)
np.testing.assert_allclose(result.force_balance, 0, atol=1e-12)
print(f"Maximum ux: {result.displacement[:, 0].max():.6f} m")
print(f"Strain energy: {result.energy:.6f} J")
print(f"Force balance: {result.force_balance} N")
write_vtu(
    "elastic_tension.vtu",
    mesh,
    point_data={"displacement_m": result.displacement, "reaction_N": result.reactions},
    cell_data={
        "stress_Pa": {"tetra4": result.stress, "triangle3": np.zeros((2, 6))},
        "strain": {"tetra4": result.strain, "triangle3": np.zeros((2, 6))},
    },
)
