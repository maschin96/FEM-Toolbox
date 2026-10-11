"""Direct and inverse d33 effect in a transversely constrained piezo block.

L=20 mm, A=10x10 mm², C33=1 GPa, e33=5 C/m², kappa33=1e-8 F/m.
100 V gives u_top=-0.5 µm and Q_top=17.5 nC; 10 N at short circuit
gives u_top=2 µm and Q_top=-50 nC. Output: piezo_d33.vtu.
"""

from itertools import permutations, product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.io import write_vtu
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric

length = 0.02
area = 0.01**2
points = np.array(list(product((0.0, 1.0), repeat=3))) * [0.01, 0.01, length]
cells = []
for order in permutations(range(3)):
    node = np.zeros(3, dtype=int)
    cell = [0]
    for axis in order:
        node[axis] = 1
        cell.append(int(node @ [4, 2, 1]))
    cells.append(cell)
mesh = Mesh(points=points, tetra4=cells, triangle3=[[1, 3, 7], [1, 5, 7]])

stiffness = np.diag([2e9, 2e9, 1e9, 5e8, 5e8, 5e8])
coupling = np.zeros((3, 6))
coupling[2, 2] = 5.0
material = PiezoelectricMaterial(stiffness, coupling, np.diag([1e-8] * 3))
dofs = piezoelectric_dofs(mesh)
nodes = np.arange(mesh.n_points)
bottom = np.flatnonzero(points[:, 2] == 0)
top = np.flatnonzero(points[:, 2] == length)
supports = [
    Dirichlet(dofs.indices("u", nodes, ["x", "y"]), np.zeros(2 * len(nodes))),
    Dirichlet(dofs.indices("u", bottom, ["z"]), np.zeros(len(bottom))),
    Dirichlet(dofs.indices("phi", bottom), np.zeros(len(bottom))),
]

direct = solve_piezoelectric(
    mesh,
    material,
    [*supports, Dirichlet(dofs.indices("phi", top), np.full(len(top), 100.0))],
)
inverse = solve_piezoelectric(
    mesh,
    material,
    [*supports, Dirichlet(dofs.indices("phi", top), np.zeros(len(top)))],
    tractions={0: [0, 0, 10.0 / area], 1: [0, 0, 10.0 / area]},
)
direct_charge = direct.electrode_charges[top].sum()
inverse_charge = inverse.electrode_charges[top].sum()
np.testing.assert_allclose(direct.displacement[top, 2], -5e-7, rtol=2e-11)
np.testing.assert_allclose(direct_charge, 17.5e-9, rtol=2e-11)
np.testing.assert_allclose(inverse.displacement[top, 2], 2e-6, rtol=2e-11)
np.testing.assert_allclose(inverse_charge, -50e-9, rtol=2e-11)
np.testing.assert_allclose(
    direct.displacement[top, 2].mean() / 100, inverse_charge / 10, rtol=2e-11
)
np.testing.assert_allclose(direct.force_balance, 0, atol=1e-11)
np.testing.assert_allclose(inverse.force_balance, 0, atol=1e-11)
print(
    f"Direct: u_top={direct.displacement[top, 2].mean() * 1e6:.3f} µm, "
    f"Q_top={direct_charge * 1e9:.3f} nC"
)
print(
    f"Inverse: u_top={inverse.displacement[top, 2].mean() * 1e6:.3f} µm, "
    f"Q_top={inverse_charge * 1e9:.3f} nC"
)
write_vtu(
    "piezo_d33.vtu",
    mesh,
    point_data={
        "direct_displacement_m": direct.displacement,
        "direct_potential_V": direct.potential,
        "inverse_displacement_m": inverse.displacement,
        "inverse_potential_V": inverse.potential,
    },
)
