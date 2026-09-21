"""Ideal capacitor: gap 2 mm, area 30x40 mm², absolute epsilon=2e-11 F/m.

100 V gives C=12 pF, Q_high=1.2 nC, U=60 nJ. Lateral D.n=0 models
an ideal parallel-plate field without fringing. Run after package installation.
"""

from itertools import permutations, product

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.io import write_vtu
from fem_toolbox.materials import DielectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.electrostatics import electric_dofs, solve_electrostatics

points = np.array(list(product((0.0, 1.0), repeat=3))) * [0.002, 0.03, 0.04]
cells = []
for order in permutations(range(3)):
    node = np.zeros(3, dtype=int)
    cell = [0]
    for axis in order:
        node[axis] = 1
        cell.append(int(node @ [4, 2, 1]))
    cells.append(cell)
mesh = Mesh(points=points, tetra4=cells)
low = np.flatnonzero(points[:, 0] == 0)
high = np.flatnonzero(points[:, 0] == 0.002)
dofs = electric_dofs(mesh)
result = solve_electrostatics(
    mesh,
    DielectricMaterial.isotropic(2e-11),
    [
        Dirichlet(dofs.indices("phi", low), np.zeros(len(low))),
        Dirichlet(dofs.indices("phi", high), np.full(len(high), 100.0)),
    ],
)
charge = result.electrode_charges[high].sum()
capacitance = charge / 100
np.testing.assert_allclose(capacitance, 12e-12, atol=2e-24, rtol=0)
np.testing.assert_allclose(result.energy, 60e-9, atol=2e-20, rtol=0)
np.testing.assert_allclose(result.charge_balance, 0, atol=2e-22, rtol=0)
print(f"Capacitance: {capacitance * 1e12:.6f} pF")
print(f"High electrode charge: {charge * 1e9:.6f} nC")
print(f"Stored energy: {result.energy * 1e9:.6f} nJ")
write_vtu(
    "parallel_plate.vtu",
    mesh,
    point_data={
        "potential_V": result.potential,
        "electrode_charge_C": result.electrode_charges,
    },
    cell_data={
        "E_V_m": {"tetra4": result.electric_field},
        "D_C_m2": {"tetra4": result.electric_displacement},
        "energy_J": {"tetra4": result.element_energy},
    },
)
