"""Run the existing independent four-mesh reference with the iterative backend."""

import runpy
from pathlib import Path


def test_iterative_heat_convergence(cube_mesh, petsc_backend):
    from fem_toolbox.solvers.petsc_backend import PetscFieldSplitBackend

    reference = Path(__file__).parents[1] / "convergence" / "test_heat_convergence.py"
    case = runpy.run_path(str(reference))["test_heat_l2_and_energy_convergence"]
    case(cube_mesh, backend=PetscFieldSplitBackend())
