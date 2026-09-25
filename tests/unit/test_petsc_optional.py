"""Core and adapter imports must not initialize optional PETSc."""

import os
import subprocess
import sys


def test_core_and_petsc_adapter_import_without_optional_dependency():
    code = """
import sys
class RejectPETSc:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "petsc4py" or fullname.startswith("petsc4py."):
            raise ModuleNotFoundError("test: PETSc unavailable")
sys.meta_path.insert(0, RejectPETSc())
from fem_toolbox.model import Model
from fem_toolbox.solvers.petsc_backend import PetscBackend
from fem_toolbox.errors import ValidationError
assert "petsc4py" not in sys.modules
try:
    PetscBackend()
except ValidationError as error:
    assert "petsc4py" in str(error)
else:
    raise AssertionError("missing optional runtime must be diagnosed")
"""
    subprocess.run(
        [sys.executable, "-c", code],
        check=True,
        env={**os.environ, "PYTHONPATH": "src"},
    )
