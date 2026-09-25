"""Optional locally; CI requires a usable PETSc runtime, never a skip."""

import os

import pytest


@pytest.fixture
def petsc_backend():
    from fem_toolbox.errors import ValidationError
    from fem_toolbox.solvers.petsc_backend import PetscBackend

    try:
        return PetscBackend()
    except ValidationError as error:
        if os.environ.get("FEM_REQUIRE_PETSC") == "1":
            pytest.fail(str(error))
        pytest.skip(f"optional PETSc runtime unavailable: {error}")
