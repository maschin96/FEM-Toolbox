"""Select Ubuntu's packaged scalar build before any PETSc import and require it."""

import argparse
import json
import os
import platform
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--scalar", choices=("real", "complex"), required=True)
args, pytest_args = parser.parse_known_args()
paths = list(
    Path("/usr/lib/petscdir").glob(
        f"petsc*/**/*-{args.scalar}/lib/python3/dist-packages"
    )
)
if len(paths) != 1:
    raise SystemExit(f"Expected one {args.scalar} PETSc runtime, found {paths}")
sys.path.insert(0, str(paths[0]))
sys.path.insert(0, str(Path("src").resolve()))
os.environ["FEM_REQUIRE_PETSC"] = "1"

import numpy as np  # noqa: E402
import petsc4py  # noqa: E402
import pytest  # noqa: E402
from petsc4py import PETSc  # noqa: E402

actual = "complex" if np.dtype(PETSc.ScalarType).kind == "c" else "real"
if actual != args.scalar:
    raise SystemExit(f"Expected {args.scalar}, got {actual}")
print(
    json.dumps(
        {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "petsc": PETSc.Sys.getVersion(),
            "petsc4py": petsc4py.__version__,
            "scalar": str(np.dtype(PETSc.ScalarType)),
            "integer": str(np.dtype(PETSc.IntType)),
            "numpy": np.__version__,
        }
    ),
    flush=True,
)
raise SystemExit(pytest.main(pytest_args or ["-q", "tests/petsc"]))
