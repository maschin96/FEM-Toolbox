"""External mesh-format adapters."""

from pathlib import Path

from fem_toolbox.io.meshio_adapter import read_gmsh, read_vtu, write_vtu
from fem_toolbox.mesh import Mesh


def export_result_bundle(
    path: str | Path, mesh: Mesh, result: object
) -> tuple[Path, Path]:
    """Export raw SI fields and a versioned JSON manifest beside a VTU file."""
    from fem_toolbox.io.result_bundle import export_result_bundle as export

    return export(path, mesh, result)


__all__ = ("export_result_bundle", "read_gmsh", "read_vtu", "write_vtu")
