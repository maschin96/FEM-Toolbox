import csv
import json
import xml.etree.ElementTree as ET

import meshio
import numpy as np
import pytest

from fem_toolbox.io.heat_series import HeatSeriesWriter, export_heat_series
from fem_toolbox.materials import ThermalMaterial
from fem_toolbox.solvers.transient_heat import solve_transient_heat


def test_paraview_fields_and_probe_csv(tmp_path, cube_mesh):
    mesh = cube_mesh(1)
    result = solve_transient_heat(
        mesh,
        ThermalMaterial.isotropic(1, density=1, specific_heat=1),
        [0, 0.5, 1],
        300,
        source=2,
        probes={"center": [0.5, 0.5, 0.5]},
        output_every=2,
    )
    path = export_heat_series(tmp_path / "series", mesh, result)
    manifest = json.loads((path / "manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["schema_version"] == 1
    datasets = ET.parse(path / "heat.pvd").findall(".//DataSet")
    assert [float(d.attrib["timestep"]) for d in datasets] == [0.0, 1.0]
    last = meshio.read(path / datasets[-1].attrib["file"])
    np.testing.assert_allclose(last.point_data["temperature_K"], 302.0)
    assert "heat_flux_W_m2" in last.cell_data
    with (path / "probes.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 3
    assert float(rows[-1]["center_K"]) == pytest.approx(302.0)
    assert len(list(csv.DictReader((path / "balance.csv").open()))) == 2
    with pytest.raises(FileExistsError):
        export_heat_series(path, mesh, result)


def test_interrupted_stream_is_explicitly_incomplete(tmp_path, cube_mesh):
    mesh = cube_mesh(1)
    writer = HeatSeriesWriter(tmp_path / "partial", mesh)

    def source(x, t):
        if t > 0.5:
            raise RuntimeError("simulated interrupted input")
        return 1.0

    with pytest.raises(RuntimeError):
        solve_transient_heat(
            mesh,
            ThermalMaterial.isotropic(1, density=1, specific_heat=1),
            [0, 0.5, 1],
            300,
            source=source,
            on_frame=writer.write,
        )
    manifest = json.loads((writer.path / "manifest.json").read_text())
    assert manifest["status"] == "incomplete"
    assert [f["time_s"] for f in manifest["frames"]] == [0.0, 0.5]
    for dataset in ET.parse(writer.path / "heat.pvd").findall(".//DataSet"):
        meshio.read(writer.path / dataset.attrib["file"])


def test_export_rejects_different_geometry(tmp_path, cube_mesh):
    from fem_toolbox.errors import ValidationError
    from fem_toolbox.mesh import Mesh

    mesh = cube_mesh(1)
    result = solve_transient_heat(
        mesh, ThermalMaterial.isotropic(1, density=1, specific_heat=1), [0, 1], 300
    )
    other = Mesh(points=mesh.points * 2, tetra4=mesh.tetra4, triangle3=mesh.triangle3)
    with pytest.raises(ValidationError, match="mesh"):
        export_heat_series(tmp_path / "wrong", other, result)
