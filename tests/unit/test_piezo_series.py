"""Round-trip full fields and probes, retaining complex phasors in ParaView."""

import csv
import json
import xml.etree.ElementTree as ET

import meshio
import numpy as np
import pytest
from numpy.testing import assert_allclose
from test_piezo_dynamic_model import case

from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.model import Model


@pytest.mark.parametrize(
    "physics", ["piezo_modes", "transient_piezoelectric", "harmonic_piezoelectric"]
)
def test_series_roundtrip_and_probes(tmp_path, physics):
    mesh, document = case(physics)
    result = Model(mesh=mesh, document=document).solve()
    path = tmp_path / "fields"
    export_piezo_series(path, result, probes={"tip": 1, "center": [0.25, 0.25, 0.25]})
    manifest = json.loads((path / "manifest.json").read_text())
    assert manifest["status"] == "complete"
    frames = ET.parse(path / "piezo.pvd").findall(".//DataSet")
    output = meshio.read(path / frames[0].attrib["file"])
    assert "strain_real" in output.cell_data
    assert "stress_real" in output.cell_data
    assert "electric_field_real" in output.cell_data
    with (path / "probes.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 2 * len(frames)
    if physics == "harmonic_piezoelectric":
        solution = result.field_result.solution[0]
        assert manifest["axis"] == "frequency_hz"
        assert_allclose(
            output.point_data["u_imag"], solution.displacement.imag, atol=1e-14
        )
        assert_allclose(
            output.point_data["u_amplitude"], abs(solution.displacement), atol=1e-14
        )
        assert_allclose(
            output.point_data["u_phase_rad"][1, 0],
            np.angle(solution.displacement[1, 0]),
            atol=1e-14,
        )
        center = next(r for r in rows if r["probe"] == "center" and r["frame"] == "0")
        assert_allclose(
            float(center["u_x_real"]), solution.displacement[1, 0].real / 4, atol=1e-14
        )
        assert_allclose(
            float(center["u_x_imag"]), solution.displacement[1, 0].imag / 4, atol=1e-14
        )
        assert_allclose(
            output.cell_data_dict["strain_imag"]["tetra"][0, 0],
            solution.displacement[1, 0].imag,
            atol=1e-14,
        )
    with pytest.raises(FileExistsError):
        export_piezo_series(path, result)


def test_invalid_probe_does_not_publish_directory(tmp_path):
    mesh, document = case()
    result = Model(mesh=mesh, document=document).solve()
    with pytest.raises(ValueError, match="outside"):
        export_piezo_series(tmp_path / "bad", result, probes={"outside": [9, 9, 9]})
    assert not (tmp_path / "bad").exists()


def test_interrupted_export_stays_incomplete(tmp_path, monkeypatch):
    import fem_toolbox.io.piezo_series as adapter

    mesh, document = case()
    result = Model(mesh=mesh, document=document).solve()
    original = adapter.write_vtu
    calls = 0

    def fail_second(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected storage failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(adapter, "write_vtu", fail_second)
    with pytest.raises(OSError, match="storage"):
        export_piezo_series(tmp_path / "partial", result)
    manifest = json.loads((tmp_path / "partial/manifest.json").read_text())
    assert manifest["status"] == "incomplete"
    assert len(manifest["frames"]) == 1


def test_transient_manifest_preserves_energy_and_field_residuals(tmp_path):
    mesh, document = case("transient_piezoelectric")
    result = Model(mesh=mesh, document=document).solve()
    export_piezo_series(tmp_path / "time", result)
    manifest = json.loads((tmp_path / "time/manifest.json").read_text())
    solution = result.field_result.solution
    assert_allclose(
        [f["total_energy_J"] for f in manifest["frames"]],
        solution.total_energy,
        atol=1e-15,
    )
    assert_allclose(
        [f["relative_residual_u"] for f in manifest["frames"]],
        solution.relative_residual_u,
        atol=1e-15,
    )
