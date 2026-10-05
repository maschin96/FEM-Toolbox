"""Incremental ParaView PVD/VTU heat series, probe CSV and balance manifest."""

from __future__ import annotations

import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from femantis.errors import ValidationError
from femantis.io.meshio_adapter import write_vtu
from femantis.mesh import Mesh
from femantis.solvers.transient_heat import (
    HeatFrame,
    TransientHeatResult,
    heat_mesh_fingerprint,
)


class HeatSeriesWriter:
    """Write into a new directory; existing paths are never overwritten.

    Pass write as solve_transient_heat(on_frame=writer.write), then call
    finish(result). Failures leave a readable series marked incomplete.
    Probe/balance CSV and final provenance are published only on finish.
    Inactive mesh nodes/triangle fields are zero-filled with validity masks.
    """

    def __init__(self, path: str | Path, mesh: Mesh) -> None:
        self.path = Path(path)
        self.mesh = mesh
        self.mesh_sha256 = heat_mesh_fingerprint(mesh)
        self.active = np.unique(mesh.cells)
        self._frames: list[dict[str, Any]] = []
        self._complete = False
        self.path.mkdir(parents=True, exist_ok=False)
        self._publish()

    def _publish(self, result: TransientHeatResult | None = None) -> None:
        root = ET.Element(
            "VTKFile", type="Collection", version="0.1", byte_order="LittleEndian"
        )
        collection = ET.SubElement(root, "Collection")
        for frame in self._frames:
            ET.SubElement(
                collection,
                "DataSet",
                timestep=repr(frame["time_s"]),
                group="",
                part="0",
                file=frame["file"],
            )
        temporary = self.path / "heat.pvd.tmp"
        ET.ElementTree(root).write(temporary, encoding="utf-8", xml_declaration=True)
        temporary.replace(self.path / "heat.pvd")
        manifest: dict[str, Any] = {
            "schema": "fem-toolbox-heat-series",
            "schema_version": 2 if self.mesh.n_tetra10 else 1,
            "mesh_sha256": self.mesh_sha256,
            "status": "complete" if result is not None else "incomplete",
            "frames": self._frames,
            "node_ids": self.active.tolist(),
            "units": {
                "time": "s",
                "temperature": "K",
                "heat_flux": "W/m^2",
                "gradient": "K/m",
                "reactions": "W",
                "residual": "W",
            },
            "initial_reactions": "not_evaluated_zero_filled",
            "inactive_values": "zero_filled_use_validity_masks",
        }
        if self.mesh.n_tetra10:
            manifest.update(
                cell_type=self.mesh.cell_type,
                field_location="cell_centroid",
                field_evaluation_points_m=self.mesh.points[self.mesh.cells[:, :4]]
                .mean(axis=1)
                .tolist(),
                cell_field_names=[
                    "gradient_K_m_at_centroid",
                    "heat_flux_W_m2_at_centroid",
                ],
            )
        if result is not None:
            manifest.update(
                provenance=result.provenance.to_dict(),
                times_s=result.times.tolist(),
                probes=list(result.probe_names),
                probe_file="probes.csv",
                balance_file="balance.csv",
            )
        temporary = self.path / "manifest.json.tmp"
        temporary.write_text(
            json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        temporary.replace(self.path / "manifest.json")

    def write(self, frame: HeatFrame) -> None:
        """Publish one complete frame before adding its PVD/manifest entry."""
        if self._complete or (
            self._frames
            and (
                frame.time <= self._frames[-1]["time_s"]
                or frame.index <= self._frames[-1]["index"]
            )
        ):
            raise ValidationError(
                "series frames must increase and remain open", field="export"
            )
        if frame.mesh_sha256 != self.mesh_sha256:
            raise ValidationError("frame belongs to a different mesh", field="export")
        if frame.temperature.shape != self.active.shape:
            raise ValidationError("frame node count differs from mesh", field="export")
        if self.mesh.n_tetra10:
            expected = self.mesh.points[self.mesh.cells[:, :4]].mean(axis=1)
            if frame.field_points is None or not np.array_equal(
                frame.field_points, expected
            ):
                raise ValidationError(
                    "Tet10 frame centroid points differ from mesh", field="export"
                )
        point_data: dict[str, Any] = {}
        for name, values in (
            ("temperature_K", frame.temperature),
            ("reactions_W", frame.reactions),
            ("residual_W", frame.residual),
        ):
            field = np.zeros(self.mesh.n_points)
            field[self.active] = values
            point_data[name] = field
        valid = np.zeros(self.mesh.n_points, dtype=np.uint8)
        valid[self.active] = 1
        point_data["active_temperature"] = valid
        cell_data: dict[str, Any] = {}
        for name, values in (
            ("gradient_K_m", frame.gradient),
            ("heat_flux_W_m2", frame.heat_flux),
        ):
            if self.mesh.n_tetra10:
                name += "_at_centroid"
            blocks = {self.mesh.cell_type: values}
            if self.mesh.n_facets:
                blocks[self.mesh.facet_type] = np.zeros((self.mesh.n_facets, 3))
            cell_data[name] = blocks
        mask = {self.mesh.cell_type: np.ones(self.mesh.n_cells, dtype=np.uint8)}
        if self.mesh.n_facets:
            mask[self.mesh.facet_type] = np.zeros(self.mesh.n_facets, dtype=np.uint8)
        cell_data["valid_volume_field"] = mask
        name = f"heat_{frame.index:06d}.vtu"
        temporary = self.path / f".{name}"
        write_vtu(temporary, self.mesh, point_data=point_data, cell_data=cell_data)
        temporary.replace(self.path / name)
        self._frames.append(
            {
                "index": frame.index,
                "time_s": frame.time,
                "file": name,
                "sha256": hashlib.sha256((self.path / name).read_bytes()).hexdigest(),
            }
        )
        self._publish()

    def finish(self, result: TransientHeatResult) -> Path:
        """Mark a matching completed run after writing all probe/balance records."""
        expected = [(frame.index, frame.time) for frame in result.frames]
        actual = [(frame["index"], frame["time_s"]) for frame in self._frames]
        if (
            self._complete
            or result.mesh_sha256 != self.mesh_sha256
            or expected != actual
            or not np.array_equal(result.node_ids, self.active)
        ):
            raise ValidationError("result does not match open series", field="export")
        with (self.path / "probes.csv").open(
            "w", encoding="utf-8", newline=""
        ) as stream:
            writer = csv.writer(stream)
            writer.writerow(["time_s", *(f"{name}_K" for name in result.probe_names)])
            for time, values in zip(
                result.times, result.probe_temperature, strict=True
            ):
                writer.writerow([time, *values])
        with (self.path / "balance.csv").open(
            "w", encoding="utf-8", newline=""
        ) as stream:
            fields = list(result.balances[0].__dataclass_fields__)
            balance_writer = csv.DictWriter(stream, fieldnames=fields)
            balance_writer.writeheader()
            balance_writer.writerows(asdict(balance) for balance in result.balances)
        self._publish(result)
        self._complete = True
        return self.path


def export_heat_series(
    path: str | Path, mesh: Mesh, result: TransientHeatResult
) -> Path:
    """Export completed fields and measurements into a new directory."""
    if result.mesh_sha256 != heat_mesh_fingerprint(mesh):
        raise ValidationError("result belongs to a different mesh", field="export")
    writer = HeatSeriesWriter(path, mesh)
    for frame in result.frames:
        writer.write(frame)
    return writer.finish(result)
