"""Mesh-bound piezo fields/probes, preserving complex frequency-domain data."""

from __future__ import annotations

import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.dynamic_result import PiezoDynamicResult
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.elements.piezoelectric import piezoelectric_response
from fem_toolbox.errors import ValidationError
from fem_toolbox.io.meshio_adapter import write_vtu
from fem_toolbox.model import ModelResult
from fem_toolbox.solvers.piezo_dynamics import PiezoelectricModes
from fem_toolbox.solvers.piezo_transient import TransientPiezoelectricResult
from fem_toolbox.solvers.transient_heat import _probes, heat_mesh_fingerprint


def _frames(
    result: PiezoDynamicResult,
) -> Iterator[tuple[float, float | None, dict[str, Any]]]:
    solution = result.solution
    if isinstance(solution, PiezoelectricModes):
        for i, hz in enumerate(solution.frequencies_hz):
            yield (
                float(i + 1),
                float(hz),
                {"u": solution.displacement[i], "phi": solution.potential[i]},
            )
    elif isinstance(solution, TransientPiezoelectricResult):
        for i, time in enumerate(solution.times):
            yield (
                float(time),
                None,
                {
                    "u": solution.displacement[i],
                    "v": solution.velocity[i],
                    "a": solution.acceleration[i],
                    "phi": solution.potential[i],
                    "mechanical_residual": solution.residual_u[i],
                    "electric_residual": solution.residual_phi[i],
                },
            )
    else:
        for state in solution:
            yield (
                state.frequency_hz,
                state.frequency_hz,
                {
                    "u": state.displacement,
                    "v": state.velocity,
                    "a": state.acceleration,
                    "phi": state.potential,
                    "mechanical_residual": state.residual_u,
                    "electric_residual": state.residual_phi,
                },
            )


def _parts(
    name: str, array: ArrayLike, harmonic: bool
) -> dict[str, NDArray[np.float64]]:
    values = np.asarray(array, dtype=np.complex128)
    parts = {name + "_real": np.asarray(values.real, dtype=np.float64)}
    if harmonic:
        amplitude = np.abs(values)
        parts.update(
            {
                name + "_imag": values.imag,
                name + "_amplitude": amplitude,
                name + "_phase_rad": np.where(amplitude > 0, np.angle(values), 0.0),
                name + "_phase_valid": (amplitude > 0).astype(float),
            }
        )
    return parts


def export_piezo_series(
    path: str | Path,
    result: ModelResult,
    *,
    probes: Mapping[str, int | ArrayLike] | None = None,
) -> Path:
    """Write a new directory with PVD/VTU fields, probe CSV and manifest.

    Requires a schema-3 ModelResult bound to its solved mesh/materials. Probe
    locations are mesh node indices or points (3,) m inside Tet4s. Interpolate
    complex u/v/a/phi before computing amplitude/phase. Cell strain/stress/E/D
    remain unsmoothed Tet4 values. Inactive nodes/facets are zero-filled with
    validity masks. Amplitude and phase are componentwise (peak, radians);
    zero amplitude has phase=0 and phase_valid=0. Real modes/transients only
    emit *_real, never a fictitious phase.

    ParaView's collection coordinate is time_s, frequency_hz or mode_index,
    explicitly distinguished in the manifest. Modal vectors are mass-normalized
    shapes, not physical amplitudes. Existing directories are never overwritten;
    an interrupted export remains marked incomplete with its completed frames.
    """
    if not isinstance(result, ModelResult) or not isinstance(
        result.field_result, PiezoDynamicResult
    ):
        raise ValidationError("expected schema-3 dynamic ModelResult", field="result")
    bundle = result.field_result
    mesh = bundle.mesh
    active = np.unique(mesh.tetra4)
    geometries = [Tet4Geometry(mesh.points[cell]) for cell in mesh.tetra4]
    names, weights = _probes(mesh, active, geometries, probes or {})
    solution = bundle.solution
    harmonic = isinstance(solution, tuple)
    axis = (
        "frequency_hz"
        if harmonic
        else "time_s"
        if isinstance(solution, TransientPiezoelectricResult)
        else "mode_index"
    )
    destination = Path(path)
    destination.mkdir(parents=True, exist_ok=False)
    manifest: dict[str, Any] = {
        "schema": "fem-toolbox-piezo-series",
        "schema_version": 1,
        "status": "incomplete",
        "axis": axis,
        "mesh_sha256": heat_mesh_fingerprint(mesh),
        "frames": [],
        "units": {
            "u": "m",
            "v": "m/s",
            "a": "m/s^2",
            "phi": "V",
            "strain": "1",
            "stress": "Pa",
            "electric_field": "V/m",
            "electric_displacement": "C/m^2",
            "mechanical_residual": "N",
            "electric_residual": "C",
            "phase": "rad",
        },
        "modal_normalization": "u.T M u = 1; shapes, not excitation amplitudes"
        if axis == "mode_index"
        else None,
        "phasor": "Re(x_hat exp(+i omega t)); peak amplitudes" if harmonic else None,
        "inactive_values": "zero-filled; use validity masks",
        "phase_zero_amplitude": "0, invalid",
        "provenance": result.provenance.to_dict(),
        "probes": {
            name: (
                int(location)
                if isinstance(location, (int, np.integer))
                else np.asarray(location).tolist()
            )
            for name, location in (probes or {}).items()
        },
    }

    def publish() -> None:
        root = ET.Element(
            "VTKFile", type="Collection", version="0.1", byte_order="LittleEndian"
        )
        collection = ET.SubElement(root, "Collection")
        for frame in manifest["frames"]:
            ET.SubElement(
                collection,
                "DataSet",
                timestep=repr(frame["axis_value"]),
                group="",
                part="0",
                file=frame["file"],
            )
        temporary = destination / "piezo.pvd.tmp"
        ET.ElementTree(root).write(temporary, encoding="utf-8", xml_declaration=True)
        temporary.replace(destination / "piezo.pvd")
        temporary = destination / "manifest.json.tmp"
        temporary.write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        temporary.replace(destination / "manifest.json")

    publish()
    with (destination / "probes.csv").open("w", newline="", encoding="utf-8") as stream:
        writer: csv.DictWriter[str] | None = None
        for index, (coordinate, hz, fields) in enumerate(_frames(bundle)):
            nodal: dict[str, NDArray[np.float64]] = {}
            valid = np.zeros(mesh.n_points)
            valid[active] = 1
            nodal["active_node"] = valid
            for name, values in fields.items():
                for key, part in _parts(name, values, harmonic).items():
                    expanded = np.zeros((mesh.n_points, *part.shape[1:]))
                    expanded[active] = part
                    nodal[key] = expanded
            cell_values: dict[str, list[Any]] = {
                name: []
                for name in (
                    "strain",
                    "stress",
                    "electric_field",
                    "electric_displacement",
                )
            }
            for cell, geometry, material in zip(
                mesh.tetra4, geometries, bundle.materials, strict=True
            ):
                positions = np.searchsorted(active, cell)
                u = np.asarray(fields["u"][positions], dtype=complex)
                phi = np.asarray(fields["phi"][positions], dtype=complex)
                real = piezoelectric_response(geometry, material, u.real, phi.real)
                imaginary = piezoelectric_response(geometry, material, u.imag, phi.imag)
                for name in cell_values:
                    cell_values[name].append(
                        getattr(real, name) + 1j * getattr(imaginary, name)
                    )
            cell_data: dict[str, dict[str, NDArray[np.float64]]] = {}
            cell_parts = {"tetra_valid": np.ones(mesh.n_tetra4)}
            for name, values in cell_values.items():
                cell_parts.update(_parts(name, np.asarray(values), harmonic))
            for name, values in cell_parts.items():
                cell_data[name] = {"tetra4": values}
                if mesh.n_triangle3:
                    cell_data[name]["triangle3"] = np.zeros(
                        (mesh.n_triangle3, *values.shape[1:])
                    )
            filename = f"frame_{index:06d}.vtu"
            write_vtu(
                destination / filename, mesh, point_data=nodal, cell_data=cell_data
            )
            for probe_index, name in enumerate(names):
                row: dict[str, Any] = {
                    "frame": index,
                    "probe": name,
                    "axis_value": coordinate,
                    "frequency_hz": hz if hz is not None else "",
                    "angular_frequency": 2 * np.pi * hz if hz is not None else "",
                }
                for field in ("u", "v", "a", "phi"):
                    if field not in fields:
                        continue
                    interpolated = weights[probe_index] @ fields[field]
                    entries: list[tuple[str, Any]]
                    if field == "phi":
                        entries = [(field, interpolated)]
                    else:
                        entries = [
                            (field + "_" + component, interpolated[j])
                            for j, component in enumerate("xyz")
                        ]
                    for key, value in entries:
                        row.update(
                            {
                                part: float(scalar)
                                for part, scalar in _parts(key, value, harmonic).items()
                            }
                        )
                if writer is None:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                writer.writerow(row)
            if not names and writer is None:
                writer = csv.DictWriter(
                    stream,
                    fieldnames=[
                        "frame",
                        "probe",
                        "axis_value",
                        "frequency_hz",
                        "angular_frequency",
                    ],
                )
                writer.writeheader()
            manifest["frames"].append(
                {
                    "index": index,
                    "axis_value": coordinate,
                    "frequency_hz": hz,
                    "angular_frequency": 2 * np.pi * hz if hz is not None else None,
                    "file": filename,
                    "sha256": hashlib.sha256(
                        (destination / filename).read_bytes()
                    ).hexdigest(),
                }
            )
            frame_metadata = manifest["frames"][-1]
            if isinstance(solution, (PiezoelectricModes, TransientPiezoelectricResult)):
                frame_metadata["relative_residual_u"] = float(
                    solution.relative_residual_u[index]
                )
                frame_metadata["relative_residual_phi"] = float(
                    solution.relative_residual_phi[index]
                )
                if isinstance(solution, TransientPiezoelectricResult):
                    frame_metadata["total_energy_J"] = float(
                        solution.total_energy[index]
                    )
            else:
                frame_metadata["relative_residual_u"] = solution[
                    index
                ].relative_residual_u
                frame_metadata["relative_residual_phi"] = solution[
                    index
                ].relative_residual_phi
            publish()
    manifest["status"] = "complete"
    manifest["probe_sha256"] = hashlib.sha256(
        (destination / "probes.csv").read_bytes()
    ).hexdigest()
    publish()
    return destination
