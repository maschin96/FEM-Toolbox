"""Versioned raw SI result export for stationary Tet4 analyses."""

from __future__ import annotations

import json
import os
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

import numpy as np
from numpy.typing import NDArray

from fem_toolbox import __version__
from fem_toolbox.errors import ValidationError
from fem_toolbox.io.meshio_adapter import write_vtu
from fem_toolbox.mesh import Mesh
from fem_toolbox.provenance import PROVENANCE_SCHEMA_VERSION, RunProvenance

_BUNDLE_VERSION = 1

_POINT: dict[str, dict[str, tuple[str, str]]] = {
    "heat": {
        "temperature": ("temperature_K", "K"),
        "reactions": ("heat_reactions_W", "W"),
        "residual": ("heat_residual_W", "W"),
    },
    "elasticity": {
        "displacement": ("displacement_m", "m"),
        "reactions": ("mechanical_reactions_N", "N"),
        "applied_loads": ("applied_force_N", "N"),
        "residual": ("mechanical_residual_N", "N"),
    },
    "electrostatics": {
        "potential": ("potential_V", "V"),
        "reactions": ("electric_reactions_C", "C"),
        "electrode_charges": ("electrode_charges_C", "C"),
        "residual": ("electric_residual_C", "C"),
    },
    "piezoelectric": {
        "displacement": ("displacement_m", "m"),
        "potential": ("potential_V", "V"),
        "residual_u": ("mechanical_residual_N", "N"),
        "residual_phi": ("electric_residual_C", "C"),
        "reactions_u": ("mechanical_reactions_N", "N"),
        "reactions_phi": ("electric_reactions_C", "C"),
        "electrode_charges": ("electrode_charges_C", "C"),
        "applied_force": ("applied_force_N", "N"),
        "applied_electric_load": ("applied_electric_load_C", "C"),
    },
}

_CELL: dict[str, dict[str, tuple[str, str]]] = {
    "heat": {
        "gradient": ("temperature_gradient_K_m", "K/m"),
        "heat_flux": ("heat_flux_W_m2", "W/m^2"),
    },
    "elasticity": {
        "strain": ("strain_1", "1"),
        "stress": ("stress_Pa", "Pa"),
        "element_energy": ("mechanical_energy_J", "J"),
    },
    "electrostatics": {
        "electric_field": ("electric_field_V_m", "V/m"),
        "electric_displacement": ("electric_displacement_C_m2", "C/m^2"),
        "element_energy": ("electric_energy_J", "J"),
    },
    "piezoelectric": {
        "strain": ("strain_1", "1"),
        "stress": ("stress_Pa", "Pa"),
        "electric_field": ("electric_field_V_m", "V/m"),
        "electric_displacement": ("electric_displacement_C_m2", "C/m^2"),
        "mechanical_energy": ("mechanical_energy_J", "J"),
        "coupling_energy": ("coupling_energy_J", "J"),
        "dielectric_enthalpy": ("dielectric_enthalpy_J", "J"),
    },
}

_SCALARS: dict[str, dict[str, str]] = {
    "heat": {
        "source_power": "source_power_W",
        "prescribed_outflow": "prescribed_outflow_W",
        "convection_outflow": "convection_outflow_W",
        "dirichlet_inflow": "dirichlet_inflow_W",
        "power_balance": "power_balance_W",
        "relative_residual": "relative_residual",
    },
    "elasticity": {"relative_residual": "relative_residual"},
    "electrostatics": {
        "volume_charge": "volume_charge_C",
        "prescribed_surface_charge": "prescribed_surface_charge_C",
        "dirichlet_outward_charge": "dirichlet_outward_charge_C",
        "charge_balance": "charge_balance_C",
        "relative_residual": "relative_residual",
    },
    "piezoelectric": {
        "charge_balance": "charge_balance_C",
        "relative_residual_u": "relative_residual_u",
        "relative_residual_phi": "relative_residual_phi",
        "mechanical_scale": "mechanical_scale_N_m",
        "electric_scale": "electric_scale_C_V",
    },
}


def _components(result: object) -> tuple[str, dict[str, Any], RunProvenance | None]:
    from fem_toolbox.model import ModelResult
    from fem_toolbox.solvers.elasticity import ElasticResult
    from fem_toolbox.solvers.electrostatics import ElectrostaticResult
    from fem_toolbox.solvers.heat import HeatResult
    from fem_toolbox.solvers.piezoelectric import PiezoelectricResult
    from fem_toolbox.solvers.thermomechanics import ThermomechanicalResult
    from fem_toolbox.solvers.thermopiezoelectric import ThermopiezoelectricResult

    provenance = None
    if isinstance(result, ModelResult):
        provenance = result.provenance
        result = result.field_result
    if isinstance(result, HeatResult):
        return "heat", {"heat": result}, provenance
    if isinstance(result, ElasticResult):
        return "elasticity", {"elasticity": result}, provenance
    if isinstance(result, ElectrostaticResult):
        return "electrostatics", {"electrostatics": result}, provenance
    if isinstance(result, PiezoelectricResult):
        return "piezoelectric", {"piezoelectric": result}, provenance
    if isinstance(result, ThermomechanicalResult):
        return (
            "thermomechanics",
            {"heat": result.heat, "elasticity": result.elasticity},
            provenance,
        )
    if isinstance(result, ThermopiezoelectricResult):
        return (
            "thermopiezoelectric",
            {"heat": result.heat, "piezoelectric": result.piezoelectric},
            provenance,
        )
    raise ValidationError("unsupported stationary result type", field="result")


def _array(value: object, rows: int, name: str) -> NDArray[np.float64]:
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise ValidationError(f"{name} must be numeric", field=name) from error
    if (
        array.ndim not in (1, 2)
        or array.shape[0] != rows
        or not np.all(np.isfinite(array))
    ):
        raise ValidationError(f"{name} must have {rows} finite rows", field=name)
    return array


def _components_for(name: str, values: NDArray[np.float64]) -> list[str]:
    if values.ndim == 1:
        return ["scalar"]
    if values.shape[1] == 6 and name.removesuffix("_at_centroid") in (
        "strain_1",
        "stress_Pa",
    ):
        return ["xx", "yy", "zz", "yz", "xz", "xy"]
    if values.shape[1] == 3:
        return ["x", "y", "z"]
    raise ValidationError(f"unsupported component shape for {name}", field=name)


def _point_field(
    mesh: Mesh, part: Any, attribute: str, name: str
) -> tuple[NDArray[np.float64], list[int]]:
    ids = np.asarray(part.node_ids, dtype=np.int64)
    if (
        ids.ndim != 1
        or len(np.unique(ids)) != len(ids)
        or np.any(ids < 0)
        or np.any(ids >= mesh.n_points)
    ):
        raise ValidationError("invalid active mesh node IDs", field="node_ids")
    values = _array(getattr(part, attribute), len(ids), name)
    padded = np.zeros((mesh.n_points, *values.shape[1:]), dtype=np.float64)
    padded[ids] = values
    return padded, ids.tolist()


def _metrics(kind: str, part: Any) -> dict[str, Any]:
    values: dict[str, Any] = {
        alias: float(getattr(part, attribute))
        for attribute, alias in _SCALARS[kind].items()
    }
    if kind in ("elasticity", "piezoelectric"):
        values["force_balance_N"] = np.asarray(part.force_balance).tolist()
        values["moment_balance_N_m"] = np.asarray(part.moment_balance).tolist()
    if kind in ("elasticity", "electrostatics"):
        values["total_energy_J"] = float(part.energy)
    if kind == "piezoelectric":
        values["total_enthalpy_J"] = float(part.total_enthalpy)
        values["electrodes"] = [
            {
                "name": state.name,
                "potential_V": state.potential,
                "charge_C": state.charge,
                "prescribed_charge_C": state.prescribed_charge,
            }
            for state in part.electrode_states
        ]
    return values


def export_result_bundle(
    path: str | Path, mesh: Mesh, result: object
) -> tuple[Path, Path]:
    """Write `<path>.vtu` and `<path>.json` with raw SI fields and metrics.

    The destination parent must exist; existing output files are never replaced.
    Fields are staged first and rolled back if final publication fails. Inactive
    mesh nodes are zero padded in VTU and listed by active ID in the manifest.
    """
    if not isinstance(mesh, Mesh):
        raise ValidationError("expected Mesh", field="mesh")
    prefix = Path(path)
    if prefix.suffix:
        raise ValidationError(
            "expected a filename prefix without extension", field="path"
        )
    vtu, manifest = prefix.with_suffix(".vtu"), prefix.with_suffix(".json")
    if vtu.exists() or manifest.exists():
        raise ValueError(f"result bundle destination exists: {prefix}")
    if not prefix.parent.is_dir():
        raise ValueError(f"result bundle parent does not exist: {prefix.parent}")
    physics, parts, provenance = _components(result)
    quadratic = bool(mesh.n_tetra10)
    evaluation_points = None
    if quadratic:
        if physics not in (
            "heat",
            "electrostatics",
            "elasticity",
            "thermomechanics",
            "piezoelectric",
        ):
            raise ValidationError(
                "only Tet10 heat/electrostatics/elasticity/thermomechanics/piezoelectric export is supported",
                field="result",
            )
        expected = mesh.points[mesh.cells[:, :4]].mean(axis=1)
        for kind, part in parts.items():
            points = getattr(part, "field_points", None)
            if points is None or not np.array_equal(points, expected):
                raise ValidationError(
                    f"Tet10 {kind} centroid evaluation points differ from mesh",
                    field="result",
                )
        evaluation_points = expected
    point_data: dict[str, NDArray[np.float64]] = {}
    cell_data: dict[str, dict[str, NDArray[np.float64]]] = {}
    point_meta: dict[str, Any] = {}
    cell_meta: dict[str, Any] = {}
    metrics: dict[str, Any] = {}
    solver_provenance: dict[str, Any] = {}
    for kind, part in parts.items():
        solver_provenance[kind] = part.provenance.to_dict()
        metrics[kind] = _metrics(kind, part)
        for attribute, (name, unit) in _POINT[kind].items():
            values, ids = _point_field(mesh, part, attribute, name)
            point_data[name] = values
            point_meta[name] = {
                "unit": unit,
                "kind": "raw",
                "shape": list(values.shape),
                "components": _components_for(name, values),
                "active_node_ids": ids,
                "inactive_padding": 0.0,
            }
        for attribute, (name, unit) in _CELL[kind].items():
            integral = quadratic and attribute in (
                "element_energy",
                "mechanical_energy",
                "coupling_energy",
                "dielectric_enthalpy",
            )
            if quadratic and not integral:
                name += "_at_centroid"
            values = _array(getattr(part, attribute), mesh.n_cells, name)
            cell_data[name] = {mesh.cell_type: values}
            if mesh.n_facets:
                cell_data[name][mesh.facet_type] = np.zeros(
                    (mesh.n_facets, *values.shape[1:])
                )
            cell_meta[name] = {
                "unit": unit,
                "kind": "cell_integral" if integral else "raw",
                "shape": list(values.shape),
                "components": _components_for(name, values),
                "cell_type": mesh.cell_type,
                **(
                    {"location": "cell_volume" if integral else "cell_centroid"}
                    if quadratic
                    else {}
                ),
                f"{mesh.facet_type}_padding": 0.0 if mesh.n_facets else None,
            }
    if provenance is None:
        ordered = {
            kind: value["input_sha256"]
            for kind, value in sorted(solver_provenance.items())
        }
        input_hash = sha256(json.dumps(ordered, sort_keys=True).encode()).hexdigest()
        provenance_data: dict[str, Any] = {
            "input_sha256": input_hash,
            "source": "component solver hashes",
        }
        model_version = None
    else:
        provenance_data = provenance.to_dict()
        model_version = provenance.solver_parameters.get("model_schema_version")
    payload: dict[str, Any] = {
        "bundle_schema_version": 2 if quadratic else _BUNDLE_VERSION,
        "toolbox_version": __version__,
        "provenance_schema_version": PROVENANCE_SCHEMA_VERSION,
        "physics": physics,
        "model_schema_version": model_version,
        "solver_status": "converged",
        "provenance": provenance_data,
        "solver_provenance": solver_provenance,
        "mesh": {
            "points": mesh.n_points,
            f"{mesh.cell_type}_cells": mesh.n_cells,
            f"{mesh.facet_type}_facets": mesh.n_facets,
            f"{mesh.cell_type}_ids": (
                mesh.cell_ids.tolist()
                if mesh.cell_ids is not None
                else list(range(mesh.n_cells))
            ),
            "volume_sets": sorted(mesh.volume_sets),
            "facet_sets": sorted(mesh.facet_sets),
        },
        "fields": {"point": point_meta, mesh.cell_type: cell_meta},
        **(
            {"field_evaluation_points_m": evaluation_points.tolist()}
            if evaluation_points is not None
            else {}
        ),
        "metrics": metrics,
        "limitations": [
            (
                "Stationary linear Tet10 result; no mesh convergence or experimental validation implied"
                if quadratic
                else "Stationary linear Tet4 result; no mesh convergence or experimental validation implied"
            ),
            (
                "Energies are cell integrals; other cell fields are raw centroid samples, not elementwise constant; no nodal averaging"
                if quadratic
                else "Cell fields are raw element values; no nodal averaging"
            ),
            *(
                ["Thermal coupling is one-way; no reciprocal or transient feedback"]
                if physics in ("thermomechanics", "thermopiezoelectric")
                else []
            ),
        ],
        "files": {"fields": vtu.name},
    }
    with TemporaryDirectory(prefix=".fem-bundle-", dir=prefix.parent) as temporary:
        staged_vtu = Path(temporary) / vtu.name
        staged_manifest = Path(temporary) / manifest.name
        write_vtu(staged_vtu, mesh, point_data=point_data, cell_data=cell_data)
        payload["files"]["fields_sha256"] = sha256(staged_vtu.read_bytes()).hexdigest()
        staged_manifest.write_text(
            json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        published = False
        try:
            os.link(staged_vtu, vtu)
            published = True
            os.link(staged_manifest, manifest)
        except OSError:
            if published:
                vtu.unlink()
            raise
    return vtu, manifest
