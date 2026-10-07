"""Run the declarative one-Tet4 thermal/piezoelectric reference case."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fem_toolbox.io import write_vtu
from fem_toolbox.model import Model
from fem_toolbox.solvers.thermopiezoelectric import ThermopiezoelectricResult


def run(output: Path) -> tuple[Path, Path]:
    """Write raw SI fields and a deterministic, example-scoped JSON report."""
    model = Model.from_file(Path(__file__).with_name("case.toml"))
    solved = model.solve()
    result = solved.field_result
    if not isinstance(result, ThermopiezoelectricResult):
        raise RuntimeError("example requires thermopiezoelectric physics")
    mesh = model.mesh
    heat = result.heat
    piezo = result.piezoelectric
    active = np.arange(mesh.n_points)
    if not np.array_equal(heat.node_ids, active) or not np.array_equal(
        piezo.node_ids, active
    ):
        raise RuntimeError("example mesh must use every node")

    output.mkdir(parents=True, exist_ok=True)
    fields = output / "fields.vtu"
    report = output / "report.json"
    write_vtu(
        fields,
        mesh,
        point_data={
            "temperature_K": heat.temperature,
            "displacement_m": piezo.displacement,
            "potential_V": piezo.potential,
        },
        cell_data={
            "heat_flux_W_m2": {"tetra4": heat.heat_flux},
            "stress_Pa": {"tetra4": piezo.stress},
            "electric_displacement_C_m2": {"tetra4": piezo.electric_displacement},
        },
    )
    payload = {
        "schema_version": 1,
        "report_kind": "thermopiezoelectric-example-v1",
        "limitations": [
            "Synthetic SI parameters; no experimental material validation",
            "One Tet4 and prescribed nodal temperatures; no mesh convergence claim",
            "One-way stationary heat transfer; no reciprocal or transient coupling",
        ],
        "model": {
            "schema_version": model.to_document()["schema_version"],
            "physics": solved.physics,
            "mesh_points": mesh.n_points,
            "tetra4_cells": mesh.n_tetra4,
        },
        "provenance": solved.provenance.to_dict(),
        "solver_provenance": {
            "heat": heat.provenance.to_dict(),
            "piezoelectric": piezo.provenance.to_dict(),
        },
        "fields": {
            "point": {
                "temperature_K": "raw absolute temperature",
                "displacement_m": "raw Cartesian displacement",
                "potential_V": "raw electric potential",
            },
            "tetra4": {
                "heat_flux_W_m2": "raw cell heat flux",
                "stress_Pa": "raw cell engineering-Voigt stress [xx, yy, zz, yz, xz, xy]",
                "electric_displacement_C_m2": "raw cell electric displacement",
            },
        },
        "metrics": {
            "free_displacement_x_m": float(piezo.displacement[1, 0]),
            "free_potential_V": float(piezo.potential[1]),
            "heat_flux_x_W_m2": float(heat.heat_flux[0, 0]),
            "total_electric_enthalpy_J": piezo.total_enthalpy,
        },
        "balances": {
            "heat_power_W": heat.power_balance,
            "force_N": piezo.force_balance.tolist(),
            "moment_N_m": piezo.moment_balance.tolist(),
            "charge_C": piezo.charge_balance,
        },
        "relative_residuals": {
            "heat": heat.relative_residual,
            "mechanical": piezo.relative_residual_u,
            "electric": piezo.relative_residual_phi,
        },
        "files": {"fields": fields.name},
    }
    report.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return fields, report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="result directory")
    args = parser.parse_args()
    fields, report = run(args.output)
    print(f"{fields}\n{report}")
