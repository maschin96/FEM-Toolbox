"""Import prepared SI files and report geometry/BC gates; no CFD solve."""

import argparse
import json
from pathlib import Path

from fem_toolbox.cfd import FluidDomain


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(__file__).with_name("domain.json")
    )
    parser.add_argument("--output", type=Path, default=Path("fluid-preflight.json"))
    args = parser.parse_args()
    domain = FluidDomain.from_file(args.input)
    mesh, walls = domain.mesh, domain.wall_quality
    report = {
        "cells": len(mesh.quad4),
        "boundary_faces": len(mesh.line2),
        "minimum_cell_area_m2": float(mesh.cell_areas_m2.min()),
        "maximum_edge_aspect_ratio": float(mesh.cell_aspect_ratio.max()),
        "maximum_nonorthogonality_deg": mesh.max_nonorthogonality_deg,
        "maximum_skewness": mesh.max_skewness,
        "minimum_wall_distance_m": float(walls.wall_distance_m.min()),
        "maximum_layer_growth": float(walls.layer_growth.max()),
        "inlet_mass_flow_kg_s": domain.inlet_mass_flow_kg_s,
        "wall_resolution": "not_checked_without_wall_shear",
        "solver_status": "not_run",
    }
    # Explicitly do not replace an existing result report.
    with args.output.open("x", encoding="utf-8") as stream:
        stream.write(
            json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"
        )
    print(json.dumps(report, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
