"""Solve the prepared SI Quad4 channel and export independent final diagnostics."""

import argparse
import csv
import json
from pathlib import Path

from fem_toolbox.cfd import FluidDomain, solve_quad_flow


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input", type=Path, default=Path(__file__).with_name("domain.json")
    )
    parser.add_argument("--output", type=Path, default=Path("quad_flow_results"))
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"output already exists: {args.output}")
    domain = FluidDomain.from_file(args.input)
    result = solve_quad_flow(domain, velocity_scale_m_s=1, length_scale_m=2)
    report = {
        key: getattr(result, key)
        for key in (
            "converged",
            "iterations",
            "simple_steps",
            "inlet_mass_flow_kg_s",
            "outlet_mass_flow_kg_s",
            "global_mass_balance_kg_s",
            "viscous_dissipation_W",
            "gradient_dissipation_W",
            "viscous_discretization_power_W",
            "convective_power_W",
            "numerical_transport_power_W",
            "pressure_boundary_power_W",
            "pressure_interpolation_power_W",
            "body_power_W",
            "viscous_boundary_power_W",
            "power_balance_W",
        )
    }
    report["global_momentum_balance_N"] = result.global_momentum_balance_N.tolist()
    report["provenance"] = result.provenance.to_dict()
    fields = {
        key: getattr(result, key).tolist()
        for key in (
            "velocity_m_s",
            "pressure_Pa",
            "mass_flux_kg_s",
            "face_velocity_m_s",
            "mass_residual_kg_m3_s",
            "momentum_residual_N_m3",
            "nonorthogonal_correction_residual_N_m3",
            "wall_traction_Pa",
        )
    }
    args.output.mkdir(parents=True, exist_ok=False)
    for name, document in (("report.json", report), ("fields.json", fields)):
        with (args.output / name).open("x", encoding="utf-8") as stream:
            stream.write(
                json.dumps(document, sort_keys=True, indent=2, allow_nan=False) + "\n"
            )
    with (args.output / "residuals.csv").open(
        "x", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "momentum_N_m3",
                "mass_kg_m3_s",
                "global_mass_kg_s",
                "global_momentum_N",
                "power_W",
            ]
        )
        writer.writerows(result.residual_history)
    print(json.dumps(report, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
