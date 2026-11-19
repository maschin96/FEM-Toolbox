"""Run a prepared stationary FSI lid case from its strict TOML and export artifacts."""

import argparse
import json
import sys
from pathlib import Path

from fem_toolbox.coupling import FsiCouplingError
from fem_toolbox.coupling.fsi_config import run_fsi_config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["laminar", "sa"])
    parser.add_argument("--config", type=Path, help="override the prepared fsi.toml")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config = args.config or Path(__file__).with_name(args.case) / "fsi.toml"
    try:
        result = run_fsi_config(config, args.output)
    except FsiCouplingError as error:
        print(
            f"FSI stopped ({error.kind}); artifacts in {args.output}", file=sys.stderr
        )
        return 2
    final = result.history[-1]
    print(
        json.dumps(
            {
                "status": "converged",
                "iterations": result.iterations,
                "lid_force_y_N": float(final.fluid_force_N[1]),
                "displacement_residual_m": final.displacement_residual_m,
                "output": str(args.output),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
