"""Run an installed M11 topology study using its declarative input file."""

import argparse
import json
from pathlib import Path

from fem_toolbox.errors import ValidationError
from fem_toolbox.topology_study import TopologyStudy

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--config", type=Path, default=Path(__file__).parent / "reference/study.toml"
)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
try:
    study = TopologyStudy.from_file(args.config)
    result = study.run(args.output)
except (OSError, ValueError, ValidationError, KeyError, TypeError) as error:
    parser.exit(2, f"topology study: {error}\n")
print(
    json.dumps(
        {
            "status": result.status,
            "compliance_J": result.final.objective,
            "volume_fraction": result.final.volume_fraction,
            "output": str(args.output),
        },
        indent=2,
    )
)
raise SystemExit(0 if result.status == "converged" else 2)
