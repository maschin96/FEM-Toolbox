"""Export a saved heat study without solving again."""

import argparse
from pathlib import Path

from fem_toolbox.study_reports import Quantity, export_study_table, read_study_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, default=Path("parameter_study_results"))
    parser.add_argument("--output", type=Path, default=Path("parameter_study_table"))
    args = parser.parse_args()
    export_study_table(
        args.study / "catalog.json",
        args.output,
        [
            Quantity(
                "temperature",
                "point",
                "temperature_K",
                "K",
                "maximum temperature over active nodes",
                reduction="max",
            ),
            Quantity(
                "power_balance",
                "metric",
                "heat.power_balance_W",
                "W",
                "signed global heat power imbalance",
                absolute_tolerance=1e-11,
            ),
        ],
    )
    for row in read_study_table(args.output)["rows"]:
        print(row["case_id"], row["status"], row["values"])


if __name__ == "__main__":
    main()
