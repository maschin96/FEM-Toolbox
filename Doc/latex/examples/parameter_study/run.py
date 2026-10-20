"""Run or resume a two-case heat study; verify an independent Tet4 reference."""

import argparse
import json
from math import isclose
from pathlib import Path

import meshio

from fem_toolbox.studies import Study


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("parameter_study_results"))
    parser.add_argument("--max-cases", type=int)
    args = parser.parse_args()
    catalog = Study.from_file(Path(__file__).with_name("study.toml")).run(
        args.output, max_cases=args.max_cases
    )
    for case in catalog["cases"]:
        print(case["case_id"], case["status"], f"attempts={case['attempts']}")
        if case["status"] != "succeeded":
            continue
        k = case["parameters"]["conductivity"]["value"]
        fields = meshio.read(args.output / case["result"]["fields"])
        # Unit-axis tetrahedron: K11=k/6, f1=Q/24, T1=300+Q/(4k).
        assert isclose(
            fields.point_data["temperature_K"][1], 300 + 2 / k, rel_tol=0, abs_tol=1e-11
        )
        bundle = json.loads((args.output / case["result"]["manifest"]).read_text())
        assert abs(bundle["metrics"]["heat"]["power_balance_W"]) < 1e-11


if __name__ == "__main__":
    main()
