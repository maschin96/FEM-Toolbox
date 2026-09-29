"""Train, verify and reload a surrogate from a real stationary FEM study."""

import argparse
import json
from pathlib import Path

import numpy as np

from fem_toolbox.studies import Study
from fem_toolbox.study_reports import Quantity, export_study_table
from fem_toolbox.surrogates import (
    DataSplit,
    PolynomialSurrogate,
    SurrogateData,
    fit_surrogate,
)

# Declared before training: numerical roundoff tolerances, not tuned from errors.
TEMPERATURE_ATOL_K = 1e-10
POWER_ATOL_W = 1e-12


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("surrogate_heat_results"))
    args = parser.parse_args()
    root = args.output
    root.mkdir(parents=True, exist_ok=False)
    Study.from_file(Path(__file__).with_name("study.toml")).run(root / "study")
    export_study_table(
        root / "study/catalog.json",
        root / "table",
        [
            Quantity(
                "temperature",
                "point",
                "temperature_K",
                "K",
                "maximum active-node temperature",
                reduction="max",
            ),
            Quantity(
                "source_power",
                "metric",
                "heat.source_power_W",
                "W",
                "integrated volumetric source",
            ),
            Quantity(
                "balance",
                "metric",
                "heat.power_balance_W",
                "W",
                "signed global power imbalance",
                absolute_tolerance=1e-11,
            ),
        ],
    )
    data = SurrogateData.from_table(
        root / "table",
        parameter_columns=["parameter:source"],
        target_columns=["temperature", "source_power"],
    )
    split = DataSplit(data.case_ids[:6], data.case_ids[6:])
    model = fit_surrogate(data, split, degree=1, ridge=0)
    # Independent unit-tetrahedron algebra: T1=300+Q/(4*k), V=1/6, k=2.
    expected = np.column_stack([300 + data.X[:, 0] / 8, data.X[:, 0] / 6])
    np.testing.assert_allclose(
        data.Y[:, 0], expected[:, 0], atol=TEMPERATURE_ATOL_K, rtol=0
    )
    np.testing.assert_allclose(data.Y[:, 1], expected[:, 1], atol=POWER_ATOL_W, rtol=0)
    model.save(root / "surrogate.json")
    loaded = PolynomialSurrogate.load(root / "surrogate.json")
    heldout = loaded.predict(data.X[6:])
    np.testing.assert_allclose(
        heldout.values[:, 0], data.Y[6:, 0], atol=TEMPERATURE_ATOL_K, rtol=0
    )
    np.testing.assert_allclose(
        heldout.values[:, 1], data.Y[6:, 1], atol=POWER_ATOL_W, rtol=0
    )
    assert not heldout.extrapolated.any()
    (root / "report.json").write_text(
        json.dumps(loaded.report, sort_keys=True, indent=2, allow_nan=False) + "\n"
    )
    print("Held-out FEM error:", loaded.report["test"])
    print("Loaded prediction Q=5 W/m³:", loaded.predict([5.0]).values.tolist())


if __name__ == "__main__":
    main()
