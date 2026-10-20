"""Synthetic input uncertainty on a one-tetrahedron stationary heat model."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose

from fem_toolbox.study_reports import Quantity
from fem_toolbox.uncertainty import (
    Normal,
    SamplingPlan,
    UncertainParameter,
    UncertaintyStudy,
    Uniform,
)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=Path("uncertainty-heat-output"))
args = parser.parse_args()
plan = SamplingPlan(
    (
        UncertainParameter(
            "source",
            "loads.source.value",
            "W/m^3",
            Uniform(4, 12),
            "Synthetic teaching law, no measurement evidence",
        ),
        UncertainParameter(
            "temperature",
            "boundary.dirichlet.0.value",
            "K",
            Normal(300, 1),
            "Synthetic teaching law, no measurement evidence",
        ),
    ),
    size=128,
    seed=64,
)
quantity = Quantity(
    "Tmax", "point", "temperature_K", "K", "maximum nodal temperature", reduction="max"
)
study = UncertaintyStudy(plan, model="model.toml", root=Path(__file__).parent)
# This checks solver roundoff against the exact discrete one-tetrahedron
# response, not continuum discretization error. The latter remains unknown.
report = study.run(
    args.output,
    quantities=[quantity],
    numerical_error_bounds={"Tmax": None},
    sensitivity_center=[8, 300],
    sensitivity_steps=[0.1, 0.01],
)
assert report["complete"]
x = plan.sample()
# Independent discrete equation for this tetrahedron: Tmax = Tboundary + Q/8.
expected = x[:, 1] + x[:, 0] / 8
values = np.array([s["values"]["Tmax"] for s in report["samples"]])
assert_allclose(values, expected, atol=1e-11, rtol=0)
assert_allclose(
    report["sensitivity"]["Tmax"]["gradient"], [1 / 8, 1], atol=2e-10, rtol=0
)
summary = report["summaries"]["Tmax"]
# Exact law: E[Tmax]=301 K, Var[Tmax]=1+Var(Q)/64=13/12 K^2.
# Predeclared six-standard-error diagnostic for the mean (not a confidence
# guarantee for this small sample or a statistical validation of a real model).
assert abs(summary["mean"] - 301) < 6 * np.sqrt((13 / 12) / plan.size)
verification = {
    "sample_count": plan.size,
    "seed": plan.seed,
    "analytic_mean_K": 301.0,
    "analytic_variance_K2": 13 / 12,
    "maximum_discrete_error_K": float(np.max(np.abs(values - expected))),
    "numerical_continuum_error": "not established",
    "experimental_validation": "not performed",
}
(args.output / "verification.json").write_text(
    json.dumps(verification, indent=2) + "\n"
)
print(json.dumps({"verification": verification, "summary": summary}, indent=2))
