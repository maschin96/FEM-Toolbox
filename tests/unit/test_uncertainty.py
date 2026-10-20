"""Independent probability laws and derivatives for M9 uncertainty propagation."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.uncertainty import (
    Normal,
    SamplingPlan,
    UncertainParameter,
    Uniform,
    local_sensitivity,
    summarize_samples,
)


def parameter(name, law):
    return UncertainParameter(name, f"loads.{name}.value", "1", law, "synthetic law")


def test_correlated_linear_normal_reference_and_reproducibility():
    # Exact: E[4+2X-Y]=8; Var=4*4+9-4*Cov(X,Y)=13.
    n = 100_000
    plan = SamplingPlan(
        (parameter("x", Normal(1, 2)), parameter("y", Normal(-2, 3))),
        size=n,
        seed=641,
        correlation=((1, 0.5), (0.5, 1)),
    )
    x = plan.sample()
    assert_allclose(x, plan.sample(), atol=0, rtol=0)
    assert not x.flags.writeable
    with pytest.raises(ValueError):
        x.setflags(write=True)
    values = 4 + 2 * x[:, 0] - x[:, 1]
    report = summarize_samples(values, unit="K", numerical_error_bound=0.01)
    assert abs(report["mean"] - 8) < 6 * np.sqrt(13 / n)
    assert abs(report["variance"] - 13) < 6 * 13 * np.sqrt(2 / (n - 1))
    assert report["numerical_error_bound"] == 0.01
    assert report["variance"] != pytest.approx(13 + 0.01**2, abs=1e-9)
    assert report["convergence"][-1]["size"] == n
    assert report["quantile_probabilities"] == [0.025, 0.5, 0.975]
    assert plan.to_document()["seed"] == 641
    assert plan.to_document()["correlation"] == [[1.0, 0.5], [0.5, 1.0]]


def test_uniform_linear_reference_and_constant_output():
    n = 50_000
    plan = SamplingPlan((parameter("x", Uniform(-1, 1)),), size=n, seed=64)
    x = plan.sample()[:, 0]
    assert np.all((x >= -1) & (x <= 1))
    report = summarize_samples(3 + 2 * x, unit="m", numerical_error_bound=None)
    assert abs(report["mean"] - 3) < 6 * np.sqrt(4 / (3 * n))
    # Var((2X)^2)=16*(1/5-1/9), independent uniform fourth moment.
    assert abs(report["variance"] - 4 / 3) < 6 * np.sqrt(16 * (1 / 5 - 1 / 9) / n)
    constant = summarize_samples(np.full(n, 7), unit="K", numerical_error_bound=0)
    assert constant["variance"] == 0
    assert constant["quantiles"] == [7, 7, 7]


def test_sensitivity_is_a_local_derivative_and_reports_steps():
    report = local_sensitivity(
        lambda x: 4 + 2 * x[0] - 3 * x[1],
        center=[1, 2],
        steps=[0.01, 0.02],
    )
    assert_allclose(report["gradient"], [2, -3], atol=1e-12, rtol=0)
    assert report["center"] == [1.0, 2.0]
    assert report["steps"] == [0.01, 0.02]
    assert report["method"] == "central_difference"


@pytest.mark.parametrize(
    "law",
    [
        lambda: Normal(0, 0),
        lambda: Normal(np.nan, 1),
        lambda: Uniform(1, 1),
        lambda: Uniform(2, 1),
    ],
)
def test_invalid_laws(law):
    with pytest.raises(ValueError):
        law()


@pytest.mark.parametrize(
    "kwargs",
    [
        {"size": 1, "seed": 1},
        {"size": 2, "seed": True},
        {"size": 2, "seed": -1},
        {"size": 2.5, "seed": 1},
    ],
)
def test_invalid_sampling_controls(kwargs):
    with pytest.raises(ValueError):
        SamplingPlan((parameter("x", Normal(0, 1)),), **kwargs)


@pytest.mark.parametrize(
    "correlation", [((1, 2), (2, 1)), ((1, 0.1), (0.2, 1)), ((2, 0), (0, 1)), ((1,),)]
)
def test_invalid_correlation(correlation):
    with pytest.raises(ValueError):
        SamplingPlan(
            (parameter("x", Normal(0, 1)), parameter("y", Normal(0, 1))),
            size=10,
            seed=1,
            correlation=correlation,
        )


def test_uniform_dependence_is_rejected_and_singular_normal_is_supported():
    with pytest.raises(ValueError, match="normal"):
        SamplingPlan(
            (parameter("x", Uniform(0, 1)), parameter("y", Normal(0, 1))),
            size=10,
            seed=1,
            correlation=((1, 0.5), (0.5, 1)),
        )
    plan = SamplingPlan(
        (parameter("x", Normal(0, 1)), parameter("y", Normal(0, 1))),
        size=10,
        seed=1,
        correlation=((1, 1), (1, 1)),
    )
    assert_allclose(plan.sample()[:, 0], plan.sample()[:, 1], atol=1e-14, rtol=0)


@pytest.mark.parametrize(
    "values", [[1], [1, np.nan], [True, False], [1 + 1j, 2 + 2j], [[1, 2], [3, 4]]]
)
def test_invalid_output_is_never_silently_dropped(values):
    with pytest.raises(ValueError):
        summarize_samples(values, unit="K", numerical_error_bound=None)


def test_sampling_detaches_correlation_and_preserves_prefix():
    parameters = (parameter("x", Normal(0, 1)), parameter("y", Normal(2, 3)))
    correlation = np.array([[1, 0.4], [0.4, 1]])
    small = SamplingPlan(parameters, size=10, seed=42, correlation=correlation)
    large = SamplingPlan(parameters, size=100, seed=42, correlation=correlation)
    correlation[0, 1] = 99
    assert small.to_document()["correlation"][0][1] == 0.4
    assert_allclose(small.sample(), large.sample()[:10], atol=0, rtol=0)
    document = small.to_document()
    document["parameters"][0]["distribution"]["mean"] = 99
    assert small.to_document()["parameters"][0]["distribution"]["mean"] == 0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"unit": "", "numerical_error_bound": None},
        {"unit": "K", "numerical_error_bound": -1},
        {"unit": "K", "numerical_error_bound": np.nan},
        {"unit": "K", "numerical_error_bound": True},
        {
            "unit": "K",
            "numerical_error_bound": None,
            "quantile_probabilities": [0.5, 0.5],
        },
        {
            "unit": "K",
            "numerical_error_bound": None,
            "quantile_probabilities": [-0.1, 0.5],
        },
    ],
)
def test_invalid_summary_contract(kwargs):
    with pytest.raises(ValueError):
        summarize_samples([1, 2, 3], **kwargs)


def test_invalid_sensitivity_contract_and_callback_failure():
    with pytest.raises(ValueError):
        local_sensitivity(lambda x: x[0], center=[1], steps=[0])
    with pytest.raises(ValueError):
        local_sensitivity(lambda x: np.nan, center=[1], steps=[0.1])
    with pytest.raises(ValueError):
        local_sensitivity(lambda x: x[0], center=[1e100], steps=[1e-100])
