"""Render existing stationary Model results; no graphical session required."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.model import Model
from fem_toolbox.post.visualization import plot_field, save_jpeg
from fem_toolbox.solvers.thermopiezoelectric import ThermopiezoelectricResult


def run(output: Path) -> None:
    model = Model.from_file(Path(__file__).with_name("case.toml"))
    result = model.solve().field_result
    assert isinstance(result, ThermopiezoelectricResult)
    output.mkdir(parents=True, exist_ok=True)
    for filename, values, association, name, unit, component, nodes, limits in (
        (
            "temperature",
            result.heat.temperature,
            "node",
            "T",
            "K",
            None,
            result.heat.node_ids,
            (300, 310),
        ),
        (
            "displacement",
            result.piezoelectric.displacement,
            "node",
            "u",
            "m",
            "magnitude",
            result.piezoelectric.node_ids,
            None,
        ),
        ("heat_flux", result.heat.heat_flux, "cell", "q", "W/m²", "x", None, (-25, 0)),
        (
            "potential",
            result.piezoelectric.potential,
            "node",
            "phi",
            "V",
            None,
            result.piezoelectric.node_ids,
            None,
        ),
        (
            "stress",
            result.piezoelectric.stress,
            "cell",
            "sigma",
            "Pa",
            "xx",
            None,
            None,
        ),
    ):
        plot = plot_field(
            model.mesh,
            values,
            association=association,
            name=name,
            unit=unit,
            component=component,
            nodes=nodes,
            limits=limits,
        )
        try:
            save_jpeg(plot, output / f"{filename}.jpg")
        finally:
            plot.close()
    # A second independently solved temperature case with the same color limits.
    document = model.to_document()
    document["boundary"]["dirichlet"][0]["value"] = [300.0, 305.0, 300.0, 300.0]
    second = Model(mesh=model.mesh, document=document).solve().field_result
    assert isinstance(second, ThermopiezoelectricResult)
    plot = plot_field(
        model.mesh,
        second.heat.temperature,
        association="node",
        nodes=second.heat.node_ids,
        name="T",
        unit="K",
        limits=(300, 310),
    )
    try:
        save_jpeg(plot, output / "temperature_second.jpg")
    finally:
        plot.close()
    np.testing.assert_allclose(result.heat.heat_flux[0, 0], -20, atol=1e-12, rtol=0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("field_images"))
    run(parser.parse_args().output)
