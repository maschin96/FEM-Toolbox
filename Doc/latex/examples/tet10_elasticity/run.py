"""Verify a quadratic mechanical patch from a Tet10 Gmsh/TOML model."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io import export_result_bundle
from fem_toolbox.model import Model


def run(output):
    model = Model.from_file(Path(__file__).resolve().parent / "model.toml")
    solved = model.solve()
    result = solved.field_result
    expected = np.zeros_like(model.mesh.points)
    expected[:, 0] = 0.01 * model.mesh.points[:, 0] ** 2
    np.testing.assert_allclose(result.displacement, expected, rtol=0, atol=2e-14)
    stress = np.zeros_like(result.stress)
    stress[:, 0] = 2 * result.field_points[:, 0]
    np.testing.assert_allclose(result.stress, stress, rtol=0, atol=2e-11)
    assert abs(result.energy - 0.02 / 3) < 2e-15
    np.testing.assert_allclose(result.force_balance, 0, rtol=0, atol=2e-12)
    np.testing.assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-12)
    output.mkdir(parents=True, exist_ok=False)
    export_result_bundle(output / "solid", model.mesh, solved)
    print(f"Verified Tet10 elasticity: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_elastic_results"))
    run(parser.parse_args().output)
