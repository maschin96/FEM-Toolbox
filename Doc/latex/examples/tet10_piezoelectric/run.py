"""Verify a Tet10 piezoelectric Gmsh/TOML model with a floating electrode."""

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
    expected[:, 0] = -0.02 * model.mesh.points[:, 0]
    np.testing.assert_allclose(result.displacement, expected, rtol=0, atol=2e-13)
    np.testing.assert_allclose(
        result.potential, model.mesh.points[:, 0], rtol=0, atol=2e-12
    )
    np.testing.assert_allclose(result.stress, 0, rtol=0, atol=2e-11)
    assert abs(result.total_enthalpy + 1.52) < 2e-12
    assert abs(result.electrode_states[1].charge - 3.04) < 2e-11
    assert abs(result.charge_balance) < 2e-11
    np.testing.assert_allclose(result.force_balance, 0, rtol=0, atol=2e-11)
    np.testing.assert_allclose(result.moment_balance, 0, rtol=0, atol=2e-11)
    output.mkdir(parents=True, exist_ok=False)
    export_result_bundle(output / "piezo", model.mesh, solved)
    print(f"Verified Tet10 piezoelectricity: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_piezo_results"))
    run(parser.parse_args().output)
