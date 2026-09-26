"""Verify a quadratic electric potential and integrated stored energy from TOML."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io import export_result_bundle
from fem_toolbox.model import Model


def run(output):
    model = Model.from_file(Path(__file__).resolve().parent / "model.toml")
    solved = model.solve()
    result = solved.field_result
    np.testing.assert_allclose(
        result.potential, model.mesh.points[:, 0] ** 2, rtol=0, atol=2e-13
    )
    expected = np.zeros_like(result.electric_field)
    expected[:, 0] = -2 * result.field_points[:, 0]
    np.testing.assert_allclose(result.electric_field, expected, rtol=0, atol=2e-12)
    np.testing.assert_allclose(
        result.electric_displacement, 2e-11 * expected, rtol=0, atol=1e-22
    )
    assert abs(result.energy - 4e-11 / 3) < 2e-24
    assert abs(result.charge_balance) < 2e-23
    output.mkdir(parents=True, exist_ok=False)
    export_result_bundle(output / "electric", model.mesh, solved)
    print(f"Verified Tet10 electrostatics: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_electric_results"))
    run(parser.parse_args().output)
