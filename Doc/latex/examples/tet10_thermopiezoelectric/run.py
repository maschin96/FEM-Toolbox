"""Verify a thermally induced open-circuit voltage on a straight Tet10 mesh."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io import export_result_bundle
from fem_toolbox.model import Model


def run(output):
    model = Model.from_file(Path(__file__).resolve().parent / "model.toml")
    solved = model.solve()
    heat, piezo = solved.field_result.heat, solved.field_result.piezoelectric
    x = model.mesh.points[:, 0]
    np.testing.assert_allclose(heat.temperature, 1 + x, rtol=0, atol=4e-12)
    np.testing.assert_allclose(piezo.potential, -0.01 / 3 * x * x, rtol=0, atol=2e-13)
    np.testing.assert_allclose(piezo.displacement, 0, rtol=0, atol=2e-13)
    np.testing.assert_allclose(piezo.electric_displacement, 0, rtol=0, atol=2e-12)
    stress = np.zeros_like(piezo.stress)
    stress[:, 0] = -76 / 75 * piezo.field_points[:, 0]
    np.testing.assert_allclose(piezo.stress, stress, rtol=0, atol=2e-11)
    assert abs(piezo.total_enthalpy - (0.005 / 3 + 0.0002 / 9)) < 2e-13
    assert abs(heat.power_balance) < 2e-11
    assert abs(piezo.charge_balance) < 2e-11
    np.testing.assert_allclose(piezo.force_balance, 0, rtol=0, atol=2e-11)
    np.testing.assert_allclose(piezo.moment_balance, 0, rtol=0, atol=2e-11)
    output.mkdir(parents=True, exist_ok=False)
    export_result_bundle(output / "coupled", model.mesh, solved)
    print(f"Verified Tet10 thermopiezoelectricity: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("tet10_thermopiezo_results")
    )
    run(parser.parse_args().output)
