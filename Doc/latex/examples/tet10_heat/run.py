"""Solve and independently verify stationary and transient Tet10 TOML examples."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io import export_result_bundle
from fem_toolbox.io.heat_series import export_heat_series
from fem_toolbox.model import Model


def run(output):
    root = Path(__file__).resolve().parent
    output.mkdir(parents=True, exist_ok=False)
    model = Model.from_file(root / "stationary.toml")
    solved = model.solve()
    result = solved.field_result
    np.testing.assert_allclose(
        result.temperature, 300 + model.mesh.points[:, 0] ** 2, rtol=0, atol=4e-12
    )
    expected_flux = np.zeros_like(result.heat_flux)
    expected_flux[:, 0] = -4 * result.field_points[:, 0]
    np.testing.assert_allclose(result.heat_flux, expected_flux, rtol=0, atol=5e-12)
    assert abs(result.power_balance) < 2e-11
    export_result_bundle(output / "stationary", model.mesh, solved)
    transient = Model.from_file(root / "transient.toml")
    history = transient.solve().field_result
    for frame in history.frames:
        np.testing.assert_allclose(
            frame.temperature, 300 + frame.time, rtol=0, atol=5e-12
        )
    assert max(abs(b.balance) for b in history.balances) < 2e-10
    export_heat_series(output / "transient", transient.mesh, history)
    print(f"Verified Tet10 heat results: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_heat_results"))
    run(parser.parse_args().output)
