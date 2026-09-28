"""Verify a complex Tet10 frequency sweep with an independent corner oscillator."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.model import Model


def run(output):
    model = Model.from_file(Path(__file__).resolve().parent / "model.toml")
    solved = model.solve()
    for state in solved.field_result.solution:
        # int N_corner=-1/120; m=rho/420; k=(C+e²/kappa)/10.
        expected = (
            -0.05
            * (1 + 2j)
            / ((10 + 4 / 3) / 10 - 6 / 420 * state.angular_frequency**2)
        )
        np.testing.assert_allclose(
            state.displacement[1, 0], expected, rtol=3e-14, atol=0
        )
        np.testing.assert_allclose(
            state.potential[1], 2 / 3 * expected, rtol=3e-14, atol=0
        )
    export_piezo_series(output, solved, probes={"interior": [0.2, 0.1, 0.1]})
    print(f"Verified Tet10 harmonic sweep: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_harmonic_results"))
    run(parser.parse_args().output)
