"""Verify a mass-normalized Tet10 corner oscillator from Gmsh/TOML."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.model import Model


def run(output):
    model = Model.from_file(Path(__file__).resolve().parent / "model.toml")
    solved = model.solve()
    modes = solved.field_result.solution
    # N1=x(2x-1): integral N1²=1/420 and integral (dN1/dx)²=1/10.
    mass = 6 / 420
    omega = np.sqrt((10 + 4 / 3) / 10 / mass)
    np.testing.assert_allclose(modes.angular_frequencies, [omega], rtol=3e-14)
    np.testing.assert_allclose(
        modes.displacement[0, 1, 0], 1 / np.sqrt(mass), rtol=3e-14
    )
    np.testing.assert_allclose(modes.potential[0, 1], 2 / 3 / np.sqrt(mass), rtol=3e-14)
    export_piezo_series(output, solved, probes={"interior": [0.2, 0.1, 0.1]})
    print(f"Verified Tet10 modal oscillator: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_modes_results"))
    run(parser.parse_args().output)
