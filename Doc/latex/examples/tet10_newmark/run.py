"""Verify free Tet10 Newmark motion, energy and algebraic potential from TOML."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.model import Model


def run(output):
    model = Model.from_file(Path(__file__).resolve().parent / "model.toml")
    solved = model.solve()
    history = solved.field_result.solution
    # Exact scalar Tet10 integrals: m=rho/420, k=(C+e²/kappa)/10.
    stiffness = (10 + 4 / 3) / 10
    omega = np.sqrt(stiffness / (6 / 420))
    phase = np.arange(len(history.times)) * 2 * np.arctan(omega * 0.01 / 2)
    np.testing.assert_allclose(
        history.displacement[:, 1, 0], 0.01 * np.cos(phase), rtol=3e-12, atol=1e-14
    )
    np.testing.assert_allclose(
        history.total_energy, 0.5 * stiffness * 0.01**2, rtol=3e-12
    )
    np.testing.assert_allclose(
        history.potential[:, 1],
        2 / 3 * history.displacement[:, 1, 0],
        rtol=3e-12,
        atol=1e-14,
    )
    export_piezo_series(output, solved, probes={"interior": [0.2, 0.1, 0.1]})
    print(f"Verified Tet10 Newmark oscillator: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("tet10_newmark_results"))
    run(parser.parse_args().output)
