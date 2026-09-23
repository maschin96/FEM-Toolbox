"""Run the synthetic scalar piezo benchmark through Model and ParaView export."""

from __future__ import annotations

import argparse
from pathlib import Path

from fem_toolbox.io.piezo_series import export_piezo_series
from fem_toolbox.model import Model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("piezo_dynamic_results"))
    output = parser.parse_args().output
    model = Model.from_file(Path(__file__).with_name("model.toml"))
    probes = {"tip": 1, "center": [0.25, 0.25, 0.25]}
    transient = model.solve()
    export_piezo_series(output / "transient", transient, probes=probes)
    document = model.to_document()
    document["physics"] = "harmonic_piezoelectric"
    document["analysis"] = {
        "frequencies": {"unit": "Hz", "values": [0, 0.2, 0.4, 0.6, 0.8, 1]},
        "load_factor": {"real": 1.0, "imag": 0.5},
    }
    harmonic = Model(mesh=model.mesh, document=document).solve()
    export_piezo_series(output / "harmonic", harmonic, probes=probes)
    document["physics"] = "piezo_modes"
    document["base"]["loads"] = {}
    document["analysis"] = {"count": 1}
    modes = Model(mesh=model.mesh, document=document).solve()
    export_piezo_series(output / "modes", modes, probes=probes)
    print(f"Transient, harmonic and modal fields/probes: {output.resolve()}")


if __name__ == "__main__":
    main()
