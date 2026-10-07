"""Run a sourced aluminium heat model with measured source and ambient step."""

from pathlib import Path

from fem_toolbox.io.heat_series import export_heat_series
from fem_toolbox.model import Model
from fem_toolbox.solvers.transient_heat import TransientHeatResult


def main() -> None:
    model = Model.from_file(Path(__file__).with_name("model.toml"))
    result = model.solve().field_result
    assert isinstance(result, TransientHeatResult)
    path = export_heat_series(Path.cwd() / "transient_heat_results", model.mesh, result)
    print(f"ParaView: {path / 'heat.pvd'}")
    print(f"Probes: {path / 'probes.csv'}")
    print(
        f"Maximum balance defect: {max(abs(b.balance) for b in result.balances):.3e} W"
    )


if __name__ == "__main__":
    main()
