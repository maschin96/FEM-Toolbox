"""Strict versioned TOML for the supported uniform stationary-flow contract."""

import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from fem_toolbox.cfd.flow import StationaryFlowResult, solve_stationary_flow
from fem_toolbox.cfd.flow_model import (
    CartesianFluidGrid,
    ChannelBoundary,
    ClosedWalls,
    IncompressibleFluid,
)
from fem_toolbox.errors import ModelSchemaError, ValidationError


def _keys(value: object, keys: set[str], path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        raise ModelSchemaError(path=path, reason=f"expected keys {sorted(keys)}")
    return value


@dataclass(frozen=True, slots=True)
class FlowProblem:
    """Validated grid, fluid and boundary data; reproducible zero-body-force run."""

    grid: CartesianFluidGrid
    fluid: IncompressibleFluid
    boundary: ChannelBoundary | ClosedWalls

    def __post_init__(self) -> None:
        if (
            not isinstance(self.grid, CartesianFluidGrid)
            or not isinstance(self.fluid, IncompressibleFluid)
            or not isinstance(self.boundary, (ChannelBoundary, ClosedWalls))
        ):
            raise ValidationError("validated flow objects required", field="model")
        if (
            isinstance(self.boundary, ChannelBoundary)
            and len(self.boundary.inlet_velocity_m_s) != self.grid.cells[1]
        ):
            raise ValidationError("inlet profile must match ny", field="inlet")

    @classmethod
    def from_document(cls, document: Mapping[str, object]) -> "FlowProblem":
        d = _keys(
            document, {"flow_schema_version", "grid", "fluid", "boundary"}, "flow"
        )
        if type(d["flow_schema_version"]) is not int or d["flow_schema_version"] != 1:
            raise ModelSchemaError(
                path="flow_schema_version", reason="expected integer 1"
            )
        g = _keys(d["grid"], {"cells", "size_m"}, "grid")
        f = _keys(d["fluid"], {"density_kg_m3", "viscosity_Pa_s"}, "fluid")
        raw = d["boundary"]
        if not isinstance(raw, Mapping):
            raise ModelSchemaError(path="boundary", reason="expected table")
        if raw.get("type") == "channel":
            b = _keys(
                raw, {"type", "inlet_velocity_m_s", "outlet_pressure_Pa"}, "boundary"
            )
            boundary: ChannelBoundary | ClosedWalls = ChannelBoundary(
                b["inlet_velocity_m_s"], b["outlet_pressure_Pa"]
            )
        elif raw.get("type") == "closed":
            b = _keys(raw, {"type", "tangential_velocities_m_s"}, "boundary")
            boundary = ClosedWalls(b["tangential_velocities_m_s"])
        else:
            raise ModelSchemaError(
                path="boundary.type", reason="expected channel or closed"
            )
        return cls(
            CartesianFluidGrid(g["cells"], g["size_m"]),
            IncompressibleFluid(f["density_kg_m3"], f["viscosity_Pa_s"]),
            boundary,
        )

    @classmethod
    def from_file(cls, path: str | Path) -> "FlowProblem":
        """Read UTF-8 TOML and diagnose unreadable files/syntax as schema errors."""
        try:
            document = tomllib.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise ModelSchemaError(path=str(path), reason=str(error)) from error
        return cls.from_document(document)

    def to_document(self) -> dict[str, Any]:
        """Detached schema-1 data; values and field names explicitly use SI."""
        boundary = (
            {
                "type": "channel",
                "inlet_velocity_m_s": self.boundary.inlet_velocity_m_s.tolist(),
                "outlet_pressure_Pa": self.boundary.outlet_pressure_Pa,
            }
            if isinstance(self.boundary, ChannelBoundary)
            else {
                "type": "closed",
                "tangential_velocities_m_s": list(
                    self.boundary.tangential_velocities_m_s
                ),
            }
        )
        return dict(
            flow_schema_version=1,
            grid=dict(cells=list(self.grid.cells), size_m=list(self.grid.size_m)),
            fluid=dict(
                density_kg_m3=self.fluid.density_kg_m3,
                viscosity_Pa_s=self.fluid.viscosity_Pa_s,
            ),
            boundary=boundary,
        )

    def solve(self) -> StationaryFlowResult:
        """Solve with documented default controls; custom controls use Python API."""
        return solve_stationary_flow(self.grid, self.fluid, self.boundary)
