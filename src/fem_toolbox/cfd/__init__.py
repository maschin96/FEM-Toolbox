"""Experimental CFD operators, separate from the stable FEM Model API."""

from fem_toolbox.cfd.channel import ChannelResult, solve_channel
from fem_toolbox.cfd.flow import StationaryFlowResult, solve_stationary_flow
from fem_toolbox.cfd.flow_config import FlowProblem
from fem_toolbox.cfd.flow_model import (
    CartesianFluidGrid,
    ChannelBoundary,
    ClosedWalls,
    IncompressibleFluid,
)
from fem_toolbox.cfd.stokes import StokesResult, solve_stokes

__all__ = [
    "CartesianFluidGrid",
    "ChannelBoundary",
    "ChannelResult",
    "ClosedWalls",
    "FlowProblem",
    "IncompressibleFluid",
    "StationaryFlowResult",
    "StokesResult",
    "solve_channel",
    "solve_stationary_flow",
    "solve_stokes",
]
