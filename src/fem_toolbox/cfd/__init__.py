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
from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import FluidMesh2D, WallQuality
from fem_toolbox.cfd.stokes import StokesResult, solve_stokes

__all__ = [
    "CartesianFluidGrid",
    "ChannelBoundary",
    "ChannelResult",
    "ClosedWalls",
    "FlowProblem",
    "FluidDomain",
    "FluidMesh2D",
    "IncompressibleFluid",
    "StationaryFlowResult",
    "StokesResult",
    "WallQuality",
    "solve_channel",
    "solve_stationary_flow",
    "solve_stokes",
]
