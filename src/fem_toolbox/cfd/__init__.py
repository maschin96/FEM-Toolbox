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
from fem_toolbox.cfd.quad_flow import (
    QuadFlowConvergenceError,
    QuadFlowResult,
    solve_quad_flow,
)
from fem_toolbox.cfd.quad_fv import FaceGradient, QuadFVGeometry, StressFlux
from fem_toolbox.cfd.stokes import StokesResult, solve_stokes

__all__ = [
    "CartesianFluidGrid",
    "ChannelBoundary",
    "ChannelResult",
    "ClosedWalls",
    "FaceGradient",
    "FlowProblem",
    "FluidDomain",
    "FluidMesh2D",
    "IncompressibleFluid",
    "QuadFVGeometry",
    "QuadFlowConvergenceError",
    "QuadFlowResult",
    "StationaryFlowResult",
    "StokesResult",
    "StressFlux",
    "WallQuality",
    "solve_channel",
    "solve_quad_flow",
    "solve_stationary_flow",
    "solve_stokes",
]
