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
from fem_toolbox.cfd.mesh_motion import (
    MeshMotionError,
    MeshMotionResult,
    move_fluid_domain,
)
from fem_toolbox.cfd.quad_flow import (
    QuadFlowConvergenceError,
    QuadFlowResult,
    solve_quad_flow,
)
from fem_toolbox.cfd.quad_fv import FaceGradient, QuadFVGeometry, StressFlux
from fem_toolbox.cfd.rans import RansConvergenceError, RansResult, solve_rans
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
    "MeshMotionError",
    "MeshMotionResult",
    "QuadFVGeometry",
    "QuadFlowConvergenceError",
    "QuadFlowResult",
    "RansConvergenceError",
    "RansResult",
    "StationaryFlowResult",
    "StokesResult",
    "StressFlux",
    "WallQuality",
    "move_fluid_domain",
    "solve_channel",
    "solve_quad_flow",
    "solve_rans",
    "solve_stationary_flow",
    "solve_stokes",
]
