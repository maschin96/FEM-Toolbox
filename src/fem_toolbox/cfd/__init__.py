"""Experimental CFD operators, separate from the stable FEM Model API."""

from fem_toolbox.cfd.channel import ChannelResult, solve_channel
from fem_toolbox.cfd.stokes import StokesResult, solve_stokes

__all__ = ["ChannelResult", "StokesResult", "solve_channel", "solve_stokes"]
