"""Experimental CFD operators, separate from the stable FEM Model API."""

from fem_toolbox.cfd.channel import ChannelResult, solve_channel

__all__ = ["ChannelResult", "solve_channel"]
