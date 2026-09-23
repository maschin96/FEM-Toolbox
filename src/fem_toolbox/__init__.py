"""Finite-element research tools for verified 3D multiphysics simulations.

Stable stationary user APIs and their validity limits are documented in the
public API contract. Internal research building blocks may evolve separately.
"""

from importlib.metadata import version

__version__ = version("fem-toolbox")

__all__ = ("__version__",)
