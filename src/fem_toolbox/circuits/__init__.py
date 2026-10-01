"""Experimental linear circuit API; isolated from the stable FEM Model API."""

from fem_toolbox.circuits.ac import ACResult, solve_ac
from fem_toolbox.circuits.dc import CircuitTopologyError, DCResult, solve_dc
from fem_toolbox.circuits.model import Circuit, Component
from fem_toolbox.circuits.ports import ElectricalPort, ThermalPort
from fem_toolbox.circuits.transient import TransientResult, solve_transient

__all__ = [
    "ACResult",
    "Circuit",
    "CircuitTopologyError",
    "Component",
    "DCResult",
    "ElectricalPort",
    "ThermalPort",
    "TransientResult",
    "solve_ac",
    "solve_dc",
    "solve_transient",
]
