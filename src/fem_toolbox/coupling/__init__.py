"""Experimental stationary planar fluid/structure interface operators."""

from fem_toolbox.coupling.fsi import (
    FsiCouplingError,
    FsiIteration,
    StationaryFsiProblem,
    StationaryFsiResult,
    solve_stationary_fsi,
)
from fem_toolbox.coupling.interface import (
    InterfaceLoads,
    InterfaceMap,
    TransferAudit,
    build_domain_interface,
    build_interface_map,
)
from fem_toolbox.coupling.structure import OneWayResult, solve_one_way

__all__ = [
    "FsiCouplingError",
    "FsiIteration",
    "InterfaceLoads",
    "InterfaceMap",
    "OneWayResult",
    "StationaryFsiProblem",
    "StationaryFsiResult",
    "TransferAudit",
    "build_domain_interface",
    "build_interface_map",
    "solve_one_way",
    "solve_stationary_fsi",
]
