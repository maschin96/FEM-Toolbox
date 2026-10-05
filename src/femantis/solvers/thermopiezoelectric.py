"""One-way stationary heat transfer into the monolithic piezoelectric system."""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from femantis.boundary import Dirichlet, Electrode
from femantis.materials import (
    PiezoelectricMaterial,
    ThermalExpansion,
    ThermalMaterial,
)
from femantis.mesh import Mesh
from femantis.physics.elasticity import VectorField
from femantis.physics.electrostatics import ChargeDensity
from femantis.physics.heat import Convection, ScalarField
from femantis.solvers.heat import HeatResult, solve_heat
from femantis.solvers.linear_backend import LinearBackend
from femantis.solvers.piezoelectric import PiezoelectricResult, solve_piezoelectric


@dataclass(frozen=True, slots=True)
class ThermopiezoelectricResult:
    """Thermal and coupled results, each with its own residuals and SI balances."""

    heat: HeatResult
    piezoelectric: PiezoelectricResult


def solve_thermopiezoelectric(
    mesh: Mesh,
    thermal_material: ThermalMaterial | Sequence[ThermalMaterial],
    piezoelectric_material: PiezoelectricMaterial | Sequence[PiezoelectricMaterial],
    expansion: ThermalExpansion | Sequence[ThermalExpansion],
    thermal_conditions: Iterable[Dirichlet],
    piezoelectric_conditions: Iterable[Dirichlet],
    *,
    thermal_backend: LinearBackend | None = None,
    piezoelectric_backend: LinearBackend | None = None,
    source: ScalarField = 0.0,
    outward_fluxes: Mapping[int, ScalarField] | None = None,
    convection: Mapping[int, Convection] | None = None,
    body_force: VectorField = (0.0, 0.0, 0.0),
    tractions: Mapping[int, VectorField] | None = None,
    point_loads: Mapping[int, ArrayLike] | None = None,
    volume_charge: ChargeDensity = 0.0,
    surface_charges: Mapping[int, ChargeDensity] | None = None,
    electrodes: Sequence[Electrode] | None = None,
    electric_thermal_strain: bool = False,
    block_scaling: bool = True,
    method: str = "direct",
) -> ThermopiezoelectricResult:
    """Solve heat first, then transfer absolute nodal K on the same Tet4/Tet10 mesh.

    Expansion is defined per cell in global engineering-Voigt order (1/K)
    with a reference temperature (K). ``electric_thermal_strain`` controls
    the optional thermal contribution to electric displacement explicitly.
    Separate optional thermal_backend/piezoelectric_backend hooks select linear
    backends (None uses SciPy). There is no reciprocal heat feedback or time dependence.
    """
    heat = solve_heat(
        mesh,
        thermal_material,
        thermal_conditions,
        backend=thermal_backend,
        source=source,
        outward_fluxes=outward_fluxes,
        convection=convection,
    )
    temperature = np.zeros(mesh.n_points)
    temperature[heat.node_ids] = heat.temperature
    piezoelectric = solve_piezoelectric(
        mesh,
        piezoelectric_material,
        piezoelectric_conditions,
        backend=piezoelectric_backend,
        body_force=body_force,
        tractions=tractions,
        point_loads=point_loads,
        volume_charge=volume_charge,
        surface_charges=surface_charges,
        electrodes=electrodes,
        expansion=expansion,
        temperature=temperature,
        electric_thermal_strain=electric_thermal_strain,
        block_scaling=block_scaling,
        method=method,
    )
    return ThermopiezoelectricResult(heat, piezoelectric)
