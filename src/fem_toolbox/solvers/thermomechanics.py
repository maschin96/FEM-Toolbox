"""Sequential stationary heat to elasticity coupling on one Tet4 mesh."""

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.elasticity import VectorField
from fem_toolbox.physics.heat import Convection, ScalarField
from fem_toolbox.solvers.elasticity import ElasticResult, solve_elasticity
from fem_toolbox.solvers.heat import HeatResult, solve_heat


@dataclass(frozen=True)
class ThermomechanicalResult:
    heat: HeatResult
    elasticity: ElasticResult


def solve_thermomechanics(
    mesh: Mesh,
    thermal_material: ThermalMaterial | Sequence[ThermalMaterial],
    elastic_material: ElasticMaterial | Sequence[ElasticMaterial],
    expansion: ThermalExpansion | Sequence[ThermalExpansion],
    thermal_conditions: Iterable[Dirichlet],
    mechanical_conditions: Iterable[Dirichlet],
    *,
    source: ScalarField = 0.0,
    outward_fluxes: Mapping[int, ScalarField] | None = None,
    convection: Mapping[int, Convection] | None = None,
    body_force: VectorField = (0.0, 0.0, 0.0),
    tractions: Mapping[int, VectorField] | None = None,
    point_loads: Mapping[int, ArrayLike] | None = None,
) -> ThermomechanicalResult:
    """Solve heat then mechanics, transferring absolute nodal K on the same mesh."""
    heat = solve_heat(
        mesh,
        thermal_material,
        thermal_conditions,
        source=source,
        outward_fluxes=outward_fluxes,
        convection=convection,
    )
    temperature = np.zeros(mesh.n_points)
    temperature[heat.node_ids] = heat.temperature
    elasticity = solve_elasticity(
        mesh,
        elastic_material,
        mechanical_conditions,
        expansion=expansion,
        temperature=temperature,
        body_force=body_force,
        tractions=tractions,
        point_loads=point_loads,
    )
    return ThermomechanicalResult(heat, elasticity)
