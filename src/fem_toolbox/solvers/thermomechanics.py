"""Sequential stationary heat to elasticity coupling on one Tet4 or straight Tet10 mesh."""

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
from fem_toolbox.solvers.linear_backend import LinearBackend


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
    thermal_backend: LinearBackend | None = None,
    mechanical_backend: LinearBackend | None = None,
    source: ScalarField = 0.0,
    outward_fluxes: Mapping[int, ScalarField] | None = None,
    convection: Mapping[int, Convection] | None = None,
    body_force: VectorField = (0.0, 0.0, 0.0),
    tractions: Mapping[int, VectorField] | None = None,
    point_loads: Mapping[int, ArrayLike] | None = None,
) -> ThermomechanicalResult:
    """Solve heat then mechanics, transferring absolute nodal K on the same mesh.

    Optional thermal_backend/mechanical_backend select independent linear
    backends; None uses SciPy. Coupling is one-way, stationary, without thermal
    feedback from deformation. Tet10 interpolates all ten nodal temperatures
    in the eigenstrain integral; mechanical energy is integrated independently.
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
    elasticity = solve_elasticity(
        mesh,
        elastic_material,
        mechanical_conditions,
        backend=mechanical_backend,
        expansion=expansion,
        temperature=temperature,
        body_force=body_force,
        tractions=tractions,
        point_loads=point_loads,
    )
    return ThermomechanicalResult(heat, elasticity)
