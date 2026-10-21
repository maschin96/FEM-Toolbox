"""One-way interface load adapter to the existing verified Tet4 mechanics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from fem_toolbox.coupling.interface import InterfaceLoads, InterfaceMap, _fail
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.model import Model
from fem_toolbox.solvers.elasticity import ElasticResult, solve_elasticity


@dataclass(frozen=True, slots=True)
class OneWayResult:
    """Conservative loads and reference-geometry static structural response."""

    loads: InterfaceLoads
    structure: ElasticResult
    maximum_displacement_ratio: float
    maximum_displacement_gradient: float


def solve_one_way(
    mapping: InterfaceMap,
    model: Model,
    traction_Pa: ArrayLike,
    *,
    external_pressure_Pa: float = 0.0,
) -> OneWayResult:
    """Apply sigma_f*n_f [nsegment,2] to a validated elasticity Model.

    Existing body forces, tractions and prescribed displacements are preserved.
    Forces are added through existing point-load DOFs; there is no new element
    formulation. Structure mesh coordinates must equal the map coordinates.
    Small-displacement and gradient bounds are both 0.02, as in Doc/71.
    """
    inputs = model.elasticity_inputs()
    mesh = inputs["mesh"]
    if mesh.points.shape != mapping.structure_points_m.shape or not np.array_equal(
        mesh.points, mapping.structure_points_m
    ):
        _fail("structural Model and interface map coordinates do not match")
    loads = mapping.loads_from_fluid_traction(
        traction_Pa, external_pressure_Pa=external_pressure_Pa
    )
    point_loads = {
        int(n): f for n, f in enumerate(loads.nodal_forces_N) if np.any(f != 0)
    }
    result = solve_elasticity(**inputs, point_loads=point_loads)
    displacement = np.zeros_like(mesh.points)
    displacement[result.node_ids] = result.displacement
    ratio = float(np.max(np.linalg.norm(displacement, axis=1)) / mapping.length_scale_m)
    gradient = max(
        float(
            np.linalg.norm(
                displacement[cell].T @ Tet4Geometry(mesh.points[cell]).shape_gradients
            )
        )
        for cell in mesh.tetra4
    )
    if ratio > 0.02 or gradient > 0.02:
        _fail(
            f"small-deformation structure gate failed: displacement/Lref={ratio}, gradient={gradient}"
        )
    if result.relative_residual > 1e-10:
        _fail(f"structure equilibrium gate failed: {result.relative_residual}")
    return OneWayResult(loads, result, ratio, gradient)
