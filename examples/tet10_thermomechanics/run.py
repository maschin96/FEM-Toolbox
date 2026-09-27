"""Stationary Tet10 heat -> mechanics with independent quadratic reference fields."""

import argparse
from pathlib import Path

import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.io import export_result_bundle, read_gmsh
from fem_toolbox.materials import ElasticMaterial, ThermalExpansion, ThermalMaterial
from fem_toolbox.solvers.elasticity import mechanical_dofs
from fem_toolbox.solvers.heat import thermal_dofs
from fem_toolbox.solvers.thermomechanics import solve_thermomechanics


def run(output):
    mesh = read_gmsh(Path(__file__).resolve().parent / "body.msh", length_scale_to_m=1)
    left = np.unique(mesh.facets[mesh.facet_sets["left"]])
    result = solve_thermomechanics(
        mesh,
        ThermalMaterial.isotropic(2),
        ElasticMaterial.isotropic(100, 0),
        ThermalExpansion([0.01, 0, 0, 0, 0, 0], 1),
        [Dirichlet(thermal_dofs(mesh).indices("T", left), np.ones(len(left)))],
        [Dirichlet(mechanical_dofs(mesh).indices("u", left), np.zeros(3 * len(left)))],
        source=-4,
        outward_fluxes={int(f): -4 for f in mesh.facet_sets["right"]},
        body_force=lambda p: [-2 + 2 * p[0], 0, 0],
        tractions={int(f): [1, 0, 0] for f in mesh.facet_sets["right"]},
    )
    np.testing.assert_allclose(
        result.heat.temperature, 1 + mesh.points[:, 0] ** 2, rtol=0, atol=4e-12
    )
    exact = np.zeros_like(mesh.points)
    exact[:, 0] = 0.01 * mesh.points[:, 0] ** 2
    np.testing.assert_allclose(
        result.elasticity.displacement, exact, rtol=0, atol=1e-13
    )
    x = result.elasticity.field_points[:, 0]
    stress = np.zeros_like(result.elasticity.stress)
    stress[:, 0] = 2 * x - x**2
    np.testing.assert_allclose(result.elasticity.stress, stress, rtol=0, atol=2e-11)
    assert abs(result.elasticity.energy - 0.008 / 3) < 2e-14
    assert abs(result.heat.power_balance) < 2e-11
    np.testing.assert_allclose(result.elasticity.force_balance, 0, rtol=0, atol=2e-12)
    np.testing.assert_allclose(result.elasticity.moment_balance, 0, rtol=0, atol=2e-12)
    output.mkdir(parents=True, exist_ok=False)
    export_result_bundle(output / "coupled", mesh, result)
    print(f"Verified Tet10 thermomechanics: {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("tet10_thermal_expansion_results")
    )
    run(parser.parse_args().output)
