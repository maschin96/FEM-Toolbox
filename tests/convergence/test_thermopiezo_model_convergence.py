"""Four regular 3D meshes against an independent layered continuum solution."""

from __future__ import annotations

from itertools import product

import numpy as np
import pytest
from numpy.testing import assert_allclose

from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from tests.reference.thermopiezo_layer_1d import (
    BODY_FORCE,
    HEAT_SOURCE,
    LOWER,
    T_BOTTOM,
    T_REFERENCE,
    T_TOP,
    TOP_CHARGE,
    UPPER,
    VOLTAGE,
    VOLUME_CHARGE,
    heat_flux,
    layer,
    reference,
    temperature,
)

LEVELS = (2, 4, 6, 8)


def test_continuum_reference_satisfies_constitutive_and_boundary_data() -> None:
    analytic = reference(top_voltage=VOLTAGE)
    assert_allclose(
        [temperature(0), temperature(1), temperature(0.5)],
        [310, 300, 303 + 2 / 3],
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        [heat_flux(0), heat_flux(0.5), heat_flux(1)],
        [70 / 3, 82 / 3, 94 / 3],
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(
        [
            analytic.displacement(0),
            analytic.displacement(1),
            analytic.potential(0),
            analytic.potential(1),
        ],
        [0, 0, 0, VOLTAGE],
        atol=1e-12,
        rtol=0,
    )
    assert_allclose(analytic.top_charge, 1.8137459449929476, atol=1e-12, rtol=0)
    assert_allclose(analytic.heat_energy, 268.5, atol=1e-11, rtol=0)
    for z in (0.1, 0.4, 0.6, 0.9):
        material = layer(z)
        u_prime, phi_prime = analytic.derivatives(z)
        mechanical = u_prime - material.expansion * (
            float(temperature(z)) - T_REFERENCE
        )
        stress = material.stiffness * mechanical + material.coupling * phi_prime
        displacement = (
            material.coupling * mechanical - material.permittivity * phi_prime
        )
        assert_allclose(stress, analytic.stress(z), atol=1e-13, rtol=0)
        assert_allclose(
            displacement, analytic.electric_displacement(z), atol=1e-13, rtol=0
        )


def _model(
    cube_mesh,
    n: int,
    *,
    voltage: float | None = None,
    charge: float | None = None,
    top_displacement: float = 0.0,
) -> Model:
    original = cube_mesh(n)
    points = original.points
    cell_mid = points[original.tetra4, 2].mean(axis=1)
    face_z = points[original.triangle3, 2]
    lower_faces = np.flatnonzero(np.all(face_z == 0, axis=1))
    upper_faces = np.flatnonzero(np.all(face_z == 1, axis=1))
    mesh = Mesh(
        points=points,
        tetra4=original.tetra4,
        triangle3=original.triangle3,
        volume_sets={
            "lower": np.flatnonzero(cell_mid < 0.5),
            "upper": np.flatnonzero(cell_mid > 0.5),
        },
        facet_sets={"bottom": lower_faces, "top": upper_faces},
    )

    def piezo(material):
        coupling = np.zeros((3, 6))
        coupling[2, 2] = material.coupling
        return {
            "stiffness": np.diag([20, 20, material.stiffness, 4, 4, 4]).tolist(),
            "stress_coupling": coupling.tolist(),
            "permittivity": np.diag([5, 5, material.permittivity]).tolist(),
            "units": {
                "stiffness": "Pa",
                "stress_coupling": "C/m^2",
                "permittivity": "F/m",
            },
        }

    def expansion(material):
        return {
            "alpha": [0, 0, material.expansion, 0, 0, 0],
            "alpha_unit": "1/K",
            "reference_temperature": T_REFERENCE,
            "reference_temperature_unit": "K",
        }

    boundary = {
        "dirichlet": [
            {"field": "T", "set": "bottom", "value": T_BOTTOM, "unit": "K"},
            {"field": "T", "set": "top", "value": T_TOP, "unit": "K"},
            {
                "field": "u",
                "nodes": list(range(mesh.n_points)),
                "components": ["x", "y"],
                "value": 0,
                "unit": "m",
            },
            {
                "field": "u",
                "set": "bottom",
                "components": ["z"],
                "value": 0,
                "unit": "m",
            },
            {
                "field": "u",
                "set": "top",
                "components": ["z"],
                "value": top_displacement,
                "unit": "m",
            },
        ],
        "electrodes": [
            {"name": "bottom", "set": "bottom", "potential": 0, "potential_unit": "V"},
            (
                {
                    "name": "top",
                    "set": "top",
                    "potential": voltage,
                    "potential_unit": "V",
                }
                if voltage is not None
                else {"name": "top", "set": "top", "charge": charge, "charge_unit": "C"}
            ),
        ],
    }
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "thermopiezoelectric",
        "materials": {
            "T": {
                "lower": {"conductivity": LOWER.conductivity, "unit": "W/(m K)"},
                "upper": {"conductivity": UPPER.conductivity, "unit": "W/(m K)"},
            },
            "piezo": {"lower": piezo(LOWER), "upper": piezo(UPPER)},
            "expansion": {"lower": expansion(LOWER), "upper": expansion(UPPER)},
        },
        "regions": {
            "T": {"lower": "lower", "upper": "upper"},
            "piezo": {"lower": "lower", "upper": "upper"},
            "expansion": {"lower": "lower", "upper": "upper"},
        },
        "boundary": boundary,
        "loads": {
            "source": {"value": HEAT_SOURCE, "unit": "W/m^3"},
            "body_force": {"value": [0, 0, BODY_FORCE], "unit": "N/m^3"},
            "volume_charge": {"value": VOLUME_CHARGE, "unit": "C/m^3"},
        },
        "solver": {"electric_thermal_strain": True},
    }
    return Model(mesh=mesh, document=document)


def _rates(errors: np.ndarray) -> np.ndarray:
    return (
        np.log(errors[:-1] / errors[1:])
        / np.log(np.asarray(LEVELS[1:]) / LEVELS[:-1])[:, None]
    )


def _tetra_quadrature() -> tuple[np.ndarray, np.ndarray]:
    roots, weights = np.polynomial.legendre.leggauss(5)
    roots, weights = (roots + 1) / 2, weights / 2
    barycentric, quadrature = [], []
    for i, j, k in product(range(5), repeat=3):
        r, s, t = roots[[i, j, k]]
        a, b, c = r, (1 - r) * s, (1 - r) * (1 - s) * t
        barycentric.append([1 - a - b - c, a, b, c])
        quadrature.append(weights[i] * weights[j] * weights[k] * (1 - r) ** 2 * (1 - s))
    return np.asarray(barycentric), np.asarray(quadrature)


@pytest.mark.parametrize("electrode", ["voltage", "charge"])
def test_coupled_model_converges_with_material_interface(
    cube_mesh, electrode: str
) -> None:
    analytic = (
        reference(top_voltage=VOLTAGE)
        if electrode == "voltage"
        else reference(top_charge=TOP_CHARGE)
    )
    errors: list[list[float]] = []
    energy_errors: list[list[float]] = []
    barycentric, quadrature = _tetra_quadrature()
    for n in LEVELS:
        model = _model(
            cube_mesh,
            n,
            voltage=VOLTAGE if electrode == "voltage" else None,
            charge=TOP_CHARGE if electrode == "charge" else None,
        )
        run = model.solve()
        heat = run.field_result.heat
        piezo = run.field_result.piezoelectric
        cells = model.mesh.tetra4
        z = barycentric @ model.mesh.points[cells, 2].T
        z = z.T
        exact_u_prime, exact_phi_prime = analytic.derivatives_array(z)
        exact_temperature_prime = -(heat_flux(0) + HEAT_SOURCE * z) / np.where(
            z < 0.5, LOWER.conductivity, UPPER.conductivity
        )
        exact = np.stack(
            [
                temperature(z),
                analytic.displacement_array(z),
                analytic.potential_array(z),
            ],
            axis=2,
        )
        numerical = np.stack(
            [
                heat.temperature[cells] @ barycentric.T,
                piezo.displacement[cells, 2] @ barycentric.T,
                piezo.potential[cells] @ barycentric.T,
            ],
            axis=2,
        )
        integration = quadrature[None, :] / n**3
        squared_l2 = np.sum(
            integration[:, :, None] * (numerical - exact) ** 2, axis=(0, 1)
        )
        squared_h1 = np.array(
            [
                np.sum(
                    integration
                    * (
                        np.sum(heat.gradient[:, None, :2] ** 2, axis=2)
                        + (heat.gradient[:, None, 2] - exact_temperature_prime) ** 2
                    )
                ),
                np.sum(
                    integration
                    * (
                        piezo.strain[:, None, 3] ** 2
                        + piezo.strain[:, None, 4] ** 2
                        + (piezo.strain[:, None, 2] - exact_u_prime) ** 2
                    )
                ),
                np.sum(
                    integration
                    * (
                        np.sum(piezo.electric_field[:, None, :2] ** 2, axis=2)
                        + (-piezo.electric_field[:, None, 2] - exact_phi_prime) ** 2
                    )
                ),
            ]
        )
        errors.append(np.sqrt(squared_l2).tolist() + np.sqrt(squared_h1).tolist())

        volume = 1 / (6 * n**3)
        cell_k = np.where(
            model.mesh.points[model.mesh.tetra4, 2].mean(axis=1) < 0.5,
            LOWER.conductivity,
            UPPER.conductivity,
        )
        heat_energy = float(np.sum(cell_k * np.sum(heat.gradient**2, axis=1) * volume))
        energy_errors.append(
            [
                abs(heat_energy - analytic.heat_energy),
                abs(piezo.total_enthalpy - analytic.enthalpy),
            ]
        )

        states = {state.name: state for state in piezo.electrode_states}
        # Integrated electrode outputs inherit the piecewise smooth O(h^2)
        # error. The 10^-3 h^2 envelope is below 0.03% at the coarsest mesh.
        assert_allclose(
            states["top"].charge, analytic.top_charge, atol=1e-3 / n**2, rtol=0
        )
        assert_allclose(
            states["top"].potential, analytic.voltage, atol=1e-3 / n**2, rtol=0
        )
        assert abs(heat.power_balance) < 1e-8
        assert_allclose(heat.source_power, HEAT_SOURCE, atol=1e-11, rtol=0)
        assert_allclose(heat.dirichlet_inflow, -HEAT_SOURCE, atol=1e-8, rtol=0)
        assert_allclose(piezo.force_balance, 0, atol=1e-8, rtol=0)
        assert_allclose(piezo.moment_balance, 0, atol=1e-8, rtol=0)
        assert_allclose(piezo.reactions_u[:, 2].sum(), -BODY_FORCE, atol=1e-8, rtol=0)
        top_nodes = model.mesh.points[piezo.node_ids, 2] == 1
        assert_allclose(
            piezo.reactions_u[top_nodes, 2].sum(),
            analytic.stress(1.0),
            atol=0.02 / n**2,
            rtol=0,
        )
        assert abs(piezo.charge_balance) < 1e-8
        assert_allclose(
            sum(state.charge for state in states.values()),
            -VOLUME_CHARGE,
            atol=1e-8,
            rtol=0,
        )
        assert heat.relative_residual < 1e-10
        assert piezo.relative_residual_u < 1e-10
        assert piezo.relative_residual_phi < 1e-10
        assert (
            run.provenance.solver_parameters["effective_materials"]["piezo"]["catalog"][
                "lower"
            ]["units"]["permittivity"]
            == "F/m"
        )

        middle = model.mesh.points[model.mesh.tetra4, 2].mean(axis=1)
        lower = np.isclose(middle, 0.5 - 0.5 / n)
        upper = np.isclose(middle, 0.5 + 0.5 / n)
        for field, exact, slope, tolerance in (
            (heat.heat_flux[:, 2], heat_flux(0.5), HEAT_SOURCE, 1e-9),
            (piezo.stress[:, 2], analytic.stress(0.5), -BODY_FORCE, 0.02 / n**2),
            (
                piezo.electric_displacement[:, 2],
                analytic.electric_displacement(0.5),
                VOLUME_CHARGE,
                0.004 / n**2,
            ),
        ):
            left_limit = field[lower].mean() + slope * 0.5 / n
            right_limit = field[upper].mean() - slope * 0.5 / n
            assert_allclose([left_limit, right_limit], exact, atol=tolerance, rtol=0)

    errors_array = np.asarray(errors)
    energy_array = np.asarray(energy_errors)
    rates = _rates(errors_array)
    energy_rates = _rates(energy_array)
    assert np.all((rates[:, :3] > 1.75) & (rates[:, :3] < 2.25))
    assert np.all((rates[:, 3:] > 0.75) & (rates[:, 3:] < 1.25))
    assert np.all((energy_rates > 1.5) & (energy_rates < 2.5))


def test_reversing_prescribed_electrode_charge_changes_coupled_solution(
    cube_mesh,
) -> None:
    correct = _model(cube_mesh, 2, charge=TOP_CHARGE).solve().field_result.piezoelectric
    reversed_charge = (
        _model(cube_mesh, 2, charge=-TOP_CHARGE).solve().field_result.piezoelectric
    )
    expected = reference(top_charge=TOP_CHARGE)
    assert_allclose(correct.electrode_states[1].charge, expected.top_charge, atol=1e-10)
    assert_allclose(reversed_charge.electrode_states[1].charge, -TOP_CHARGE, atol=1e-10)
    assert abs(reversed_charge.electrode_states[1].potential - expected.voltage) > 0.1


def test_conductor_and_mechanical_port_reciprocity(cube_mesh) -> None:
    n = 4
    delta = 0.1

    def ports(*, displacement: float, voltage: float) -> tuple[float, float]:
        model = _model(cube_mesh, n, voltage=voltage, top_displacement=displacement)
        result = model.solve().field_result.piezoelectric
        top = model.mesh.points[result.node_ids, 2] == 1
        top_force = float(result.reactions_u[top, 2].sum())
        top_charge = next(
            state.charge for state in result.electrode_states if state.name == "top"
        )
        return top_force, top_charge

    force_high, _ = ports(displacement=0, voltage=VOLTAGE + delta)
    force_low, _ = ports(displacement=0, voltage=VOLTAGE - delta)
    _, charge_high = ports(displacement=delta, voltage=VOLTAGE)
    _, charge_low = ports(displacement=-delta, voltage=VOLTAGE)
    reciprocal_force = (force_high - force_low) / (2 * delta)
    reciprocal_charge = -(charge_high - charge_low) / (2 * delta)
    assert_allclose(reciprocal_force, reciprocal_charge, atol=2e-10, rtol=0)


def test_short_circuit_and_open_circuit_have_distinct_limits(cube_mesh) -> None:
    n = 4
    short = _model(cube_mesh, n, voltage=0).solve().field_result.piezoelectric
    opened = _model(cube_mesh, n, charge=0).solve().field_result.piezoelectric
    short_reference = reference(top_voltage=0)
    open_reference = reference(top_charge=0)
    short_top = next(state for state in short.electrode_states if state.name == "top")
    open_top = next(state for state in opened.electrode_states if state.name == "top")
    assert_allclose(short_top.potential, 0, atol=1e-12, rtol=0)
    assert_allclose(open_top.charge, 0, atol=1e-10, rtol=0)
    assert_allclose(
        short_top.charge, short_reference.top_charge, atol=1e-3 / n**2, rtol=0
    )
    assert_allclose(
        open_top.potential, open_reference.voltage, atol=1e-3 / n**2, rtol=0
    )
    assert abs(short_top.charge) > 0.05
    assert abs(open_top.potential) > 0.02
