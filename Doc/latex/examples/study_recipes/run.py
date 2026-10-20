"""Reproduce a stationary recipe and its lifecycle in a fresh output directory."""

import argparse
import copy
import json
import platform
import time
import tomllib
from pathlib import Path
from shutil import copytree

import meshio
import numpy as np

import fem_toolbox
from fem_toolbox.io import read_gmsh
from fem_toolbox.model import Model
from fem_toolbox.studies import Study
from fem_toolbox.study_reports import Quantity, export_study_table, read_study_table


def quantities(recipe):
    if recipe == "heat":
        return [
            Quantity(
                "temperature",
                "point",
                "temperature_K",
                "K",
                "temperature at node 1",
                indices=(1,),
            ),
            Quantity("source", "metric", "heat.source_power_W", "W", "total source"),
            Quantity(
                "inflow",
                "metric",
                "heat.dirichlet_inflow_W",
                "W",
                "signed Dirichlet inflow",
            ),
            Quantity(
                "balance",
                "metric",
                "heat.power_balance_W",
                "W",
                "power imbalance",
                absolute_tolerance=1e-11,
            ),
        ]
    result = [
        Quantity(
            "displacement",
            "point",
            "displacement_m",
            "m",
            "mean top displacement",
            component="z",
            indices=(1, 3, 5, 7),
            reduction="mean",
        ),
        Quantity(
            "charge",
            "metric",
            "piezoelectric.electrodes.1.charge_C",
            "C",
            "driven electrode charge",
        ),
        Quantity(
            "force_balance",
            "metric",
            "piezoelectric.force_balance_N",
            "N",
            "force imbalance norm",
            component="norm",
            absolute_tolerance=1e-10,
        ),
        Quantity(
            "charge_balance",
            "metric",
            "piezoelectric.charge_balance_C",
            "C",
            "charge imbalance",
            absolute_tolerance=1e-18,
        ),
    ]
    for name, field in [
        ("mechanical", "mechanical_energy_J"),
        ("coupling", "coupling_energy_J"),
        ("dielectric", "dielectric_enthalpy_J"),
    ]:
        result.append(
            Quantity(
                name,
                "cell",
                field,
                "J",
                "integral over specimen",
                reduction="sum",
                indices=tuple(range(6)),
            )
        )
    return result


def run(recipe, output):
    output.mkdir(parents=True, exist_ok=False)
    inputs = output / "inputs"
    copytree(Path(__file__).with_name(recipe), inputs)
    manifest = tomllib.loads((inputs / "study.toml").read_text())
    runner = Study.from_file(inputs / "study.toml")
    started = time.perf_counter()
    partial = runner.run(output / "runs", max_cases=1)
    assert partial["cases"][0]["status"] == "succeeded"
    assert all(c["status"] == "pending" for c in partial["cases"][1:])
    export_study_table(
        output / "runs/catalog.json", output / "partial", quantities(recipe)
    )
    catalog = runner.run(output / "runs")
    assert all(c["status"] == "succeeded" for c in catalog["cases"]), catalog
    resumed = runner.run(output / "runs")
    assert catalog == resumed
    elapsed = time.perf_counter() - started
    export_study_table(
        output / "runs/catalog.json", output / "table", quantities(recipe)
    )
    rows = read_study_table(output / "table")["rows"]
    checks = []
    for case, row in zip(catalog["cases"], rows, strict=True):
        assert not row["missing"], row["missing"]
        assert row["balance_status"] == "passed", row
        # Construct separate direct inputs from the explicit recipe parameters.
        doc = tomllib.loads((inputs / "model.toml").read_text())
        p = case["parameters"]
        if recipe == "heat":
            k = p["conductivity"]["value"]
            doc["materials"]["T"]["thermal"]["conductivity"] = k
            expected = {"temperature": 300 + 2 / k, "source": 8 / 6, "inflow": -8 / 6}
            tolerances = {"temperature": 1e-11, "source": 1e-11, "inflow": 1e-11}
        else:
            geometry = p["geometry"]["value"]
            coupling = p["coupling"]["value"]
            doc["mesh"]["path"] = geometry
            doc["materials"]["piezo"]["solid"]["stress_coupling"] = coupling
            length = {"short.msh": 0.02, "long.msh": 0.04}[geometry]
            e = coupling[2][2]
            expected = {
                "displacement": -e / 1e9 * 100,
                "charge": 1e-4 / length * (1e-8 + e * e / 1e9) * 100,
                "mechanical": 0.5 * 1e-4 / length * e * e / 1e9 * 100**2,
                "coupling": -1e-4 / length * e * e / 1e9 * 100**2,
                "dielectric": -0.5 * 1e-4 / length * 1e-8 * 100**2,
            }
            tolerances = {
                "displacement": 1e-16,
                "charge": 1e-18,
                "mechanical": 1e-16,
                "coupling": 1e-16,
                "dielectric": 1e-16,
            }
        mesh = read_gmsh(inputs / doc["mesh"]["path"], length_scale_to_m=1.0)
        direct = Model(mesh=mesh, document=doc).solve().field_result
        fields = meshio.read(output / "runs" / case["result"]["fields"])
        if recipe == "heat":
            np.testing.assert_allclose(
                fields.point_data["temperature_K"],
                direct.temperature,
                rtol=2e-11,
                atol=1e-11,
            )
            direct_values = {
                "temperature": float(direct.temperature[1]),
                "source": direct.source_power,
                "inflow": direct.dirichlet_inflow,
                "balance": direct.power_balance,
            }
        else:
            np.testing.assert_allclose(
                fields.point_data["displacement_m"],
                direct.displacement,
                rtol=2e-11,
                atol=1e-16,
            )
            np.testing.assert_allclose(
                fields.point_data["potential_V"],
                direct.potential,
                rtol=2e-11,
                atol=1e-11,
            )
            direct_values = {
                "displacement": float(direct.displacement[[1, 3, 5, 7], 2].mean()),
                "charge": direct.electrode_states[1].charge,
                "mechanical": float(direct.mechanical_energy.sum()),
                "coupling": float(direct.coupling_energy.sum()),
                "dielectric": float(direct.dielectric_enthalpy.sum()),
                "force_balance": float(np.linalg.norm(direct.force_balance)),
                "charge_balance": direct.charge_balance,
            }
        for name, value in expected.items():
            np.testing.assert_allclose(
                row["values"][name], value, rtol=2e-11, atol=tolerances[name]
            )
        for name, value in direct_values.items():
            np.testing.assert_allclose(
                row["values"][name], value, rtol=2e-11, atol=tolerances.get(name, 1e-18)
            )
        checks.append(
            {
                "case_id": case["case_id"],
                "expected": expected,
                "actual": row["values"],
                "direct": direct_values,
                "absolute_tolerances": tolerances,
                "relative_tolerance": 2e-11,
            }
        )
    # One physically invalid material variant remains visible next to valid cases.
    failure = copy.deepcopy(manifest)
    if recipe == "heat":
        failure["parameters"][0]["values"] = [2.0, -1.0]
    else:
        failure["parameters"] = [
            {
                "name": "permittivity",
                "target": "materials.piezo.solid.permittivity",
                "unit": "F/m",
                "values": [np.eye(3).tolist(), (-np.eye(3)).tolist()],
            }
        ]
    failed = Study(failure, root=inputs).run(output / "failure")
    assert [c["status"] for c in failed["cases"]] == ["succeeded", "failed"]
    export_study_table(
        output / "failure/catalog.json", output / "failure_table", quantities(recipe)
    )
    # Baseline change retains IDs but changes fingerprints and recomputes each case.
    model_path = inputs / "model.toml"
    before = model_path.read_text()
    changed = (
        before.replace("value = 8.0", "value = 16.0")
        if recipe == "heat"
        else before.replace("potential = 100.0", "potential = 50.0")
    )
    assert changed != before
    model_path.write_text(changed)
    invalidated = runner.run(output / "runs")
    assert all(
        a["case_id"] == b["case_id"]
        and a["fingerprints"]["input"] != b["fingerprints"]["input"]
        and b["attempts"] == a["attempts"] + 1
        and b["status"] == "succeeded"
        for a, b in zip(catalog["cases"], invalidated["cases"], strict=True)
    )
    model_path.write_text(before)
    runner.run(output / "runs")
    report = {
        "recipe": recipe,
        "toolbox_version": fem_toolbox.__version__,
        "package_path": fem_toolbox.__file__,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "study_seconds": elapsed,
        "cases": checks,
        "resume_attempts_unchanged": True,
        "input_change_invalidated": True,
        "failed_cases": 1,
        "failure_error": failed["cases"][1].get("error"),
        "direct_and_analytic_checks_passed": True,
    }
    (output / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("recipe", choices=("heat", "piezo"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.recipe, args.output)
