"""Hand-calculated DC references, physical balances and topology diagnostics."""

import pytest


def solve(text, **kwargs):
    from fem_toolbox.circuits import Circuit, solve_dc

    return solve_dc(Circuit.from_netlist(text), **kwargs)


def balances(result):
    # Small O(1-10 V), mA and mW references: absolute limits are far above
    # float64 roundoff but far below the physical signals. No exact equality.
    assert max(abs(v) for v in result.kcl_residual_A.values()) < 1e-12
    assert abs(result.power_balance_W) < 1e-12
    assert result.relative_residual < 1e-12


def test_voltage_divider_currents_power_and_immutable_output():
    r = solve("V1 in 0 12\nR1 in out 1k\nR2 out 0 2k")
    assert r.node_potentials_V == pytest.approx({"0": 0, "in": 12, "out": 8}, abs=1e-12)
    assert r.branch_currents_A == pytest.approx(
        {"V1": -0.004, "R1": 0.004, "R2": 0.004}, abs=1e-14
    )
    assert r.branch_powers_W == pytest.approx(
        {"V1": -0.048, "R1": 0.016, "R2": 0.032}, abs=1e-13
    )
    balances(r)
    with pytest.raises(TypeError):
        r.node_potentials_V["out"] = 999
    with pytest.raises(TypeError):
        r.branch_currents_A["R1"] = 999


def test_current_source_orientation_and_parallel_resistors():
    r = solve("I1 0 a 3m\nR1 a 0 4k\nR2 a 0 4k")
    assert r.node_potentials_V["a"] == pytest.approx(6, abs=1e-12)
    assert r.branch_currents_A["I1"] == pytest.approx(0.003, abs=1e-14)
    assert r.branch_currents_A["R1"] == pytest.approx(0.0015, abs=1e-14)
    assert r.branch_powers_W["I1"] == pytest.approx(-0.018, abs=1e-13)
    balances(r)


def test_reversed_voltage_source_and_negative_current_source():
    r = solve("V1 0 a 5\nR1 a 0 1k\nI1 a 0 -1m")
    assert r.node_potentials_V["a"] == pytest.approx(-5, abs=1e-12)
    assert r.branch_currents_A["R1"] == pytest.approx(-0.005, abs=1e-14)
    assert r.branch_currents_A["V1"] == pytest.approx(-0.006, abs=1e-14)
    balances(r)


def test_floating_voltage_source_between_grounded_nodes():
    # v_a-v_b=3, v_a/1000+v_b/2000=0 -> v_a=1, v_b=-2.
    r = solve("V1 a b 3\nR1 a 0 1k\nR2 b 0 2k")
    assert r.node_potentials_V["a"] == pytest.approx(1, abs=1e-12)
    assert r.node_potentials_V["b"] == pytest.approx(-2, abs=1e-12)
    assert r.branch_currents_A["V1"] == pytest.approx(-0.001, abs=1e-14)
    balances(r)


def test_dc_capacitor_open_inductor_short():
    r = solve("V1 in 0 10\nR1 in a 1k\nL1 a b 2m\nR2 b 0 1k\nC1 a 0 1u")
    assert r.node_potentials_V["a"] == pytest.approx(5, abs=1e-12)
    assert r.node_potentials_V["b"] == pytest.approx(5, abs=1e-12)
    assert r.branch_currents_A["L1"] == pytest.approx(0.005, abs=1e-14)
    assert r.branch_currents_A["C1"] == pytest.approx(0, abs=1e-14)
    assert r.branch_voltages_V["L1"] == pytest.approx(0, abs=1e-12)
    balances(r)


def test_balanced_bridge_against_independent_branch_reduction():
    r = solve("V1 in 0 6\nR1 in a 1k\nR2 a 0 2k\nR3 in b 2k\nR4 b 0 4k\nR5 a b 5k")
    assert r.node_potentials_V["a"] == pytest.approx(4, abs=1e-12)
    assert r.node_potentials_V["b"] == pytest.approx(4, abs=1e-12)
    assert r.branch_currents_A["R5"] == pytest.approx(0, abs=1e-14)
    assert r.branch_currents_A["V1"] == pytest.approx(-0.003, abs=1e-14)
    balances(r)


def test_no_source_zero_solution_and_direct_voltage_constraint():
    for text in ("R1 a 0 10\nC1 a 0 1u", "V1 a 0 0"):
        r = solve(text)
        assert r.node_potentials_V["a"] == pytest.approx(0, abs=1e-12)
        balances(r)


@pytest.mark.parametrize(
    ("text", "kind"),
    [
        ("R1 a 0 1k\nR2 b c 1k", "floating"),
        ("I1 a 0 1m", "floating"),
        ("C1 a 0 1u", "floating"),
        ("R1 a 0 1k\nC1 a b 1u", "floating"),
        ("V1 a 0 1\nV2 a 0 2", "inconsistent_voltage_loop"),
        ("V1 a 0 1\nV2 a 0 1", "redundant_voltage_loop"),
        ("V1 a 0 1\nV2 b a 2\nV3 b 0 4", "inconsistent_voltage_loop"),
        ("V1 a 0 1\nV2 b a 2\nV3 b 0 3", "redundant_voltage_loop"),
        ("V1 a 0 1\nL1 a 0 1m", "inconsistent_voltage_loop"),
        ("L1 a 0 1m\nL2 a 0 2m", "redundant_voltage_loop"),
    ],
)
def test_invalid_dc_topology_before_factorization(text, kind, monkeypatch):
    from fem_toolbox.circuits import CircuitTopologyError, dc

    def forbidden(*args, **kwargs):
        raise AssertionError("topology must be rejected before numerical solve")

    monkeypatch.setattr(dc, "splu", forbidden)
    with pytest.raises(CircuitTopologyError) as caught:
        solve(text)
    assert caught.value.kind == kind


def test_deterministic_provenance_across_order_comments_and_suffixes():
    a = solve("V1 a 0 5\nR1 a 0 1k")
    b = solve("* same circuit\nr1 A 0 1000\nv1 A 0 5\n.end")
    assert a.provenance.to_dict() == b.provenance.to_dict()
    changed = solve("V1 a 0 5\nR1 a 0 2k")
    assert a.provenance.input_sha256 != changed.provenance.input_sha256


@pytest.mark.parametrize("value", [True, -1, float("nan"), float("inf")])
def test_bad_balance_tolerance(value):
    from fem_toolbox.errors import ValidationError

    with pytest.raises(ValidationError):
        solve("R1 a 0 1", current_tolerance_A=value)


def test_unrepresentable_conductance_fails_explicitly():
    from fem_toolbox.errors import ValidationError

    with pytest.raises(ValidationError):
        solve("V1 a 0 1\nR1 a 0 1e-320")


def test_large_inconsistent_voltage_loop_has_concrete_diagnosis():
    from fem_toolbox.circuits import CircuitTopologyError

    with pytest.raises(CircuitTopologyError) as caught:
        solve("V1 a 0 1e308\nV2 b 0 -1e308\nV3 a b 1")
    assert caught.value.kind == "inconsistent_voltage_loop"


def test_unbalanced_bridge_against_hand_solved_two_node_equations():
    # 3*a-b=10, -a+2.5*b=10 -> a=70/13, b=80/13 V.
    r = solve("V1 in 0 10\nR1 in a 1k\nR2 a 0 1k\nR3 in b 1k\nR4 b 0 2k\nR5 a b 1k")
    assert r.node_potentials_V["a"] == pytest.approx(70 / 13, abs=1e-12)
    assert r.node_potentials_V["b"] == pytest.approx(80 / 13, abs=1e-12)
    assert r.branch_currents_A["R5"] == pytest.approx(-1 / 1300, abs=1e-14)
    balances(r)


@pytest.mark.parametrize("resistance", [1e-12, 1e12])
def test_equilibration_across_conductance_scales(resistance):
    r = solve(f"V1 a 0 1\nR1 a 0 {resistance}")
    assert r.node_potentials_V["a"] == pytest.approx(1, abs=1e-12)
    assert r.branch_currents_A["R1"] == pytest.approx(1 / resistance, rel=1e-12)
    assert r.branch_currents_A["V1"] == pytest.approx(-1 / resistance, rel=1e-12)
    assert abs(r.power_balance_W) < 1e-12
