"""Strict reproducible input: geometry, SI fluid and named boundary mode."""

import pytest

from fem_toolbox.errors import ValidationError


def document():
    return dict(
        flow_schema_version=1,
        grid=dict(cells=[4, 3], size_m=[2, 1]),
        fluid=dict(density_kg_m3=1, viscosity_Pa_s=1),
        boundary=dict(
            type="channel", inlet_velocity_m_s=[0.1, 0.2, 0.1], outlet_pressure_Pa=5
        ),
    )


def load(doc):
    from fem_toolbox.cfd import FlowProblem

    return FlowProblem.from_document(doc)


def test_roundtrip_and_provenance():
    p = load(document())
    r = p.solve()
    q = load(p.to_document())
    assert q.solve().provenance == r.provenance
    d = p.to_document()
    d["grid"]["cells"][0] = 99
    assert p.grid.cells == (4, 3)


def test_toml_from_file(tmp_path):
    path = tmp_path / "flow.toml"
    path.write_text("""flow_schema_version=1
[grid]
cells=[4,3]
size_m=[2,1]
[fluid]
density_kg_m3=1
viscosity_Pa_s=1
[boundary]
type="channel"
inlet_velocity_m_s=[0.1,0.2,0.1]
outlet_pressure_Pa=5
""")
    from fem_toolbox.cfd import FlowProblem

    assert FlowProblem.from_file(path).to_document() == document()


@pytest.mark.parametrize(
    "kind", ["unknown", "version", "boundary", "missing", "syntax", "profile"]
)
def test_invalid_document(kind, tmp_path):
    d = document()
    if kind == "unknown":
        d["grid"]["bogus"] = 1
    elif kind == "version":
        d["flow_schema_version"] = True
    elif kind == "boundary":
        d["boundary"]["type"] = "outlet-only"
    elif kind == "missing":
        del d["fluid"]["density_kg_m3"]
    elif kind == "profile":
        d["boundary"]["inlet_velocity_m_s"] = [1, 2]
    with pytest.raises(ValidationError):
        if kind == "syntax":
            from fem_toolbox.cfd import FlowProblem

            path = tmp_path / "bad.toml"
            path.write_text("[")
            FlowProblem.from_file(path)
        else:
            load(d)
