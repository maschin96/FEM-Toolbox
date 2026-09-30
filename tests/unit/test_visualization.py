from pathlib import Path

import numpy as np
import pytest

from fem_toolbox.mesh import Mesh
from fem_toolbox.post.visualization import plot_field, save_jpeg

pytest.importorskip("matplotlib", reason="requires visualization extra")


def mesh():
    return Mesh(
        points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1], [9, 9, 9]],
        tetra4=[[0, 1, 2, 3], [0, 2, 1, 4]],
    )


def test_surface_and_cell_colors(tmp_path: Path):
    plot = plot_field(
        mesh(), [2, 8], association="cell", name="flux", unit="W/m²", limits=(0, 10)
    )
    assert len(plot.faces) == 6
    assert not any(set(face) == {0, 1, 2} for face in plot.faces)
    assert plot.face_values.tolist() == [2, 8, 2, 8, 2, 8]
    assert plot.artist.norm.vmin == 0
    assert plot.artist.norm.vmax == 10
    assert plot.axes.get_xlabel() == "x [m]"
    path = tmp_path / "field.jpeg"
    save_jpeg(plot, path, dpi=100)
    from PIL import Image

    with Image.open(path) as image:
        assert image.format == "JPEG"
        assert image.size == (800, 600)
    with pytest.raises(FileExistsError):
        save_jpeg(plot, path)
    plot.close()


def test_nodal_reduction_and_vectors():
    m = mesh()
    nodes = np.arange(5)
    p = plot_field(
        m,
        np.tile([3, 4, 0], (5, 1)),
        association="node",
        nodes=nodes,
        component="magnitude",
        name="u",
        unit="m",
    )
    np.testing.assert_allclose(p.face_values, 5, atol=1e-14, rtol=0)
    p.close()
    p = plot_field(
        m, m.points[nodes, 0], association="node", nodes=nodes, name="T", unit="K"
    )
    np.testing.assert_allclose(
        p.face_values, m.points[p.faces, 0].mean(axis=1), atol=1e-14, rtol=0
    )
    p.close()


@pytest.mark.parametrize(
    "component,index",
    [("xx", 0), ("yy", 1), ("zz", 2), ("yz", 3), ("xz", 4), ("xy", 5)],
)
def test_stress(component, index):
    p = plot_field(
        mesh(),
        np.tile(np.arange(6), (2, 1)),
        association="cell",
        name="stress",
        unit="Pa",
        component=component,
    )
    assert np.all(p.face_values == index)
    p.close()


@pytest.mark.parametrize(
    "values,kwargs",
    [
        ([1], {}),
        ([np.nan, 2], {}),
        ([1 + 1j, 2], {}),
        ([1, 2], {"limits": (2, 1)}),
        ([[1, 2, 3], [1, 2, 3]], {"component": "xx"}),
        ([1, 2], {"association": "bad"}),
    ],
)
def test_invalid(values, kwargs):
    options = dict(association="cell", name="T", unit="K")
    options.update(kwargs)
    with pytest.raises(ValueError):
        plot_field(mesh(), values, **options)


def test_missing_optional_dependency(monkeypatch):
    import builtins

    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.startswith(("matplotlib", "mpl_toolkits")):
            raise ImportError("blocked")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(ImportError, match=r"fem-toolbox\[visualization\]"):
        plot_field(mesh(), [1, 2], association="cell", name="T", unit="K")


def test_active_nodes_and_export_errors(tmp_path):
    with pytest.raises(ValueError, match="row count"):
        plot_field(mesh(), np.zeros(6), association="node", name="T", unit="K")
    with pytest.raises(ValueError, match="exactly once"):
        plot_field(
            mesh(),
            np.zeros(5),
            nodes=[0, 1, 2, 3, 3],
            association="node",
            name="T",
            unit="K",
        )
    p = plot_field(mesh(), [1, 1], association="cell", name="phi", unit="V")
    assert p.artist.norm.vmin < 1 < p.artist.norm.vmax
    with pytest.raises(ValueError, match="extension"):
        save_jpeg(p, tmp_path / "a.png")
    with pytest.raises(FileNotFoundError):
        save_jpeg(p, tmp_path / "missing" / "a.jpg")
    save_jpeg(p, tmp_path / "a.jpg")
    save_jpeg(p, tmp_path / "a.jpg", overwrite=True)
    p.close()


@pytest.mark.parametrize("component,index", [("x", 0), ("y", 1), ("z", 2)])
def test_vector_components(component, index):
    p = plot_field(
        mesh(),
        [[3, 4, 0], [3, 4, 0]],
        association="cell",
        component=component,
        name="E",
        unit="V/m",
        limits=(1, 2),
    )
    assert np.all(p.face_values == [3, 4, 0][index])
    assert p.figure.axes[1].get_ylabel() == f"E ({component}) [V/m]"
    p.close()


def test_no_mutation_or_registered_figures():
    from matplotlib import pyplot

    m = mesh()
    values = np.arange(5.0)
    before = values.copy()
    figures = pyplot.get_fignums()
    for _ in range(3):
        p = plot_field(m, values, association="node", name="T", unit="K")
        p.figure.canvas.draw()
        p.close()
        assert not p.figure.axes
    assert pyplot.get_fignums() == figures
    np.testing.assert_array_equal(values, before)
    assert not m.points.flags.writeable


def test_tet10_rejected():
    points = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)]
    points = np.vstack([points, [(points[a] + points[b]) / 2 for a, b in edges]])
    m = Mesh(points=points, tetra10=[list(range(10))])
    with pytest.raises(ValueError, match="Tet10"):
        plot_field(m, [1], association="cell", name="T", unit="K")
