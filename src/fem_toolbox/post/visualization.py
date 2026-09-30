"""Optional, headless Tet4 exterior field rendering (SI coordinates)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.mesh import Mesh, inspect_boundary_topology

if TYPE_CHECKING:
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection


@dataclass
class FieldPlot:
    """Owned figure; call close after saving. Face arrays describe raw colors."""

    figure: Figure
    axes: Axes3D
    artist: Poly3DCollection
    faces: NDArray[np.int64]
    face_values: NDArray[np.float64]

    def close(self) -> None:
        """Release figure artists; figures are never registered with pyplot."""
        self.figure.clear()


def plot_field(
    mesh: Mesh,
    values: ArrayLike,
    *,
    association: str,
    name: str,
    unit: str,
    nodes: ArrayLike | None = None,
    component: str | None = None,
    limits: tuple[float, float] | None = None,
    cmap: str = "viridis",
    view: tuple[float, float] = (25, -60),
    figsize: tuple[float, float] = (8, 6),
) -> FieldPlot:
    """Render real scalar, xyz vector or Voigt stress on Tet4 exterior.

    Node rows follow explicitly supplied active ``nodes`` (default: sorted
    volume nodes); cell rows follow tetra4. Every active node is required.
    Scalar nodal colors are arithmetic face means *after* component/magnitude
    reduction. Cell colors use the owning volume, without smoothing.
    Limits are finite and strictly increasing; clipped colors have colorbar
    extensions. Constant automatic ranges expand by max(1, abs(value))*1e-6.
    Angles are elevation/azimuth in degrees; size is in inches. No show(),
    pyplot registration, style or backend changes occur.
    """
    if mesh.tetra10.shape[0]:
        raise ValueError("Tet10 field visualization is unsupported; use ParaView")
    if association not in ("node", "cell"):
        raise ValueError("association must be node or cell")
    if not name.strip() or not unit.strip():
        raise ValueError("name and SI unit must be nonempty")
    raw = np.asarray(values)
    if np.iscomplexobj(raw):
        raise ValueError("complex fields are unsupported")
    data = np.asarray(raw, dtype=float)
    if data.ndim not in (1, 2) or not np.all(np.isfinite(data)):
        raise ValueError("field must be finite scalar/vector rows")
    active = np.unique(mesh.tetra4)
    if not len(active):
        raise ValueError("nonempty Tet4 volume required")
    if association == "node":
        selected = active if nodes is None else np.asarray(nodes)
        if (
            selected.ndim != 1
            or selected.dtype.kind not in "iu"
            or not np.array_equal(np.sort(selected), active)
        ):
            raise ValueError("nodes must contain each active volume node exactly once")
        expected = len(selected)
    else:
        if nodes is not None:
            raise ValueError("nodes only applies to node fields")
        selected = active
        expected = len(mesh.tetra4)
    if len(data) != expected:
        raise ValueError(
            "field row count does not match association; no export padding"
        )
    label = name
    if data.ndim == 1:
        if component is not None:
            raise ValueError("scalar field has no component")
        scalar = data
    else:
        components = {3: ("x", "y", "z"), 6: ("xx", "yy", "zz", "yz", "xz", "xy")}
        names = components.get(data.shape[1], ())
        if component == "magnitude" and data.shape[1] == 3:
            scalar = np.linalg.norm(data, axis=1)
        elif component in names:
            scalar = data[:, names.index(component)]
        else:
            raise ValueError(
                "choose xyz/magnitude for vectors or xx,yy,zz,yz,xz,xy for Voigt"
            )
        label += f" ({component})"
    exterior = inspect_boundary_topology(mesh).exterior_faces
    faces = np.array([face.nodes for face in exterior], dtype=np.int64)
    if association == "node":
        lookup = np.full(mesh.n_points, -1, dtype=int)
        lookup[selected] = np.arange(len(selected))
        colors = scalar[lookup[faces]].mean(axis=1)
        reduction = "nodal face mean"
    else:
        colors = scalar[[face.cell_index for face in exterior]]
        reduction = "owner cell value"
    if not np.all(np.isfinite(colors)):
        raise ValueError("field reduction overflowed")
    low, high = (float(colors.min()), float(colors.max())) if limits is None else limits
    if limits is None and low == high:
        width = max(1, abs(low)) * 1e-6
        low, high = low - width, high + width
    if not np.isfinite([low, high]).all() or low >= high:
        raise ValueError("limits must be finite and increasing")
    if (
        len(figsize) != 2
        or not np.isfinite(figsize).all()
        or min(figsize) <= 0
        or len(view) != 2
        or not np.isfinite(view).all()
    ):
        raise ValueError("finite view and positive figsize required")
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.colors import Normalize
        from matplotlib.figure import Figure
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as error:
        raise ImportError(
            "Install visualization: pip install 'fem-toolbox[visualization]'"
        ) from error
    figure = Figure(figsize=figsize)
    FigureCanvasAgg(figure)
    axes = figure.add_subplot(111, projection="3d")
    artist = Poly3DCollection(mesh.points[faces], cmap=cmap, norm=Normalize(low, high))
    artist.set_array(colors)
    axes.add_collection3d(artist)
    points = mesh.points[active]
    for axis, setter in enumerate((axes.set_xlim, axes.set_ylim, axes.set_zlim)):
        setter(float(points[:, axis].min()), float(points[:, axis].max()))
    axes.set_box_aspect(np.ptp(points, axis=0))
    axes.set(
        xlabel="x [m]", ylabel="y [m]", zlabel="z [m]", title=f"{label}: {reduction}"
    )
    axes.view_init(elev=view[0], azim=view[1])
    under, over = bool(np.any(colors < low)), bool(np.any(colors > high))
    extend = (
        "both" if under and over else "min" if under else "max" if over else "neither"
    )
    figure.colorbar(
        artist, ax=axes, shrink=0.65, pad=0.12, label=f"{label} [{unit}]", extend=extend
    )
    return FieldPlot(figure, axes, artist, faces, colors)


def save_jpeg(
    plot: FieldPlot, path: str | Path, *, dpi: float = 150, overwrite: bool = False
) -> None:
    """Save JPEG in existing directory; exclusive creation unless overwrite=True.

    Pixel size is figure inches times DPI. OS write errors propagate with path.
    The caller retains ownership and must close the plot.
    """
    target = Path(path)
    if target.suffix.lower() not in (".jpg", ".jpeg"):
        raise ValueError("JPEG requires .jpg or .jpeg extension")
    if not np.isfinite(dpi) or dpi <= 0:
        raise ValueError("dpi must be finite and positive")
    with target.open("wb" if overwrite else "xb") as stream:
        plot.figure.savefig(stream, format="jpeg", dpi=dpi)
