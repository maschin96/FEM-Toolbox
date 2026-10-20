"""Deterministic conforming edge-star bisection for affine Tet4 meshes."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.assembly import indices
from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.errors import ValidationError
from fem_toolbox.mesh.model import Mesh
from fem_toolbox.mesh.topology import inspect_boundary_topology


@dataclass(frozen=True)
class Refinement:
    """Child-to-input row maps; external IDs are deliberately not reused."""

    mesh: Mesh
    cell_parents: NDArray[np.int64]
    facet_parents: NDArray[np.int64]

    def __post_init__(self) -> None:
        for name in ("cell_parents", "facet_parents"):
            value = np.asarray(getattr(self, name))
            if value.ndim != 1 or value.dtype.kind not in "iu":
                raise ValidationError("invalid refinement parent array", field=name)
            value = np.asarray(value, dtype=np.int64)
            object.__setattr__(
                self, name, np.frombuffer(value.tobytes(), dtype=np.int64)
            )


def validate_tet4(mesh: Mesh) -> None:
    """Require nonempty, complete affine Tet4 topology and nondegenerate cells."""
    if not mesh.n_tetra4 or mesh.n_tetra10:
        raise ValidationError("requires a nonempty Tet4 mesh", field="mesh")
    inspect_boundary_topology(mesh, require_complete=True)
    for cell in mesh.tetra4:
        Tet4Geometry(mesh.points[cell])


def refine_tet4(
    mesh: Mesh, marked: ArrayLike, *, max_points: int | None = None
) -> Refinement:
    """Bisect each marked cell's longest edge, including its entire edge star.

    Edges are chosen on the INPUT mesh, ties by sorted node pair, processed in
    sorted order. Every incident volume and boundary simplex is split at the same
    midpoint, so no hanging nodes arise. Closure may refine unmarked cells.
    Cell/facet tags and overlapping named sets are inherited via parent rows.
    Geometry and material interfaces are preserved exactly for straight elements.
    No smoothing, coarsening, curved surfaces or uniform shape-quality bound.
    max_points limits active thermal DOFs; exceeded budget leaves input untouched.
    External IDs are removed: child identity is given by the returned parent maps.
    """
    validate_tet4(mesh)
    raw = np.asarray(marked)
    rows = indices(raw if raw.size else np.empty(0, dtype=np.int64), mesh.n_tetra4)
    edges: set[tuple[int, int]] = set()
    for row in rows:
        candidates = sorted(
            tuple(sorted(pair)) for pair in combinations(mesh.tetra4[row], 2)
        )
        edge = max(
            candidates,
            key=lambda e: float(np.linalg.norm(mesh.points[e[0]] - mesh.points[e[1]])),
        )
        edges.add((int(edge[0]), int(edge[1])))
    active_count = len(np.unique(mesh.tetra4))
    if max_points is not None and active_count + len(edges) > max_points:
        raise ValidationError("refinement exceeds DOF budget", field="budget")
    if not edges:
        return Refinement(mesh, np.arange(mesh.n_tetra4), np.arange(mesh.n_triangle3))
    points = mesh.points.tolist()
    cells = mesh.tetra4.tolist()
    facets = mesh.triangle3.tolist()
    cell_parents = list(range(mesh.n_tetra4))
    facet_parents = list(range(mesh.n_triangle3))

    def split(
        rows: list[list[int]], parents: list[int], edge: tuple[int, int], mid: int
    ) -> tuple[list[list[int]], list[int]]:
        children: list[list[int]] = []
        origins: list[int] = []
        a, b = edge
        for row, parent in zip(rows, parents, strict=True):
            if a in row and b in row:
                left, right = row.copy(), row.copy()
                left[left.index(b)] = mid
                right[right.index(a)] = mid
                children.extend((left, right))
                origins.extend((parent, parent))
            else:
                children.append(row)
                origins.append(parent)
        return children, origins

    for edge in sorted(edges):
        mid = len(points)
        points.append((mesh.points[edge[0]] / 2 + mesh.points[edge[1]] / 2).tolist())
        cells, cell_parents = split(cells, cell_parents, edge, mid)
        facets, facet_parents = split(facets, facet_parents, edge, mid)
    cp, fp = np.array(cell_parents), np.array(facet_parents)
    refined = Mesh(
        points=points,
        tetra4=cells,
        triangle3=facets,
        cell_tags=None if mesh.cell_tags is None else mesh.cell_tags[cp],
        facet_tags=None if mesh.facet_tags is None else mesh.facet_tags[fp],
        volume_sets={
            name: np.flatnonzero(np.isin(cp, rows))
            for name, rows in mesh.volume_sets.items()
        },
        facet_sets={
            name: np.flatnonzero(np.isin(fp, rows))
            for name, rows in mesh.facet_sets.items()
        },
    )
    validate_tet4(refined)
    return Refinement(refined, cp, fp)


def validate_refinement(original: Mesh, result: Refinement) -> None:
    """Reject invalid parentage, domain changes and lost physical groups.

    Children must lie in their parent simplex; summed child measures must equal
    each parent measure. Together with complete conforming boundary topology this
    checks the edge-bisection adapter contract, not arbitrary mesh intersection.
    """
    mesh = result.mesh
    validate_tet4(mesh)
    for kind, parents, before, after, tags, sets in (
        (
            "cell",
            result.cell_parents,
            original.tetra4,
            mesh.tetra4,
            "cell_tags",
            "volume_sets",
        ),
        (
            "facet",
            result.facet_parents,
            original.triangle3,
            mesh.triangle3,
            "facet_tags",
            "facet_sets",
        ),
    ):
        if (
            parents.shape != (len(after),)
            or parents.dtype.kind not in "iu"
            or np.any((parents < 0) | (parents >= len(before)))
        ):
            raise ValidationError("invalid refinement parent map", field=kind)
        old_tags, new_tags = getattr(original, tags), getattr(mesh, tags)
        if (old_tags is None) != (new_tags is None) or (
            old_tags is not None and not np.array_equal(new_tags, old_tags[parents])
        ):
            raise ValidationError("refinement changed tags", field=tags)
        old_sets, new_sets = getattr(original, sets), getattr(mesh, sets)
        if old_sets.keys() != new_sets.keys() or any(
            not np.array_equal(
                np.sort(new_sets[name]),
                np.flatnonzero(np.isin(parents, old_sets[name])),
            )
            for name in old_sets
        ):
            raise ValidationError("refinement changed named groups", field=sets)
        measures = np.zeros(len(before))
        expected = np.zeros(len(before))
        for row, nodes in enumerate(before):
            points = original.points[nodes]
            edges = (points[1:] - points[0]).T
            expected[row] = (
                abs(np.linalg.det(edges)) / 6
                if kind == "cell"
                else np.linalg.norm(np.cross(edges[:, 0], edges[:, 1])) / 2
            )
            # Precompute inverse/pseudoinverse once per parent.
            inverse = np.linalg.pinv(edges)
            for child in np.flatnonzero(parents == row):
                coords = mesh.points[after[child]]
                bary = (coords - points[0]) @ inverse.T
                reconstruction = bary @ edges.T + points[0]
                scale = max(float(np.linalg.norm(edges)), np.finfo(float).tiny)
                if (
                    np.any(bary < -1e-10)
                    or np.any(bary.sum(axis=1) > 1 + 1e-10)
                    or np.max(np.abs(reconstruction - coords)) > 1e-10 * scale
                ):
                    raise ValidationError("child lies outside its parent", field=kind)
                differences = (coords[1:] - coords[0]).T
                measures[row] += (
                    abs(np.linalg.det(differences)) / 6
                    if kind == "cell"
                    else np.linalg.norm(np.cross(differences[:, 0], differences[:, 1]))
                    / 2
                )
        if not np.allclose(measures, expected, rtol=1e-10, atol=0):
            raise ValidationError("refinement changed parent measure", field=kind)
