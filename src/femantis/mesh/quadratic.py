"""Quadratic simplex node ordering and conforming edge incidence (no operators)."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from femantis.errors import ConnectivityError, MeshValidationError

TET10_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
TRI6_EDGES = ((0, 1), (1, 2), (2, 0))


def edge_nodes(nodes: NDArray[np.int64]) -> dict[tuple[int, ...], int]:
    """Map endpoint IDs to midpoint IDs in VTK/meshio local ordering."""
    edges = TET10_EDGES if len(nodes) == 10 else TRI6_EDGES
    start = 4 if len(nodes) == 10 else 3
    return {
        tuple(sorted((int(nodes[a]), int(nodes[b])))): int(nodes[start + i])
        for i, (a, b) in enumerate(edges)
    }


def validate_edge_incidence(
    cells: NDArray[np.int64], ids: NDArray[np.int64] | None
) -> None:
    """Reject hanging/duplicated edge DOFs and a midpoint reused as a vertex."""
    known: dict[tuple[int, ...], int] = {}
    reverse: dict[int, tuple[int, ...]] = {}
    corners = set(int(n) for n in cells[:, :4].flat)
    for index, nodes in enumerate(cells):
        for edge, midpoint in edge_nodes(nodes).items():
            if (
                midpoint in corners
                or known.get(edge, midpoint) != midpoint
                or reverse.get(midpoint, edge) != edge
            ):
                raise ConnectivityError(
                    entity_kind="tetra10",
                    entity_id=index if ids is None else int(ids[index]),
                    reason=f"nonconforming edge {edge} with midpoint {midpoint}",
                )
            known[edge], reverse[midpoint] = midpoint, edge


def require_straight(points: NDArray[np.float64], nodes: NDArray[np.int64]) -> None:
    """Require mid-edge coordinates within 1e-10 of each local edge length.

    Compare relative coordinates so absolute translation does not relax the
    tolerance. Curved geometry needs separate quadrature and is unsupported.
    """
    for (a, b), midpoint in edge_nodes(nodes).items():
        delta = points[b] - points[a]
        error = points[midpoint] - points[a] - 0.5 * delta
        if np.linalg.norm(error) > 1e-10 * np.linalg.norm(delta):
            raise MeshValidationError(
                field="geometry",
                reason=f"curved or incorrectly ordered quadratic edge {(a, b)}; expected midpoint",
            )


def reordered_tet10(
    nodes: NDArray[np.int64], corners: NDArray[np.int64]
) -> NDArray[np.int64]:
    """Transport all six edge nodes with a corner permutation."""
    edges = edge_nodes(nodes)
    return np.array(
        [
            *corners,
            *(
                edges[tuple(sorted((int(corners[a]), int(corners[b]))))]
                for a, b in TET10_EDGES
            ),
        ],
        dtype=np.int64,
    )
