"""Conforming midpoint enrichment of independent test meshes."""

import numpy as np

from fem_toolbox.mesh import Mesh


def elevate(mesh):
    points = mesh.points.tolist()
    edges = {}
    cells = []
    for cell in mesh.tetra4:
        mids = []
        for a, b in ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3)):
            edge = tuple(sorted((int(cell[a]), int(cell[b]))))
            if edge not in edges:
                edges[edge] = len(points)
                points.append(
                    ((mesh.points[edge[0]] + mesh.points[edge[1]]) / 2).tolist()
                )
            mids.append(edges[edge])
        cells.append([*cell, *mids])
    facets = [
        [
            *face,
            *(
                edges[tuple(sorted((int(face[a]), int(face[b]))))]
                for a, b in ((0, 1), (1, 2), (2, 0))
            ),
        ]
        for face in mesh.triangle3
    ]
    return Mesh(
        points=np.array(points),
        tetra10=np.array(cells),
        triangle6=np.array(facets),
        volume_sets={"body": np.arange(mesh.n_cells)},
        facet_sets={
            name: np.flatnonzero(np.all(mesh.points[mesh.triangle3, 0] == x, axis=1))
            for name, x in (("left", 0), ("right", 1))
        },
    )
