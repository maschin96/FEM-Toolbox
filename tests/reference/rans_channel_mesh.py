"""Authored wall-mesh family for M12 verification; not a product mesh generator."""

import numpy as np

from fem_toolbox.cfd import FluidDomain, FluidMesh2D


def channel_domain(
    cells_per_half,
    *,
    beta=4.5,
    half_height_m=1.0,
    width_m=1.0,
    length_m=1.0,
    density_kg_m3=1.0,
    kinematic_viscosity_m2_s=1 / 587.19,
    pressure_drive_N_m3=1.0,
):
    nx = 6
    n = int(cells_per_half)
    if n < 4 or n != cells_per_half:
        raise ValueError("at least four cells per half required")
    half = np.expm1(beta * np.linspace(0, 1, n + 1)) / np.expm1(beta)
    y = half_height_m * np.r_[half, 2 - half[-2::-1]]
    ny = 2 * n
    points = np.array([(length_m * i / nx, yy) for yy in y for i in range(nx + 1)])

    def node(i, j):
        return j * (nx + 1) + i

    cells = [
        [node(i, j), node(i + 1, j), node(i + 1, j + 1), node(i, j + 1)]
        for j in range(ny)
        for i in range(nx)
    ]
    edges = [[node(i, 0), node(i + 1, 0)] for i in range(nx)]
    edges += [[node(nx, j), node(nx, j + 1)] for j in range(ny)]
    edges += [[node(i + 1, ny), node(i, ny)] for i in range(nx)]
    edges += [[node(0, j + 1), node(0, j)] for j in range(ny)]
    groups = {
        "bottom": np.arange(nx),
        "right": np.arange(nx, nx + ny),
        "top": np.arange(nx + ny, 2 * nx + ny),
        "left": np.arange(2 * nx + ny, 2 * (nx + ny)),
    }
    mesh = FluidMesh2D(points, cells, edges, groups)
    drop = pressure_drive_N_m3 * length_m
    return FluidDomain.from_document(
        {
            "fluid_domain_schema_version": 1,
            "mesh": "mesh.json",
            "width_m": width_m,
            "fluid": {
                "density_kg_m3": density_kg_m3,
                "viscosity_Pa_s": density_kg_m3 * kinematic_viscosity_m2_s,
            },
            "boundaries": {
                "bottom": {"type": "wall"},
                "top": {"type": "wall"},
                "left": {
                    "type": "periodic",
                    "partner": "right",
                    "translation_m": [length_m, 0],
                    "pressure_drop_Pa": drop,
                },
                "right": {
                    "type": "periodic",
                    "partner": "left",
                    "translation_m": [-length_m, 0],
                    "pressure_drop_Pa": -drop,
                },
            },
            "pressure_reference_cell": 0,
        },
        mesh=mesh,
    )
