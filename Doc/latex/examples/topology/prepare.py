"""Deterministic six-tetrahedron-per-brick reference, without a mesher dependency."""

from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path

import numpy as np


def make_inputs(output, divisions=2, radius=0.3):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    nx, ny, nz = 2 * divisions, divisions, divisions
    points = np.array(
        list(
            product(
                np.linspace(0, 2, nx + 1),
                np.linspace(0, 1, ny + 1),
                np.linspace(0, 1, nz + 1),
            )
        )
    )

    def node(p):
        return (p[0] * (ny + 1) + p[1]) * (nz + 1) + p[2]

    cells = []
    for base in product(range(nx), range(ny), range(nz)):
        for order in permutations(range(3)):
            position = np.array(base)
            cell = [node(position)]
            for axis in order:
                position = position.copy()
                position[axis] += 1
                cell.append(node(position))
            cells.append(cell)
    cells = np.array(cells)
    counts = Counter(
        tuple(sorted(face)) for cell in cells for face in combinations(cell, 3)
    )
    faces = sorted(face for face, count in counts.items() if count == 1)
    lines = [
        "$MeshFormat",
        "2.2 0 8",
        "$EndMeshFormat",
        "$PhysicalNames",
        "5",
        '3 1 "design"',
        '3 2 "passive"',
        '2 3 "clamp"',
        '2 4 "end"',
        '2 5 "side"',
        "$EndPhysicalNames",
        "$Nodes",
        str(len(points)),
    ]
    lines += [
        f"{i + 1} {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}" for i, p in enumerate(points)
    ]
    lines += ["$EndNodes", "$Elements", str(len(cells) + len(faces))]
    for i, cell in enumerate(cells):
        x = points[cell, 0]
        tag = 2 if np.count_nonzero(x == 0) >= 3 or np.count_nonzero(x == 2) >= 3 else 1
        lines.append(
            f"{i + 1} 4 2 {tag} {tag} " + " ".join(str(int(n) + 1) for n in cell)
        )
    for j, face in enumerate(faces):
        x = points[list(face), 0]
        tag = 3 if np.all(x == 0) else 4 if np.all(x == 2) else 5
        lines.append(
            f"{len(cells) + j + 1} 2 2 {tag} {tag} "
            + " ".join(str(n + 1) for n in face)
        )
    lines += ["$EndElements"]
    (output / "mesh.msh").write_text("\n".join(lines) + "\n")
    for name, force in [("y", "[0.0, -1.0, 0.0]"), ("z", "[0.0, 0.0, -1.0]")]:
        (output / f"{name}.toml").write_text(f"""schema_version = 1
units = "SI"
physics = "elasticity"
[mesh]
path = "mesh.msh"
length_scale_to_m = 1.0
[materials.u.solid]
young = 1000.0
poisson = 0.3
unit = "Pa"
[regions.u]
design = "solid"
passive = "solid"
[[boundary.dirichlet]]
field = "u"
set = "clamp"
value = 0.0
unit = "m"
[[boundary.traction]]
set = "end"
value = {force}
unit = "N/m^2"
""")
    (output / "study.toml").write_text(f"""schema_version = 1
models = ["y.toml", "z.toml"]
weights = [0.5, 0.5]
[design]
volume_fraction = 0.6
filter_radius = {radius}
penalty = 3.0
minimum_stiffness = 0.001
solid_sets = ["passive"]
[optimizer]
max_iterations = 100
max_evaluations = 300
tolerance = 1e-7
""")
    return output / "study.toml"


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--divisions", type=int, default=2)
    parser.add_argument("--radius", type=float, default=0.3)
    args = parser.parse_args()
    if args.divisions < 1 or args.radius <= 0:
        parser.error("positive divisions and radius required")
    print(make_inputs(args.output, args.divisions, args.radius))
