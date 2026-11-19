"""Authored elastically guided lid benchmark for stationary FSI (#92).

Doc/71 fixes the benchmark: a parallel channel of height H0, length L and
width b carries a fixed volume flow Q from a developed inlet profile to a
pressure outlet. Its upper wall is a guided lid on a linear spring k_s.

The lid is a stiff Tet4 plate on a soft nu=0 spring layer. With x/z guided at
every node and the top fixed, the spring layer is uniaxial, so its stiffness
k_s = E_s b L / h_s holds exactly. The plate modulus is 1e9 E_s, so a beam on
an elastic foundation with beta*L < 1e-2 distributes the linear pressure
almost uniformly. This is an authored test structure, not a product mesh
generator. No production CFD/FSI routine is used for the analytic values.
"""

from collections import Counter
from dataclasses import dataclass
from itertools import permutations, product

import numpy as np
from scipy.optimize import brentq

from fem_toolbox.cfd import FluidDomain, FluidMesh2D
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model

STIFFNESS_RATIO = 1e9


def box_tetra(xs, ys, zs):
    """Kuhn subdivision of a tensor grid; returns points, Tet4 and boundary."""
    xs, ys, zs = (np.asarray(v, float) for v in (xs, ys, zs))
    shape = (len(xs), len(ys), len(zs))
    points = np.array([(x, y, z) for x in xs for y in ys for z in zs])

    def index(p):
        return (p[0] * shape[1] + p[1]) * shape[2] + p[2]

    cells = []
    for base in product(*(range(n - 1) for n in shape)):
        for order in permutations(range(3)):
            p = np.array(base)
            cell = [index(p)]
            for axis in order:
                p = p.copy()
                p[axis] += 1
                cell.append(index(p))
            cells.append(cell)
    cells = np.array(cells)
    faces = Counter(
        tuple(sorted(np.delete(cell, opposite)))
        for cell in cells
        for opposite in range(4)
    )
    boundary = np.array([face for face, count in faces.items() if count == 1])
    return points, cells, boundary


def lid_structure(
    *,
    nx,
    height_m,
    length_m=1.0,
    width_m=1.0,
    spring_stiffness_N_m,
    lid_thickness_m=0.1,
    spring_height_m=0.1,
    layers=1,
    nz=1,
):
    """Guided lid over y=height_m; returns the Model and the wet facet set name."""
    xs = np.linspace(0, length_m, nx + 1)
    lid_y = height_m + lid_thickness_m * np.linspace(0, 1, layers + 1)
    spring_y = lid_y[-1] + spring_height_m * np.linspace(0, 1, layers + 1)[1:]
    ys = np.r_[lid_y, spring_y]
    zs = np.linspace(0, width_m, nz + 1)
    points, cells, boundary = box_tetra(xs, ys, zs)
    centroid_y = points[cells, 1].mean(axis=1)
    lid = np.flatnonzero(centroid_y < lid_y[-1])
    spring = np.flatnonzero(centroid_y > lid_y[-1])
    wet = np.flatnonzero(np.all(points[boundary, 1] == height_m, axis=1))
    mesh = Mesh(
        points=points,
        tetra4=cells,
        triangle3=boundary,
        facet_sets={"wet": wet},
        volume_sets={"lid": lid, "spring": spring},
    )
    spring_young = spring_stiffness_N_m * spring_height_m / (length_m * width_m)
    top = np.flatnonzero(points[:, 1] == ys[-1]).tolist()
    document = {
        "schema_version": 1,
        "units": "SI",
        "physics": "elasticity",
        "mesh": {"path": "lid.msh", "length_scale_to_m": 1.0},
        "materials": {
            "u": {
                "lid": {
                    "young": STIFFNESS_RATIO * spring_young,
                    "poisson": 0.0,
                    "unit": "Pa",
                },
                "spring": {"young": spring_young, "poisson": 0.0, "unit": "Pa"},
            }
        },
        "regions": {"u": {"lid": "lid", "spring": "spring"}},
        "boundary": {
            "dirichlet": [
                {
                    "field": "u",
                    "nodes": list(range(len(points))),
                    "components": ["x", "z"],
                    "value": 0.0,
                    "unit": "m",
                },
                {
                    "field": "u",
                    "nodes": top,
                    "components": ["y"],
                    "value": 0.0,
                    "unit": "m",
                },
            ]
        },
    }
    return Model(mesh=mesh, document=document), "wet"


def lid_channel(
    *,
    nx,
    y_nodes_m,
    inlet_profile,
    length_m=1.0,
    width_m=1.0,
    density_kg_m3=1.0,
    viscosity_Pa_s=1.0,
):
    """Inlet/outlet channel; inlet_profile(y_lo, y_hi) -> (U[nface,2], a[nface])."""
    y = np.asarray(y_nodes_m, float)
    ny = len(y) - 1
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
    # Inlet faces follow the left group order (top to bottom).
    velocity, nu_tilde = inlet_profile(y[:-1][::-1], y[1:][::-1])
    return FluidDomain.from_document(
        {
            "fluid_domain_schema_version": 1,
            "mesh": "lid.json",
            "width_m": width_m,
            "fluid": {
                "density_kg_m3": density_kg_m3,
                "viscosity_Pa_s": viscosity_Pa_s,
            },
            "boundaries": {
                "bottom": {"type": "wall"},
                "top": {"type": "wall", "fsi": True},
                "left": {
                    "type": "inlet",
                    "velocity_m_s": np.asarray(velocity).tolist(),
                    "nu_tilde_m2_s": np.asarray(nu_tilde).tolist(),
                },
                "right": {"type": "pressure_outlet", "pressure_Pa": 0.0},
            },
        },
        mesh=mesh,
    )


@dataclass(frozen=True)
class LaminarLid:
    """Doc/71 laminar lid reference: Poiseuille pressure on a rigid lid."""

    flow_m3_s: float = 1e-3
    height_m: float = 0.1
    length_m: float = 1.0
    width_m: float = 1.0
    viscosity_Pa_s: float = 1.0
    spring_stiffness_N_m: float = 1e4

    def pressure_drop_Pa(self, height_m):
        return (
            12
            * self.viscosity_Pa_s
            * self.length_m
            * self.flow_m3_s
            / (self.width_m * height_m**3)
        )

    def lid_force_N(self, height_m):
        # Linear pressure from Delta p to the zero outlet: mean Delta p / 2.
        return self.width_m * self.length_m * self.pressure_drop_Pa(height_m) / 2

    @property
    def one_way_displacement_m(self):
        return self.lid_force_N(self.height_m) / self.spring_stiffness_N_m

    @property
    def coupled_displacement_m(self):
        """Positive root of d (H0+d)^3 = 6 mu L^2 Q / k_s."""
        rhs = (
            6
            * self.viscosity_Pa_s
            * self.length_m**2
            * self.flow_m3_s
            / (self.width_m * self.spring_stiffness_N_m)
        )
        return brentq(
            lambda d: d * (self.height_m + d) ** 3 - rhs,
            0,
            self.one_way_displacement_m,
            xtol=1e-16,
            rtol=1e-15,
        )

    def inlet_profile(self, y_lo, y_hi):
        """Face-averaged parabola, so the discrete inlet carries exactly Q."""
        h = self.height_m
        s0, s1 = np.asarray(y_lo) / h, np.asarray(y_hi) / h
        mean = ((s1**2 - s0**2) / 2 - (s1**3 - s0**3) / 3) / (s1 - s0)
        u = 6 * self.flow_m3_s / (self.width_m * h) * mean
        return np.column_stack([u, np.zeros_like(u)]), np.zeros_like(u)


def graded_wall_nodes(cells_per_half, height_m, beta=4.5):
    """Symmetric tanh-free exponential wall grading as in rans_channel_mesh."""
    half = np.expm1(beta * np.linspace(0, 1, cells_per_half + 1)) / np.expm1(beta)
    return height_m / 2 * np.r_[half, 2 - half[-2::-1]]


@dataclass(frozen=True)
class SaLid:
    """Doc/71 SA lid: fixed Q from the independent 1D SA BVP at H0.

    nu=(0.05 m)(0.01 m/s)/587.19, f_x=0.002 N/m^3 at H0 gives Re_tau=587.19.
    For a developed channel at fixed Q and nu, Q/(b nu)=2 Re_tau U_b+(Re_tau)
    fixes Re_tau independently of H, so Delta p(H)=Delta p0 (H0/H)^3 and the
    lid root has the laminar form. pressure_drop_Pa(H) still solves the BVP
    at H for the fixed Q, so this similarity is checked, not assumed.
    """

    height_m: float = 0.1
    length_m: float = 1.0
    width_m: float = 1.0
    density_kg_m3: float = 1.0
    kinematic_viscosity_m2_s: float = 0.05 * 0.01 / 587.19
    pressure_drive_N_m3: float = 0.002
    spring_stiffness_N_m: float = 1.0

    def _reference(self, height_m, drive):
        from tests.reference.sa_channel_1d import solve_channel_reference

        return solve_channel_reference(
            half_height_m=height_m / 2,
            density_kg_m3=self.density_kg_m3,
            kinematic_viscosity_m2_s=self.kinematic_viscosity_m2_s,
            pressure_drive_N_m3=drive,
        )

    @property
    def base(self):
        return self._reference(self.height_m, self.pressure_drive_N_m3)

    @property
    def flow_m3_s(self):
        return self.base.flow_per_span_m2_s * self.width_m

    def pressure_drop_Pa(self, height_m, flow_m3_s=None):
        """Independent BVP: drive giving the fixed Q at height H, times L."""
        q = self.flow_m3_s if flow_m3_s is None else flow_m3_s
        guess = self.pressure_drive_N_m3 * (self.height_m / height_m) ** 3
        drive = brentq(
            lambda f: (
                self._reference(height_m, f).flow_per_span_m2_s * self.width_m - q
            ),
            0.9 * guess,
            1.1 * guess,
            xtol=1e-14 * guess,
            rtol=1e-12,
        )
        return drive * self.length_m

    def lid_force_N(self, height_m):
        return self.width_m * self.length_m * self.pressure_drop_Pa(height_m) / 2

    @property
    def one_way_displacement_m(self):
        f0 = self.width_m * self.length_m * self.pressure_drive_N_m3 * self.length_m
        return f0 / 2 / self.spring_stiffness_N_m

    @property
    def coupled_displacement_m(self):
        h0, d1 = self.height_m, self.one_way_displacement_m
        return brentq(
            lambda d: d * (h0 + d) ** 3 - d1 * h0**3, 0, d1, xtol=1e-16, rtol=1e-15
        )

    def inlet_profile(self, y_lo, y_hi, samples=64):
        """Face averages of the BVP U and nu_tilde at H0 (Gauss-Legendre)."""
        base = self.base
        t, w = np.polynomial.legendre.leggauss(samples)
        lo, hi = np.asarray(y_lo), np.asarray(y_hi)
        y = 0.5 * (hi + lo)[:, None] + 0.5 * (hi - lo)[:, None] * t
        u, a = base.evaluate(y.ravel())
        u, a = u.reshape(y.shape) @ w / 2, a.reshape(y.shape) @ w / 2
        return np.column_stack([u, np.zeros_like(u)]), a


def _gmsh_structure(path, *, nx, height_m, lid_thickness_m=0.1, spring_height_m=0.1):
    """Gmsh 2.2: volumes lid/spring, disjoint facets wet/top/side."""
    xs = np.linspace(0, 1.0, nx + 1)
    ys = height_m + np.array([0.0, lid_thickness_m, lid_thickness_m + spring_height_m])
    points, cells, boundary = box_tetra(xs, ys, [0.0, 1.0])
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$PhysicalNames", "5"]
    lines += ['3 1 "lid"', '3 2 "spring"', '2 3 "wet"', '2 4 "top"', '2 5 "side"']
    lines += ["$EndPhysicalNames", "$Nodes", str(len(points))]
    lines += [
        f"{i + 1} {p[0]:.17g} {p[1]:.17g} {p[2]:.17g}" for i, p in enumerate(points)
    ]
    lines += ["$EndNodes", "$Elements", str(len(cells) + len(boundary))]
    for i, cell in enumerate(cells):
        tag = 1 if points[cell, 1].mean() < ys[1] else 2
        lines.append(f"{i + 1} 4 2 {tag} {tag} " + " ".join(str(n + 1) for n in cell))
    for j, face in enumerate(sorted(map(tuple, boundary))):
        y = points[list(face), 1]
        tag = 3 if np.all(y == ys[0]) else 4 if np.all(y == ys[-1]) else 5
        lines.append(
            f"{len(cells) + j + 1} 2 2 {tag} {tag} "
            + " ".join(str(n + 1) for n in face)
        )
    lines.append("$EndElements")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _structure_toml(spring_stiffness_N_m, spring_height_m=0.1):
    spring = spring_stiffness_N_m * spring_height_m
    lid = STIFFNESS_RATIO * spring
    guide = "".join(
        f"""[[boundary.dirichlet]]
field = "u"
set = "{name}"
components = {components}
value = 0.0
unit = "m"
"""
        for name, components in (
            ("wet", '["x", "z"]'),
            ("side", '["x", "z"]'),
            ("top", '["x", "y", "z"]'),
        )
    )
    return f"""schema_version = 1
units = "SI"
physics = "elasticity"
[mesh]
path = "structure.msh"
length_scale_to_m = 1.0
[materials.u.lid]
young = {lid!r}
poisson = 0.0
unit = "Pa"
[materials.u.spring]
young = {spring!r}
poisson = 0.0
unit = "Pa"
[regions.u]
lid = "lid"
spring = "spring"
{guide}"""


EXAMPLES = {
    # Fast installation-level cases: the coarsest #92 mesh pair of each lid.
    "laminar": {"reference": LaminarLid, "nx": 10, "ny": 4, "uref": 1.0},
    "sa": {"reference": SaLid, "nx": 4, "cells_per_half": 32, "uref": 0.2},
}


def write_example(directory, case):
    """Write the prepared, versioned example inputs of one lid case."""
    import json
    from pathlib import Path

    from fem_toolbox.cfd.mesh_motion import _domain_document

    spec = EXAMPLES[case]
    reference = spec["reference"]()
    target = Path(directory)
    target.mkdir(parents=True, exist_ok=False)
    if case == "laminar":
        y = np.linspace(0, reference.height_m, spec["ny"] + 1)
        viscosity = reference.viscosity_Pa_s
    else:
        y = graded_wall_nodes(spec["cells_per_half"], reference.height_m)
        viscosity = reference.density_kg_m3 * reference.kinematic_viscosity_m2_s
    domain = lid_channel(
        nx=spec["nx"],
        y_nodes_m=y,
        inlet_profile=reference.inlet_profile,
        viscosity_Pa_s=viscosity,
    )
    document = _domain_document(domain)
    document["mesh"] = "mesh.json"

    def dump(name, value):
        (target / name).write_text(
            json.dumps(value, indent=1, allow_nan=False) + "\n", encoding="utf-8"
        )

    dump("mesh.json", domain.mesh.to_document())
    dump("domain.json", document)
    _gmsh_structure(
        target / "structure.msh", nx=spec["nx"], height_m=reference.height_m
    )
    (target / "structure.toml").write_text(
        _structure_toml(reference.spring_stiffness_N_m), encoding="utf-8"
    )
    fluid = 'domain = "domain.json"\nmodel = "laminar"\n'
    if case == "sa":
        nu = reference.kinematic_viscosity_m2_s
        start = np.maximum(
            reference.base.evaluate(domain.mesh.cell_centres_m[:, 1])[1], 3 * nu
        )
        dump("nu_tilde.json", {"nu_tilde_m2_s": start.tolist()})
        fluid = (
            'domain = "domain.json"\nmodel = "sa-noft2"\n'
            'initial_nu_tilde = "nu_tilde.json"\n'
            "[fluid.options]\nsa_relaxation = 0.1\n"
        )
    (target / "fsi.toml").write_text(
        f"""fsi_schema_version = 1
[fluid]
{fluid}[structure]
model = "structure.toml"
facet_set = "wet"
[coupling]
wall = "top"
length_scale_m = {reference.height_m!r}
velocity_scale_m_s = {spec["uref"]!r}
external_pressure_Pa = 0.0
relaxation = 1.0
max_coupling_iterations = 200
divergence_window = 10
span_origin_m = 0.0
preserve_inlet_flow = true
""",
        encoding="utf-8",
    )
    return target / "fsi.toml"
