"""Stationary harmonic mesh displacement on a fixed Quad4 reference domain.

This solves two scalar Laplace problems, never a physical time/ALE equation.
All displacements are absolute SI metres relative to the supplied reference.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Any, NoReturn

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import sparse  # type: ignore[import-untyped]
from scipy.sparse.linalg import splu  # type: ignore[import-untyped]

from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import FluidMesh2D, _float_array, _freeze_float
from fem_toolbox.errors import ValidationError


def _hash(document: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            document, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def _domain_document(domain: FluidDomain) -> dict[str, Any]:
    document: dict[str, Any] = {
        "fluid_domain_schema_version": 1,
        "mesh": "mesh.json",
        "width_m": domain.width_m,
        "fluid": {
            "density_kg_m3": domain.fluid.density_kg_m3,
            "viscosity_Pa_s": domain.fluid.viscosity_Pa_s,
        },
        "boundaries": {
            name: {
                k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in bc.items()
            }
            for name, bc in domain.boundaries.items()
        },
    }
    if domain.pressure_reference_cell is not None:
        document["pressure_reference_cell"] = domain.pressure_reference_cell
    return document


class MeshMotionError(ValidationError):
    """Rejected stationary motion with stage, reference identity and last field."""

    def __init__(
        self,
        message: str,
        *,
        stage: str,
        reference_mesh_sha256: str,
        displacement_m: NDArray[np.float64],
    ) -> None:
        self.stage = stage
        self.reference_mesh_sha256 = reference_mesh_sha256
        self.displacement_m = _freeze_float(displacement_m)
        super().__init__(f"mesh motion {stage}: {message}", field="mesh_motion")


@dataclass(frozen=True, slots=True)
class MeshMotionResult:
    """New immutable domain and absolute harmonic displacement [npoint,2] m.

    Laplace residuals contain free reduced DOFs (periodic pairs share a DOF).
    Boundary group/connectivity identities are retained. Hashes distinguish
    reference/current meshes, full domains and prescribed interface traces.
    """

    domain: FluidDomain
    displacement_m: NDArray[np.float64]
    laplace_residual_m: NDArray[np.float64]
    maximum_relative_laplace_residual: float
    minimum_area_ratio: float
    reference_mesh_sha256: str
    mesh_sha256: str
    reference_domain_sha256: str
    domain_sha256: str
    interface_sha256: str
    input_sha256: str
    settings: Mapping[str, Any]


def _laplace(mesh: FluidMesh2D) -> Any:
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    for cell in mesh.quad4:
        xy = mesh.points[cell] - mesh.points[cell[0]]
        local = np.zeros((4, 4))
        for xi in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
            for eta in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
                derivative = 0.25 * np.array(
                    [
                        [-(1 - eta), -(1 - xi)],
                        [1 - eta, -(1 + xi)],
                        [1 + eta, 1 + xi],
                        [-(1 + eta), 1 - xi],
                    ]
                )
                jacobian = xy.T @ derivative
                gradient = derivative @ np.linalg.inv(jacobian)
                local += np.linalg.det(jacobian) * gradient @ gradient.T
        rows.extend(np.repeat(cell, 4))
        columns.extend(np.tile(cell, 4))
        values.extend(local.ravel())
    return sparse.csr_matrix(
        (values, (rows, columns)), shape=(len(mesh.points), len(mesh.points))
    )


def move_fluid_domain(
    reference: FluidDomain,
    interface_displacements_m: Mapping[str, ArrayLike],
    *,
    length_scale_m: float,
    preserve_inlet_flow: bool = True,
    max_nodes: int = 50000,
) -> MeshMotionResult:
    """Harmonic absolute motion with prescribed FSI endpoint traces.

    Each explicitly marked FSI wall needs [nsegment,2 endpoints,2] metres,
    in its line2 endpoint order. Other walls stay fixed. Inlet/outlet chains
    interpolate their two wall-corner displacements; other open boundaries
    have natural zero normal gradient. Explicit periodic node pairs are tied.
    Fixed inlet volume flow is preserved by scaling each original velocity
    profile against its new signed face areas; original SA inlet values remain.

    Reject max|d|/Lref>0.02, inverted/poor-quality geometry, area ratios<0.5,
    conflicting corners, noncongruent periodic motion and resource exhaustion.
    The new FluidDomain rebuilds wall distances, growth and all BC preflight.
    Final y+ still requires a subsequent fluid solve, not a mesh-only estimate.
    """
    mesh = reference.mesh
    source_hash = _hash(mesh.to_document())
    displacement = np.zeros_like(mesh.points)

    def fail(stage: str, message: str) -> NoReturn:
        raise MeshMotionError(
            message,
            stage=stage,
            reference_mesh_sha256=source_hash,
            displacement_m=displacement,
        )

    try:
        scale_array = _float_array(length_scale_m, "length_scale_m")
    except ValidationError as error:
        fail("input", str(error))
    if scale_array.ndim != 0 or scale_array <= 0:
        fail("input", "positive scalar length_scale_m required")
    scale = float(scale_array)
    if type(preserve_inlet_flow) is not bool:
        fail("input", "preserve_inlet_flow must be boolean")
    if type(max_nodes) is not int or max_nodes <= 0:
        fail("input", "positive integer node budget required")
    if len(mesh.points) > max_nodes:
        fail("resource", f"node budget {max_nodes} exhausted")
    expected = {
        name
        for name, bc in reference.boundaries.items()
        if bc["type"] == "wall" and bc["fsi"]
    }
    if (
        not isinstance(interface_displacements_m, Mapping)
        or set(interface_displacements_m) != expected
        or not expected
    ):
        fail("input", f"interface traces required for exactly {sorted(expected)}")
    prescribed: dict[int, NDArray[np.float64]] = {}
    canonical = {}
    tolerance = 1e-12 * scale

    def assign(node: int, value: NDArray[np.float64], label: str) -> None:
        if node in prescribed and np.linalg.norm(prescribed[node] - value) > tolerance:
            fail("boundary", f"conflicting corner/node {node}: {label}")
        prescribed[node] = value.copy()
        displacement[node] = value

    for name in sorted(expected):
        edges = mesh.line2[mesh.boundary_sets[name]]
        try:
            data = _float_array(interface_displacements_m[name], f"interface.{name}")
        except ValidationError as error:
            fail("input", str(error))
        if data.shape != (len(edges), 2, 2):
            fail("input", f"interface {name}: expected ({len(edges)},2,2)")
        canonical[name] = data.tolist()
        for edge, values_at_edge in zip(edges, data, strict=True):
            for node, value in zip(edge, values_at_edge, strict=True):
                assign(int(node), value, name)
    if np.max(np.linalg.norm(displacement, axis=1)) > 0.02 * scale:
        fail("displacement", "maximum displacement exceeds 0.02 Lref")
    for name, bc in reference.boundaries.items():
        if bc["type"] == "wall" and not bc["fsi"]:
            for node in np.unique(mesh.line2[mesh.boundary_sets[name]]):
                assign(int(node), np.zeros(2), name)
    for name, bc in reference.boundaries.items():
        if bc["type"] not in ("inlet", "pressure_outlet"):
            continue
        edges = mesh.line2[mesh.boundary_sets[name]]
        neighbours: dict[int, list[int]] = {}
        for a, b in edges:
            neighbours.setdefault(int(a), []).append(int(b))
            neighbours.setdefault(int(b), []).append(int(a))
        endpoints = sorted(
            n for n, adjacent in neighbours.items() if len(adjacent) == 1
        )
        if len(endpoints) != 2 or any(
            len(adjacent) > 2 for adjacent in neighbours.values()
        ):
            fail("boundary", f"{name}: inlet/outlet must be one simple open chain")
        if any(n not in prescribed for n in endpoints):
            fail(
                "boundary", f"{name}: both endpoint wall corners need displacement data"
            )
        chain, previous, current = [], -1, endpoints[0]
        while True:
            chain.append(current)
            following = [n for n in neighbours[current] if n != previous]
            if not following:
                break
            previous, current = current, following[0]
            if current in chain:
                fail("boundary", f"{name}: loop in inlet/outlet chain")
        if len(chain) != len(neighbours):
            fail("boundary", f"{name}: disconnected inlet/outlet chain")
        arc = np.r_[
            0.0, np.cumsum(np.linalg.norm(np.diff(mesh.points[chain], axis=0), axis=1))
        ]
        first, last = prescribed[chain[0]].copy(), prescribed[chain[-1]].copy()
        for node, fraction in zip(chain, arc / arc[-1], strict=True):
            assign(node, (1 - fraction) * first + fraction * last, name)

    parent = np.arange(len(mesh.points))

    def root(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = int(parent[node])
        return node

    for face, partner in reference.periodic_partners.items():
        if face > partner:
            continue
        a, b = mesh.line2[face], mesh.line2[partner]
        translation = mesh.points[b].mean(axis=0) - mesh.points[a].mean(axis=0)
        for node in a:
            distance = np.linalg.norm(
                mesh.points[b] - (mesh.points[node] + translation), axis=1
            )
            if distance.min() > 1e-10 * scale:
                fail("periodic", "periodic endpoint pairing is not congruent")
            other = int(b[int(np.argmin(distance))])
            ra, rb = root(int(node)), root(other)
            parent[ra] = rb
    roots = np.array([root(n) for n in range(len(mesh.points))])
    _, reduced_ids = np.unique(roots, return_inverse=True)
    count = int(reduced_ids.max()) + 1
    reduction = sparse.csr_matrix(
        (np.ones(len(roots)), (np.arange(len(roots)), reduced_ids)),
        shape=(len(roots), count),
    )
    reduced_values = np.zeros((count, 2))
    fixed = np.zeros(count, dtype=bool)
    for node, value in prescribed.items():
        row = reduced_ids[node]
        if fixed[row] and np.linalg.norm(reduced_values[row] - value) > tolerance:
            fail("periodic", f"inconsistent periodic displacement at node {node}")
        fixed[row], reduced_values[row] = True, value
    if not np.any(fixed):
        fail("boundary", "Laplace motion has no displacement anchor")
    matrix = (reduction.T @ _laplace(mesh) @ reduction).tocsr()
    free = ~fixed
    try:
        if np.any(free):
            rhs = -(matrix[free][:, fixed] @ reduced_values[fixed])
            reduced_values[free] = splu(matrix[free][:, free].tocsc()).solve(rhs)
    except RuntimeError as error:
        fail("linear", str(error))
    displacement = np.asarray(reduction @ reduced_values)
    if not np.all(np.isfinite(displacement)):
        fail("linear", "nonfinite harmonic displacement")
    if np.max(np.linalg.norm(displacement, axis=1)) > 0.02 * scale:
        fail("displacement", "maximum displacement exceeds 0.02 Lref")
    residual = np.asarray(matrix @ reduced_values)[free]
    magnitude = np.asarray(abs(matrix) @ abs(reduced_values))[free]
    relative = float(
        np.max(abs(residual), initial=0)
        / max(float(np.max(magnitude, initial=0)), 1e-14 * scale)
    )
    if relative > 1e-12:
        fail("linear", f"Laplace backward-error gate failed: {relative}")
    try:
        moved = FluidMesh2D(
            mesh.points + displacement,
            mesh.quad4,
            mesh.line2,
            mesh.boundary_sets,
            mesh.volume_sets,
        )
    except ValidationError as error:
        fail("geometry", str(error))
    minimum_ratio = float(np.min(moved.cell_areas_m2 / mesh.cell_areas_m2))
    if minimum_ratio < 0.5:
        fail("area", f"minimum cell area ratio {minimum_ratio} below 0.5")
    document = _domain_document(reference)
    if preserve_inlet_flow:
        for name, bc in reference.boundaries.items():
            if bc["type"] != "inlet":
                continue
            ids = mesh.boundary_sets[name]
            old_area = mesh.face_lengths_m[mesh.boundary_face_indices[ids]]
            new_area = moved.face_lengths_m[moved.boundary_face_indices[ids]]
            old_flux = float(
                np.einsum("ij,ij->i", bc["velocity_m_s"], mesh.boundary_normals[ids])
                @ old_area
            )
            new_flux = float(
                np.einsum("ij,ij->i", bc["velocity_m_s"], moved.boundary_normals[ids])
                @ new_area
            )
            if (
                new_flux == 0
                or old_flux / new_flux <= 0
                or not np.isfinite(old_flux / new_flux)
            ):
                fail(
                    "boundary",
                    f"{name}: fixed-flow inlet profile has no admissible scaling",
                )
            document["boundaries"][name]["velocity_m_s"] = (
                bc["velocity_m_s"] * (old_flux / new_flux)
            ).tolist()
    try:
        domain = FluidDomain.from_document(document, mesh=moved)
    except ValidationError as error:
        fail("quality", str(error))
    reference_domain_hash = _hash(
        {"mesh": mesh.to_document(), "domain": _domain_document(reference)}
    )
    interface_hash = _hash(canonical)
    settings = {
        "length_scale_m": scale,
        "preserve_inlet_flow": preserve_inlet_flow,
        "max_nodes": max_nodes,
        "maximum_displacement_ratio": 0.02,
        "minimum_area_ratio": 0.5,
        "maximum_relative_laplace_residual": 1e-12,
        "corner_tolerance_m": tolerance,
        "periodic_tolerance_m": 1e-10 * scale,
        "solver": "Quad4-Laplace-2x2-Gauss-SparseLU",
    }
    return MeshMotionResult(
        domain,
        _freeze_float(displacement),
        _freeze_float(residual),
        relative,
        minimum_ratio,
        source_hash,
        _hash(moved.to_document()),
        reference_domain_hash,
        _hash({"mesh": moved.to_document(), "domain": document}),
        interface_hash,
        _hash(
            {
                "reference_domain_sha256": reference_domain_hash,
                "interface_sha256": interface_hash,
                "settings": settings,
            }
        ),
        MappingProxyType(settings),
    )
