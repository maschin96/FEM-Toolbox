"""Topology and constraint preflight for declarative stationary models."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Literal, cast

import numpy as np

from fem_toolbox.boundary import Electrode
from fem_toolbox.errors import ModelSchemaError
from fem_toolbox.mesh import Mesh
from fem_toolbox.physics.heat import Convection
from fem_toolbox.solvers.elasticity import mechanical_dofs
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs


@dataclass(frozen=True, slots=True)
class PrescribedDof:
    """One resolved Dirichlet value with a path back to the model document."""

    field: Literal["T", "u", "phi"]
    dof: int
    node: int
    value: float
    path: str
    selection: str | None


def _components(mesh: Mesh) -> tuple[np.ndarray, ...]:
    """Return Tet4 node-connected components in sorted node order."""
    active = np.unique(mesh.tetra4)
    parent = {int(node): int(node) for node in active}

    def root(node: int) -> int:
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for cell in mesh.tetra4:
        for a, b in combinations((int(n) for n in cell), 2):
            parent[root(b)] = root(a)
    grouped: dict[int, list[int]] = {}
    for node in active:
        grouped.setdefault(root(int(node)), []).append(int(node))
    return tuple(
        np.asarray(nodes, dtype=np.int64)
        for nodes in sorted(grouped.values(), key=lambda group: group[0])
    )


def _mechanical_rank(
    mesh: Mesh, nodes: np.ndarray, fixed: set[int], coupled: bool
) -> int:
    """Use the same normalized six rigid modes and tolerance as the solver."""
    dofs = piezoelectric_dofs(mesh) if coupled else mechanical_dofs(mesh)
    xyz = mesh.points[nodes]
    xyz = xyz - xyz.mean(axis=0)
    xyz /= np.max(np.linalg.norm(xyz, axis=1))
    rigid = np.zeros((len(nodes) * 3, 6))
    for i, (x, y, z) in enumerate(xyz):
        rigid[3 * i : 3 * i + 3, :3] = np.eye(3)
        rigid[3 * i : 3 * i + 3, 3:] = [[0, z, -y], [-z, 0, x], [y, -x, 0]]
    constrained = np.isin(dofs.indices("u", nodes), list(fixed))
    return int(np.linalg.matrix_rank(rigid[constrained], tol=1e-10))


def validate_model_constraints(
    mesh: Mesh,
    physics: str,
    prescribed: tuple[PrescribedDof, ...],
    convection: dict[int, Convection],
    electrodes: tuple[Electrode, ...],
    surface_charges: dict[int, float],
) -> None:
    """Reject ambiguous values and missing references before solver assembly.

    Structural checks use the same component and rigid-mode semantics as the
    direct solvers. Load-dependent compatibility remains a solver diagnostic.
    """
    components = _components(mesh)
    component_of = {
        int(node): i for i, nodes in enumerate(components) for node in nodes
    }
    assigned: dict[tuple[str, int], PrescribedDof] = {}
    by_field: dict[str, list[PrescribedDof]] = {"T": [], "u": [], "phi": []}
    for item in prescribed:
        prior = assigned.get((item.field, item.dof))
        if prior is not None and prior.value != item.value:
            selection = f", set {item.selection!r}" if item.selection else ""
            raise ModelSchemaError(
                path=item.path,
                reason=(
                    f"field {item.field}, component {component_of[item.node]}, "
                    f"node {item.node}, DOF {item.dof}{selection}: value {item.value} "
                    f"conflicts with {prior.value} from {prior.path}"
                ),
            )
        assigned[(item.field, item.dof)] = item
        by_field[item.field].append(item)

    fixed_phi_nodes = {item.node for item in by_field["phi"]}
    electrode_nodes: dict[int, str] = {}
    names: set[str] = set()
    for i, electrode in enumerate(electrodes):
        path = f"boundary.electrodes[{i}]"
        if electrode.name in names:
            raise ModelSchemaError(
                path=f"{path}.name",
                reason=f"duplicate electrode name {electrode.name!r}",
            )
        names.add(electrode.name)
        facet_set = cast(str, electrode.facets)
        facets = mesh.facet_sets[facet_set]
        collision = set(int(face) for face in facets) & set(surface_charges)
        if collision:
            face = min(collision)
            raise ModelSchemaError(
                path=path,
                reason=f"electrode {electrode.name!r} on set {facet_set!r} "
                f"conflicts with prescribed surface charge on facet {face}",
            )
        nodes = np.unique(mesh.triangle3[facets])
        for node in nodes:
            node_id = int(node)
            if node_id in electrode_nodes:
                raise ModelSchemaError(
                    path=path,
                    reason=(
                        f"field phi, component {component_of[node_id]}, node {node_id}: "
                        f"electrode {electrode.name!r} on set {facet_set!r} "
                        f"overlaps {electrode_nodes[node_id]!r}"
                    ),
                )
            electrode_nodes[node_id] = electrode.name
            user_values = [item for item in by_field["phi"] if item.node == node_id]
            if electrode.charge is not None and user_values:
                raise ModelSchemaError(
                    path=path,
                    reason=f"field phi, component {component_of[node_id]}, node {node_id}: "
                    f"floating electrode {electrode.name!r} on set {facet_set!r} "
                    "conflicts with Dirichlet potential",
                )
            if electrode.potential is not None:
                for user in user_values:
                    if user.value != electrode.potential:
                        raise ModelSchemaError(
                            path=path,
                            reason=f"field phi, component {component_of[node_id]}, node {node_id}, "
                            f"DOF {user.dof}: electrode {electrode.name!r} on set "
                            f"{facet_set!r} potential "
                            f"conflicts with {user.path}",
                        )
                fixed_phi_nodes.add(node_id)

    thermal = physics in ("heat", "thermopiezoelectric")
    mechanical = physics in ("elasticity", "piezoelectric", "thermopiezoelectric")
    electric = physics in ("electrostatics", "piezoelectric", "thermopiezoelectric")
    fixed_t = {item.node for item in by_field["T"]}
    fixed_u = {item.dof for item in by_field["u"]}
    for component, nodes in enumerate(components):
        node_set = set(int(node) for node in nodes)
        node_ids = nodes.tolist()
        if thermal:
            anchored = bool(fixed_t & node_set)
            for face, condition in convection.items():
                if condition.coefficient > 0 and node_set.intersection(
                    int(node) for node in mesh.triangle3[face]
                ):
                    anchored = True
            if not anchored:
                raise ModelSchemaError(
                    path="boundary.T",
                    reason=f"field T, component {component}, nodes {node_ids}: "
                    "missing temperature reference or positive convection",
                )
        if mechanical:
            rank = _mechanical_rank(
                mesh,
                nodes,
                fixed_u,
                physics in ("piezoelectric", "thermopiezoelectric"),
            )
            if rank < 6:
                raise ModelSchemaError(
                    path="boundary.u",
                    reason=f"field u, component {component}, nodes {node_ids}: "
                    f"unrestrained rigid modes (support rank {rank}/6)",
                )
        if electric and not fixed_phi_nodes.intersection(node_set):
            raise ModelSchemaError(
                path="boundary.phi",
                reason=f"field phi, component {component}, nodes {node_ids}: "
                "missing potential reference; floating charge does not fix gauge",
            )
