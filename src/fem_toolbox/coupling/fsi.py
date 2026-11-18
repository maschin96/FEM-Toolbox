"""Stationary partitioned Dirichlet-Neumann fluid/structure iteration (#91).

Doc/71 section 5: d^k -> harmonic fluid geometry -> full fluid solve ->
conservative loads on the current overlap -> linear structure on the
unchanged reference mesh -> r_d = d_solved - d^k -> d^(k+1) = d^k + omega r_d.
Iteration numbers are no physical time; there is no wall velocity or ALE flux.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from itertools import pairwise
from types import MappingProxyType
from typing import Any, NoReturn, cast

import numpy as np
from numpy.typing import ArrayLike, NDArray

from fem_toolbox.cfd.fluid_domain import FluidDomain
from fem_toolbox.cfd.fluid_mesh import _freeze_float
from fem_toolbox.cfd.mesh_motion import (
    MeshMotionError,
    MeshMotionResult,
    move_fluid_domain,
)
from fem_toolbox.cfd.quad_flow import QuadFlowResult, solve_quad_flow
from fem_toolbox.cfd.rans import RansResult, solve_rans
from fem_toolbox.coupling.interface import (
    InterfaceLoads,
    InterfaceMap,
    build_domain_interface,
)
from fem_toolbox.coupling.structure import _check_coordinates, _solve_structure
from fem_toolbox.errors import ValidationError
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.solvers.elasticity import ElasticResult

__all__ = [
    "FAILURE_KINDS",
    "FLUID_MODELS",
    "FSI_SCHEMA_VERSION",
    "FsiCouplingError",
    "FsiIteration",
    "StationaryFsiProblem",
    "StationaryFsiResult",
    "solve_stationary_fsi",
]

FSI_SCHEMA_VERSION = 1
FLUID_MODELS = ("laminar", "sa-noft2")
FAILURE_KINDS = (
    "input",
    "interface",
    "mesh",
    "fluid",
    "structure",
    "divergence",
    "budget",
)
# Scales, solver choice and the SA start field are coupling inputs, never options.
_RESERVED_OPTIONS = frozenset(
    {"velocity_scale_m_s", "length_scale_m", "initial_nu_tilde_m2_s", "model"}
)

FloatArray = NDArray[np.float64]


def _hash(document: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            document, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def _array_hash(*arrays: ArrayLike) -> str:
    digest = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        digest.update(str((value.dtype.str, value.shape)).encode())
        digest.update(value.tobytes())
    return digest.hexdigest()


def _input_error(name: str, reason: str) -> NoReturn:
    raise ValidationError(f"stationary FSI input {name}: {reason}", field=name)


@dataclass(frozen=True, slots=True)
class FsiIteration:
    """One complete coupling step on geometry d^k.

    displacement_residual_m is max nodal |d_solved-d^k| (unrelaxed); load_change_N
    is max nodal |F^k-F^(k-1)| (inf at k=0). Forces are sums over structure
    nodes; transfer errors are the raw conservation audit of the current map.
    """

    iteration: int
    displacement_residual_m: float
    relative_displacement_residual: float
    load_change_N: float
    fluid_force_N: FloatArray
    external_pressure_force_N: FloatArray
    structure_reaction_N: FloatArray
    transfer_force_error_N: FloatArray
    transfer_moment_error_N_m: FloatArray
    transfer_relative_error: float
    maximum_displacement_ratio: float
    maximum_displacement_gradient: float
    minimum_area_ratio: float
    fluid_iterations: int


class FsiCouplingError(ValidationError):
    """Controlled stop with kind, complete history and last iterate d^k [nnode,3]."""

    def __init__(
        self,
        kind: str,
        reason: str,
        history: list[FsiIteration],
        displacement_m: FloatArray,
    ) -> None:
        self.kind = kind
        self.history = tuple(history)
        self.displacement_m = _freeze_float(displacement_m)
        last = history[-1].relative_displacement_residual if history else None
        super().__init__(
            f"stationary FSI stopped ({kind}): {reason}; completed iterations={len(history)}, last relative displacement residual={last}",
            field="fsi_coupling",
        )


@dataclass(frozen=True, slots=True)
class StationaryFsiProblem:
    """Explicit inputs of one stationary bidirectional coupling.

    fluid_domain is the unchanged reference domain; wall_name must be an
    fsi=true wall whose segments are covered by structure facet_set. The
    structure is a validated static elasticity Model on the reference mesh;
    its Dirichlet data, body forces and tractions are kept. Pressure loads use
    p - external_pressure_Pa. fluid_options pass unchanged to solve_quad_flow
    (laminar) or solve_rans (sa-noft2), except the reserved scale/model keys.
    Default relaxation 0.2 and budget 200 follow Doc/71; no adaptive omega.
    """

    fluid_domain: FluidDomain
    wall_name: str
    structure: Model
    facet_set: str
    length_scale_m: float
    velocity_scale_m_s: float
    fluid_model: str = "laminar"
    initial_nu_tilde_m2_s: ArrayLike | None = None
    fluid_options: Mapping[str, Any] = field(default_factory=dict)
    external_pressure_Pa: float = 0.0
    relaxation: float = 0.2
    max_coupling_iterations: int = 200
    divergence_window: int = 10
    span_origin_m: float = 0.0
    preserve_inlet_flow: bool = True


@dataclass(frozen=True, slots=True)
class StationaryFsiResult:
    """Accepted stationary state with the full iteration table.

    geometry_displacement_m is d^k used for the final fluid geometry; the
    structure_displacement_m solves the final loads, max difference <=1e-8 Lref.
    Fluid, loads and structure were all freshly computed on that final
    geometry; no stale loads are reported. fsi_schema_version is in settings.
    """

    fluid: QuadFlowResult | RansResult
    mesh_motion: MeshMotionResult
    interface: InterfaceMap
    loads: InterfaceLoads
    structure: ElasticResult
    geometry_displacement_m: FloatArray
    structure_displacement_m: FloatArray
    iterations: int
    history: tuple[FsiIteration, ...]
    converged: bool
    settings: Mapping[str, Any]
    input_sha256: str


def _diverging(residuals: list[float], window: int) -> bool:
    """True after `window` consecutive strict residual increases."""
    if len(residuals) <= window:
        return False
    tail = residuals[-window - 1 :]
    return all(b > a for a, b in pairwise(tail))


def _validated(problem: StationaryFsiProblem) -> dict[str, Any]:
    # Runtime guards for untyped callers.
    if not isinstance(cast(object, problem), StationaryFsiProblem):
        _input_error("problem", "StationaryFsiProblem required")
    if not isinstance(cast(object, problem.fluid_domain), FluidDomain):
        _input_error("fluid_domain", "validated FluidDomain required")
    if not isinstance(cast(object, problem.structure), Model):
        _input_error("structure", "validated elasticity Model required")
    bc = problem.fluid_domain.boundaries.get(problem.wall_name)
    if bc is None or bc["type"] != "wall" or not bc.get("fsi", False):
        _input_error("wall_name", "named explicitly FSI-marked wall required")
    fsi_walls = sorted(
        name
        for name, value in problem.fluid_domain.boundaries.items()
        if value["type"] == "wall" and value.get("fsi", False)
    )
    if fsi_walls != [problem.wall_name]:
        _input_error(
            "wall_name", f"exactly one coupled FSI wall supported: {fsi_walls}"
        )
    scalars: dict[str, float] = {}
    for name in (
        "length_scale_m",
        "velocity_scale_m_s",
        "external_pressure_Pa",
        "relaxation",
        "span_origin_m",
    ):
        value = getattr(problem, name)
        if isinstance(value, bool) or not isinstance(value, int | float):
            _input_error(name, "finite real SI scalar required")
        if not np.isfinite(value):
            _input_error(name, "finite real SI scalar required")
        scalars[name] = float(value)
    if min(scalars["length_scale_m"], scalars["velocity_scale_m_s"]) <= 0:
        _input_error("scales", "positive SI length and velocity scales required")
    if not 0 < scalars["relaxation"] <= 1:
        _input_error("relaxation", "constant relaxation in (0,1] required")
    for name in ("max_coupling_iterations", "divergence_window"):
        value = getattr(problem, name)
        if type(value) is not int or value < 1:
            _input_error(name, "positive integer budget required")
    if type(cast(object, problem.preserve_inlet_flow)) is not bool:
        _input_error("preserve_inlet_flow", "explicit boolean required")
    if problem.fluid_model not in FLUID_MODELS:
        _input_error("fluid_model", f"one of {FLUID_MODELS} required")
    options = dict(problem.fluid_options)
    reserved = sorted(_RESERVED_OPTIONS & options.keys())
    if reserved:
        _input_error("fluid_options", f"reserved coupling inputs {reserved}")
    nu_tilde: list[float] | None = None
    if problem.fluid_model == "sa-noft2":
        if problem.initial_nu_tilde_m2_s is None:
            _input_error("initial_nu_tilde_m2_s", "explicit SA start field required")
        nu_tilde = np.asarray(problem.initial_nu_tilde_m2_s, float).tolist()
    elif problem.initial_nu_tilde_m2_s is not None:
        _input_error("initial_nu_tilde_m2_s", "only valid for sa-noft2")
    try:
        option_document = json.loads(
            json.dumps(options, allow_nan=False, default=_json)
        )
    except (TypeError, ValueError) as error:
        _input_error("fluid_options", f"JSON-representable values required: {error}")
    return {
        "fsi_schema_version": FSI_SCHEMA_VERSION,
        "wall_name": problem.wall_name,
        "facet_set": problem.facet_set,
        "fluid_model": problem.fluid_model,
        "fluid_options": option_document,
        "initial_nu_tilde_m2_s": nu_tilde,
        "max_coupling_iterations": problem.max_coupling_iterations,
        "divergence_window": problem.divergence_window,
        "preserve_inlet_flow": problem.preserve_inlet_flow,
        "displacement_gate_relative": 1e-8,
        "load_gate_relative": 1e-6,
        **scalars,
    }


def _json(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"unsupported option value {type(value).__name__}")


def _structure_hash(model: Model, facet_set: str) -> str:
    mesh = model.elasticity_inputs()["mesh"]
    rows = mesh.facet_sets.get(facet_set)
    if rows is None:
        _input_error("facet_set", "named structural Triangle3 facet set required")
    return _hash(
        {
            "document": model.to_document(),
            "mesh": _array_hash(mesh.points, mesh.tetra4, mesh.triangle3, rows),
        }
    )


def solve_stationary_fsi(problem: StationaryFsiProblem) -> StationaryFsiResult:
    """Run the Doc/71 partitioned stationary coupling with constant relaxation.

    Stops successfully only when max|r_d| <= 1e-8 Lref AND the nodal load
    change max|F^k-F^(k-1)| <= 1e-6 rho U_ref^2 (b Lref), with every fluid,
    mesh, transfer and structure gate passed on the same final geometry. A
    small relaxed update alone is never success. Failures raise
    FsiCouplingError with kind in FAILURE_KINDS and the full history.
    """
    settings = _validated(problem)
    domain = problem.fluid_domain
    lref = settings["length_scale_m"]
    uref = settings["velocity_scale_m_s"]
    omega = settings["relaxation"]
    model = problem.structure
    reference_mesh: Mesh = model.elasticity_inputs()["mesh"]
    structure_sha = _structure_hash(model, problem.facet_set)
    area_ref = domain.width_m * lref
    displacement_gate = 1e-8 * lref
    load_gate = 1e-6 * domain.fluid.density_kg_m3 * uref**2 * area_ref
    history: list[FsiIteration] = []
    d = np.zeros_like(reference_mesh.points)

    def fail(kind: str, reason: str) -> FsiCouplingError:
        return FsiCouplingError(kind, reason, history, d)

    try:
        reference_map = build_domain_interface(
            domain,
            wall_name=problem.wall_name,
            structure=reference_mesh,
            facet_set=problem.facet_set,
            length_scale_m=lref,
            span_origin_m=settings["span_origin_m"],
        )
        _check_coordinates(reference_map, model)
    except ValidationError as error:
        raise fail("interface", str(error)) from error
    rows = domain.mesh.boundary_sets[problem.wall_name]
    wall_lookup = {
        int(row): i for i, row in enumerate(domain.wall_quality.boundary_indices)
    }
    wall_rows = np.array([wall_lookup[int(row)] for row in rows], dtype=np.int64)
    options = dict(problem.fluid_options)
    previous_forces: FloatArray | None = None
    input_sha: str | None = None

    for k in range(settings["max_coupling_iterations"]):
        try:
            trace = reference_map.segment_endpoint_displacements(d)
        except ValidationError as error:
            raise fail("interface", str(error)) from error
        try:
            motion = move_fluid_domain(
                domain,
                {problem.wall_name: trace},
                length_scale_m=lref,
                preserve_inlet_flow=settings["preserve_inlet_flow"],
            )
        except MeshMotionError as error:
            raise fail("mesh", str(error)) from error
        if input_sha is None:
            input_sha = _hash(
                {
                    "settings": settings,
                    "reference_domain_sha256": motion.reference_domain_sha256,
                    "structure_sha256": structure_sha,
                }
            )
        try:
            fluid: QuadFlowResult | RansResult
            if problem.fluid_model == "laminar":
                fluid = solve_quad_flow(
                    motion.domain,
                    velocity_scale_m_s=uref,
                    length_scale_m=lref,
                    **options,
                )
                flow = fluid
            else:
                fluid = solve_rans(
                    motion.domain,
                    initial_nu_tilde_m2_s=problem.initial_nu_tilde_m2_s,  # type: ignore[arg-type]
                    velocity_scale_m_s=uref,
                    length_scale_m=lref,
                    **options,
                )
                flow = fluid.flow
        except ValidationError as error:
            raise fail("fluid", str(error)) from error
        traction = np.asarray(flow.wall_traction_Pa)[wall_rows]
        current_mesh = Mesh(
            points=reference_mesh.points + d,
            tetra4=reference_mesh.tetra4,
            triangle3=reference_mesh.triangle3,
            facet_sets=reference_mesh.facet_sets,
        )
        try:
            current_map = build_domain_interface(
                motion.domain,
                wall_name=problem.wall_name,
                structure=current_mesh,
                facet_set=problem.facet_set,
                length_scale_m=lref,
                span_origin_m=settings["span_origin_m"],
            )
            loads = current_map.loads_from_fluid_traction(
                traction, external_pressure_Pa=settings["external_pressure_Pa"]
            )
        except ValidationError as error:
            raise fail("interface", str(error)) from error
        forces = np.asarray(loads.nodal_forces_N)
        try:
            structure, solved, ratio, gradient = _solve_structure(model, forces, lref)
        except ValidationError as error:
            raise fail("structure", str(error)) from error
        residual = solved - d
        r_max = float(np.max(np.linalg.norm(residual, axis=1)))
        load_change = (
            float("inf")
            if previous_forces is None
            else float(np.max(np.linalg.norm(forces - previous_forces, axis=1)))
        )
        history.append(
            FsiIteration(
                iteration=k,
                displacement_residual_m=r_max,
                relative_displacement_residual=r_max / lref,
                load_change_N=load_change,
                fluid_force_N=_freeze_float(np.asarray(loads.fluid_force_N)),
                external_pressure_force_N=_freeze_float(
                    np.asarray(loads.external_pressure_force_N)
                ),
                structure_reaction_N=_freeze_float(structure.reactions.sum(axis=0)),
                transfer_force_error_N=loads.audit.force_error_N,
                transfer_moment_error_N_m=loads.audit.moment_error_N_m,
                transfer_relative_error=loads.audit.maximum_relative_error,
                maximum_displacement_ratio=ratio,
                maximum_displacement_gradient=gradient,
                minimum_area_ratio=motion.minimum_area_ratio,
                fluid_iterations=fluid.iterations,
            )
        )
        if not np.all(np.isfinite(residual)) or not np.isfinite(r_max):
            raise fail("divergence", "nonfinite displacement residual")
        if r_max <= displacement_gate and load_change <= load_gate:
            assert input_sha is not None
            return StationaryFsiResult(
                fluid=fluid,
                mesh_motion=motion,
                interface=current_map,
                loads=loads,
                structure=structure,
                geometry_displacement_m=_freeze_float(d),
                structure_displacement_m=_freeze_float(solved),
                iterations=k + 1,
                history=tuple(history),
                converged=True,
                settings=MappingProxyType(
                    {
                        **settings,
                        "displacement_gate_m": displacement_gate,
                        "load_gate_N": load_gate,
                        "reference_area_m2": area_ref,
                    }
                ),
                input_sha256=input_sha,
            )
        if _diverging(
            [row.displacement_residual_m for row in history],
            settings["divergence_window"],
        ):
            raise fail(
                "divergence",
                f"unrelaxed residual increased in {settings['divergence_window']} consecutive iterations",
            )
        previous_forces = forces
        d = d + omega * residual
    raise fail(
        "budget",
        f"coupling budget {settings['max_coupling_iterations']} exhausted before both gates",
    )
