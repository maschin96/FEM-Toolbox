"""Experimental mesh/BC/material preflight, deliberately without a CFD solver."""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any

import numpy as np

from fem_toolbox.cfd.channel import _scalar
from fem_toolbox.cfd.flow_model import IncompressibleFluid
from fem_toolbox.cfd.fluid_mesh import (
    FluidMesh2D,
    WallQuality,
    _fail,
    _float_array,
    _freeze_float,
    _keys,
)
from fem_toolbox.errors import ModelSchemaError


@dataclass(frozen=True, slots=True, init=False)
class FluidDomain:
    """Immutable SI mesh/material/BC preflight, not a RANS solve.

    Boundary indices refer to input line2 rows. Inlet velocity=(2,) or
    (nselection,2) m/s; nu_tilde scalar or (nselection,) m²/s. Periodic
    translations and pressure drops are paired with opposite signs.
    A single connected mesh needs one outlet group or one pressure gauge.
    Width is explicit for conversion of planar fluxes to kg/s.
    """

    mesh: FluidMesh2D
    fluid: IncompressibleFluid
    boundaries: Mapping[str, Mapping[str, Any]]
    wall_quality: WallQuality
    inlet_mass_flow_kg_s: float
    pressure_reference_cell: int | None
    width_m: float
    periodic_partners: Mapping[int, int]

    @classmethod
    def from_document(
        cls, document: Mapping[str, Any], *, mesh: FluidMesh2D
    ) -> FluidDomain:
        _keys(
            document,
            {"fluid_domain_schema_version", "mesh", "fluid", "boundaries", "width_m"},
            {"pressure_reference_cell"},
            "fluid_domain",
        )
        if (
            type(document["fluid_domain_schema_version"]) is not int
            or document["fluid_domain_schema_version"] != 1
        ):
            _fail("fluid_domain_schema_version", "expected integer 1")
        mesh = _validated_mesh(mesh)
        _mesh_path(document["mesh"])
        f = document["fluid"]
        _keys(f, {"density_kg_m3", "viscosity_Pa_s"}, set(), "fluid")
        fluid = IncompressibleFluid(f["density_kg_m3"], f["viscosity_Pa_s"])
        width = _scalar(document["width_m"], "width_m")
        if width <= 0:
            _fail("width_m", "positive SI width required")
        raw = document["boundaries"]
        if not isinstance(raw, Mapping) or not raw:
            _fail("boundaries", "nonempty named boundary table required")
        covered = np.zeros(len(mesh.line2), dtype=np.int64)
        boundaries: dict[str, Mapping[str, Any]] = {}
        walls = []
        outlets = 0
        inlet_flow = 0.0
        normals = mesh.boundary_normals
        for name, bc in raw.items():
            path = f"boundaries.{name}"
            if name not in mesh.boundary_sets:
                _fail(path, "unknown boundary group")
            ids = mesh.boundary_sets[name]
            covered[ids] += 1
            if not isinstance(bc, Mapping):
                _fail(path, "expected boundary table")
            kind = bc.get("type")
            data: dict[str, Any] = {"type": kind}
            if kind == "wall":
                _keys(bc, {"type"}, {"fsi"}, path)
                marker = bc.get("fsi", False)
                if type(marker) is not bool:
                    _fail(f"{path}.fsi", "expected boolean")
                data["fsi"] = marker
                walls.append(name)
            elif kind == "symmetry":
                _keys(bc, {"type"}, set(), path)
            elif kind == "pressure_outlet":
                _keys(bc, {"type", "pressure_Pa"}, set(), path)
                data["pressure_Pa"] = _scalar(bc["pressure_Pa"], f"{path}.pressure_Pa")
                outlets += 1
            elif kind == "inlet":
                _keys(bc, {"type", "velocity_m_s", "nu_tilde_m2_s"}, set(), path)
                velocity = _float_array(bc["velocity_m_s"], f"{path}.velocity_m_s")
                if velocity.shape == (2,):
                    velocity = np.tile(velocity, (len(ids), 1))
                if velocity.shape != (len(ids), 2):
                    _fail(f"{path}.velocity_m_s", "expected (2,) or (nface,2)")
                a = _float_array(bc["nu_tilde_m2_s"], f"{path}.nu_tilde_m2_s")
                if a.ndim == 0:
                    a = np.full(len(ids), float(a))
                if a.shape != (len(ids),) or np.any(a < 0):
                    _fail(
                        f"{path}.nu_tilde_m2_s",
                        "finite nonnegative scalar or (nface,) values required",
                    )
                outward = np.einsum("ij,ij->i", velocity, normals[ids])
                if np.any(outward > 0) or not np.any(outward < 0):
                    _fail(
                        f"{path}.velocity_m_s",
                        "inward normal flux required; no outward inlet face",
                    )
                lengths = mesh.face_lengths_m[mesh.boundary_face_indices[ids]]
                inlet_flow -= float(fluid.density_kg_m3 * width * (outward @ lengths))
                data.update(
                    velocity_m_s=_freeze_float(velocity), nu_tilde_m2_s=_freeze_float(a)
                )
            elif kind == "periodic":
                _keys(
                    bc,
                    {"type", "partner", "translation_m", "pressure_drop_Pa"},
                    set(),
                    path,
                )
                if not isinstance(bc["partner"], str) or bc["partner"] == name:
                    _fail(f"{path}.partner", "a different named partner required")
                translation = _float_array(bc["translation_m"], f"{path}.translation_m")
                if translation.shape != (2,) or not np.any(translation):
                    _fail(f"{path}.translation_m", "nonzero translation (2,) required")
                data.update(
                    partner=bc["partner"],
                    translation_m=_freeze_float(translation),
                    pressure_drop_Pa=_scalar(
                        bc["pressure_drop_Pa"], f"{path}.pressure_drop_Pa"
                    ),
                )
            else:
                _fail(
                    f"{path}.type",
                    "expected wall/inlet/pressure_outlet/symmetry/periodic",
                )
            boundaries[name] = MappingProxyType(data)
        if np.any(covered != 1):
            i = int(np.flatnonzero(covered != 1)[0])
            _fail(
                f"line2[{i}]", f"exactly one physical BC required, found {covered[i]}"
            )
        if not walls:
            _fail("boundaries", "wall assignment required for wall-resolved SA")
        gauge = document.get("pressure_reference_cell")
        if outlets:
            if outlets != 1 or gauge is not None:
                _fail(
                    "pressure_reference_cell",
                    "one pressure outlet group and no additional gauge required",
                )
        elif type(gauge) is not int or not 0 <= gauge < len(mesh.quad4):
            _fail(
                "pressure_reference_cell",
                "one valid cell gauge required without pressure outlet",
            )
        if inlet_flow > 0 and not outlets:
            _fail("boundaries", "nonzero inlet mass flux requires a pressure outlet")
        if not np.isfinite(inlet_flow):
            _fail("boundaries", "inlet mass flux exceeds float64 range")
        partners = _periodic_pairs(mesh, boundaries)
        obj = object.__new__(cls)
        for key, value in dict(
            mesh=mesh,
            fluid=fluid,
            boundaries=MappingProxyType(boundaries),
            wall_quality=mesh.wall_quality(walls),
            inlet_mass_flow_kg_s=inlet_flow,
            pressure_reference_cell=gauge,
            width_m=width,
            periodic_partners=MappingProxyType(partners),
        ).items():
            object.__setattr__(obj, key, value)
        return obj

    @classmethod
    def from_file(cls, path: str | Path) -> FluidDomain:
        """Read schema-1 JSON domain and mesh confined to the same input tree."""
        source = Path(path).resolve()
        try:
            d = json.loads(source.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise ModelSchemaError(path=str(source), reason=str(error)) from error
        if not isinstance(d, dict) or "mesh" not in d:
            _fail("mesh", "missing mesh reference")
        relative = _mesh_path(d["mesh"])
        target = (source.parent / relative).resolve()
        if not target.is_relative_to(source.parent):
            _fail("mesh", "mesh symlink escapes domain input directory")
        return cls.from_document(d, mesh=FluidMesh2D.from_file(target))


def _mesh_path(value: Any) -> Path:
    if not isinstance(value, str) or not value or "\x00" in value:
        _fail("mesh", "nonempty relative JSON mesh path required")
    p = Path(value)
    if p.is_absolute() or ".." in p.parts:
        _fail("mesh", "mesh path must stay within input directory")
    return p


def _periodic_pairs(
    mesh: FluidMesh2D, boundaries: Mapping[str, Mapping[str, Any]]
) -> dict[int, int]:
    pairs = {}
    scale = float(np.ptp(mesh.points, axis=0).max())
    tolerance = 1e-10 * scale
    for name, bc in boundaries.items():
        if bc["type"] != "periodic":
            continue
        partner = boundaries.get(bc["partner"])
        path = f"boundaries.{name}.partner"
        if (
            partner is None
            or partner["type"] != "periodic"
            or partner["partner"] != name
        ):
            _fail(path, "reciprocal periodic group required")
        if (
            not np.allclose(
                bc["translation_m"], -partner["translation_m"], atol=tolerance, rtol=0
            )
            or bc["pressure_drop_Pa"] != -partner["pressure_drop_Pa"]
        ):
            _fail(path, "opposite translations and pressure drops required")
        own, other = mesh.boundary_sets[name], mesh.boundary_sets[bc["partner"]]
        if len(own) != len(other):
            _fail(path, "congruent periodic face counts required")
        matched = set()
        for i in own:
            target = mesh.points[mesh.line2[i]] + bc["translation_m"]
            candidates = []
            for j in other:
                points = mesh.points[mesh.line2[j]]
                if np.allclose(target, points, atol=tolerance, rtol=0) or np.allclose(
                    target, points[::-1], atol=tolerance, rtol=0
                ):
                    candidates.append(int(j))
            if len(candidates) != 1 or candidates[0] in matched:
                _fail(path, f"face {i}: unique translated match required")
            j = candidates[0]
            if not np.allclose(
                mesh.boundary_normals[i], -mesh.boundary_normals[j], atol=1e-10, rtol=0
            ):
                _fail(path, "opposing outward periodic normals required")
            matched.add(j)
            pairs[int(i)] = j
    return pairs


def _validated_mesh(value: object) -> FluidMesh2D:
    if not isinstance(value, FluidMesh2D):
        _fail("mesh", "validated FluidMesh2D required")
    return value
