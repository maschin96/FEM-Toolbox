"""Mesh-bound dynamic model result; keeps material/field export associations."""

from dataclasses import dataclass

from femantis.materials import PiezoelectricMaterial
from femantis.mesh import Mesh
from femantis.provenance import RunProvenance
from femantis.solvers.piezo_dynamics import PiezoelectricModes
from femantis.solvers.piezo_harmonic import HarmonicPiezoelectricResult
from femantis.solvers.piezo_transient import TransientPiezoelectricResult


@dataclass(frozen=True, eq=False)
class PiezoDynamicResult:
    """Dynamic solution bound to the exact mesh and effective cell materials.

    solution is a modal result, transient history, or ordered tuple of harmonic
    states. Provenance includes the resolved schema-3 document and solver hashes.
    """

    mesh: Mesh
    materials: tuple[PiezoelectricMaterial, ...]
    solution: (
        PiezoelectricModes
        | TransientPiezoelectricResult
        | tuple[HarmonicPiezoelectricResult, ...]
    )
    provenance: RunProvenance
