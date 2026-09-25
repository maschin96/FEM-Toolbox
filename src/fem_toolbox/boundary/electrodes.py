"""Conducting equipotential boundaries for quasistatic electric fields."""

from dataclasses import dataclass

from numpy.typing import ArrayLike


@dataclass(frozen=True)
class Electrode:
    """Named exterior Triangle3/Triangle6 set with fixed potential (V) or total charge (C).

    ``facets`` is a named mesh facet set or active surface row indices. Exactly one
    of ``potential`` and ``charge`` must be supplied. A floating open-circuit
    electrode has ``charge=0``; no potential reference is inferred.
    """

    name: str
    facets: str | ArrayLike
    potential: float | None = None
    charge: float | None = None


@dataclass(frozen=True)
class ElectrodeState:
    """Solved common electrode potential (V) and adjacent conductor charge (C)."""

    name: str
    potential: float
    charge: float
    prescribed_charge: float | None
