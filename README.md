# FEM Toolbox

Eine wissenschaftlich ausgerichtete Python-Toolbox für lineare, dreidimensionale
Finite-Elemente-Simulationen mit gekoppelter

- Festkörpermechanik,
- Wärmeleitung,
- Elektrostatik und
- Piezoelektrizität.

Der erste Elementtyp ist der lineare Vierknoten-Tetraeder **Tet4**. Das Projekt
beginnt bewusst mit kleinen Deformationen, linearen Materialien und
quasistatischen Feldgleichungen. Diese Einschränkung schafft eine überprüfbare
Basis, auf der transiente, dynamische und nichtlineare Modelle später aufgebaut
werden können.

> Projektstatus: Forschungs- und Implementierungsplanung. Es existiert noch kein
> freigegebener numerischer Solver.

## Wissenschaftlicher Anspruch

Die Toolbox soll nicht nur Ergebnisse erzeugen, sondern deren Vertrauenswürdigkeit
nachweisbar machen. Jede Physik wird deshalb zunächst separat verifiziert. Erst
danach werden Kopplungsterme aktiviert. Dimensionsanalyse, Energie- und
Reziprozitätsprüfungen, Patch-Tests, analytische Vergleichsfälle und
Netzkonvergenz sind Teil der Definition of Done.

## Geplanter Funktionsumfang der ersten Version

| Gebiet | Modell | Primärvariable |
|---|---|---|
| Mechanik | kleine Verzerrungen, linear isotrop/anisotrop, statisch | Verschiebung `u = (u_x,u_y,u_z)` |
| Wärme | lineare stationäre Wärmeleitung | Temperatur `T` |
| Elektrostatik | lineares Dielektrikum, freie Ladung | Potential `phi` |
| Piezoelektrizität | lineare Stress-Charge-Form | `u`, `phi` |
| Thermomechanik | thermische Dehnung, zunächst einseitig gekoppelt | `T -> u` |

Jeder Knoten besitzt im vollständig belegten Modell fünf Freiheitsgrade:
`[u_x, u_y, u_z, T, phi]`. Die Implementierung soll jedoch feldweise DOF-Tabellen
verwenden, damit Teilmodelle keine unbenutzten Freiheitsgrade tragen.

## Dokumentation

Der Einstiegspunkt ist [Doc/README.md](Doc/README.md). Besonders wichtig sind:

- [Doc/01_SCOPE_AND_REQUIREMENTS.md](Doc/01_SCOPE_AND_REQUIREMENTS.md) — Ziele,
  Nichtziele und Abnahmekriterien
- [Doc/02_MATHEMATICAL_MODEL.md](Doc/02_MATHEMATICAL_MODEL.md) — Kontinuumsmodell,
  Konventionen und schwache Form
- [Doc/03_TET4_FORMULATION.md](Doc/03_TET4_FORMULATION.md) — diskrete
  Tet4-Gleichungen
- [Doc/04_SOFTWARE_ARCHITECTURE.md](Doc/04_SOFTWARE_ARCHITECTURE.md) — Python-
  Architektur und Datenmodell
- [Doc/05_VERIFICATION_AND_VALIDATION.md](Doc/05_VERIFICATION_AND_VALIDATION.md) —
  V&V-Programm
- [Doc/06_ROADMAP.md](Doc/06_ROADMAP.md) — implementierbare Arbeitspakete
- [Doc/07_RISKS_AND_DECISIONS.md](Doc/07_RISKS_AND_DECISIONS.md) — Risiken und
  Architekturentscheidungen

## Geplante Nutzung

Die angestrebte öffentliche API soll deklarativ und reproduzierbar sein:

```python
from fem_toolbox import Model

model = Model.from_file("cantilever.yaml")
result = model.solve()
result.write("cantilever.vtu")
```

Die Eingabedatei beschreibt Netz, Materialien, Gebiete, Randbedingungen,
Lasten, Analyseart und Solveroptionen. Intern bleiben Netz- und Physikobjekte
Python-APIs, damit wissenschaftliche Experimente ohne Dateiformat möglich sind.

## Vorgesehener Technologie-Stack

- Python 3.11 oder neuer
- NumPy für Elementalgebra und Feldspeicher
- SciPy für Sparse-Matrizen und Referenzsolver
- meshio für Gmsh-/VTK-Ein- und Ausgabe
- pytest für Verifikationstests
- optional später PETSc/petsc4py für große gekoppelte Systeme

Die Kernformulierung soll NumPy/SciPy als einzige zwingende numerische
Abhängigkeiten haben. Meshing und Visualisierung bleiben austauschbare Adapter.

## Entwicklung

Die Projektstruktur wird in Meilenstein M0 angelegt. Danach liefert jeder
Meilenstein ein ausführbares Referenzproblem und dokumentierte
Verifikationsdaten. Eine Kopplung gilt erst als implementiert, wenn beide
Einzelfelder zuvor ihre Tests bestehen.

## Quellenbasis

Die Planung stützt sich auf etablierte Kontinuums- und FEM-Formulierungen sowie
offizielle Referenzdokumentation. Die annotierte Liste steht in
[Doc/08_REFERENCES.md](Doc/08_REFERENCES.md).
