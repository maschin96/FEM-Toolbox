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

> Projektstatus: `v0.1.0` stellt das installierbare Projektfundament bereit. Auf
> `main` entsteht im Entwicklungszyklus `0.2.0.dev0` die Netz- und
> Tet4-Geometriebasis. Es existiert noch kein freigegebener numerischer Solver.

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
- [Doc/09_RELEASE_PROCESS.md](Doc/09_RELEASE_PROCESS.md) — Versionszyklus,
  Release-Gates und Tagging-Regeln

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

## Lokale Installation

Voraussetzung ist Python 3.11 oder neuer. Der aktuelle Entwicklungsstand wird
aus dem Repository inklusive Smoke-Test installiert:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m mypy src/fem_toolbox
```

Das Paket stellt die installierte Version, mathematische und SI-Konventionen,
spezifische Validierungsfehler sowie das deterministische Provenienzmodell
`RunProvenance` bereit. Der M1-Entwicklungsstand ergänzt ein validiertes
Mesh-Datenmodell. Numerische Solver werden in den folgenden Milestones
testgetrieben ergänzt.

## Mesh-Datenmodell

`Mesh` kopiert seine Eingaben in unveränderliche NumPy-Arrays. Konnektivitäten
und Auswahlsätze verwenden ausschließlich nullbasierte interne Indizes; IDs aus
externen Netzformaten bleiben getrennte Metadaten.

```python
import numpy as np

from fem_toolbox.mesh import Mesh

mesh = Mesh(
    points=np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    ),
    tetra4=np.array([[0, 1, 2, 3]]),
    tetra4_ids=np.array([101]),
    cell_tags=np.array([7]),
    volume_sets={"ceramic": np.array([0]), "active": np.array([0])},
)
```

Verschiedene benannte Sätze dürfen sich überschneiden. Doppelte Indizes im
selben Satz, ungültige Formen, nichtganzzahlige Konnektivitäten, Indexüberläufe
und topologisch doppelte Entitäten werden beim Erzeugen des Objekts abgelehnt.

## Gmsh-Import und VTU-Export

Die `meshio`-Adapter lesen lineare Tetraeder und Dreiecksflächen aus Gmsh 4.1
oder 2.2. Physische Namen werden getrennt als Volumen- und Flächensätze
übernommen. Das erhält auch mehrere physische Gruppen pro Gmsh-Entität. Der
einzelne numerische Kompatibilitätstag entspricht dagegen dem ersten von Gmsh
gemeldeten Tag. Andere Zelltypen führen zu einem `UnsupportedCellTypeError`.

Gmsh speichert keine verbindliche Längeneinheit. Der Faktor
`length_scale_to_m` macht die Umrechnung in die intern vorausgesetzten Meter
explizit; ohne Angabe werden Eingabekoordinaten als Meter interpretiert.

```python
import numpy as np

from fem_toolbox.io import read_gmsh, write_vtu

mesh = read_gmsh("specimen.msh", length_scale_to_m=1.0e-3)  # Millimeter -> Meter
write_vtu(
    "result.vtu",
    mesh,
    point_data={"temperature_K": np.full(mesh.n_points, 293.15)},
)
```

`write_vtu` schreibt Punktdaten und Zellfelder. Ein Zellfeld muss für jeden im
Netz vorhandenen Zelltyp ein Array unter `tetra4` beziehungsweise `triangle3`
enthalten. Benannte Auswahlsätze werden als voneinander unabhängige
`fem_toolbox:set:<name>`-Masken exportiert und bleiben daher auch bei
Überlappung erhalten.

## Tet4-Geometrie

`Tet4Geometry` berechnet das signierte und positive Volumen, konstante
Formfunktionsgradienten, eine positive lokale Orientierung, alle äußeren
Flächennormalen sowie skalierungsunabhängige Qualitätsindikatoren:

```python
from fem_toolbox.elements import Tet4Geometry

geometry = Tet4Geometry(mesh.points[mesh.tetra4[0]], entity_id=101)
print(geometry.volume, geometry.quality.scaled_volume)
```

Entartete oder extrem abgeflachte Tetraeder werden vor einem linearen Lösevorgang
mit ihrer Entitäts-ID abgelehnt. Die Standardgrenze und ihre mathematische
Definition stehen in [Doc/03_TET4_FORMULATION.md](Doc/03_TET4_FORMULATION.md).

GitHub Actions führt Tests auf Python 3.11, 3.12 und 3.13 aus. Linting,
Formatprüfung und Typprüfung laufen auf Python 3.13. Der Workflow besitzt nur
Leserechte auf Repository-Inhalte; verwendete Actions sind auf vollständige
Commit-SHAs fixiert.

## Entwicklung

Die Projektstruktur und Qualitätsinfrastruktur wurden in Meilenstein M0
angelegt. Die folgenden Meilensteine liefern jeweils ein ausführbares
Referenzproblem und dokumentierte Verifikationsdaten. Eine Kopplung gilt erst
als implementiert, wenn beide Einzelfelder zuvor ihre Tests bestehen.

Die Umsetzung wird in [GitHub-Milestones und Issues](https://github.com/maschin96/FEM-Toolbox/milestones)
verfolgt. Jeder abgeschlossene Milestone erzeugt genau eine freigegebene Version
und einen annotierten Git-Tag. Der Planungsstand ist als `v0.0.1` markiert, das
Projektfundament als `v0.1.0`; der vollständige stationäre Workflow ist für
`v1.0.0` vorgesehen.

## Quellenbasis

Die Planung stützt sich auf etablierte Kontinuums- und FEM-Formulierungen sowie
offizielle Referenzdokumentation. Die annotierte Liste steht in
[Doc/08_REFERENCES.md](Doc/08_REFERENCES.md).
