# FEM Toolbox

Eine wissenschaftlich ausgerichtete Python-Toolbox im Aufbau für lineare,
dreidimensionale Finite-Elemente-Simulationen. Ziel sind verifizierte Modelle für
Festkörpermechanik, Wärmeleitung, Elektrostatik und deren thermomechanische und
piezoelektrische Kopplungen auf Basis linearer Vierknoten-Tetraeder (**Tet4**).

**Aktueller Stand: v1.0.0 — vollständiger stationärer Workflow.**
Verfügbar sind das Netzmodell, Tet4-Geometrie, drei stationäre Einzelfelder,
die einseitige Kopplung `T → u`, der monolithische piezoelektrische
`u–φ`-Solver mit festen und schwebenden Elektroden sowie `T → (u, φ)`.
Die versionierte TOML-`Model`-API ergänzt Materialregionen, Preflight und
Ergebnisbündel mit Provenienz. Umfang und Migration stehen in den
[Release Notes](Doc/releases/v1.0.0.md). Die
[Gültigkeitsübersicht](Doc/26_VALIDITY_AND_LIMITS.md) beschreibt die
verifizierten Annahmen und numerischen Grenzen.

[Installation](#installation) · [Schnellstart](#schnellstart) ·
[Dokumentation](#dokumentation) · [Entwicklung](#entwicklung) ·
[Roadmap](#roadmap)

## Verfügbare Funktionen

| Bereich | Umfang in v1.0.0 |
|---|---|
| Netzmodell | Unveränderliche NumPy-Arrays für 3D-Punkte, Tet4-Zellen und Triangle3-Flächen; validierte Indizes, externe IDs und Tags |
| Auswahlsätze | Benannte Volumen- und Flächensätze, auch mit überlappenden Mitgliedschaften |
| Tet4-Geometrie | Volumen, Formfunktionen, Gradienten, positive lokale Orientierung, Flächen und äußere Normalen |
| Netzqualität | Skalierungsunabhängige Qualitätsindikatoren und Erkennung entarteter Tetraeder |
| Ein- und Ausgabe | Gmsh-2.2/4.1-Import und VTU-Export über `meshio`, einschließlich Feldern und Auswahlsätzen |
| Projektbasis | Mathematische und SI-Konventionen, spezifische Validierungsfehler und deterministische Laufprovenienz mit `RunProvenance` |
| Elastizität | Kleine lineare Tet4-Verformungen, anisotrope Materialien, Volumen-/Flächen-/Punktlasten sowie Spannungs-, Energie- und Reaktionsauswertung |
| Wärme | Stationäre Tet4-Wärmeleitung mit anisotroper Leitfähigkeit, Quellen, äußerem Wärmestrom, Konvektion und Leistungsbilanz |
| Elektrostatik | Stationäres Tet4-Potential mit anisotroper Permittivität, freien Ladungen, Gauge-Diagnostik, E-/D-Feldern und Ladungsbilanz |
| Thermomechanik | Sequenzielle stationäre Wärme-zu-Mechanik-Kopplung mit anisotroper thermischer Expansion und regionsweisen Referenztemperaturen |
| Piezoelektrizität | Stress-Charge-Materialdaten, Konversion und Achsrotation; gekoppelte Tet4-Operatoren und monolithische Sparse-LU-Lösung mit feldweiser Skalierung, Reaktionen, Ladungs- und Energiebilanzen |
| Thermopiezoelektrizität | Einseitiger stationärer Ablauf `T → (u, φ)` mit regionsweiser thermischer Eigenverzerrung und getrennten Wärme-, Kraft- und Ladungsbilanzen in v1.0.0 |

Ein vollständiger deklarativer Lauf mit Gmsh-Eingabe, VTU-Feldern und
Provenienzbericht steht unter [examples/thermopiezoelectric](examples/thermopiezoelectric/README.md).

## Installation

Die vollständige [Installationsanleitung](INSTALL.md) beschreibt die geprüften
Plattformen, CPython 3.11–3.13, virtuelle Umgebungen, Projektabhängigkeiten,
die optionale Gmsh-Installation, einen ersten Solverlauf und Fehlerbehebung.
Zum Einlesen der mitgelieferten `.msh`-Dateien ist Gmsh nicht erforderlich.

```sh
git clone https://github.com/maschin96/FEM-Toolbox.git
cd FEM-Toolbox
python -m venv .venv
```

Virtuelle Umgebung aktivieren — unter macOS/Linux:

```sh
source .venv/bin/activate
```

Unter Windows mit PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Anschließend das Paket samt Entwicklungswerkzeugen installieren:

```sh
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
python -m pytest
```

Die editierbare Installation verwendet den aktuellen Repository-Stand. Die
[Installationsanleitung](INSTALL.md) enthält zusätzlich eine Prüfung in einer
sauberen Umgebung und den Gmsh-Aufruf für ein eigenes Testnetz.

## Schnellstart

Dieses vollständige Beispiel erzeugt einen Tetraeder mit Kanten entlang der
Koordinatenachsen von jeweils einem Meter, berechnet sein Volumen und exportiert
das Netz. Nach der Installation lässt es sich als Python-Skript ausführen:

```python
import numpy as np

from fem_toolbox.elements import Tet4Geometry
from fem_toolbox.io import write_vtu
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
    volume_sets={"specimen": np.array([0])},
)

geometry = Tet4Geometry(mesh.points[mesh.tetra4[0]], entity_id=101)
print(f"Volumen: {geometry.volume:.6f} m³")
write_vtu("tetrahedron.vtu", mesh)
```

Erwartete Ausgabe: `Volumen: 0.166667 m³`. Die Datei `tetrahedron.vtu` wird im
aktuellen Arbeitsverzeichnis erstellt und kann beispielsweise in ParaView
geöffnet werden.

### Gmsh-Netz importieren

Aus dem Repository-Root kann das mitgelieferte Referenznetz direkt eingelesen
und mit einem beispielhaften Punktfeld exportiert werden:

```python
import numpy as np

from fem_toolbox.io import read_gmsh, write_vtu

mesh = read_gmsh("tests/fixtures/overlapping_groups.msh", length_scale_to_m=1.0)
write_vtu(
    "mesh_with_temperature.vtu",
    mesh,
    point_data={"temperature_K": np.full(mesh.n_points, 293.15)},
)
```

Das Temperaturfeld ist vorgegeben und dient ausschließlich als Exportbeispiel.
Gmsh legt keine verbindliche Längeneinheit fest: `length_scale_to_m=1.0`
interpretiert Koordinaten als Meter, `1.0e-3` rechnet Millimeter in Meter um.

### Konventionen und Grenzen

- Interne Konnektivitäten und Auswahlsätze verwenden **nullbasierte Indizes**.
  Externe Entitäts-IDs bleiben getrennte Metadaten.
- `Mesh` kopiert seine Eingaben und prüft unter anderem Array-Formen,
  Indexgrenzen und topologisch doppelte Entitäten. Geometrische Entartung wird
  durch `Tet4Geometry` und beim Gmsh-Import geprüft.
- Der Import unterstützt ausschließlich lineare Tetraeder und Dreiecksflächen;
  andere Zelltypen führen zu einem `UnsupportedCellTypeError`.
- Mehrfache physische Gmsh-Gruppen bleiben als benannte Auswahlsätze erhalten.
  Der einzelne numerische Tag bildet nur die erste gemeldete Mitgliedschaft ab.
- VTU-Zellfelder benötigen für jeden vorhandenen Zelltyp ein Array unter
  `tetra4` beziehungsweise `triangle3`. Auswahlsätze werden als separate
  `fem_toolbox:set:<name>`-Masken exportiert.
- Sehr flache oder entartete Tetraeder werden abgelehnt. Definition und
  Standardgrenze der Qualitätsprüfung stehen in der
  [Tet4-Formulierung](Doc/03_TET4_FORMULATION.md).

## Wissenschaftlicher Anspruch

Jedes physikalische Einzelfeld soll vor seiner Kopplung separat verifiziert
werden. Dazu gehören analytische Vergleichsfälle, Patch-Tests, Energie- und
Reziprozitätsprüfungen, globale Bilanzen und Netzkonvergenz.

Für v1.0.0 umfasst der Nachweis Patch-, Bilanz-, analytische Benchmark- und
vierstufige Konvergenztests für Mechanik, Wärme, Elektrostatik,
Thermomechanik und Piezoelektrizität. Direkter und inverser Piezoeffekt
werden zusätzlich gegen eine unabhängige 1D-FE-Implementierung geprüft.
**Verifikation** prüft die korrekte Umsetzung des mathematischen Modells;
**Validierung** prüft dessen Eignung gegenüber geeigneten Referenz- oder
experimentellen Daten und ist für diese Freigabe noch nicht beansprucht. Das
[V&V-Konzept](Doc/05_VERIFICATION_AND_VALIDATION.md) legt die Anforderungen fest.

## Roadmap

Der erste stationäre Workflow beschränkt sich auf kleine Deformationen und
lineare Materialien. Vorgesehen sind Verschiebung `u`, Temperatur `T` und
elektrisches Potential `phi` mit feldweiser Freiheitsgradverwaltung.

| Meilenstein | Zielversion | Inhalt | Stand |
|---|---|---|---|
| M0 | v0.1.0 | Paket, Konventionen, Provenienz und Qualitätsinfrastruktur | Freigegeben |
| M1 | v0.2.0 | Netz und Tet4-Geometrie | Freigegeben |
| M2 | v0.3.0 | Stationäre Elastizität, Wärmeleitung und Elektrostatik | Freigegeben |
| M3 | v0.4.0 | Sequenzielle Thermomechanik `T → u` | Freigegeben |
| M4 | v0.5.0 | Monolithische Piezoelektrizität `u–phi` | Freigegeben |
| M5 | v1.0.0 | Elektroden, Modellformat und vollständiger stationärer Workflow | Freigegeben |
| M6 | v1.1.0 | Transiente Wärmeleitung und lineare Dynamik | Geplant |
| M7 | v1.2.0 | PETSc-Backend, Tet10 und Skalierung | Geplant |
| M8 | v1.3.0* | Reproduzierbare Parameterstudien | Geplant |
| M9 | v1.4.0* | Experimentelle Validierung und Unsicherheit | Geplant |
| M10 | v1.5.0* | Kontrollierte adaptive Netzstudien | Geplant |
| M11 | Offen | Topologieoptimierung | Grober Horizont |
| M12 | Offen | RANS-CFD mit Fluid-Struktur-Kopplung | Grober Horizont |

*Vorläufige Zielversionen ohne Terminzusage. Produktvision, Bestandsanalyse,
Abhängigkeiten und Abnahmegates stehen in der
[zweiten Planungswelle](Doc/13_PRODUCT_VISION.md).

Die deklarative API `Model.from_file(...).solve()` ist für Wärme, Elastizität,
Elektrostatik, Piezoelektrizität und den Ablauf `T → (u, φ)` verfügbar. Das [TOML-Schema](Doc/21_MODEL_SCHEMA.md)
enthält ein ausführbares Beispiel und die Format-Migrationsregeln. Die
[API-Referenz](Doc/29_PUBLIC_API_CONTRACT.md) beschreibt öffentliche Python-
Einstiegspunkte und den Kompatibilitätsvertrag für v1.0. Weitere
Arbeitspakete und Abnahmekriterien stehen in der
[Roadmap](Doc/06_ROADMAP.md); die Umsetzung wird über
[GitHub-Milestones](https://github.com/maschin96/FEM-Toolbox/milestones) verfolgt.

## Dokumentation

Der [Dokumentationsindex](Doc/README.md) führt durch die fachlichen und
technischen Entwurfsdokumente.

| Thema | Dokument |
|---|---|
| Ziele, Umfang und Nichtziele | [Anforderungen](Doc/01_SCOPE_AND_REQUIREMENTS.md) |
| SI-Einheiten, Vorzeichen und schwache Formen | [Mathematisches Modell](Doc/02_MATHEMATICAL_MODEL.md) |
| Elementgleichungen und Geometrie | [Tet4-Formulierung](Doc/03_TET4_FORMULATION.md) |
| Modulgrenzen, Datenmodelle und Ziel-API | [Softwarearchitektur](Doc/04_SOFTWARE_ARCHITECTURE.md) |
| Wissenschaftliche Nachweise | [Verifikation und Validierung](Doc/05_VERIFICATION_AND_VALIDATION.md) |
| Stationärer v1.0-Gültigkeitsbereich und Grenzen | [Nutzungsgrenzen](Doc/26_VALIDITY_AND_LIMITS.md) |
| Stationäre thermische Ausdehnung | [Thermomechanik](Doc/14_THERMOMECHANICS.md) |
| Fachliche und technische Entscheidungen | [Risiken und Entscheidungen](Doc/07_RISKS_AND_DECISIONS.md) |
| Literatur und Referenzdokumentation | [Quellen](Doc/08_REFERENCES.md) |
| Versionierung und Freigaben | [Release-Prozess](Doc/09_RELEASE_PROCESS.md) |

## Entwicklung

Der Python-Code liegt unter `src/fem_toolbox/`. Die [Teststruktur](tests/README.md)
unterscheidet Unit-, Patch-, Konvergenz- und Regressionstests;
[examples/](examples/README.md) enthält drei Einzelfeldbeispiele und ein
reproduzierbares Thermomechanik-Beispiel.

Aus dem Repository-Root mit aktivierter Entwicklungsumgebung prüfen:

```sh
python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m mypy src/fem_toolbox
git diff --check
```

GitHub Actions führt die Tests auf Python 3.11, 3.12 und 3.13 aus; Linting,
Formatierung und Typprüfung laufen auf Python 3.13.

Beiträge folgen den [Entwicklungsrichtlinien](AGENTS.md): zuerst ein passendes
Issue, anschließend ein fokussierter Branch und Pull Request. Ausführbare
Änderungen entstehen testgetrieben; reine Dokumentationsänderungen werden auf
Inhalt, Links und Konsistenz geprüft. Jeder abgeschlossene Meilenstein erhält
nach bestandenem Release-Gate eine eigene Version und einen annotierten Git-Tag.

## Lizenz und Weitergabe

Der aktuelle Quellstand steht unter der [MIT-Lizenz](LICENSE) und darf unter
Beibehaltung des Copyright- und Lizenzhinweises verwendet, verändert und
weitergegeben werden. Die früher veröffentlichten Tags und Artefakte erhalten
durch diese Änderung nicht nachträglich neue Lizenzmetadaten. Laufzeitpakete
wie NumPy, SciPy und meshio behalten ihre eigenen Bedingungen; die
[Bestandsaufnahme](Doc/27_DISTRIBUTION_INVENTORY.md) nennt die geprüften
Drittbestandteile und den Umfang der Paketarchive. Wheel und sdist werden als GitHub-Release-Artefakte bereitgestellt. Ein
öffentlicher PyPI-Upload ist damit nicht verbunden.
