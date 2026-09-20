# Softwarearchitektur

Status: Entwurf

## 1. Grundsätze

Die Architektur trennt mathematische Verantwortlichkeiten. Elementroutinen
kennen weder Dateiformate noch globale Solver. Randbedingungen verändern keine
Materialobjekte. Postprocessing rekonstruiert Größen aus einer unveränderlichen
Lösung und dokumentiert, ob Werte roh, gemittelt oder projiziert sind.

Zentrale Entwurfsregeln:

- explizite Datenflüsse statt globalem Zustand
- unveränderliche, validierte Materialparameter
- feldweise DOF-Verwaltung
- kleine Referenzimplementierung vor Optimierung
- Sparse-Assembly über COO-Tripletts, danach CSR/CSC
- Abhängigkeiten auf Infrastrukturgrenzen beschränken

## 2. Zielstruktur

```text
FEM-Toolbox/
├── README.md
├── Doc/
├── pyproject.toml
├── src/fem_toolbox/
│   ├── mesh/            # Netz, Tags, Gmsh-Adapter, Qualitätsmaße
│   ├── fields/          # Feld- und DOF-Definitionen
│   ├── materials/       # Elastik, Wärme, Dielektrik, Piezo
│   ├── elements/        # Tet4-Geometrie und Elementoperatoren
│   ├── physics/         # schwache Formen und Lastbeiträge
│   ├── assembly/        # globale Sparse-Assembly
│   ├── boundary/        # Dirichlet-, Neumann- und Elektrodenmodelle
│   ├── solvers/         # lineare und später transiente Solver
│   ├── post/            # abgeleitete Felder, Bilanzen, Export
│   ├── io/              # Modellkonfiguration und Provenienz
│   └── model.py         # öffentliche Orchestrierungs-API
├── tests/
│   ├── unit/
│   ├── patch/
│   ├── convergence/
│   └── regression/
└── examples/
```

## 3. Kerndatenmodelle

### Mesh

```text
points:          float[nnode, 3]
tetra4:          int[nelem, 4]
triangle3:       int[nface, 3]
cell_tags:       int[nelem]
facet_tags:      int[nface]
```

Indizes sind intern nullbasiert. Externe IDs werden beim Import getrennt
gespeichert. Tags adressieren Gebiete; überlappende semantische Gruppen werden
nicht durch stilles Überschreiben eines einzelnen Tags modelliert, sondern durch
benannte Auswahlsätze.

### Field und DofMap

Ein `Field` enthält Namen, Komponenten, Trägerknoten und physikalische Dimension.
Eine `DofMap` ordnet lokale Feld-DOFs globalen Indizes zu. Blockinformationen
bleiben erhalten, damit der Solver später Field-Split-Präconditioner bilden kann.

### Material

Materialklassen akzeptieren SI-Werte und liefern konstitutive Matrizen in der
projektweiten Konvention. Jede Instanz validiert:

- Tensorformen und Symmetrie,
- positive Definitheit der uncoupled Energieterme,
- orthonormale Materialorientierung,
- endliche Werte und plausible Größenordnungen.

Die gekoppelte Stabilität wird zusätzlich über die elektrische Enthalpie und
geeignete Schur-Komplemente getestet.

### BoundaryCondition

Dirichlet-Bedingungen werden als `(field, components, selection, value)`
repräsentiert. Werte dürfen konstant oder räumlich/zeitlich auswertbare Callables
sein. Widersprüchliche Zuweisungen an denselben DOF führen zu einem Fehler.

## 4. Assembly-Pipeline

```text
Modell validieren
  -> aktive Felder und DOF-Tabellen erzeugen
  -> Elemente/Gebiete durchlaufen
  -> lokale Blöcke und Lasten berechnen
  -> COO-Tripletts sammeln
  -> globale CSR/CSC-Matrix bilden
  -> Neumann-Beiträge integrieren
  -> Dirichlet-Bedingungen symmetrieerhaltend eliminieren
  -> skalieren und lösen
  -> Residuen/Reaktionen auf dem ursprünglichen System auswerten
  -> abgeleitete Größen und Bilanzen berechnen
```

Die ursprüngliche Matrix und rechte Seite bleiben logisch verfügbar, damit
Reaktionskräfte und Ladungen korrekt aus dem unmodifizierten Residuum folgen.

## 5. Lineare Solver

Die erste Implementierung verwendet SciPy:

- `spsolve` als Referenz für kleine und mittlere Probleme,
- CG nur für nachweislich symmetrisch-positive Einzelfeldsysteme,
- MINRES optional für das symmetrisch-indefinite piezoelektrische System,
- GMRES nur, wenn eine spätere Formulierung tatsächlich unsymmetrisch ist.

Eine dimensionsbasierte diagonale beziehungsweise blockweise Skalierung ist bei
gekoppelten Systemen Pflicht: mechanische und elektrische Koeffizienten können
viele Größenordnungen auseinanderliegen. Abbruchkriterien werden zusätzlich in
feldweisen physikalischen Residuen ausgewertet.

Für große Modelle ist ein optionales petsc4py-Backend mit FieldSplit- und
Schur-Komplement-Präconditionierung vorgesehen. Dieses Backend darf die
Element- oder Material-API nicht verändern.

## 6. Eingabe und Ausgabe

- Gmsh `.msh` als primäres externes Netzformat
- VTU/PVD für ParaView-Ausgabe
- YAML oder TOML für reproduzierbare Modellbeschreibungen
- NumPy-Arrays und Python-Objekte als direkte Forschungs-API

Der Importer prüft Zelltypen, Einheitenmetadaten, Tags und Orientierung. Wegen
der uneinheitlichen Behandlung mehrfacher physischer Gmsh-Gruppen wird das
Mapping in eigenen Tests mit repräsentativen Netzen abgesichert.

## 7. Öffentliche API

Die API soll drei Ebenen bieten:

1. `Model.from_file(...)` für reproduzierbare Fälle,
2. `Model(...)` für Python-Skripte,
3. reine Element- und Assemblyfunktionen für Forschung und Tests.

Keine dieser Ebenen soll versteckte globale Solveroptionen besitzen.

## 8. Abhängigkeitsstrategie

Zwingend: NumPy, SciPy. Standardadapter: meshio. Entwicklungsabhängigkeiten:
pytest, ruff, mypy. SymPy darf Testdaten für Manufactured Solutions erzeugen,
aber keine Laufzeitabhängigkeit des Solvers sein. Visualisierung erfolgt extern
über ParaView oder optionales PyVista.

## 9. Ausführbare Konventionen und Provenienz

`fem_toolbox.conventions` ist die ausführbare Quelle für räumliche Komponenten,
Voigt-Reihenfolge, technische Schubverzerrungen, Vorzeichen des elektrischen
Feldes und interne SI-Einheiten. Änderungen müssen gleichzeitig mit
`02_MATHEMATICAL_MODEL.md` und den zugehörigen Vertragstests erfolgen.

Jeder Solverlauf erhält eine `RunProvenance` mit Schema-Version,
Toolbox-Version, SHA-256 der kanonischen Modelleingabe, Solverkennung und allen
ergebniswirksamen Solverparametern. Die Parameter sind tief kopiert,
JSON-kompatibel und frei von nichtendlichen Gleitkommawerten. Zeitstempel und
lokale Dateipfade gehören nicht in den deterministischen Kern; sie können in
einer späteren äußeren Ausführungshülle ergänzt werden.
