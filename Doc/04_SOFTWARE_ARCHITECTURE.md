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

Die ausführbare Referenz ist `fem_toolbox.mesh.Mesh`. Koordinaten werden als
unveränderliches `float64`-Array, Konnektivitäten, IDs, Tags und Auswahlsätze als
unveränderliche `int64`-Arrays gespeichert. Der Konstruktor kopiert Eingaben und
prüft Formen, endliche Koordinaten, ganzzahlige Indizes, Indexgrenzen, eindeutige
externe IDs sowie wiederholte Knoten und topologisch doppelte Entitäten. Eine
Fehlermeldung für eine Konnektivitätszeile nennt deren externe ID, falls diese
vorliegt, andernfalls ihren nullbasierten internen Index.

`volume_sets` referenziert Tet4-Zeilen, `facet_sets` Triangle3-Zeilen. Indizes
innerhalb eines Satzes sind eindeutig; verschiedene benannte Sätze dürfen sich
bewusst überschneiden. Geometrische Qualitätsmerkmale wie Volumen, Orientierung
und Entartung gehören zu den Tet4-Operatoren des folgenden Arbeitspakets und
nicht in dieses reine Topologiemodell.

Vor Randintegralen prüft `fem_toolbox.mesh.inspect_boundary_topology(mesh)`
die Triangle3-Zeilen gegen die Tet4-Facetteninzidenz. Jede eingetragene
Dreiecksfläche muss genau einem Tet4 gehören. Der Bericht liefert je
kanonischer Facette die Eigentümer-IDs und je Randfläche Volumenelement,
externe IDs, nach außen orientierte Knoten,
Flächeninhalt und Normale; die Zuordnung `triangle_faces` folgt der
Triangle3-Zeilenfolge. Innere, fremde und nichtmannigfaltige Flächen werden
mit Entitäts-ID und physischer Gruppe diagnostiziert. Topologisch doppelte
Triangle3-Zeilen lehnt bereits `Mesh` ab. Nicht aufgeführte äußere Flächen
sind für Teilrandnetze zulässig; `require_complete=True` verlangt eine
vollständige Triangle3-Abdeckung und meldet fehlende Facetten. Überlappende
Sätze bleiben zulässig. Erst die spätere `Model`-Validierung bewertet, ob
die auf denselben Flächen verlangten Randbedingungen fachlich kollidieren.

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

`assign_region_materials` ist der Preflight für benannte `volume_sets`: Es
ordnet jedem aktiven Tet4 pro Feld genau ein typisiertes Material zu und
liefert die zellweise Herkunft für Ergebnis und Provenienz. Identische
Überlappungen sind zulässig, widersprüchliche Zuordnungen, Lücken und
ungenutzte Definitionen nicht. Der `Model`-Adapter (#25) ruft ihn
vor der Assembly auf. Vertrag und Zweimaterialnachweis stehen in
[Doc/19](19_MATERIAL_REGIONS.md).

### BoundaryCondition

Dirichlet-Bedingungen werden als `(field, components, selection, value)`
repräsentiert. Werte dürfen konstant oder räumlich/zeitlich auswertbare Callables
sein. Widersprüchliche Zuweisungen an denselben DOF führen zu einem Fehler.
Der piezoelektrische Referenzsolver unterstützt zusätzlich feste und
schwebende `Electrode`-Randbedingungen über benannte Triangle3-Sätze oder
Zeilenindizes. Schwebende Elektroden werden auf je einen elektrischen DOF
kondensiert; der Gesamtladungsvektor folgt der negativen elektrischen
Residualkonvention. [Doc/20](20_ELECTRODES.md) beschreibt API und Nachweis.

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
- GMRES/FGMRES auch für unsymmetrisch präkonditionierte Systeme, insbesondere
  die approximative Schur-Blockfaktorisierung aus ADR-026.

Eine dimensionsbasierte diagonale beziehungsweise blockweise Skalierung ist bei
gekoppelten Systemen Pflicht: mechanische und elektrische Koeffizienten können
viele Größenordnungen auseinanderliegen. Abbruchkriterien werden zusätzlich in
feldweisen physikalischen Residuen ausgewertet.

Der implementierte piezoelektrische Referenzsolver verwendet eine kongruente
Skalierung der beiden Feldblöcke nach ihrer medianen positiven Diagonale und
SciPy Sparse LU. CG ist für diesen indefiniten Operator gesperrt; MINRES bleibt
geplant, bis ein geeigneter symmetrischer Präconditioner verifiziert ist.
Details und Nachweise stehen in [Doc/17](17_PIEZOELECTRIC_SOLVER.md).

Die interne Backend-Grenze aus #34a erhält freie unskalierte Operatoren,
Feldpartitionen nach Kondensation und die positive Diagonalskalierung.
`ScipyBackend.factorize(...).solve(rhs)` liefert ein unveränderliches Ergebnis
mit `SolveReport`; wiederverwendete Faktoren besitzen ihre eigenen Daten.
Der [Backend-Vertrag](37_LINEAR_BACKEND.md) beschreibt reelle/komplexe Systeme,
M6-Teilprobleme, Fehler und die Abgrenzung zum Eigenwertsolver (ADR-024).

Der [optionale PETSc-Adapter (#34b)](38_PETSC_BACKEND.md) implementiert
einprozessiges Direktlösen mit reellen und komplexen Scalar-Konfigurationen
(ADR-025). [FieldSplit/Schur mit FGMRES (#34c)](39_PETSC_FIELDSPLIT.md)
ergänzt eine iterative Variante mit LU-Teilblöcken und expliziter
Nichtkonvergenzdiagnose (ADR-026). Dieses Backend darf die
Element- oder Material-API nicht verändern.

## 6. Eingabe und Ausgabe

- Gmsh `.msh` als primäres externes Netzformat
- VTU/PVD für ParaView-Ausgabe
- YAML oder TOML für reproduzierbare Modellbeschreibungen
- NumPy-Arrays und Python-Objekte als direkte Forschungs-API

Der Importer prüft Zelltypen, Einheitenmetadaten, Tags und Orientierung. Wegen
der uneinheitlichen Behandlung mehrfacher physischer Gmsh-Gruppen wird das
Mapping in eigenen Tests mit repräsentativen Netzen abgesichert.

Der ausführbare Adapter `fem_toolbox.io.read_gmsh` akzeptiert ausschließlich
`tetra` (Tet4) und `triangle` (Triangle3), normiert die Orientierung der
Volumenelemente und verlangt einen expliziten Skalierungsfaktor, wenn die
Gmsh-Koordinaten nicht in Metern vorliegen. `meshio` liefert bei Gmsh 4.1 alle
physischen Mitgliedschaften als Cell-Sets; der Adapter vereinigt diese mit den
aus Physical-Tags ableitbaren Sätzen. Der VTU-Adapter schreibt jeden Satz als
eigene binäre Zellmaske. Der einzelne `cell_tags`- beziehungsweise
`facet_tags`-Wert bleibt eine Kompatibilitätsdarstellung und ersetzt keine
überlappenden Auswahlsätze.

## 7. Öffentliche API

Die API bietet drei Ebenen:

1. `Model.from_file(...)` für reproduzierbare Fälle,
2. `Model(...)` für Python-Skripte,
3. reine Element- und Assemblyfunktionen für Forschung und Tests.

Keine dieser Ebenen soll versteckte globale Solveroptionen besitzen.
Der Stand des versionierten TOML-Adapters und seine Grenzen stehen in
[Doc/21](21_MODEL_SCHEMA.md).
Der stationäre Kombinationsmodus löst Wärme zuerst und übergibt absolute
Knotentemperaturen an den monolithischen Piezo-Solver; seine beiden Ergebnisse
und feldweisen Bilanzen bleiben getrennt ([Doc/22](22_THERMOPIEZOELECTRIC.md)).
Der vorgesehene stabile v1.0-Umfang, Formen, Einheiten, Fehler und Migration
der direkten Solver- und Model-API stehen in [Doc/29](29_PUBLIC_API_CONTRACT.md).
Die modellweite Prüfung von Randwertkonflikten und notwendigen Referenzen
läuft vor der Assembly ([Doc/23](23_MODEL_PREFLIGHT.md)).

## 8. Abhängigkeitsstrategie

Zwingend: NumPy, SciPy. Standardadapter: meshio ab Version 5.3.5 (MIT-Lizenz);
ältere 5.3-Versionen scheitern mit NumPy 2 beim Import
([Verträglichkeitsnachweis](28_DEPENDENCY_COMPATIBILITY.md)).
Entwicklungsabhängigkeiten: pytest, ruff, mypy. Die Entwicklungsumgebung begrenzt
NumPy derzeit auf Versionen
vor 2.5, weil die NumPy-2.5-Typstubs mindestens Python 3.12 voraussetzen, während
mypy den unterstützten Python-3.11-Vertrag prüft; die Laufzeitabhängigkeit bleibt
davon unberührt. SymPy darf Testdaten für Manufactured Solutions erzeugen, aber
keine Laufzeitabhängigkeit des Solvers sein. Visualisierung erfolgt extern über
ParaView oder optionales PyVista.

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

## 10. Algebraische Referenzbasis (Issue #11)

`fem_toolbox.fields.Field(name, components, nodes, unit)` beschreibt ausschließlich
aktive Trägerknoten. `unit` ist verpflichtende SI-Metainformation, kein automatischer
Einheitenkonverter. `DofMap` nummeriert Felder in Eingabereihenfolge und innerhalb
eines Feldes zuerst Knoten, dann Komponenten; für das volle Modell sind die Felder
als `[u, T, phi]` zu übergeben. Knoten- und Komponentenreihenfolgen bleiben erhalten.
`block(name)` liefert den zusammenhängenden Feldblock für spätere Blocksolver.
Leere Träger erzeugen keine DOFs. Die Zugehörigkeit der Träger zum konkreten Netz
muss der aufrufende Modelladapter prüfen.

`fem_toolbox.assembly.assemble(size, matrices, vectors, format="csr")` summiert
`MatrixBlock(rows, columns, values)` über COO-Tripletts zu CSR oder CSC. Separate
Zeilen- und Spaltenindizes erlauben rechteckige Kopplungsblöcke. Wiederholte Indizes
werden auch in `VectorBlock(dofs, values)` addiert. Formen, endliche reelle Werte
und Indexgrenzen werden geprüft. `LinearSystem` besitzt eigene Kopien; seine
Matrix- und Lastzugriffe liefern Kopien zum Schutz des Originalsystems.

`fem_toolbox.boundary.Dirichlet(dofs, values)` arbeitet auf aufgelösten globalen
DOFs. Räumliche oder zeitliche Funktionen werden vorab im jeweiligen Adapter
ausgewertet. `eliminate_dirichlet` setzt mit der Freiselektion `P` und dem Vektor
vorgegebener Werte `g` das System `P K P + I_fixed`, `P (f - K g) + g` auf.
Identische doppelte Vorgaben sind zulässig; unterschiedliche Werte werden ohne
implizite, einheitenabhängige Toleranz abgelehnt. Das ursprüngliche System bleibt
für `original.residual(solution) = K solution - f` erhalten. An freien DOFs steht
das Gleichgewichtsresiduum, an gebundenen DOFs die algebraische Reaktion in der
jeweiligen Gleichungskonvention.

Die Protokolle `fem_toolbox.physics.VolumeLoad` und `SurfaceLoad` definieren
`integrate(mesh, cell/facet, dofs) -> VectorBlock`. Zell- bzw. Flächenindizes sind
nullbasiert; die konkrete Physik verantwortet Quadratur, äußere Orientierung und
Lastvorzeichen. Die Infrastruktur legt keine Wärme- oder Ladungskonvention neu
fest. Konkrete Integratoren und Feldsolver folgen in #12–#14.

Die algebraischen Tests verwenden unabhängig aufgestellte Mini-Matrizen und ein
Zwei-Feder-System mit vorgeschriebener Verschiebung 2 und Endlast 3: Lösung
`[2, 5, 8]`, Reaktion `[-3, 0, 0]`. Ganzzahlige Zuordnungen und exakt darstellbare
Assembly-Summen werden exakt verglichen; für Lösung und Residuum gilt bei diesen
gut konditionierten Größen der Ordnung eins `atol=1e-14, rtol=0` als
Rundungstoleranz. Physikalische Patch- und Netzkonvergenztests sind erst mit den
konkreten Elementoperatoren anwendbar.


### M7-Erweiterung: quadratische Netztopologie

#35a ergänzt `Mesh` um Tet10/Triangle6 und generische Zugriffe auf die aktive
Zellfamilie. Netz-I/O, Randtopologie und Materialregionen verwenden diese Zugriffe;
die Tet4-Operatoren bleiben unverändert. Der thermische Teil von #35d integriert
Tet10 in stationäre und transiente Wärmesolver, Model, Preflight und Export.
Der [elektrische Teil](46_TET10_ELECTROSTATICS.md) ergänzt stationäre Elektrostatik
mit Schwerpunktfeldern und separat integrierter Energie.
[Statische Elastizität](47_TET10_ELASTICITY.md) ergänzt diesen Vertrag für
Spannungen, Dehnungen und Verzerrungsenergie. Die
[stationäre Thermomechanik](48_TET10_THERMOMECHANICS.md) integriert
punktabhängige Eigendehnungen und Energie und überträgt T → u über die Python-API.
[Stationäre Tet10-Piezoelektrizität](49_TET10_PIEZOELECTRIC.md) ergänzt
den monolithischen Solver, Elektroden, Model und Export.
[Stationäre Tet10-Thermopiezoelektrizität](50_TET10_THERMOPIEZOELECTRIC.md)
ergänzt punktabhängige Eigenverzerrungen und den Ablauf T → (u, φ).
Tet10-Modalanalyse, Newmark und Harmonik sind ergänzt.
Der punktabhängige Ergebnisvertrag
ist in [Dokument 45](45_TET10_HEAT_INTEGRATION.md) beschrieben. Knotenordnung, Kompatibilität und
Geometriegrenzen sind in [Dokument 40](40_TET10_TOPOLOGY.md) festgelegt.

## Tet10-Modalanalyse (#35d)

Konsistente Masse, offene/kurzgeschlossene Elektroden, Model und Modenexport
sind integriert. [Umfang und Verifikation](51_TET10_MODES.md).
Gesamt-Abnahme-/Kostenvergleiche bleiben offen.

## Tet10-Newmark (#35d)

Zeitintegration mit konsistenter Masse, bewegten Rändern, Last-/Spannungsereignissen,
Model und Zeitserienexport ist integriert. [Nachweise und Grenzen](52_TET10_NEWMARK.md).
Gesamt-Abnahme-/Kostenvergleiche bleiben offen.

## Tet10-Harmonik (#35d)

Komplexe Frequenzantworten einschließlich DC, Elektroden, Model und P2-Export
sind integriert. [Nachweise und Grenzen](53_TET10_HARMONIC.md).
Gesamt-Abnahme und Genauigkeits-/Kostenvergleiche bleiben offen.

## M8: Feldbildadapter

Der optionale Adapter `post.visualization` verwendet Matplotlib/Agg und
vorhandene Tet4-Randtopologie. Solver bleiben unabhängig; [Vertrag](62_FIELD_VISUALIZATION.md).

## Experimenteller CFD-Einstieg

`cfd/channel.py` enthält einen eigenständigen konservativen Kanaloperator.
Keine Abhängigkeit von Festkörper-Elementoperatoren oder der Model-API.
[Vertrag und Grenzen](68_EXPERIMENTAL_CFD_CORE.md).

## Experimentelle Druck-Geschwindigkeits-Kopplung

`cfd/stokes.py` besitzt den uniformen 2D-MAC-Operator, das monolithische
Sparse-Sattelpunktsystem und getrennte Feld-/Bilanzrekonstruktion.
API und Grenzen sind in [Doc/69](69_CFD_PRESSURE_COUPLING.md) dokumentiert.

## Vollständiger stationärer CFD-Referenzkern (#86)

`cfd/flow_model.py` besitzt validierte uniforme Netze, SI-Fluid und Randtypen.
`mac.py` assembliert Diffusion/Divergenz/Gradient, konservative Transportflüsse
und Oseen-Schritte. `flow.py` besitzt die nichtlineare Iteration und finale
feldweise Rekonstruktion. `flow_config.py` ist der strikte TOML-Adapter.
Die experimentelle Python-API bleibt außerhalb der stabilen FEM-Model-API.
[Vertrag, Grenzen und Referenzen](70_STATIONARY_FLOW_CORE.md).

## M12: RANS-/FSI-Entwicklungsarchitektur (#84, ADR-049)

Die geplante körperangepasste Quad4-Finite-Volumen-API bleibt getrennt vom
uniformen MAC-Referenzkern und vom stabilen FEM-Model-Schema. Zellzentrierte
RANS-Felder verwenden konservierte Face-Massenströme, SIMPLE und Rhie-Chow;
#85 liefert Import/Preflight; #226 verifiziert den neuen Netzoperator.
SA-noft2, InterfaceMap, harmonische MeshMotion und partitionierte stationäre
Kopplung besitzen getrennte Zuständigkeiten. Die Struktur verwendet weiterhin
den verifizierten Kleinverformungssolver auf ihrem Referenznetz.
[Ziel-API, Datenformen, Adapter, Provenienz und Gates](71_M12_RANS_FSI_CONTRACT.md)
sind beschlossen, noch nicht verfügbar. NumPy/SciPy/meshio bleiben die
Abhängigkeitsbasis; keine externe CFD-/Kopplungsbibliothek wird eingebettet.

## Fluidnetz-/Randadapter (#85)

`cfd/fluid_mesh.py` besitzt immutable Quad4-Geometrie, Face-Eigentümer,
Wandsegmentabstände und Netz-/y+-Gates. `cfd/fluid_domain.py` validiert die
SI-Domäne, vollständige BC-Zuordnung, Einlassfluss und periodische Paarung.
JSON-I/O ist eine getrennte Factory-Grenze; Geometrie/Preflight rechnen ohne
Dateizugriff. Kein neuer Solver, keine neue Bibliothek, keine FSI-Netzbewegung.
[API, Formate und Grenzen](72_FLUID_MESH_ADAPTER.md).

## Quad-FV-Grundoperatoren (#229 / #226a)

`cfd/quad_fv.py` besitzt SI-FV-Geometrie, periodische gemeinsame Faces,
gewichtete Least-Squares-Gradienten, nichtorthogonale Facekorrektur, volle
variable Spannung und kompakte Druck-/Rhie-Chow-Operatoren. Reine Feldoperationen
und Sparse-Assembly ohne I/O oder PDE-Orchestrierung. #230 ergänzt SIMPLE,
BC-Rekonstruktion, Transport und vollständige finale Solverdiagnostik;
#226 bleibt bis zur Abnahme beider Arbeitspakete offen.
[API, Vorzeichen und Operatornachweise](73_QUAD_FV_OPERATORS.md).

## Stationärer Quad-FV-Solver (#230 / #226b)

`cfd/quad_system.py` stellt die affine konservative Rekonstruktion und volle
Spannungs-/Transportblöcke bereit. `cfd/quad_flow.py` orchestriert Picard,
Krylov-beschleunigte SIMPLE-Korrekturen und frische SI-Feld-/Bilanzprüfung.
Eine konsistente kubische Zell-/Face-Rekonstruktion vermeidet Randübergangsfehler;
der kompakte Anteil der inneren Normalableitung bleibt erhalten.
Die Solver-API verlangt Uref/Lref und liefert immutable Felder, Faceflüsse,
Wandtraktion und sämtliche diskreten Energieanteile mit Provenienz.
FluidDomain bleibt Eingabeadapter, der MAC-Kern unabhängige Referenz.
[Algorithmus und Nachweise](74_QUAD_FV_SIMPLE.md).
