# M11: Topologie-Studien und Geometrie-Reanalyse

19.10.2026 · Freigabestand 1.6.0 · Veröffentlichungsnachweis in [#83](https://github.com/maschin96/FEM-Toolbox/issues/83).
Die Themenpriorität M11 wurde ausdrücklich als nächster Ausbau beauftragt.
Deshalb ist nach v1.5.0 das nächste Minor-Ziel v1.6.0; die ältere v1.9.0-Angabe
in Planungswelle 3 war vorläufig. Andere Milestones werden nicht als freigegeben
vorgezogen. [Mathematik und vorab festgelegte Gates](83_M11_TOPOLOGY_CONTRACT.md).

Zum Einstieg: [Was ist SIMP-Optimierung?](86_SIMP_EINFUEHRUNG.md)

## Vorbereitetes Beispiel

```sh
python -m pip install -e '.[dev]'
python -I examples/topology/run.py --output results/topology
```

Die Eingabe `examples/topology/reference/study.toml` referenziert zwei normale
Model-TOML-Dateien und ein Gmsh-Netz. Ein explizites `--config` wählt eine andere
Studie. Ausgabeordner müssen neu sein; es gibt keinen Resume-/Überschreibmodus.
Die Quaderabmessungen sind 2x1x1 m, E=1000 Pa, nu=0.3; die Stirnfläche x=0
ist vollständig fest, x=2 trägt in getrennten Fällen je 1 N in -y und -z.
Die positive Gewichtung 0.5/0.5 wird nicht automatisch normiert. Es handelt
sich um einen synthetischen Methodenbenchmark mit bewusst kleinen SI-Werten.

```toml
schema_version = 1
models = ["y.toml", "z.toml"]
weights = [0.5, 0.5]
[design]
volume_fraction = 0.6
filter_radius = 0.3
penalty = 3.0
minimum_stiffness = 0.001
solid_sets = ["passive"]
[optimizer]
max_iterations = 100
max_evaluations = 300
tolerance = 1e-7
```

Relative Model-/Netzpfade müssen innerhalb des Studienverzeichnisses bleiben.
Unbekannte Schlüssel werden abgelehnt. `design` erlaubt außerdem
`fixed_solid`, `fixed_void` (Zellindexlisten), `solid_sets`, `void_sets`
(Namen vorhandener Volumengruppen). Feste Voll-/Leergebiete dürfen sich nicht
überschneiden. Alle Modelle müssen dasselbe Tet4-Netz und dieselben effektiven
Materialtensoren besitzen. Die unveränderte Model-API übernimmt Einheiten,
Material- und Lastprüfung. Nichtnull-Verschiebungen und nichtstationäre Physik
werden abgelehnt. Nur konstante, entwurfsunabhängige Kräfte sind unterstützt.

## Python-Oberfläche

- `fem_toolbox.topology.LoadCase(name, conditions, body_force, tractions, point_loads)`
  kopiert die Eingaben; DOFs stammen aus `mechanical_dofs(mesh)`. Körperkraft
  hat N/m³, Facettentraktion N/m², Knotenlast N; keine Callables.
- `TopologyProblem(mesh, material, cases, filter_radius=..., weights=...,
  fixed_solid=..., fixed_void=..., penalty=3, minimum_stiffness=1e-3)`
  akzeptiert ein Material oder einen Tensor pro Zelle. `from_models` verwendet
  den additiven `Model.elasticity_inputs()`-Snapshot; der Snapshot ändert Model
  und Model-Dokument nicht. Es gibt keinen zweiten Material-/Lastparser.
- `physical_density(x)` und `evaluate(x)` verwenden Variablen in Eingabezellfolge
  `(ncells,)`, einschließlich der festen 0/1-Werte. `DensityEvaluation` enthält
  unveränderliche Entwurfs-/physische Dichte, C [J], v [1], dC/dx [J], dv/dx [1]
  sowie getrennte `ElasticResult`-Zustände je Lastfall.
- `optimize_topology(problem, volume_fraction, initial=None, options=...)` findet
  einen zulässigen Start oder prüft den expliziten Start. `OptimizationOptions`
  validiert Budgets und dimensionslose Toleranz. `OptimizationResult` enthält
  beste zulässige Dichte, Startziel, Auswertungszahl, Verlauf und Endstatus.
- `TopologyStudy.from_file(path).run(new_output)` ist der deklarative Weg.
  `read_topology_result(directory)` prüft Hashes und Eingabeprovenienz vor Lesen.
- `reanalyze_geometry(problem, evaluation, threshold, connect=False)` baut eine
  diskrete Vollmaterialgeometrie. `connect=True` ist der getrennte, protokollierte
  Pfad mit Verbindungszellen. `export_reanalyses(...)` führt beide Wege für drei
  Schwellen aus und exportiert alle gültigen Varianten, auch solche mit No-go-Gate.

`converged` bezeichnet einen lokalen SLSQP-Stopp, keine globale Optimalität.
`iteration_limit`, `evaluation_limit` und `solver_failed` bewahren den besten
zulässigen Zustand. Ein ungültiger oder nicht lösbarer Ausgangszustand erzeugt
ValidationError statt eines erfundenen Ergebnisses. Der Beispielaufruf liefert
Exitcode 0 nur bei konvergierter Optimierung und mindestens einem akzeptierten
Geometrieweg; Konfigurations-/Solver-/Geometriefehler ergeben Exitcode 2.

## Export und Wiederholung

`input.json` enthält effektive SI-Punkte, Zellen, Facetten, Tensoren,
Randbedingungen, alle Lasten, Gewichte, Designgebiete, Filter und Optimierer
sowie Hashes der TOML-/Netzdateien. `result.json` (eigenes Schema 1) enthält
Dichten, C/v, Verlauf, Zwischenmaterial, Bilanzen, Einheiten und RunProvenance.
`SHA256.json` bindet beide JSON-Dateien und `design.vtu`. Identische Eingaben
und Umgebung erzeugen identische JSON-Records; Laufzeiten sind nicht Teil
der deterministischen Solverprovenienz. Seeds sind null, der Algorithmus ist
nicht stochastisch. Historische Studien-/Model-Schemata bleiben unverändert.

`geometry/reanalysis.json` trennt reine Schwellen und Verbindungsvarianten,
führt hinzugefügte Originalzellzeilen, C-Verhältnis, Volumen, Residuum,
Kraft-/Momentenbilanz, kritische Zellen und Artefakthashes auf.
Geometrie-VTU enthält das Netz und maximale von-Mises-Zellspannung über die
Lastfälle. Die Dreieckswerte dieses Feldes sind Null-Platzhalter; nur Tet4
tragen Spannungswerte. Für jedes gültige Netz werden ausschließlich vorhandene
Tet4 und ihre vollständige äußere Triangle3-Hülle exportiert. Neu exponierte
Flächen sind lastfrei. Entfernte Lastflächen/-knoten werden diagnostiziert.
Körperkräfte auf entferntem Volumen sind kein unverändert reproduzierbarer
Lastfall und werden bei der Geometrie-Reanalyse abgelehnt.

## Grenzen und Fehlersuche

Eine SIMP-Leerzelle besitzt Reststeifigkeit und ist noch kein Loch. Die diskrete
Geometrie benötigt eine eigene Rechnung. Reine Schwellen können Teile trennen;
der optionale Verbindungsweg fügt ausdrücklich Material hinzu. Seine
0.9-Geometrievolumengrenze ersetzt nicht das 0.6-Optimierungsbudget.
Spannungsspitzen an idealen Klemmen/Lasten sind nicht netzkonvergent geprüft;
keine Festigkeits-, Fertigungs- oder experimentelle Freigabe. Der Filterradius
ist eine Regularisierungslänge, keine garantierte Mindestwandstärke.

Ein unerfüllbares Budget verlangt eine geänderte Problemdefinition, keine
stille Erhöhung des zulässigen Volumens. Starrkörpermoden und schlecht
konditionierte Zustände werden vom Mechaniksolver gemeldet. Budgetstopps
müssen als solche betrachtet werden. M11 verwendet SciPy in einem Prozess;
keine MPI-/GPU-/Cluster- oder allgemeine CAD-Zusage.

[Messmatrix und Freigabegrenzen](85_M11_VERIFICATION.md) dokumentieren auch
fehlgeschlagene Erstläufe und den tatsächlichen Review-/Releasezustand.
