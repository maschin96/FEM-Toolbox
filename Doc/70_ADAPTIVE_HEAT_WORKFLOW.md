# Adaptiver stationärer Wärme-Workflow (#68/#69)

M10-Entwicklungskandidat `1.5.0.dev1`, Zielrelease v1.5.0. Die finale Freigabe
#70 bleibt offen; keine abgeschlossene main-/Tag-/Release-Abnahme wird behauptet.

## Öffentliche API

```python
from fem_toolbox.adaptive import AdaptiveOptions, HeatProblem, solve_adaptive_heat
from fem_toolbox.materials import ThermalMaterial

problem = HeatProblem(
    ThermalMaterial.isotropic(1.0),
    {int(i): 300.0 for i in mesh.facet_sets["wall"]},
    source=10.0,
)
run = solve_adaptive_heat(
    mesh,
    problem,
    AdaptiveOptions(
        tolerance=0.1,
        max_steps=8,
        max_dofs=20000,
        fraction=0.5,
    ),
)
run.export("results/adaptive_heat.json")
print(run.status, run.diagnostic)
```

`mesh` ist ein vollständiges konformes Tet4-Netz in Metern. `wall` bezeichnet
Randfacettenzeilen, keine Knoten oder externen IDs. Je zusammenhängender
Wärmeregion sind geeignete Dirichlet-/Robin-Anker erforderlich. `material` ist ein
`ThermalMaterial` oder eine geordnete Folge je Eingabezelle. Dirichlet-Facetten
tragen skalare K-Werte oder räumliche Callables; Flüsse und Robin-Facetten verwenden
die bestehenden Wärmeverträge. Überlappung Dirichlet/Fluss/Robin wird abgelehnt,
gemeinsame Dirichlet-Knoten dürfen keine widersprüchlichen Werte besitzen.

`AdaptiveOptions` und `HeatProblem` sind eingefrorene Datensätze; Mapping-Eingaben
werden kopiert, Materialfolgen in Tupel überführt. Callables müssen selbst
reproduzierbar sein. Toleranz/Budgets/Stagnationspriorität und die numerische
Verfeinerungsstrategie stehen in [ADR-055](69_ADAPTIVE_REMESHING_CONTRACT.md).
Ungültige Optionen werfen `ValidationError`; erwartete Netz-/Solver-/Adapterfehler
liefern `error` und erhalten bereits abgeschlossene Stufen. Ein Fehler ist
kein erfolgreicher physikalischer Nachweis. `run.stages` ist ein Tupel von
`AdaptiveStage`, mit unveränderlichen Mesh-/Rohfeld-Arrays und Elternzeilen.

`AdaptiveHeatResult.to_dict()`/`.export(path)` schreiben striktes Schema-1-JSON:
`kind=fem-toolbox-adaptive-heat`, Endstatus, Diagnose, Optionen und Stufen. Jede
Stufe enthält vollständiges Netz samt Tags/Gruppen, SHA-256, DOFs, Markierung,
Elternzuordnung, eta/Einzelindikatoren, Temperatur K, Gradienten K/m, Flüsse W/m²,
Reaktionen/Residuen W, Quellen-/Randwärmeströme W und vorhandene Solverprovenienz.
Damit sind Netze und Ergebnisse rekonstruierbar. Zur erneuten Ausführung räumlicher
Funktionen bleibt das originale Problemskript erforderlich. Export überschreibt
einen ausdrücklich angegebenen Pfad; Wiederaufnahme einer begonnenen Schleife
und ein TOML-Adaptivitätsmodell sind kein Umfang dieses Vertrags.

`mesh.refinement.refine_tet4(mesh, marked, max_points=...)` ist der direkte
Verfeinerungseinstieg. `Refinement` enthält Mesh und unveränderliche Kind-Eltern-
Zeilen. Seine Hilfsvalidatoren sind Forschungs-/interne Verträge.

## Verifikation und Nutzengrenze

Die Unit-Tests prüfen Mehrgebietserhalt, Randtopologie, deterministische lokale
Verfeinerung, Budget, ungültige Eltern, Toleranz, Stagnation, Schrittbudget,
Solver-/Adapterfehler, gemischte Fluss-/Robin-/Dirichlet-Ränder und eine analytische
Wärmeleitfähigkeitsgrenze mit kontinuierlichem Fluss. Identische Läufe liefern
identische JSON-Records; Zeitmessungen stehen getrennt im Benchmark.

`examples/adaptive_heat/run_adaptive.py` verwendet unabhängige analytische
Referenzen auf dem Einheitswürfel (k=1 W/(m K)):

- glatt: T=300+sum_j x_j(1-x_j), Q=6 W/m³;
- lokalisiert: T=300+exp(-80(x-0.5)²), Q=(160-25600(x-0.5)²)exp(-80(x-0.5)²).

Alle Randwerte sind exakt vorgeschrieben. Wahrer Energiefehler ist
sqrt(integral |grad(T-T_h)|² dV). Gauss-Duffy-Integration Ordnung 6/7
prüft den Auswertungsfehler (relative Parität 2e-4). Diese Fehlerauswertung
verwendet weder den Residualindikator noch die Lastquadratur des Solvers.
Die Tet4-Quellenintegration bleibt Grad 2; besonders die Gaußschicht enthält
auch Diskretisierungs-/Lastquadraturfehler, die nicht verschwiegen werden.

Sieben Stufen je Strategie werden berichtet. Der glatte adaptive Fehler sinkt
von 0.25 auf 0.07753 sqrt(W K); uniform auf 0.04724. Beim lokalisierten Fall
sinkt er adaptiv von 1.96263 auf 1.34228. Die nach Pilotexploration und vor dem
finalen Nachweis eingefrorene Grenze **1.5 sqrt(W K)** wird bei 514 adaptiven
gegenüber 729 uniformen DOFs erstmals unterschritten (29.5 % weniger DOFs).
Kein unabhängiges, blind vor dem Pilotversuch gewähltes Ziel wird behauptet.
Die Regression friert diese Grenze und den Nutzenvergleich ein.
Uniform bedeutet dieselbe konforme Kantenstrategie mit Markierungsanteil 1;
keine Qualität/Geometrieänderung verfälscht den Vergleich.

Dies ist kein pauschaler Laufzeitvorteil. Kantenstern-Abschluss, Validierung und
Python-Verfeinerung haben messbaren Aufwand; der glatte Fall illustriert Grenzen.
Es gibt keine garantierte Effektivität, Formregularität, gekoppelte oder transiente
Adaptivität, Glättung/Koarsening oder dynamischen Feldtransfer.

### Gemessene Umgebung und Kosten

Die versionierte [Vergleichstabelle](evidence/m10/comparison.json) enthält alle
28 Fehler-/DOF-/Indikatorzeilen; [Umgebung und Artefakthashes](evidence/m10/environment.json)
bezeichnen den tatsächlichen Lauf. Laufzeiten umfassen Solve, Estimate, Mark,
Remesh und Vertragsvalidierung, jedoch nicht die unabhängige Fehlerintegration
oder den JSON-Export. Sie sind keine isolierten Solverzeiten.
