# Performance-Baseline (Issue #37a)

Entwicklerwerkzeug für die unveränderten SciPy-/Tet4-Referenzsolver, kein
installierter öffentlicher API-Einstieg. Aus dem Repository-Root:

```sh
PYTHONPATH=src .venv/bin/python benchmarks/profile_baseline.py \
  --sizes 2 4 6 --steps 5 --repeats 3 --output results/baseline
```

Der Ausgabeordner muss neu sein. Jeder Fall, jede Größe und jede Wiederholung
läuft in einem frischen Python-Prozess. Der Koordinator setzt die sechs im
Skript aufgeführten BLAS-/OpenMP-Threadvariablen vor dem Import auf 1.
Dies ist eine Konfiguration, keine Messung der aktiven nativen Threadzahl.
`--worker heat|electrode|transient --size N` ist der interne Einzellauf;
er übernimmt seine Threadumgebung unverändert. Alle Größen müssen positiv sein.
Bei nur einer Wiederholung ist die Stichprobenstandardabweichung `null`.

`baseline.json` enthält alle Rohmessungen und Median, Minimum, Maximum sowie
Stichprobenstandardabweichung der uninstrumentierten Solverlaufzeit.
Pro Lauf entstehen `report.json`, ein mit `pstats` lesbares `solve.prof`,
`solution.vtu` und beim transienten Fall die vollständige `series/`.
Die Berichte enthalten Git-Basiscommit/Dirty-Status, deterministischen
Solver-Eingabehash, Hardware, Betriebssystem, Python-/Bibliotheksversionen,
NumPy-Buildkonfiguration, angeforderte Threadlimits, DOFs und Nichtnullzahlen.
Eine gesperrte macOS-CPU-Abfrage wird als `cpu_probe_error` ausgewiesen.
Unter Windows ohne `resource` ist der RSS-Wert ausdrücklich `null`.

## Drei kontrollierte Fälle

Einheitswürfel mit `N³` kartesischen Zellen, je sechs konformen Tet4;
`(N+1)³` Knoten, vollständige äußere Triangle3-Flächen. Alle Angaben in SI.
Die synthetischen Materialdaten dienen der Verifikation, nicht der Darstellung
eines realen Werkstoffs.

- Wärme: `k=2`, links `T=300`, rechts `T=301`, übrige Flächen isoliert.
  Exakt `T=300+x`, `q=(-2,0,0)`.
- Piezo-Elektrode: `C=10 I₆`, `e₁₁=0.5`, alle anderen Kopplungen null,
  `κ=2 I₃`. Querbewegungen fest, `u_x=0` links, Zugspannung 2 rechts.
  Linke Elektrode geerdet, rechte schwebend mit Gesamtladung null.
  Aus `D_x=.5 ε_xx - 2 φ'=0` und `σ_xx=10 ε_xx + .5 φ'=2` folgen
  `u_x=2x/10.125`, `φ=x/20.25`, `D=0`. Innere mechanische und elektrische
  DOFs bleiben frei; die Floating-Elektrode durchläuft die reale Kondensation.
- Transiente Wärme: isolierter Würfel, `k=2`, `ρ=3`, `c=4`, `Q=12`,
  `T(0)=300`. Fünf konstante Schritte bis 1 s liefern exakt `T(t)=300+t`.
  Alle Schritte werden gerechnet und alle sechs Zustände exportiert.

Vorab festgelegte absolute Feld-/Flusstoleranz: `1e-9` in der jeweiligen
SI-Einheit. Bei Einheitsgeometrie, gut konditionierten synthetischen Materialien
und Feldern der Größenordnung 1 (Temperaturoffset 300) deckt dies akkumulierte
Double-Rundung ab; es ist keine Fehlertoleranz für allgemeine Modelle.
Wärme-/Kraft-/Ladungsbilanz wird absolut mit `1e-9 W/N/C`, freie relative
Residuen mit `1e-9` geprüft. Produktionssolver prüfen zusätzlich ihre eigenen
skalierten Bilanzen. Der instrumentierte Lauf muss das primäre Feld bitgleich
reproduzieren und erneut die analytischen Feld-/Bilanzprüfungen bestehen.
Es gibt keine Zeitgrenze als Unit-Test-Gate.

## Messgrenzen

1. `wall_s`: erster uninstrumentierter Solveraufruf nach Import, Netz- und
   Problemaufbau; enthält Geometrie, Assembly, Lösung, Nachverarbeitung und
   Provenienz, aber weder Netzgenerierung noch Export. Kein Warm-up.
2. Separater zweiter Aufruf mit cProfile und Funktionsuhren. `phases_s` ordnet
   nur die explizit in `Measurements.install` benannten Funktionsgrenzen zu:
   Geometriekonstruktion, Element-/Lastoperatoren und Sparse-Assembly,
   Dirichlet-Elimination/Elektrodenkondensation, Sparse-LU und Rückwärtseinsetzen.
   Verschachtelte gemessene Aufrufe werden nicht doppelt gezählt, sondern der
   äußeren Kategorie zugeschlagen. Prüfungen, Matrixkopien, Skalierung,
   Schleifen, Nachverarbeitung und Provenienz verbleiben teilweise im separat
   ausgewiesenen Rest `unclassified_profiled_s`. Die Kategorien sind keine
   vollständige Zerlegung fachlicher Solverphasen. cProfile verändert diese
   Zeiten; sie dürfen nicht direkt mit `wall_s` verrechnet werden.
3. Exportzeit wird anschließend ohne cProfile separat gemessen: primäres
   VTU-Feld, bei Transienten zusätzlich PVD/VTU, CSV und Manifest. Kein fsync,
   daher Betriebssystemcache-Einfluss. Das Feld heißt bewusst generisch `field`;
   seine Einheit ist K für Wärme, V für Piezo. Die Zeitserie hat SI-Feldnamen.
4. Dritter Solveraufruf ausschließlich für `tracemalloc`-Peak. Dieser erfasst
   beobachtete Python-/NumPy-Allokationen, nicht zuverlässig alle nativen
   Bibliothekspuffer. `process_peak_rss_bytes` ist zusätzlich das OS-Hochwasser
   des gesamten Workerprozesses seit Start: inklusive Imports, Profiling,
   Export und gehaltener Ergebnisse vorheriger Durchläufe. Es ist weder ein
   inkrementeller Solverbedarf noch ein prognostizierter Großmodellbedarf.

`factorized_systems` beschreibt die tatsächlich an LU übergebenen skalierten
freien Matrizen nach Randbedingungen und Elektrodenkondensation; `nnz` zählt
gespeicherte Sparse-Einträge. `dofs` im Hauptbericht zählt ursprüngliche DOFs.
Bei einem vollständig vorgeschriebenen System (z. B. Wärme mit N=1) gibt es
keine Faktorisierung; Liste und Solveanzahl sind leer/0, entsprechende Phasen
fehlen. Der transiente Fall protokolliert eine Matrix pro Zeitschritt.

Diese kleinen affine/uniformen Fälle messen Problemgrößenskalierung bei
identischem analytischem Fehlergate. Sie ersetzen keine nichttrivialen
Konvergenzbenchmarks, kein MPI-, PETSc-, Tet10- oder Großmodell-Memory-Gate.
Ergebnisse sind hardware- und lastabhängig. Der [Baseline-Bericht](../Doc/36_PERFORMANCE_BASELINE.md)
nennt die gemessene Umgebung, Streuung und Folgerungen für #37b.
