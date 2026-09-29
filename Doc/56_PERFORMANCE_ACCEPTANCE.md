# Performance-Gesamtabnahme (#37)

Stand: 28.09.2026, M7-Basis `a780649b41ea994641a9664d193d50660a726429`.
Die Implementierung aus PR #179, #187 und #198 ist integriert; die
[Tet10-Gesamtabnahme](55_TET10_ACCEPTANCE.md) aus PR #199 ergänzt die Voraussetzung
für den abschließenden Vergleich. Dieser Bericht gleicht die Kriterien aus #37
mit den unveränderten Messarchiven und dem aktuellen Code ab. Keine erneute
Messkampagne, keine neue Optimierung, keine Releasefreigabe von v1.2.0.

## Abnahmekriterien

| Kriterium | Nachweis | Ergebnis |
|---|---|---|
| 37a: reproduzierbare Baseline, getrennte Phasen, Speicher und Systemgrößen | [Baseline](36_PERFORMANCE_BASELINE.md), [Messvertrag](../benchmarks/README.md), 27 Läufe | Tet4-Wärme, Elektrodenmodell und fünf transiente Schritte, je drei Größen/Wiederholungen; Geometrie, Assembly, Randbehandlung, Faktor/Solve, Export, DOFs, nnz und Umgebung dokumentiert |
| 37b: profilgestützte Optimierung mit vorab gewähltem Ziel | [Flächengeometrie](44_GEOMETRY_PERFORMANCE.md), 30 Vorher-/30 Nachher-Läufe | Mindestens 15 % weniger Medianzeit im Elektrodenfall N=6/10; höchstens 5 % Regression in Wärmefällen: alle sechs Konfigurationen erfüllen ihre Gates |
| 37b: Kopien, COO und Wiederverwendung prüfen | Entscheidung und Cache-Vertrag unten | Keine unbelegte zusätzliche Optimierung; Zustandswiederverwendung auf unveränderliche Geometrie beschränkt |
| 37c: Tet4/Tet10 × SciPy/PETSc bei Fehlergrenzen und mehreren Größen | [Vergleich](54_ACCURACY_COST_COMPARISON.md), 72 Läufe | Zwei glatte Fälle, drei Größen, drei Wiederholungen, Fehler-/Zeit-/RSS-Budgets und explizit nicht erreichte Ziele |
| Abnahme: Parität, Verbesserung, Regressionen, Speichergrenzen, Referenz | Rohdatenabgleich unten, Geometrie-/Benchmarktests, analytische Referenzen | Gemessene Ziele erfüllt; keine Hardware-Zeitgrenze in CI, keine MPI- oder Großmodellzusage |

## Entscheidung über Optimierung und Zustand

Die getrennten Profile identifizieren die Tet4-Flächenkonstruktion als
bedeutenden Aufwand. Der implementierte Schritt berechnet Flächen erst bei
Bedarf und bewahrt sie je Geometrieobjekt. Die explizite Kreuzproduktformel
bleibt unverändert und wird gegen analytische Geometrie, Orientierung und
Permutationen geprüft. Volumen-/Lastoperatoren und Diskretisierung ändern sich
nicht. Die Gesamtzahl der Geometrieobjekte wird dadurch nicht generell reduziert.

Die weiteren Kandidaten wurden abgegrenzt:

- [Sparse-Assembly](../src/fem_toolbox/assembly/__init__.py) baut weiterhin COO
  aus Einträgen auf und konvertiert in das Ziel-Sparseformat. Matrix-/RHS-Kopien
  bewahren den bestehenden Besitzvertrag. Die Profile weisen Kopierkosten nicht
  separat vollständig aus; ein sicherer Gewinn durch Entfernen dieser Kopien
  ist damit nicht belegt. Tracing und RSS bleiben verschiedene Messgrößen.
- [Transiente Wärme](../src/fem_toolbox/solvers/transient_heat.py) bildet weiterhin
  das zum jeweiligen Zeitschritt gehörende System und löst es neu. Die Baseline
  zeigt fünf Faktorisierungen für fünf Schritte, aber auf diesen Netzen keinen
  dominierenden Faktorisierungsaufwand. Es wurde kein Faktor-/Operatorcache
  eingeführt. Eine spätere Wiederverwendung braucht eigene Tests für geänderte
  Zeitschritte, Operatoren und Randwertstrukturen.
- Der tatsächliche Cache besitzt kopierte, validierte Geometriekoordinaten;
  Zeit und Randbedingungen sind keine Cache-Eingaben. Der
  [deterministische Geometrietest](../tests/unit/test_tet4_geometry.py)
  `test_face_work_is_deferred_and_reused_per_geometry` prüft verzögerten Aufbau,
  Wiederverwendung, getrennte Instanzen, mutierte ursprüngliche Eingaben und
  unveränderliche Ergebnisse. Eine zeitabhängige Invalidierung wäre hier ohne
  entsprechenden Zustand kein fachlicher Test.

Keine zusätzliche Vektorisierung oder Parallelisierung wurde aus diesen
Messungen abgeleitet. Es gibt keine neue Leistungszusage für unveränderte Pfade.

## Erneuter Rohdatenabgleich

Am 28.09.2026 wurden die Archive unabhängig mit Python-Standardbibliothek
(`json`, `statistics`, `math`, `hashlib`) nachgerechnet, ohne die Benchmark-
Aggregationsfunktion aufzurufen. Alle vier Archive bleiben unverändert:

| Archiv | Läufe | SHA-256 |
|---|---:|---|
| [37a](../benchmarks/baselines/m7-37a-macos-arm64.json) | 27 | `c722bbfb3340b8124951585457bb6e72bf997c00261fc4350a2586df126ffcec` |
| [37b vorher](../benchmarks/baselines/m7-37b-before-macos-arm64.json) | 30 | `45da4cde681b820fd1d6ae637b33c8590fca01f0eb9702556bf1e055f4c82c9a` |
| [37b nachher](../benchmarks/baselines/m7-37b-after-macos-arm64.json) | 30 | `2ab67cae111c7f17e390e0f32cf4f7b488ad8efbee55a61afa64cc3c150e838f` |
| [37c](../benchmarks/baselines/m7-37c-linux-arm64.json) | 72 | `a7080296e64963be12d7b5c51772fa737dc323360ae412792694f3d0ebc8c2bd` |

Median, Minimum, Maximum und Stichprobenstandardabweichung aller Gruppen stimmen
mit den Einzelmessungen überein (relative Rechentoleranz 1e-14). Im Vorher-/Nachher-
Vergleich stimmen fallweise Eingabehash, Systemgrößen, nnz, Solveanzahl,
Primärfeldfehler, instrumentierte Parität und Bilanz-/Residualkennwerte exakt
überein. Die seinerzeitige Bytegleichheit der 90 VTU-Ausgaben ist im
[Messprotokoll](44_GEOMETRY_PERFORMANCE.md) dokumentiert; sie wurde bei diesem
Archivabgleich nicht erneut aus den nicht eingecheckten VTU-Dateien geprüft.

Neu aus den Rohzeiten berechnete Medianreduktion:

| Fall | N=6 | N=10 |
|---|---:|---:|
| Wärme | 46,2 % | 46,8 % |
| Elektrode | 20,5 % | 18,7 % |
| Transiente Wärme | 31,8 % | 36,6 % |

Für #37c wurden zusätzlich Wiederholungszahl, maximaler H1-Fehler, maximales RSS
und sämtliche 24 Genauigkeits- sowie 24 Kostenbudget-Auswahlen erneut berechnet.
Alle stimmen einschließlich der `null`-Einträge überein. Höchstes beobachtetes
RSS: 283,23 MiB; alle gemessenen Fälle erfüllen das 512-MiB-Vergleichsbudget.
Dies ist keine allgemeine Speicherobergrenze der Toolbox.

## Aussagegrenzen

Die 37a-Baseline war teilweise von gleichzeitig laufenden Tests beeinflusst;
sie dient der Engpasssuche. Der 37b-Vergleich wurde deshalb separat ohne parallele
Tests erhoben, jedoch blockweise und ohne exklusive CPU-Reservierung. Die
Streuung ist keine statistische Garantie. Hardware, Bibliotheken, Messgrenzen
und Threadanforderungen stehen jeweils im zugehörigen Bericht. macOS-Zeiten
aus #37a/b werden nicht mit Linux-Containerzeiten aus #37c verrechnet.

Im untersuchten Raster bleibt Tet4 bei groben Fehlergrenzen günstiger. Tet10
erreicht kleinere Fehler, aber Tet4 könnte mit zusätzlichen Netzgrößen weitere
Zeitbudgets ausnutzen. Ein nicht erreichtes Fehlerziel erlaubt kein extrapoliertes
Speedup-Verhältnis. PETSc zeigt in den kleinen direkt gelösten Vergleichsfällen
keinen allgemeinen Laufzeitvorteil. FieldSplit ist funktional verifiziert, hier
aber nicht als Performancevergleich gemessen. Keine MPI-, Adaptivitäts-,
asymptotische Skalierungs- oder experimentelle Validierungsbehauptung.

## Prüfungen und Integrationsstatus

Dieser Schnitt ist reine Dokumentation; kein künstlicher TDD-Red-Test und keine
neuen Zeitgates sind erforderlich. Die TDD-Nachweise der ausführbaren Änderungen
stehen in den Teilberichten. Die unveränderten Geometrie- und Benchmarktests
werden erneut ausgeführt; CI prüft zusätzlich den gesamten aktuellen PR-Stand.

```sh
PYTHONPATH=src python -m pytest -q tests/unit/test_tet4_geometry.py tests/unit/test_mesh_topology.py tests/regression/test_performance_baseline.py tests/regression/test_accuracy_cost.py
python -m ruff check .
python -m ruff format --check .
python -m mypy src/fem_toolbox
git diff --check
```

Lokal bestanden 36 gezielte Geometrie-/Topologie-/Benchmarktests sowie Ruff,
Formatprüfung, Mypy und `git diff --check`. Alle lokalen Linkziele der
geänderten Dokumente wurden auf Existenz geprüft.

Ein unabhängiges Review dieses Gesamtabgleichs liegt nicht vor. #37 kann nach
geprüftem Merge dieses Berichts in `codex/milestone-m7` explizit geschlossen
werden. Die Veröffentlichung mit Supportmatrix, Release-Artefakten, main-
Integration und Tag bleibt Gegenstand von #38.
