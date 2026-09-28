# Bedarfsgerechte Tet4-Flächengeometrie (#37b)

Status: integrierter Teilbeitrag zu #37. Der [Vergleich #37c](54_ACCURACY_COST_COMPARISON.md)
ist ebenfalls integriert; [Gesamtabgleich und Abnahme](56_PERFORMANCE_ACCEPTANCE.md).
Die folgenden Mess- und Testangaben dokumentieren den damaligen Optimierungsschritt.

## Vorab festgelegte Entscheidung

Die [37a-Profile](36_PERFORMANCE_BASELINE.md) zeigen `_build_faces` und
Kreuzprodukte als bedeutende Geometriekosten. Volumenoperatoren benötigen
keine Flächen, dennoch berechnet der Konstruktor vier Flächen pro Element.
Im Elektrodenmodell kommen zusätzliche Geometriekonstruktionen für den
Topologie-Preflight hinzu. Der gewählte Schritt verschiebt ausschließlich die
Flächenkonstruktion auf den ersten tatsächlichen Zugriff und bewahrt deren
Ergebnis je Geometrieobjekt. Keine globale Memoisierung, keine Änderung der
Kreuzproduktformel, keine Änderung von Solver- oder Lastoperatoren.

Vor Implementierung festgelegtes Gate: gegenüber dem M7-Stand `02853db`
mindestens 15 % weniger Median-End-to-End-Zeit im Elektrodenfall für N=6
**und** N=10; höchstens 5 % Regression in stationärer und transienter Wärme
bei beiden Größen. Je fünf frische Prozesse pro Fall/Größe/Stand, fünf
transiente Schritte, sechs angeforderte native Threadlimits auf eins.
Messungen laufen ohne parallele Tests oder weitere Benchmarkprozesse.
Streuung und Speicher werden berichtet; keine Zeitgrenze als CI-Test.

Die numerische Parität wird anhand der bestehenden analytischen Feld-/Bilanzgates
und zusätzlich durch Vergleich der exportierten primären Felder geprüft.
Die unveränderte explizite Kreuzprodukt-Flächenroutine und bestehende
analytische/permutationsinvariante Geometrietests bleiben die Referenz.

## Cache-Vertrag und Abgrenzung

Die Eingabekoordinaten werden weiterhin sofort kopiert und validiert;
Entartung und unzulässige Geometrie werden schon im Konstruktor abgelehnt.
Flächennormalen, Flächen und Schwerpunkte bleiben unveränderlich. Wiederholter
Zugriff auf dasselbe Objekt verwendet dieselben Flächen. Neue Geometrieobjekte
haben eigenen Zustand, auch nach Mutation der ursprünglich übergebenen Arrays.
Nebenläufiger Erstzugriff darf die deterministische Berechnung wiederholen;
kein Exactly-once-Vertrag über Threads wird eingeführt.

Es gibt keinen Cache für zeitabhängige Operatoren, Randwertstrukturen oder
Faktorisierungen. Daher ist hier keine Invalidierung anhand von Zeitschritten
oder Randwertänderungen erforderlich. Die in #37a betrachtete
Faktorisierungswiederverwendung war kein dominanter Engpass und wird in diesem
Schritt nicht eingebaut. Matrixkopien und COO-Speicher bleiben unverändert;
daraus werden keine weiteren Leistungszusagen abgeleitet.

## Reproduktion und Herkunft

Je Stand aus dem Repository-Root (frisches Ausgabeverzeichnis erforderlich):

```sh
PYTHONPATH=src .venv/bin/python benchmarks/profile_baseline.py \
  --sizes 6 10 --steps 5 --repeats 5 --output /tmp/fem37b-reference
# Nach Anwendung des dokumentierten Tet4Geometry-Diffs:
PYTHONPATH=src .venv/bin/python benchmarks/profile_baseline.py \
  --sizes 6 10 --steps 5 --repeats 5 --output /tmp/fem37b-candidate
```

Beide Messungen verwenden denselben Git-Basiscommit `02853db753e2b7a1b0614631b2a710088ee7f5e4`.
Die Arbeitskopie ist als dirty ausgewiesen: vor der Änderung wegen neuer Tests
und des Messplans, danach zusätzlich wegen des exakt beschriebenen
Produktionsdiffs. Während eines Messblocks bleibt der Produktionscode unverändert.
SHA-256 des Benchmarkskripts ist in beiden Fällen
`746eb54ec420bbcbde836890e8ff072b1d6fdfd271b815ff7f3f57fdf24fbbaf`;
der einzige geänderte Produktionspfad ist `src/fem_toolbox/elements/tet4.py`:

- vorher: `88e79d2c864c84a3e64c134420a9afbd0994317d25052b9c483049140ba5283e`
- nachher: `83264c4d33178347bfdf09a14b5e16a3da94fb26d7099cf98e1b821dbcf1edf8`

Die Erhebung erfolgt blockweise (zuerst Referenz, dann Kandidat), nicht
randomisiert oder unter exklusiver CPU-Reservierung. Dadurch bleiben
zeitabhängige Hintergrundlast und thermische Einflüsse mögliche Störgrößen.
Die Streuung ist eine Stichprobenbeschreibung, kein Konfidenzintervall.
Tracing-Peak und Prozess-RSS-Hochwasser haben die unterschiedlichen, im
[Messvertrag](../benchmarks/README.md) beschriebenen Bedeutungen. Keine
Aussage über maximale Modellgrößen, MPI oder verteilte Skalierung.

## Messergebnis (26.09.2026)

Alle vorab festgelegten Laufzeitgates sind erreicht. Werte in Millisekunden,
Speicher in dezimalen MB; Zeit ist der Median aus fünf Wiederholungen,
Speicher das Maximum aus fünf Wiederholungen.

| Fall | N | Vorher ms | Nachher ms | Reduktion | Tracing-Peak MB vorher → nachher | RSS-Peak MB vorher → nachher |
|---|---:|---:|---:|---:|---:|---:|
| heat | 6 | 165.30 | 88.90 | 46.2 % | 8.84 → 5.27 | 95.34 → 83.17 |
| heat | 10 | 781.32 | 415.42 | 46.8 % | 46.14 → 29.29 | 192.89 → 142.67 |
| electrode | 6 | 434.38 | 345.17 | 20.5 % | 54.59 → 51.01 | 214.94 → 207.90 |
| electrode | 10 | 1995.21 | 1621.11 | 18.7 % | 266.29 → 249.46 | 913.92 → 879.97 |
| transient | 6 | 204.60 | 139.63 | 31.8 % | 9.02 → 5.49 | 99.70 → 88.98 |
| transient | 10 | 970.19 | 614.90 | 36.6 % | 45.87 → 29.26 | 224.40 → 173.38 |

Unveränderte Rohdaten mit allen 60 Messläufen, Min/Max, Stichprobenstandardabweichung,
Phasen, DOFs, nnz und vollständiger Umgebung:
[Referenz](../benchmarks/baselines/m7-37b-before-macos-arm64.json),
[Kandidat](../benchmarks/baselines/m7-37b-after-macos-arm64.json).
Alle 90 erzeugten VTU-Dateien (einschließlich aller transienten Zustände) sind
zwischen Referenz und Kandidat bytegleich. Eingabehashes, faktorisierte Systeme,
analytische Fehler sowie Residual-/Bilanzkennwerte sind ebenfalls identisch.

Gemessen auf demselben macOS-27.0-arm64-Rechner mit zehn logischen CPUs,
CPython 3.13.7, NumPy 2.4.6/Accelerate, SciPy 1.18.1 und meshio 5.3.5.
Die Abfrage des konkreten CPU-Modells war gesperrt (`cpu_probe_error` in den
Rohdaten); native Threadzahlen wurden nicht gemessen, sondern nur angefordert.
Größter absoluter Primärfeldfehler: `4,15e-12 K/V`; größte gemessene
Bilanzabweichung: `4,13e-11 W/N/C` (verschiedene Fälle/Einheiten, kein
zusammengefasster dimensionsloser Fehler).

| Fall | N | Vorher Min–Max ms | Nachher Min–Max ms | Stichproben-Std. vorher/nachher ms |
|---|---:|---:|---:|---:|
| heat | 6 | 163.16–169.24 | 87.62–103.23 | 2.25 / 6.51 |
| heat | 10 | 771.60–791.44 | 412.71–419.03 | 8.31 / 3.11 |
| electrode | 6 | 429.15–441.62 | 340.35–355.95 | 4.53 / 5.84 |
| electrode | 10 | 1864.59–2038.10 | 1597.07–1649.40 | 66.33 / 20.29 |
| transient | 6 | 204.35–208.75 | 135.88–141.27 | 1.92 / 2.50 |
| transient | 10 | 947.57–997.77 | 613.62–652.95 | 18.80 / 17.35 |

Die getrennten Profile bestätigen die eingesparte Arbeit, nicht nur die Laufzeit:

| Fall N=10 | Flächenkonstruktionen vorher | nachher |
|---|---:|---:|
| heat | 6000 | 0 |
| electrode | 12000 | 6000 |
| transient | 6000 | 0 |

Geometriekonstruktionen selbst werden weiterhin ausgeführt. Im Elektrodenfall
berechnet der separate Randtopologiepfad weiterhin seine benötigten Flächen;
Wiederverwendung ganzer Geometrieobjekte über Solvergrenzen ist keine Behauptung
dieses Schritts. In den sechs gemessenen Konfigurationen wurden weder
Laufzeit- noch Speicherregressionen beobachtet. Unabhängige neue Systeme und
Lastprofile können andere Ergebnisse haben.

## Tests und Integrationsstand

Vor dem Refactoring: 499 Tests bestanden, 43 optionale PETSc-Tests übersprungen.
Der neue deterministische Arbeitsmengentest scheiterte vor der Änderung mit
`1 != 0` Flächenkonstruktionen beim reinen Volumenzugriff. Nach der Änderung:
57 gezielte Geometrie-/Topologietests und 500 Gesamttests bestanden;
43 optionale PETSc-Tests lokal übersprungen. Der neue Test prüft zusätzlich
Wiederverwendung, getrennte Geometrieinstanzen, Mutationen der ursprünglichen
Eingabe und unveränderliche Resultate. Vorhandene Tests prüfen Entartung
weiterhin direkt am Konstruktor, Normalen und alle Knotenpermutationen.

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

Die numerischen Gleichungen und Diskretisierungen sind unverändert; bestehende
Patch-, Bilanz- und Konvergenztests sind Teil der Gesamtsuite. Kein neues
Hardware-Zeitgate, keine abgeschwächten numerischen Toleranzen. Kein unabhängiger
Reviewer beteiligt; keine Reviewfreigabe behauptet. Der spätere Tet10-/PETSc-Vergleich und die Gesamtabnahme sind oben verlinkt.
