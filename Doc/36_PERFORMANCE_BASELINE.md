# M7-Einstieg: gemessene SciPy-/Tet4-Baseline

Stand: 25.09.2026. Teilbeitrag **37a** zu [Issue #37](https://github.com/maschin96/FEM-Toolbox/issues/37).
Status: lokal gemessen und gegen analytische Referenzen geprüft.
#37b ist inzwischen als [gemessene bedarfsgerechte Flächenauswertung](44_GEOMETRY_PERFORMANCE.md)
im Issue-Branch umgesetzt. #37c (Tet10-/PETSc-Vergleich) und der Abschluss von
#37 bleiben offen. Die folgenden Zahlen bleiben die unveränderte historische Baseline.

## Reproduktion und Messumgebung

[Werkzeug, Referenzgleichungen und Messvertrag](../benchmarks/README.md),
[unveränderte Rohdaten aller 27 Läufe](../benchmarks/baselines/m7-37a-macos-arm64.json).
Befehl aus dem Repository-Root:

```sh
PYTHONPATH=src .venv/bin/python benchmarks/profile_baseline.py \
  --sizes 2 4 6 --steps 5 --repeats 3 --output results/baseline
```

Gemessen auf Apple M4, 10 logischen CPUs, macOS 27.0 arm64, CPython 3.13.7,
NumPy 2.4.6 (Accelerate), SciPy 1.18.1 und meshio 5.3.5. Jeder Messlauf startet
in einem eigenen Prozess; sechs BLAS-/OpenMP-Umgebungsvariablen sind auf 1
gesetzt. Die tatsächlich aktive native Threadzahl wurde nicht gemessen.
Keine exklusive CPU-Reservierung oder Kontrolle anderer Systemlast; drei
Wiederholungen zeigen Streuung, kein Konfidenzintervall. Die lokale Testsuite
lief während eines Teils dieser Erhebung ebenfalls: Ausreißer bleiben in den
Rohdaten sichtbar. Ein Optimierungsvergleich muss unter vergleichbarer,
möglichst ruhiger Last neu gemessen werden.

Produktionscode: `e4817931e5ecfbbec85debd88454abba052ac4d3` (v1.1.0).
Die Arbeitskopie war wegen des neuen Benchmarks/Tests dirty; `src/` unverändert.
SHA-256 des gemessenen, formatierten Benchmarkskripts:
`7c593472f61db66922c0c0741f9359943274cb8f9e874c4ed3067ed2600caa90`.
Die Rohdaten enthalten diesen Basiscommit, Eingabehashes und vollständige
Versions-/Buildangaben. Es wurden keine Produktionspfade optimiert.

## Uninstrumentierte Laufzeit

N bezeichnet Unterteilungen je Würfelkante; Zellen = 6 N³.
Zeit ist der erste vollständige Solveraufruf nach Netz-/Problemaufbau ohne
Profiler, Speichertracing oder Export. Transient: fünf Schritte, sechs Zustände.

| Fall | N | Median ms | Min–Max ms | Stichproben-Std. ms |
|---|---:|---:|---:|---:|
| heat | 2 | 9.73 | 9.00–10.29 | 0.65 |
| heat | 4 | 50.45 | 48.72–52.31 | 1.79 |
| heat | 6 | 163.11 | 160.87–225.74 | 36.82 |
| electrode | 2 | 19.20 | 19.01–21.49 | 1.38 |
| electrode | 4 | 131.46 | 126.92–142.73 | 8.14 |
| electrode | 6 | 437.93 | 418.22–464.10 | 23.02 |
| transient | 2 | 17.35 | 14.89–20.12 | 2.62 |
| transient | 4 | 73.12 | 72.63–103.65 | 17.77 |
| transient | 6 | 216.39 | 212.63–230.64 | 9.50 |

## Systemgrößen und Speicher

Ursprüngliche DOFs vor Randbedingungen, tatsächlich faktorisierte freie DOFs
und gespeicherte nnz der skalierten LU-Eingangsmatrix sind getrennt angegeben.
Für Wärme N=6 enthalten die Sparse-Daten auch Rundungsreste kleiner Einträge;
nnz ist keine aus der Gittertopologie abgeleitete theoretische Bandbreite.
Transient sind alle fünf Faktorisierungen gleich groß.
Speicher jeweils Maximum aus drei Wiederholungen; MB = 10⁶ Bytes.

| Fall | N | DOFs ursprünglich | DOFs LU | nnz LU | Tracing-Peak MB | Prozess-RSS-Peak MB |
|---|---:|---:|---:|---:|---:|---:|
| heat | 2 | 27 | 9 | 33 | 0.314 | 68.09 |
| heat | 4 | 125 | 75 | 415 | 2.450 | 75.55 |
| heat | 6 | 343 | 245 | 2053 | 8.844 | 94.52 |
| electrode | 2 | 108 | 28 | 208 | 1.372 | 71.97 |
| electrode | 4 | 500 | 176 | 1536 | 13.759 | 100.94 |
| electrode | 6 | 1372 | 540 | 8968 | 54.586 | 221.15 |
| transient | 2 | 27 | 27 | 223 | 0.372 | 69.39 |
| transient | 4 | 125 | 125 | 1333 | 2.590 | 76.76 |
| transient | 6 | 343 | 343 | 4051 | 9.016 | 98.04 |

Tracing läuft in einem separaten dritten Solveraufruf. Prozess-RSS umfasst
auch Imports, vorherige Durchläufe, Profiler und gehaltene Ergebnisse; beide
Spalten sind unterschiedliche Messgrößen. Sie erlauben keine Hochrechnung
einer maximal lösbaren Modellgröße oder eines reinen LU-Speicherbedarfs.

## Profilingbefund und nächste Optimierungsentscheidung

Median der separat instrumentierten Funktionszeiten bei N=6 in ms:

| Gemessene Kategorie | Wärme | Elektrode | Transient |
|---|---:|---:|---:|
| Geometriekonstruktion | 247.80 | 243.42 | 247.94 |
| Element-/Lastoperatoren und Sparse-Assembly | 41.17 | 132.74 | 98.53 |
| Dirichlet/Elektrodenkondensation | 0.60 | 262.12 | 4.89 |
| LU-Faktorisierung | 0.45 | 2.11 | 4.15 |
| LU-Solve | 0.02 | 0.05 | 0.12 |
| Nicht zugeordneter Solverrest | 48.73 | 134.21 | 57.94 |
| Separater Export ohne Profiler | 1.82 | 1.86 | 20.68 |

cProfile vergrößert besonders Python-lastige Phasen; diese Werte sind nicht
Summanden der uninstrumentierten Laufzeit. Medianwerte verschiedener Spalten
müssen ebenfalls nicht den Median der Gesamtsumme ergeben. Der Rest enthält
unter anderem Schleifen, Indexierung, Matrixkopien, Skalierung, Prüfungen,
Feldnachverarbeitung und Provenienz. Der Messvertrag nennt die exakten Grenzen.

Im Elektrodenfall zeigt das Profil 2592 Tet4-Geometriekonstruktionen für
1296 Elemente: Die Topologieprüfung innerhalb der Elektroden-Randbehandlung
baut Geometrie erneut auf. `Tet4Geometry._build_faces` und NumPy-Kreuzprodukte
gehören bei allen drei Fällen zu den größten eigenen Funktionszeiten.
Daher ist die Vermeidung wiederholter Geometrie-/Flächenauswertung ein
begründeter erster **Prüfkandidat für #37b**. Vor Umsetzung dort separat als
Ziel festlegen: mindestens 15 % weniger Median-End-to-End-Zeit im
Elektrodenfall N=6 und N=10 bei fünf frischen Wiederholungen pro Stand,
unveränderter numerischer Parität und höchstens 5 % Regression in den beiden
Wärmefällen. Diese Schwellen sind ein zukünftiges Vergleichsgate, kein hier
bereits erfüllter Leistungsnachweis und kein CI-Zeitlimit.

Fünf konstante transiente Schritte faktorisieren fünfmal. Wiederverwendung
bleibt fachlich interessant, ist auf diesen kleinen Netzen zeitlich aber
kein dominanter Engpass. Vor ihrer Umsetzung sind geänderte Zeitschritte,
Randwertstrukturen und Operatoren als Invalidierungsfälle zu testen.
Kein MPI-/PETSc-Nutzen oder allgemeines asymptotisches Skalierungsgesetz lässt
sich aus dieser Baseline ableiten. Der Vergleich bei gleichem nichttrivialem
Diskretisierungsfehler folgt erst in #37c nach #34/#35.

## Verifikation und Reviewstatus

Alle 27 Läufe erfüllen die vorab im Werkzeug festgelegten analytischen
Feld-, Fluss- und Bilanzgates. Größte absolute Abweichung des primären Felds:
2.05e-12 K/V; größte geprüfte Bilanzabweichung: 1.78e-11 W/N/C (Maxima über
unterschiedliche Fälle, kein dimensionsloser Gesamtfehler). Instrumentierung
reproduziert das jeweilige primäre Feld bitgleich.

TDD: vier neue Tests zunächst wegen des fehlenden Benchmark-Einstiegs rot;
nach Implementation grün. Erweiterung um Bilanzbericht und vollständige
Zeitserie zunächst dreimal `KeyError: verification`, danach grün. VTU-Felder
werden mit meshio unabhängig eingelesen und gegen die analytischen Felder
geprüft. Ungültige Größe wird abgelehnt. Kein künstliches Performance-Zeitgate.

Lokale Prüfung:

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

400 Tests bestanden, Ruff/Formatprüfung und mypy (50 Quelldateien) bestanden.
Der lokale direkte pytest-Aufruf ohne `PYTHONPATH=src` scheitert weiterhin am
Paketimport; dies ist ein Quelltest und kein neuer Artefaktinstallationsnachweis.
Bestehende Patch-/Konvergenztests sind in der Gesamtsuite enthalten; zusätzliche
Konvergenznachweise sind ohne Änderung der Numerik nicht erforderlich.
Kein unabhängiger Reviewer beteiligt; keine Reviewfreigabe behauptet.
GitHub-CI und Integration werden am PR verfolgt und sind von diesem lokalen
Nachweis getrennt. Das übergeordnete Issue bleibt offen.
