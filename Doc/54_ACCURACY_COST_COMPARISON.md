# Genauigkeits- und Kostenvergleich Tet4/Tet10 × SciPy/PETSc (#37c, #35)

## Vorab festgelegter Messvertrag

Keine Produktionsoptimierung in diesem Schritt. Unveränderte M7-Solver werden
mit `benchmarks/accuracy_cost.py` verglichen. Zwei synthetische glatte Fälle
auf dem Einheitswürfel, jeweils exakte Dirichlet-Daten auf dem gesamten Rand:

- Wärmegradient: k=1 W/(m K), T=1+x⁴+y⁴/2+z⁴/4 K,
  Q=-12(x²+y²/2+z²/4) W/m³.
- Biegung: E=100 Pa, ν=0, a=0,01 m⁻², u=(-3ax²z,0,ax³),
  εxx=-6axz, alle anderen Verzerrungen null, σxx=-6Eaxz,
  b=(6Eaz,0,0) N/m³. Schubfreie kubische Biegekinematik mit unabhängig
  hergeleiteter Volumenlast; kein freier Kragträger oder experimenteller Versuch.

Alle Koordinaten werden numerisch in SI eingesetzt. N=2/4/6 kartesische
Teilungen je Richtung, sechs Tetraeder je Würfel, geradkantige Tet10 aus
konformer Kantenhalbierung. Kein Inkompressibilitäts- oder Singularitätsfall.

Fehler werden unabhängig von Produktionsgeometrie/-quadratur mit expliziter
P1/P2-Basis und Duffy-transformierter 6×6×6-Gauss-Legendre-Quadratur integriert.
Berichtet werden relative L2-Norm und relativer H1-Halbnormfehler
`||grad(u_h-u)||/||grad(u)||` (für Wärme entsprechend T); bei Vektoren die
Frobeniusnorm des vollständigen Gradienten, keine mechanische Energienorm.
Die Polynome werden damit exakt bis auf Rundung integriert. Unabhängiger
Affintest prüft L2-Fehler √3 und verschwindenden Gradientenfehler für beide Ordnungen.

Vor Messbeginn festgelegt:

- Relative H1-Ziele 0,5 / 0,25 / 0,125. Pro Fall, Backend und Ordnung wird die
  schnellste gemessene Konfiguration gewählt, die das Ziel erfüllt. Keine
  Interpolation oder Extrapolation; nicht erreichte Ziele bleiben `null`.
- Kostenbudgets 1 / 5 / 20 s und maximal 512 MiB Prozess-RSS. Auswahl der
  kleinsten gemessenen H1-Abweichung unter beiden Grenzen, getrennt nach
  Fall/Backend/Ordnung. Die Budgets sind Vergleichspunkte, keine CI-Gates.
- Je drei Wiederholungen in frischen Prozessen, strikt seriell. Backend und
  Elementordnung sind innerhalb jeder Größe/Wiederholung verschachtelt.
  Kein paralleler Test- oder Benchmarkprozess während der Messung.
- SciPy-Direktbackend und nativer PETSc-Direktbackend (real) aus derselben
  installierten Toolbox/Containerumgebung. Kein FieldSplit-Leistungsvergleich.
- Vor dem Workerstart sechs BLAS/OpenMP-Umgebungsvariablen auf 1 setzen;
  angeforderte Limits, keine Messung tatsächlich aktiver Threads.
- Backend-/Wiederholungsparität am identischen Netz: rtol=1e-8,
  atol=1e-10 K bzw. 1e-11 m. Physikalische Residuen und Wärme-/Kraft-/Momenten-
  bilanzen <1e-9 in den jeweiligen SI-Einheiten. Keine Referenztoleranzänderung.

## Messgrenzen und Reproduktion

`wall_s` umfasst den ersten Solverlauf einschließlich Randdaten-/Materialaufbau,
Assembly, Faktorisierung, Solve, Nachverarbeitung und Provenienz. Netzaufbau,
Import und Backendkonstruktion (auch PETSc-Initialisierung) liegen davor;
Fehlerquadratur und NPY-Ausgabe danach. Keine Exportzeitmessung in diesem Schritt.
Die Backendhülle zeichnet tatsächlich übergebene freie Matrixgröße und gespeicherte
nnz sowie Faktorisierungs-/Solvezeit auf. Diese Unterzeiten liegen innerhalb der
Gesamtzeit; Matrixstatistik und Messhülle sind darin enthalten. Keine cProfile-
Instrumentation und keine behauptete vollständige Phasenzerlegung.

RSS ist das Prozesshoch seit Start, direkt nach dem Solver und vor der
Fehlerauswertung gemessen: inklusive Import, Mesh, PETSc-Laufzeit und Ergebnis.
Es ist weder inkrementeller Solverbedarf noch zuverlässige Großmodellprognose.
Der Vorher-RSS-Wert wird separat mitgeführt. Ohne OS-Unterstützung ist RSS null;
ein Speicherbudget gilt dann als nicht nachgewiesen.

```sh
docker build -f ci/petsc.Dockerfile -t fem-toolbox-petsc:test .
docker run --rm -v "$PWD:/work:ro" -v /tmp/fem37c:/results \
  fem-toolbox-petsc:test python benchmarks/accuracy_cost.py \
  --sizes 2 4 6 --repeats 3 --petsc-scalar real \
  --revision 98ca442fe8954dcac15316cd1367f2b7819a21ab --output /results/measured
```

Neues Ausgabeverzeichnis erforderlich. `comparison.json` enthält sämtliche
72 Rohmessungen, Median/Minimum/Maximum/Stichprobenstandardabweichung,
Fehler-/Budgetauswahl und Herkunft. `samples.json` sichert Zwischenstände;
je Worker bleiben Feld-NPY und Einzelbericht zur Paritätsinspektion erhalten.
`--revision` ist eine explizite Herkunftsangabe; zusätzlich werden der tatsächlich
importierte Produktionscode, Benchmark und Meshgenerator gehasht. Das installierte
Paket wird nicht durch den lokalen src-Pfad ersetzt.

## Ergebnisse vom 28.09.2026

[Unveränderte Rohdaten](../benchmarks/baselines/m7-37c-linux-arm64.json): 72 Läufe,
24 Konfigurationen mit je drei Wiederholungen. Host Apple M4; Docker Linux
aarch64, 10 vCPUs, 8 217 059 328 Byte RAM. Python 3.12.3, NumPy 1.26.4,
SciPy 1.14.1, PETSc 3.19.6 (float64, preonly/LU). Image-ID
`sha256:e59a76579629aa941a2e3413524aa8013698b9a03c3e9ecbfde86fcdb8d1c9d3`.
Produktionsrevision und Codehashes stehen im Archiv; Produktions- und
Benchmarkhash wurden mit dem lokalen Stand abgeglichen.

Zeit: Median ± Stichprobenstandardabweichung in Sekunden, kein Konfidenzintervall.
RSS: Maximum in MiB. Min/Max, L2-Fehler und Matrixstatistik stehen im Archiv.

| Fall | Element | N | H1-Fehler | SciPy s | PETSc s | RSS SciPy/PETSc MiB |
|---|---|---:|---:|---:|---:|---:|
| heat | Tet4 | 2 | 0.477011 | 0.006025 ± 0.000265 | 0.006863 ± 0.000485 | 62.6/98.1 |
| heat | Tet4 | 4 | 0.251733 | 0.037473 ± 0.001535 | 0.039548 ± 0.001104 | 66.1/99.3 |
| heat | Tet4 | 6 | 0.169466 | 0.131161 ± 0.003109 | 0.131159 ± 0.003420 | 70.6/103.6 |
| heat | Tet10 | 2 | 0.082703 | 0.141465 ± 0.003056 | 0.139266 ± 0.001952 | 65.4/98.5 |
| heat | Tet10 | 4 | 0.021164 | 1.113278 ± 0.003938 | 1.119828 ± 0.006320 | 71.7/104.5 |
| heat | Tet10 | 6 | 0.009450 | 3.737478 ± 0.050003 | 3.812684 ± 0.041557 | 92.6/123.7 |
| bending | Tet4 | 2 | 0.489192 | 0.008788 ± 0.000314 | 0.009713 ± 0.000138 | 65.4/98.3 |
| bending | Tet4 | 4 | 0.247417 | 0.056110 ± 0.000196 | 0.056631 ± 0.000709 | 68.8/104.9 |
| bending | Tet4 | 6 | 0.165290 | 0.189921 ± 0.003439 | 0.192631 ± 0.003482 | 91.7/125.8 |
| bending | Tet10 | 2 | 0.049670 | 0.337928 ± 0.005267 | 0.339701 ± 0.009279 | 67.6/102.8 |
| bending | Tet10 | 4 | 0.012417 | 2.663061 ± 0.021537 | 2.682314 ± 0.004847 | 112.4/146.9 |
| bending | Tet10 | 6 | 0.005519 | 8.982302 ± 0.013708 | 9.172280 ± 0.162828 | 247.3/283.2 |

### Gleiche Genauigkeitsgrenze

Für beide Backends gilt dieselbe Auswahl:

| H1-Ziel | Wärme Tet4/Tet10 N | Biegung Tet4/Tet10 N |
|---|---|---|
| 50 % | 2 / 2 | 2 / 2 |
| 25 % | 6 / 2 | 4 / 2 |
| 12,5 % | nicht erreicht / 2 | nicht erreicht / 2 |

Bei 50 % und 25 % ist Tet4 im untersuchten Raster schneller. Bei 12,5 %
erfüllt nur Tet10 das Ziel; ein Tet4-Speedup-Verhältnis wäre daher unbelegt.
Tet10 unterschreitet die Grenzen teilweise deutlich: gleicher zulässiger Fehler
bedeutet hier keine künstlich angeglichene tatsächliche Abweichung.

### Gleiches Kostenbudget

Alle gemessenen Konfigurationen bleiben unter 512 MiB. Für beide Backends
ergibt sich bei 1/5/20 s dieselbe Auswahl:

| Fall/Element | 1 s: N / H1 | 5 s: N / H1 | 20 s: N / H1 |
|---|---|---|---|
| Wärme Tet4 | 6 / 0,16947 | 6 / 0,16947 | 6 / 0,16947 |
| Wärme Tet10 | 2 / 0,08270 | 6 / 0,00945 | 6 / 0,00945 |
| Biegung Tet4 | 6 / 0,16529 | 6 / 0,16529 | 6 / 0,16529 |
| Biegung Tet10 | 2 / 0,04967 | 4 / 0,01242 | 6 / 0,00552 |

Tet10 liefert innerhalb dieser Budgets die kleineren Fehler. Das erschöpfte
Netzraster begrenzt die Aussage: insbesondere Tet4 könnte das Budget mit
größeren, hier nicht gemessenen Netzen weiter ausnutzen.

### Folgerung und Grenzen

Zwischen N=4 und N=6 ist die beobachtete H1-Konvergenzordnung ungefähr 1
für Tet4 und 2 für Tet10, passend zu glatten Lösungen und den Ansatzräumen.
Die Kosten wachsen bei Tet10 deutlich; grobe Genauigkeitsziele rechtfertigen
diesen Aufwand nicht automatisch. SciPy und PETSc liefern innerhalb der
festgelegten Toleranzen dieselben Felder. PETSc zeigt hier keinen allgemeinen
Laufzeitvorteil und benötigt einschließlich Laufzeitbibliotheken mehr RSS.
Keine Aussage über MPI, iterative Vorkonditionierung, FieldSplit, große Netze,
Spitzenlasten oder allgemeine Hardware-Skalierung. Drei Wiederholungen belegen
Streuung dieser Läufe, keine statistische Überlegenheit eines Backends.

## Verifikation

TDD: Der neue Norm-/Budget-/Worker-Test scheiterte zuerst am fehlenden
Benchmarkmodul; nach Implementierung bestanden alle sieben gezielten Tests.
Native Tests vergleichen Wärme und Biegung für Tet4/Tet10 mit SciPy und PETSc,
auch im komplexen PETSc-Build. Laufzeiten sind ausdrücklich keine CI-Gates.
Die Rohdaten belegen zusätzlich alle Fehler-, Residuen-, Bilanz- und
Paritätsprüfungen der 72 Messläufe. Dies ist Code-Verifikation, keine
experimentelle Validierung.
