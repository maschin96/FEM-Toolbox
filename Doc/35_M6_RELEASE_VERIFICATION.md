# M6 / v1.1.0 — Verifikation, Migration und Grenzen

Release-Gate: [#33](https://github.com/maschin96/FEM-Toolbox/issues/33).
Die Pflichtissues #29–#31 sind in den M6-Branch integriert. #32 bleibt auf
Nutzerwunsch zurückgestellt (ADR-021), kein Entropie-/Thermoelastic-Damping-Gate.
Die finale CI, der `main`-Commit und der annotierte Tag werden im Release-Issue
festgehalten. Fehlende unabhängige Reviewer werden ausdrücklich ausgewiesen;
Integration durch den Eigentümer ersetzt keine behauptete Fachbegutachtung.

## Umfang und unveränderte Verträge

| Bereich | Freigegebene Ergänzung | Nachweis |
|---|---|---|
| Material/Zeit | Positive rho/c-Daten, konsistente/gelumpte Tet4-Speicherung, konstante quellenbelegte thermische Presets, Step/Ramp/CSV | [Speicheroperatoren](31_TET4_STORAGE.md), [Materialbibliothek](32_MATERIAL_LIBRARY.md) |
| Wärme | Backward Euler, isolierte Körper, Quellen/Flüsse/Robin/Dirichlet, Speicher-/Randbilanz | [Dokument 33](33_TRANSIENT_HEAT.md), Unit-/Raum-/Zeitkonvergenztests |
| Piezo | Mechanische Masse; algebraisches Potential; Open-/Short-Circuit, Moden, Newmark, komplexe Harmonik | [Dokument 34](34_PIEZO_DYNAMICS.md), skalare und Stabreferenzen |
| Model/Export | Schemata 2/3, PVD/VTU, Messpunkte, komplexe Komponenten, Provenienz und Abbruchkennzeichnung | Model-/CSV-Tests und unabhängiger XML-/meshio-/CSV-Roundtrip |
| Rückwärtskompatibilität | Stationäre Python-API, Model-v1 und Ergebnisbündel-v1 | Unveränderte alte Python-/TOML-Beispiele und bestehende API-/Regressionstests |

Neue Rechnungen nennen Paketversion `1.1.0` in ihrer Provenienz. Alte
Schema-1-Dokumente benötigen **keine** Formatmigration. Transiente Wärme wählt
explizit `schema_version=2, physics="transient_heat"`. Piezo-Moden, Zeit- und
Frequenzanalyse wählen Schema 3 mit `base` als Schema-1-Piezomodell.
Schema 3 liefert `ModelResult.field_result` als netz-/materialgebundenes
`PiezoDynamicResult`; dessen `solution` enthält Moden, Zeithistorie oder
harmonische Einzelzustände. Die stationären Ergebnistypen ändern sich nicht.

Für neue Modelle jeweils die passenden Serienexporter verwenden; das
stationäre `export_result_bundle` wird nicht in ein Zeitformat umgedeutet.
CSV-Import legt Einheiten, Interpolation und Extrapolation offen fest. Gemessene
Daten verleihen einer Rechnung allein noch keine experimentelle Validierung.

## Wissenschaftliche Abnahme

Die Testtoleranzen stammen aus den jeweiligen analytischen Referenzen und
wurden für das Release nicht erweitert. Kleine algebraische Residuen ersetzen
keine Raum-/Zeitverfeinerung. Alle folgenden Parameter sind synthetische
Verifikationsfälle, keine gemessenen Werkstoffdaten.

**Wärme:** Gleichmäßige Abkühlung des regulären Tetraeders mit Rate 1/s bis
`t=1 s`, 10/20/40/80 Schritte: Fehler 0,176638 / 0,090100 / 0,045512 /
0,022873 K, Ordnungen 0,971 / 0,985 / 0,993 (Fenster 0,9–1,1).
Separater Raumtest auf `n=2,4,6,8` nutzt eine quadratische Raumlösung mit
linearer Zeitabhängigkeit, für die Backward Euler zeitlich exakt ist:
L2-Ordnung 2, H1-Ordnung 1 (Fenster 1,8–2,2 bzw. 0,85–1,15).
Speicher-/Randbilanzen, Mehrmaterialfall und stationärer Grenzfall sind
zusätzlich geprüft. SI-Toleranzen und Quellen: Dokument 33.

**Moden:** Leitfähiger d33-Stab nach Dokument 34, `L=0,02 m`,
`C=1e9 Pa`, `e=5 C/m²`, `kappa=1e-8 F/m`, `rho=7500 kg/m³`.
Offene Referenz 8539,125638 Hz; Kurzschlussreferenz 4891,834344 Hz.
Dies sind verschiedene elektrische Randbedingungen, kein zeitlicher Rastereffekt.

| Unterteilungen n | Offen [Hz] | relativer Raumfehler | Kurzschluss [Hz] | relativer Raumfehler |
|---:|---:|---:|---:|---:|
| 1 | 8864,686086 | 3,81257 % | 4961,508097 | 1,42429 % |
| 2 | 8682,232035 | 1,67589 % | 4919,014028 | 0,55561 % |
| 4 | 8586,099845 | 0,55011 % | 4900,631000 | 0,17982 % |

Beide Fehler fallen und liegen auf dem feinsten Netz unter der vorgegebenen
2-%-Grenze. Maximales beobachtetes feldweises Residuum dieser Läufe < 1e-15;
Solvergrenze 1e-9. Gesamtladung der offenen Elektrode schließt innerhalb der
in Dokument 34 begründeten Kapazitäts-/Spannungsskala.

**Newmark:** Ungetriebener Ein-DOF-Oszillator aus
`tests/unit/test_piezo_newmark.py`, feste Raumdiskretisierung, `t=1 s`:

| Schritte | Verschiebungsfehler [m] | maximale relative Energiedrift |
|---:|---:|---:|
| 40 | 3,981798e-5 | 8,89e-16 |
| 80 | 9,973854e-6 | 8,89e-16 |
| 160 | 2,494675e-6 | 2,89e-15 |

Fehlerquotienten ≈4; gefordertes Fenster 3,8–4,2. Energiegrenze 3e-12 relativ;
maximale beobachtete Feldresiduen <2,3e-16. Quadratische Lagerbewegung prüft
getrennt die vorgeschriebenen Ableitungen und `M_fc*a_c`. Ereignistests prüfen,
dass Last-/Spannungssprünge keine vorzeitige Bewegung erzeugen.

**Harmonik:** Explizit `Re(x_hat exp(+i omega t))`, Spitzenwerte, keine Dämpfung.
DC-Grenze, komplexe Kraft/Spannung, Lagerträgheit und Resonanzdiagnose sind
geprüft. Bei `omega=1,3 rad/s` und passenden harmonischen Anfangswerten sind
die maximalen Zeit-/Frequenzabweichungen 7,564716e-8 m (100 Schritte) und
1,891609e-8 m (200 Schritte). Raumdiskretisierung und Anregungsfrequenz bleiben
identisch; der Quotient ≈4 misst allein den Zeitfehler. Ein Frequenzsweep
findet keine zwischen seinen Stützstellen liegenden Resonanzen automatisch.

Reproduktion: `python -m pytest tests/convergence tests/unit/test_piezo_newmark.py
tests/unit/test_piezo_harmonic.py`. Die Zahlen wurden am 23.09.2026 auf
macOS/ARM64, CPython 3.13.7, NumPy 2.4.6, SciPy 1.18.1 ermittelt. Rundungswerte
können plattformabhängig variieren; maßgeblich sind die begründeten Testgrenzen.

## Installation und Release-Gates

Die vollständige Suite enthält 396 Tests, einschließlich aller vorhandenen
stationären Regressionen und Konvergenztests. Ruff, Formatcheck, mypy und
`git diff --check` gehören zum Gate. Der Versionstest wurde zunächst gegen
`1.1.0.dev0` rot und nach Installation von `1.1.0` grün ausgeführt.

Die CI hat 17 Jobs: drei vollständige Python-Testsuiten, sechs Kombinationen
minimaler/aktueller Abhängigkeiten, Qualitätsprüfung, Artefaktbau und sechs
Wheel-/sdist-Installationen. Die Artefaktjobs haben keinen Quellbaum und laufen
mit `python -I`. Sie führen unveränderte stationäre Beispiele (`steady_heat`,
`elastic_tension`, `parallel_plate`, `thermal_expansion`, `piezo_d33`), das
stationäre TOML-Modell sowie transiente Wärme und alle drei Piezo-Analysearten
aus. Der stationäre TOML-Fall wird zusätzlich gegen seine analytischen
Verschiebungs-/Potential-/Flussreferenzen geprüft.

Paket und Archive führen MIT als Lizenz, CPython `>=3.11,<3.14` sowie
NumPy/SciPy/meshio als externe Abhängigkeiten. Es werden keine neuen
Laufzeitabhängigkeiten oder fremden Bibliotheken gebündelt. Wheel, sdist und
SHA-256-Liste werden am GitHub Release veröffentlicht; kein PyPI-Upload.

## Grenzen

Konstante lineare Stoffdaten, kleine Verformungen und gerade Tet4 bleiben
Voraussetzung. Keine thermomechanische Rückwirkung, Temperaturabhängigkeit,
Phasenwechsel, Materialhysterese, Dämpfung im Piezo-Pfad oder elektrischen
Trägheitsterme. Keine freie Starrkörper-Modalanalyse: Lager/Gauge sind explizit.
Die dichte Piezo-Referenz ist auf 2048 aktive Gesamt-DOFs begrenzt; Zeitdaten
auf 2.000.000 skalare Zustandsplätze. Wärme begrenzt die geschätzte numerische
Historie separat (Standard 256 MiB); das ist keine Prozess-RAM-Garantie.

Newmark-Lagerverschiebung und -geschwindigkeit müssen kontinuierlich sein;
allgemeine Callbacks müssen konsistente Ableitungen und Ereignisse liefern.
Für angeregte Mechanik wird keine unberechnete Arbeitsbilanz behauptet.
Wärme mit konsistenter Kapazität garantiert keine diskrete Monotonie.
Parameterpresets sind typische Werte mit Quellen, keine universellen Stoffgesetze.
Keine experimentelle Validierung, PETSc-/MPI-Skalierungs- oder Windows-CI-Zusage.
