# Experimenteller Schaltungskern: Verifikationsnachweis

Stand 01.10.2026, #114–#119. Eigenständiger linearer RLC-Kern auf
`codex/116-experimental-spice-core`. Implementiert und lokal gegen
analytische Referenzen verifiziert; keine experimentelle Validierung,
keine M13-Freigabe und keine implementierte FEM-Kopplung.

| Issue | Umsetzung und überprüfbarer Nachweis |
| --- | --- |
| #114 | [Mathematik, Portvorzeichen und Kopplungsentscheidung](66_CIRCUIT_ANALYSES_AND_PORTS.md); unabhängige Piezo- und Wärme-Referenzparameter vor der späteren Kopplung festgelegt. |
| #115 | Unveränderliche Komponenten/Ports, begrenzte Netzliste, striktes TOML-Schema 1; unbekannte Schlüssel, Einheiten, doppelte Ports, unbekannte Knoten/Elemente und inkompatible Zuordnungen werden abgelehnt. |
| #116 | Sparse-DC-MNA, Spannungsteiler, Brücke und Skalierungsfälle; unabhängige KCL-/Leistungsprüfung, Topologiefehler, Q=Cv, Li und Speicherenergie. |
| #117 | Backward Euler, explizite Anfangswerte, deterministische Quellen und Ereignisse; RC/RL-Referenzen, passive RLC-Energie und zeitliche Konvergenz. |
| #118 | Komplexe Peakphasoren, RC/RL/RLC-Impedanzen, Wirk-/Blindleistung, DC-Grenzfall und Resonanzfehler; unveränderliche Ergebnisse und Provenienz. |
| #119 | Analyseübergreifende Referenzen, Ladungs-/Flussbilanz, Speicherenergie, drei Zeitraster, Ressourcen-/Wertebereichsfehler und reproduzierbarer RLC-Beispiellauf. |

Die Tests stehen in `tests/unit/test_circuit_*.py`,
`tests/convergence/test_circuit_time_convergence.py` und
`tests/regression/test_circuit_analysis_example.py`.
Insgesamt 125 gezielte Schaltungstests, einschließlich Beispielregression.

## Unabhängige Referenzen und Toleranzen

RC/RL mit Zeitkonstante 1 s werden gegen `1-exp(-t)` verglichen.
Die Raster 40/80/160 Schritte auf [0,2] s zeigen erste Ordnung
(beobachtete Ordnung zwischen 0.9 und 1.1); der feinste Maximalfehler
liegt unter 0.003 V beziehungsweise A. Speicherenergie ist unabhängig
`C*v²/2` beziehungsweise `L*i²/2`. Die passive RLC-Referenz prüft
Energieabnahme und Konvergenz auf drei Rastern.

Der [Beispiellauf](../examples/circuit_analyses/README.md) verwendet
R=1 Ohm, L=1 H, C=1 F. Die analytische Sprungantwort hat bei dt=0.01 s
ein vorher festgelegtes Fehlerbudget 0.004 V. Bei omega=1 rad/s hebt sich
die Reaktanz auf: I=1 A, V_C=-i V, P_R=0.5 W für 1 V Peak.
Absolute Phasor-/Leistungsgrenzen 1e-12 und Energiegrenzen 1e-12 J
prüfen Rundung in diesen kleinen float64-Referenzen; Zeitfehler werden
separat geprüft. Tests setzen relative/absolute Toleranzen ausdrücklich.

## Beobachteter Red-/Green-Nachweis

- 15 Konfigurations-/Porttests scheiterten vor dem fehlenden Loader;
  nach Implementierung bestanden sie.
- 21 AC- und 19 Transiententests scheiterten an den fehlenden Analyse-APIs;
  nach Implementierung bestanden sie.
- Der DC-Speichertest scheiterte an fehlenden Ergebnisgrößen und wurde grün.
- Drei zusätzliche Regressionen reproduzierten Speicherenergie-Überlauf,
  ungültige thermische Komponentenlisten und zu große AC-Ganzzahlquellen;
  nach Fehlerkorrektur liefern sie kontrollierte Validierungsfehler.
- Der historische erste DC-/Netzlisten-Nachweis steht in
  [Doc/65](65_EXPERIMENTAL_CIRCUIT_CORE.md).

## Ausgeführte Prüfungen

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

Gesamtsuite: **821 bestanden, 83 übersprungen**, 95.51 s.
Ruff, Formatprüfung und mypy bestanden. Übersprungene optionale Prüfungen
gelten nicht als bestanden. Der Beispieltest führt das Beispiel mit neuem
Ausgabeziel aus und prüft CSV/JSON einschließlich Provenienz.

## Grenzen und Integrationsstatus

Keine nichtlinearen Bauteile, adaptiven Schritte oder allgemeinen
höherindexigen DAE-Topologien. AC lehnt singuläre beziehungsweise gemäß
konfigurierbarer skalierter Pivotgrenze unauflösbare Resonanzen ab.
Die Speichergrenze betrifft numerische Ergebnisdaten, nicht den gesamten
Prozessspeicher. Portmetadaten erzeugen keine FEM-Lasten; gekoppelte
Abnahmetests bleiben #120–#128 vorbehalten.

Der lokale Nachweis ersetzt kein unabhängiges Review und keine GitHub-CI.
Ein unabhängiger Reviewer steht in dieser Sitzung nicht zur Verfügung.
Issues bleiben bis zur geprüften Integration in `codex/milestone-m13` offen.
Der Release-Umfang von M13 ist damit nicht abgeschlossen.
