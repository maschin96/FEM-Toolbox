# M11: numerische Nachweise und offenes Release-Gate (#83)

19.10.2026 · 1.6.0.dev1 · keine finale Veröffentlichung.

## Ergebnis der vorab definierten Matrix

Die [Gates aus #76](83_M11_TOPOLOGY_CONTRACT.md) wurden für drei Netze und
zwei Radien ausgeführt. Alle sechs Fälle konvergieren innerhalb des unveränderten
100-Iterations-/300-Auswertungsbudgets, halten das Dichtevolumen 0.6 ein und
reduzieren C um mindestens 1%. Je Fall besteht mindestens ein explizit
gekennzeichneter Geometrieweg. Messdatei: [vollständige Matrix](evidence/m11/verification.json).

| Tet4 | r [m] | C_start [J] | C_final [J] | Auswertungen | Graumaß | Geometrievolumen* |
|---:|---:|---:|---:|---:|---:|---:|
| 96 | 0.3 | 0.11127684 | 0.04067543 | 17 | 0.2661 | 0.6875 |
| 96 | 0.6 | 0.11758127 | 0.06868462 | 17 | 0.6126 | 0.6771 |
| 324 | 0.3 | 0.13580981 | 0.05957262 | 22 | 0.5027 | 0.6265 |
| 324 | 0.6 | 0.13253531 | 0.08951997 | 28 | 0.8116 | 0.7901 |
| 768 | 0.3 | 0.14666601 | 0.06628759 | 29 | 0.5857 | 0.7057 |
| 768 | 0.6 | 0.14364792 | 0.09605180 | 29 | 0.8216 | 0.7786 |

*Die kleinste akzeptierte Variante aus den drei reinen und drei explizit
verbundenen Schwellen. Das zusätzliche Material ersetzt keine Einhaltung des
0.6-Dichtebudgets. Reine getrennte Schwellen und zu große Geometrien bleiben
als No-go in der Messdatei enthalten. Insbesondere r=0.6 benötigt bei feineren
Netzen einen dokumentierten Verbindungsweg bei Schwelle 0.45.

Die Wiederholung des vollständigen kleinen Mehrlastfalls liefert relative
C-Abweichung 0.0. Die separate Einzellastvariante konvergiert und erfüllt die
selben C-/Volumengates. Der unabhängige Richtungsgradient bei h=1e-5 hat
relative Abweichung 1.662e-8; die Volumenableitung weicht absolut 5.556e-12 ab.
Weitere vollständige komponentenweise Schrittweitenstudien für fast leere,
mittlere und fast volle Felder, feste Gebiete und verzerrte Netze sind Unit-Tests.
Ein separater Filtertest integriert die Volumina per Determinante und berechnet
alle Filtergewichte unabhängig von der Implementierung.

Zustandsresiduen sind in der Matrix kleiner als 1.2e-16, Kräftebilanzen kleiner
als 6e-14 N und relative Arbeit-/Energieabweichungen kleiner als 5e-14.
Geometrien besitzen eigene Solver-/Kraft-/Momentennachweise und rohe kritische
von-Mises-Zellen. Keine Festigkeitsgrenze und keine Spannungsnetzakzeptanz
werden aus diesen Zellspitzen abgeleitet.

## Was die Netz-/Filterstudie zeigt

Die Parameter sind identisch; die an Last und Klemme festgehaltenen
flächennachbarten Zellen werden bei jedem Netz neu diskretisiert. Die
Passivregion ist damit netzdefiniert und ihre physische Ausdehnung ändert sich.
Die Tabelle ist ein Nachweis des kombinierten Diskretisierungs-/Filtereinflusses,
keine reine h-Konvergenz zu einem einzigen kontinuierlichen Entwurfsgebiet.
Größerer Radius lässt deutlich mehr Zwischenmaterial zurück. Unterschiede der
lokalen Topologie oder von C beweisen keinen Fehler und keine globale
Optimalität. Eine garantierte Mindestwandstärke ist weiterhin nicht belegt.

## Fehlgeschlagene Erstläufe bleiben sichtbar

Der [erste Lauf ohne Variablenskalierung](evidence/m11/unscaled_no_go.json)
erreichte auf vier feineren Fällen das Iterationslimit. Dieser Stand stammt
vom #81-Arbeitsstand vor Koordinatenskalierung und Geometrieergänzung; seine
Paketkennung war noch 1.5.0, keine M11-Veröffentlichung. Er ist keine positive
M11-Abnahme. Die vollständigen lokalen Rohausgaben beider Läufe liegen unter
`results/m11/`; die versionierten Messungen enthalten die Ergebnis-Hashes.
Die Änderung z=x/sqrt(n_free) konditioniert das unveränderte Problem. Weder
ftol noch Budget oder Akzeptanzgrenzen wurden gelockert. Verschiedene lokale
Minima beim nichtkonvexen Problem werden nicht als identische Optima ausgegeben.

## Reproduktion und Herkunft

```sh
python -I examples/topology/verify.py --output results/topology-verification
```

Der Lauf verwendet CPython 3.13.5, NumPy 2.4.6 und SciPy 1.18.1 unter macOS ARM64.
Die Zustandsrechnungen bleiben im vorhandenen SciPy-Mechaniksolver.
`OPENBLAS_NUM_THREADS=1` wurde für die Messmatrix gesetzt. Laufkosten stehen
pro Fall in Sekunden und Designauswertungen; sie umfassen Prozessstart,
Optimierung, Export und Geometrieprüfung und sind kein allgemeiner Speedup.
Quellhashes und Nachweisdateihashes stehen in [environment.json](evidence/m11/environment.json).
Historische M9-V&V-Manifeste bleiben unverändert; M11 hat sein eigenes Nachweisarchiv.

## Freigabezustand

- #76–#82: implementierte und lokal numerisch geprüfte Arbeitspakete im M11-Branch.
- Lokale Gesamtsuite: **868 bestanden, 83 optionale Tests übersprungen**.
  Ruff und Formatprüfung bestanden; striktes mypy über 74 Quelldateien bestanden.
- Wheel und sdist separat frisch installiert: jeweils korrekter isolierter Import,
  Version 1.6.0.dev1, NumPy 2.5.4/SciPy 1.18.1, gleicher C-Wert und erfolgreiche
  Optimierungs-/Geometriegates. [Artefaktbelege](evidence/m11/artifacts.json).
- Fünf Handbücher, drei negative Dokumentfälle und Quellenarchiv-Nachbau bestanden.
  Neue Theorie-/Anwendungs-/API-/V&V-Seiten und Titel visuell geprüft; vorhandene
  Layoutwarnungen werden nicht zu einer pauschalen historischen Neuabnahme erklärt.
- Aktuelle PR-/Integrations-CI bleibt am tatsächlich gepushten Kandidaten zu prüfen.
- Unabhängiges menschliches Fachreview liegt nicht vor; unabhängige
  Finitedifferenzen/analytische Invarianten sind keine menschliche Freigabe.
- Integration nach main, grüne Zielcommit-CI, finale Paketversion, annotierter
  Tag und veröffentlichte Artefakte sind noch offen. #83 bleibt bis dahin offen.

Kein Tag, kein Release und kein automatisches Schließen von Issues wird aus
einem lokalen positiven Lauf abgeleitet. Keine experimentelle oder industrielle
Bauteilqualifikation, kein CAD-/Fertigungs- oder universeller Skalierungsnachweis.
