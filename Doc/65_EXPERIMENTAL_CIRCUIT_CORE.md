# Experimenteller Schaltungskern

Status: Arbeitsvertrag für den Experimental-Branch; keine M13-Freigabe.
Refs #114, #115, #116 und #119. Der erste Implementierungsschnitt ist ein
eigenständiger linearer DC-Referenzkern. Die vollständige Port-/Kopplungsentscheidung
aus #114 und die übrigen M13-Abnahmen bleiben offen.

## Umfang und Eingabe

`fem_toolbox.circuits` trennt unveränderliches Schaltungsmodell, Netzlistenadapter
und DC-Solver. Keine neue Abhängigkeit: NumPy/SciPy und bestehende Provenienz.
Experimentelle API: `Component`, `Circuit`, `DCResult`, `solve_dc`.

Jede Elementzeile enthält genau `Name Plusknoten Minusknoten Wert`.
Der erste Buchstabe des Namens bestimmt R, C, L, V oder I. R/C/L müssen positiv
sein; unabhängige Spannungs-/Stromquellen dürfen null oder negativ sein.
Namen beginnen mit dem Elementbuchstaben und enthalten nur ASCII-Buchstaben,
Ziffern und Unterstriche. Knoten enthalten diese Zeichen ebenfalls. Namen und
Knoten sind ohne Beachtung der Groß-/Kleinschreibung definiert; das Modell
normalisiert Elementnamen auf Großbuchstaben und Knoten auf Kleinbuchstaben.
Einzige Masse ist `0`; `gnd` wird nicht implizit umbenannt. Gleiche Anschlüsse,
doppelte Namen, leere Schaltungen und fehlende Masse sind Fehler.

Werte sind reelle Dezimalzahlen mit optionalem Exponenten und genau einem
optionalen Skalierungssuffix: `t`, `g`, `meg`, `k`, `m`, `u`, `n`, `p`, `f`.
Suffixe ignorieren Groß-/Kleinschreibung: **M bedeutet milli**, MEG bedeutet mega.
Ohm, F, H, V und A sind durch den Elementtyp festgelegt. Einheitensuffixe,
Ausdrücke und beliebige nachgestellte Zeichen werden abgelehnt.
Leerzeilen und volle Kommentarzeilen mit `*` sind erlaubt, `.end` ist optional.
Nach `.end` sind nur Leer-/Kommentarzeilen erlaubt. Keine Titelzeile, Inline-
Kommentare, `.include`, `.model`, `.param`, `.tran`, `.ac`, gesteuerten Quellen,
Halbleiter, Unterprogramme oder ausführbaren Ausdrücke. Kein beliebiger
SPICE-Dialektimport.

## DC-Gleichungen und Vorzeichen

Für jeden Zweig gilt `v = potential[plus] - potential[minus]` in V.
Der Zweigstrom `i` in A fließt von plus nach minus; `p = v*i` in W ist positive
aufgenommene Leistung (passive Vorzeichenkonvention). Eine liefernde Quelle
hat negative Leistung. Die Summe aller Zweigleistungen ist null.

R: `i = v/R`; I: vorgegebener Zweigstrom; V: vorgegebene Zweigspannung.
C ist im DC-Grenzfall offen (`i=0`), L ist kurzgeschlossen (`v=0`) mit eigenem
Stromfreiheitsgrad. Kapazitäts-/Induktivitätswerte bleiben im Modell und in
der Provenienz erhalten. Anfangswerte und gespeicherte Energie gehören noch
nicht zu diesem DC-Ergebnisvertrag.

Die modifizierte Knotenanalyse verwendet Knotenpotentiale ohne Masse und
Zweigströme der V-/L-Elemente:

```text
[G  B] [v] = [-A_I i_I]
[Bᵀ 0] [j]   [e_VL   ]
```

Jede Inzidenzspalte hat +1 am Plus- und -1 am Minusknoten; der Masseeintrag
entfällt. `G = A_R diag(1/R) A_Rᵀ`. Sortierte Namen erzeugen deterministische
DOF-Reihenfolgen. Sparse LU mit Zeilen-/Spaltenäquilibrierung dient als
Referenz. Kein CG auf dem indefiniten MNA-System.

DC-Preflight prüft vor der Faktorisierung Erdung über R/V/L (C/I verankern kein
Potential) und Zyklen idealer Spannungszweige. Widersprüchliche Zyklen haben
keine Lösung; konsistente redundante Zyklen haben unbestimmte ideale
Zweigströme und werden ebenfalls abgelehnt. Keine künstliche Leckleitfähigkeit.
Nicht darstellbare Matrizen oder Ergebnisse ergeben konkrete Fehler.

## Ergebnisse und Verifikation

Das Ergebnis enthält unveränderliche Zuordnungen für Potentiale einschließlich
Masse, Zweigspannungen, Zweigströme, aufgenommene Leistungen und KCL-Reste an
allen Knoten, außerdem Gesamtleistung und normiertes algebraisches Residuum.
`RunProvenance` hasht den kanonischen, sortierten SI-Schaltungsinput samt
Schaltungsschema und nennt Solver/NumPy/SciPy sowie Bilanzschwellen.

Handrechnungen: 12-V-Teiler mit 1 kOhm/2 kOhm ergibt 8 V und 4 mA, Quelle
-48 mW, Widerstände zusammen +48 mW. Eine 3-mA-Quelle von Masse zum Knoten
mit 2 kOhm gegen Masse ergibt 6 V und 18 mW Widerstandsleistung. Gemischte
Quellen, umgedrehte Anschlüsse, Brückennetz, DC-RLC sowie fehlerhafte Topologien
ergänzen die Referenzen. Toleranzen werden vorab in den Tests begründet.
Keine Raum-/Zeitkonvergenz für diesen algebraischen DC-Schnitt; sie wird für
Transienten und FEM-Kopplung in eigenen Issues erforderlich.

## Weitere Schritte

#117: Zeitintegration und Anfangswerte; #118: harmonische Analyse;
#119: vollständiger unabhängiger DC-/RC-/RL-/AC-Gesamtnachweis. #120 ff.:
FEM-/thermische Ports und Kopplung. Diese Entscheidungen und Referenzfälle
werden vor ihrer jeweiligen Implementierung ergänzt. Der DC-Arbeitsvertrag
schließt #114 oder #115 einschließlich Portzuordnung nicht vollständig ab.

## Nutzung und lokale Prüfungen

```python
from fem_toolbox.circuits import Circuit, solve_dc

circuit = Circuit.from_netlist("V1 in 0 12\nR1 in out 1k\nR2 out 0 2k")
result = solve_dc(circuit)
print(result.node_potentials_V["out"])  # 8 V
```

[Ausführbares Beispiel mit Netzliste](../examples/circuit_dc/README.md).
`solve_dc` prüft mit Standardschwellen `relative_tolerance=1e-10`,
`current_tolerance_A=1e-12`, `power_tolerance_W=1e-12` und
`voltage_tolerance_V=1e-12`. Absolute Schwellen gelten in der genannten
SI-Einheit; die relative Schwelle skaliert mit lokalem Stromdurchsatz,
absoluter Gesamtzweigleistung beziehungsweise Idealzweig-Potentialskala.
Andere Schwellen sind explizite Aufrufparameter und Teil der Provenienz.
Diese Standards sind keine experimentellen Abnahmeschwellen.

TDD-Nachweis vom 01.10.2026, CPython 3.13.7:

- Modell/Netzliste zuerst: 34 Tests scheiterten am fehlenden Schaltungsmodul;
  nach Implementierung bestanden dieselben 34 Tests.
- DC zuerst: 23 weitere Tests scheiterten an der fehlenden Solver-API;
  nach Implementierung bestanden alle 57 Tests.
- Numerische Randfälle: 3 zusätzliche Tests zeigten die Fehler bei
  Suffix-/Exponentenbereich, Quellenunterlauf und großer Quellenschleife.
  Nach Korrektur bestanden sämtliche 63 Schaltungstests einschließlich der
  weiteren Brücken-/Skalierungsreferenzen.

Befehle aus dem Repository-Root:

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q tests/unit/test_circuit_model.py tests/unit/test_circuit_dc.py
PYTHONPATH=src .venv/bin/python examples/circuit_dc/run.py
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

Die lokale Umgebung benötigt `PYTHONPATH=src`, weil keine editierbare
Paketinstallation vorhanden ist. Gezielte Tests, Beispiel und alle genannten
Qualitätsprüfungen bestanden. Gesamtsuite mit
`PYTHONPATH=src .venv/bin/python -m pytest -q`: **759 bestanden, 83 übersprungen**
in 94.65 s. Übersprungene optionale Prüfungen gelten nicht als bestanden.
Keine GitHub-CI, keine saubere Artefaktinstallation und kein unabhängiges
Review dieses Experimental-Standes behauptet. Kein Issue wird geschlossen.
