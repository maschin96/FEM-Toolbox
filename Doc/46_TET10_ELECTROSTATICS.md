# Tet10-Elektrostatik durch Model, Solver und Export (#35d)

**Aktueller Gesamtstand:** [Tet10-Abnahme](55_TET10_ACCEPTANCE.md).
Die folgenden Status- und Restumfangsangaben dokumentieren den jeweiligen
Implementierungsschnitt; später integrierte Analysen und Vergleiche sind im
Gesamtbericht zusammengeführt.

Dieser Integrationsschnitt erweitert die bereits vorhandene Tet10-Wärme um
stationäre Elektrostatik. Statische Mechanik ist inzwischen in
[Dokument 47](47_TET10_ELASTICITY.md) integriert. #35 bleibt für Piezo, Kopplungen und die
verbleibenden Vergleichsgates offen. Keine M7-Releasefreigabe oder experimentelle
Validierungsbehauptung.

Der nachfolgende [stationäre Piezoschnitt](49_TET10_PIEZOELECTRIC.md)
ergänzt inzwischen Tet10-Piezo; die obige Abgrenzung beschreibt den damaligen Stand.

## Umfang und Konventionen

`solve_electrostatics`, `electric_dofs` und Model-Schema 1 unterstützen
homogene geradkantige Tet10-Netze mit mittigen Kantenknoten. Zehn Knoten tragen
Potential-DOFs; benannte Triangle6-Dirichletflächen erfassen alle sechs Knoten.
Materialregionen, anisotrope absolute Permittivität, Volumenladung und äußere
D-Flüsse verwenden die bisherigen SI-Verträge. Gemischte Elementfamilien und
gekrümmte Geometrien werden weiterhin abgelehnt.

Die Gleichung bleibt `-K_phi phi = surface(D.n) - volume(rho)`. `E=-grad(phi)`,
`D=kappa E` und die Reaktionen folgen dem ursprünglichen negativen elektrischen
System. Eine Reaktion ist der unbekannte äußere D-Fluss an einem festgehaltenen
DOF; die Ladung des angrenzenden Leiters ist ihr Negatives. Jede zusammenhängende
Komponente benötigt eine explizite Potentialreferenz. Model-Preflight und direkter
Solver prüfen dies einschließlich der Kantenknoten; es gibt kein automatisches
Pinning. Floating-Elektroden bleiben beim piezoelektrischen Solver und sind kein
neuer Bestandteil des eigenständigen Elektrostatikpfads.

Der interne Dispatch verwendet die verifizierten [lokalen Operatoren](42_TET10_OPERATORS.md).
Tet10-Volumenregeln haben Grad 5, Triangle6-Flächenregeln Grad 4: bei quadratischem
N sind Raumladungen bis Grad 3 und Randladungen bis Grad 2 exakt integriert.
Komplexere räumliche Lasten benötigen einen eigenen Quadraturkonvergenznachweis.
Die vorhandene Backend-Grenze der direkten Python-API unterstützt SciPy und
optional PETSc. Model erhält keinen zusätzlichen Backend-Konfigurationsschalter.

## Feldauswertung und Energie

`ElectrostaticResult` behält seine bisherigen Arrays und ergänzt optional
`field_points` mit physischen Auswerteorten in Metern:

- Tet4: `None`, E/D sind elementweise konstant.
- Tet10: unveränderliche `(ncells, 3)`-Schwerpunktkoordinaten; E/D sind rohe
  Auswertungen an diesen Punkten, keine Knotenprojektion und kein konstantes Feld.
- `element_energy` bleibt ein `(ncells,)`-Array positiver gespeicherter Energie
  in Joule. Für Tet10 wird `1/2 integral grad(phi).kappa.grad(phi)` mit den
  punktabhängigen Gradienten über jede Zelle integriert. Ein Schwerpunktwert
  allein liefert dafür keine exakte Energie.

Konstante Potentialoffsets werden vor der Ableitung subtrahiert. Die Gesamtenergie
ist die Summe der Zellintegrale; elektrische Enthalpie und gespeicherte Energie
werden nicht verwechselt.

Das Ergebnisbündel verwendet für Tet10 weiterhin Schema 2, wie der Wärmepfad:
E-/D-Namen enden auf `_at_centroid`, mit `kind="raw"`, `location="cell_centroid"`
und `field_evaluation_points_m`. `electric_energy_J` hat **keinen**
Schwerpunktsuffix, sondern `kind="cell_integral"` und `location="cell_volume"`.
Nodale Felder bleiben unverändert. Fehlende oder zum Netz unpassende Auswerteorte
werden abgelehnt. Tet4 verwendet weiter Schema 1 und seine bisherigen Feldnamen.

## Nachweise und Toleranzen

Lokal: `PYTHONPATH=src .venv/bin/python -m pytest -q` mit 525 bestandenen und
47 optional übersprungenen Tests; Ruff, Formatprüfung, Mypy und
`git diff --check` bestanden. Native PETSc- und Artefaktprüfungen erfolgen in CI.

TDD: drei initiale Integrationstests scheiterten vor Umsetzung an fehlenden
Tet10-DOFs, Model-Sperre und Tet4-Beschränkung des Solvers. Nach Dispatch und
Exportanbindung sind sie grün. Das ausführbare Beispiel wurde ebenfalls durch
einen zunächst wegen fehlender Datei scheiternden Test eingeführt.

- Quadratischer Patch `phi=x²+y²/2+z²/4` mit anisotroper Permittivität und
  Raumladung -9e-11 C/m³: Potential, E/D, Reaktionen, Ladungsbilanz und analytische
  Energie aus Würfelmomenten `integral x_i²=1/3`, `integral x_i*x_j=1/4`.
- Triangle6-Randlast über Model, portable Dokumentreproduktion und unabhängiges
  VTU-Einlesen prüfen den vollständigen Modell- und Exportpfad.
- Plattenkondensator mit κ=2e-11 F/m und φ=x: W=1e-11 J und korrekte Leiter-/D-Fluss-
  Vorzeichen, auch nach lokaler Orientierungsumkehr und mit inaktivem Netzknoten.
- Materialgrenze κ=1e-11/8e-11 F/m: kontinuierliches Potential und normaler
  D-Fluss, achtfacher Gradientensprung und unabhängige Energie 2,8125e-12 J.
- Fehlende Referenzen, gekrümmte Kanten, fehlende/falsche Exportpunkte und
  unveränderliche Auswertekoordinaten werden geprüft.
- Native PETSc-Paritätstests laufen in den vorhandenen real/complex-CI-Jobs,
  jeweils direkt und mit FieldSplit. Ohne PETSc werden sie lokal übersprungen.
- Das [Gmsh-/TOML-Beispiel](../examples/tet10_electrostatics/README.md) wird auch
  aus frischen Wheel-/sdist-Installationen auf Python 3.11–3.13 ausgeführt.

Für kleine Referenzfälle: Potentialfehler bis 3e-13 V, E bis 2e-12 V/m,
D bis 1e-22 C/m², Ladungsbilanz bis 2e-23 C und analytische Energie bis 2e-24 J.
Die Toleranzen berücksichtigen Double-Rundung und Summation bei Permittivitäten
der Größenordnung 1e-11 F/m. Backend-Parität ist separat mit 1e-10 V,
1e-9 V/m und 1e-21 J/C begrenzt. Diese Grenzen sind keine allgemeinen
Diskretisierungsfehlerschranken.

Vier-Netze-Manufactured-Solution: `phi=x⁴+y⁴/2+z⁴/4`, κ=2e-11 F/m und
ρ=-24e-11 (x²+y²/2+z²/4). Die Fehlernormen verwenden unabhängig formulierte
baryzentrische P2-Basis/Ableitungen und eine überintegrierte Duffy-Produktregel
mit sechs Gauss-Punkten pro Richtung. Auch jedes Energieintegral wird gegen
diese unabhängige Auswertung geprüft (`rtol=2e-12`, `atol=1e-25 J`).

| Würfelteilung n | L2-Fehler | H1-Seminormfehler |
|---|---:|---:|
| 1 | 0,08030158 | 0,52014538 |
| 2 | 0,01114274 | 0,14324531 |
| 3 | 0,00334559 | 0,06476286 |
| 4 | 0,00141754 | 0,03665694 |

Die letzten beiden Intervalle liefern L2-Ordnungen 2,967/2,985 und H1-Ordnungen
1,958/1,978; vorgegebene Akzeptanzfenster sind 2,7–3,3 beziehungsweise 1,8–2,2.
