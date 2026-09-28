# Stationäre Tet10-Piezoelektrizität (#35d)

Dieser Integrationsschnitt von #35 verbindet die vorhandenen quadratischen
Elementoperatoren mit `solve_piezoelectric`, Model/Preflight, linearen Backends
und Ergebnisexport. Unterstützt sind geradkantige Tet10 mit Triangle6-Rändern,
konstante anisotrope Materialdaten pro Zelle sowie feste und schwebende Elektroden.
Thermische Piezo-Eigendehnungen, Thermopiezo und Piezo-Dynamik bleiben für Tet10
explizit abgelehnt. Die bestehende Tet4-Thermopiezofunktion bleibt erhalten.

## Diskretisierung und API

`piezoelectric_dofs(mesh)` erzeugt feldweise u- und φ-Freiheitsgrade auf allen
aktiven Zellknoten. Ungenutzte Netzpunkte erhalten keine Freiheitsgrade.
Die symmetrisch-indefinite Blockmatrix behält die Konvention
`[Kuu Kuphi; Kuphi.T -Kphi]` mit E=-grad(φ), σ=Cε-eᵀE und D=eε+κE.
Volumen-/Flächenlasten verwenden die quadratischen Integrationsroutinen.
Elektroden fassen alle sechs Triangle6-Knoten zusammen; das vorhandene
Reduktionsverfahren und seine Konflikt-/Gauge-Prüfungen bleiben maßgeblich.
Materialregionen, Blockskalierung und der optionale `backend` sind nutzbar.

`PiezoelectricResult.field_points` ist für Tet10 ein unveränderliches
(ncells,3)-Array der Zellschwerpunkte in Metern, für Tet4 weiterhin `None`.
Dehnung, Spannung, E und D sind rohe Schwerpunktwerte. Mechanische Energie,
Kopplungsenergie und negative dielektrische Enthalpie sind dagegen Integrale:

```text
Wm = integral 1/2 epsilon^T C epsilon dV
Wc = integral -epsilon^T e^T E dV
He = integral -1/2 E^T kappa E dV
```

Die bestehende Grad-5-Regel integriert die quadratischen Energiedichten exakt.
Translation und Potentialoffset werden vor der Gradientenbildung subtrahiert.
Die Gesamtenthalpie ist kein positiver Energiemaßstab. Schema 2 kennzeichnet
alle drei Energiefelder als `cell_integral`; andere Zellfelder tragen den Suffix
`_at_centroid`. Der Export prüft die Auswertekoordinaten gegen das Netz.
Tet4 verwendet unverändert Schema 1. Die Provenienz-Hashgrundlage enthält die
aktive Zellfamilie und die zugehörige Konnektivität.

## Analytische Nachweise und TDD

Die ersten drei neuen Tests scheiterten vor der Implementierung an fehlenden
Tet10-Freiheitsgraden. Der Model-Test scheiterte anschließend an der bisherigen
Tet10-Piezosperre; der Beispieltest an der noch fehlenden Datei. Nach Integration
bestehen die jeweiligen Tests. Vorhandene thermische und dynamische Sperren
werden separat geprüft; nur die überholte statische Piezosperre wurde ersetzt.

- Quadratischer gekoppelter Patch: u=(0,01x²,0,0) m, φ=x² V im Einheitswürfel,
  E=100 Pa, ν=0, eₓₓₓ=2 C/m², κ=3 I F/m. Damit σxx=6x Pa und Dₓ=-5,96x C/m²,
  b=(-6,0,0) N/m³, ρ=-5,96 C/m³. Links u=φ=0, rechts t=(6,0,0) Pa und
  D·n=-5,96 C/m²; andere Flächen lastfrei. Wm=0,02/3 J, Wc=0,08/3 J, He=-2 J.
- Feste bzw. schwebende Elektrode: φ=x, uₓ=-0,02x, σ=0, Qrechts=3,04 C,
  Gesamtenthalpie -1,52 J; alle Mittelknoten werden mitgeprüft.
- Anisotroper affiner Patch mit voller Kopplungsmatrix, Potentialoffset,
  Translation, umgekehrter Zell-/Flächenorientierung und unbenutztem Punkt.
- Zwei Materialregionen mit um Faktor 2 skalierten C/e/κ: stückweise lineares
  Potential, stetige Traktion und D·n, auch ohne Blockskalierung.
- Punktlast-Reziprozität, Bilanzen, Materiallisten, fehlende Referenz/Lagerung,
  gekrümmte Geometrie und Konflikt Elektroden-/Flächenladung.
- Model-Rundreise und gemeinsamer Export; ausführbares
  [Gmsh/TOML-Beispiel](../examples/tet10_piezoelectric/README.md).
- Native PETSc-Parität (direkt/FieldSplit, real/complex) für den quadratischen
  Patch und eine schwebende Elektrode; Wheel-/sdist-CI führt das Beispiel aus.

Toleranzen sind absolute Rundungsgrenzen für diese kleinen SI-Testfälle:
u 2e-13 m, φ 2e-12 V, σ/E/D und Bilanzen 2e-11 in der jeweiligen SI-Einheit,
Enthalpien 2e-12 J. Beim vollen anisotropen Spannungstensor gilt 2e-10 Pa wegen
der größeren Steifigkeit und Translationssubtraktion. PETSc-Parität verwendet
rtol=2e-8/atol=2e-10, der Elektrodenfall bis 2e-9 wegen iterativer Lösung.
Keine Diskretisierungsfehlergarantie und keine experimentelle Validierung.

Lokal geprüft: `PYTHONPATH=src .venv/bin/python -m pytest -q`
mit 554 bestandenen und 53 optional übersprungenen Tests.
`ruff check .`, `ruff format --check .`, `mypy src/fem_toolbox` und
`git diff --check` bestanden. Native PETSc- und Artefaktprüfungen laufen in CI.

## Vier-Netze-Konvergenz

P=x⁴+y⁴/2+z⁴/4, u=(0,01P,0,0), φ=P; E=100 Pa, ν=0,25, eₓₓₓ=2 C/m²,
κ=3 I F/m. Exakte Dirichletwerte auf allen Außenflächen und unabhängig aus den
starken Gleichungen hergeleitete Quellen:

```text
b = (-38,4x²-2,4y²-1,2z², 0, 0)
rho = -35,76x²-18y²-9z²
```

Die Fehlernormen verwenden eine unabhängig formulierte P2-Basis und
Duffy-Produktquadratur mit sechs Gauss-Punkten pro Richtung. Die unabhängige
Tensorenergie prüft jedes Zellintegral (rtol=2e-12, atol=1e-16 J).

| Würfelteilung | L2(u) | H1(u) | L2(φ) | H1(φ) |
|---|---:|---:|---:|---:|
| 1 | 8,0134e-4 | 5,20768e-3 | 8,03014e-2 | 5,20145e-1 |
| 2 | 1,14835e-4 | 1,46417e-3 | 1,11431e-2 | 1,43245e-1 |
| 3 | 3,43590e-5 | 6,60758e-4 | 3,34567e-3 | 6,47629e-2 |
| 4 | 1,44641e-5 | 3,72287e-4 | 1,41756e-3 | 3,66569e-2 |

Letzte Ordnungen: 3,007 / 1,994 / 2,985 / 1,978. Beide letzten Intervalle
müssen L2-Ordnungen 2,7–3,3 und H1-Ordnungen 1,8–2,2 einhalten. Kraft-, Momenten-
und Ladungsbilanzen werden auf jedem Netz geprüft. #35 bleibt für Thermopiezo,
Dynamik und die ausstehenden Gesamt-Abnahme-/Kostenvergleiche offen.
