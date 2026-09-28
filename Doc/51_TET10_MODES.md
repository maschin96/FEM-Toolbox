# Tet10-Modalanalyse (#35d)

Geradkantige Tet10 unterstützen `solve_piezoelectric_modes` und `Model`-Schema 3
mit `physics="piezo_modes"`. Konsistente mechanische Masse wird quadraturbasiert
integriert; das elektrische Feld besitzt keine Masse. Die bestehende elektrische
Schur-Elimination, homogene Lagerungen und feste/schwebende Triangle6-Elektroden
werden verwendet. Dichte kann global oder zellweise angegeben werden.

Der Referenzsolver bleibt auf 2048 aktive Freiheitsgrade vor Reduktion begrenzt.
Kein Mass-Lumping, kein Großmodell-Eigensolver, keine gekrümmten Elemente.
Tet10-Newmark und -Harmonik bleiben ausdrücklich gesperrt. #35 ist damit noch
nicht abgeschlossen; auch die übergreifenden Kosten-/Abnahmevergleiche fehlen.

## Numerik und Export

Der mathematisch symmetrische dielektrische Tet10-Block wird vor dem SPD-Solve
mit `(K + K.T)/2` symmetrisiert: Quadraturauslöschung erzeugt sonst asymmetrische
Rundungsreste an theoretisch verschwindenden Einträgen. Die strenge Backend-
Symmetrieprüfung bleibt erhalten. Physikalische Residuen werden weiterhin am
ursprünglichen assemblierten System gemessen.

`export_piezo_series` verwendet für Tet10 Manifest-Schema 2. Spannung, Dehnung,
E und D sind rohe Schwerpunktwerte (`*_at_centroid_real/imag`); Koordinaten und
`projection="none"` stehen im Manifest. Messpunkte interpolieren u und φ mit
allen zehn P2-Formfunktionen. Tet4 behält Schema 1 und seine Feldnamen.
Moden sind massennormierte Formen; ihre Amplituden sind keine reale Anregung.

## Unabhängige Verifikation

Für den Einheits-Rechtstetraeder gilt für die freie Eckfunktion N=x(2x−1):
∫N² dV=1/420 und ∫(∂N/∂x)² dV=1/10. Bei ρ=6, Cxx=10, exxx=2 und κxx=3
ist m=6/420 und ω²=(10+4/3)/(10m) für freies Eckpotential; bei φ=0 entfällt
4/3. Geprüft werden Frequenz, Massennormierung, elektrisches Feld, zwei
orthogonale Moden, umgekehrte Orientierung und unbenutzte Punkte. Relative
Toleranzen 3–4e-14 decken Rundung dieser kleinen analytischen Systeme ab.

Ein längs schwingender Stab (L=0,02 m, A=1e-4 m², ρ=7500 kg/m³,
Czz=1e9 Pa, ezzz=5 C/m², κ=1e-8 F/m) hat CD=Czz+ezzz²/κ.
Mit unterem festen Ende und blockierter Querverschiebung gilt
ω=λ/L sqrt(CD/ρ). Offen ist λ=π/2; kurzgeschlossen löst λ cot λ=ezzz²/(κ CD).
Die obere offene Elektrode hat Q=0; alle ihre Knoten sind equipotentiell.

| Axiale Teilungen | relativer Frequenzfehler offen | kurzgeschlossen |
|---:|---:|---:|
| 1 | 9,25560e-4 | 3,99287e-4 |
| 2 | 1,02502e-4 | 3,72655e-5 |
| 4 | 7,68848e-6 | 2,70504e-6 |
| 8 | 5,03172e-7 | 1,75838e-7 |

Die letzten Ordnungen sind 3,93 bzw. 3,94. Dies ist axiale Verfeinerung eines
Longitudinalfalls, kein allgemeiner dreidimensionaler H1/L2-Nachweis.
Feldweise relative Residuen bleiben unter 1e-9; offene Ladung wird relativ zur
Kapazität κA/L mal Modalspannung mit 1e-10 geprüft. PETSc direkt/FieldSplit wird
in real/complex gegen SciPy verglichen (Frequenzen rtol=2e-9, Residuen <1e-8).
Diese Nachweise sind Code-Verifikation, keine experimentelle Validierung.

## Reproduktion und TDD

Die Eckmodus-Tests scheiterten vor der Implementierung an der Tet4-Sperre,
der Model-Test an der Tet10-Sperre. Der Beispieltest scheiterte zunächst am
fehlenden Skript. Nach Umsetzung bestehen diese Tests einschließlich Model-
Rundreise, Schwerpunktfeldern und unabhängiger P2-Messpunktinterpolation.

```sh
PYTHONPATH=src python -m pytest tests/unit/test_tet10_modes.py tests/convergence/test_tet10_modal.py
PYTHONPATH=src python examples/tet10_modes/run.py
```

Das Gmsh/TOML-Beispiel prüft denselben analytischen Eckoszillator und exportiert
VTU/PVD, Manifest und Messpunkt-CSV. Die Artefakt-CI führt es aus installiertem
Wheel und sdist aus.
