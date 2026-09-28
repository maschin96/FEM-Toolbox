# Stationäre Tet10-Piezoelektrizität mit schwebender Elektrode

```sh
python examples/tet10_piezoelectric/run.py --output /tmp/tet10-piezo-demo
```

Das Zielverzeichnis muss neu sein. Das mitgelieferte Gmsh-2.2-Netz enthält
sechs geradkantige Tet10, zwölf Triangle6 und 27 Knoten im Einheitswürfel.
`Model.from_file` liest Materialregion, mechanische Lagerung und Elektroden
aus `model.toml`. Keine Gmsh-Installation erforderlich.

Synthetisches Material: E=100 Pa, ν=0, κ=3 I F/m, einziges eₓₓₓ=2 C/m².
Bei u=0 und φ=0 auf x=0 sowie einer schwebenden Elektrode mit Q=3,04 C auf
x=1 ist die unabhängige 1D-Lösung φ=x V, uₓ=-0,02x m, σ=0 und Dₓ=-3,04 C/m².
Alle sechs Knoten jeder Elektrodenfläche sind equipotentiell; das Vorzeichen
von Q ist das der angrenzenden Leiterladung, also entgegen D·n.

Das Skript prüft Felder, Elektrodenladung, Kraft-/Momenten-/Ladungsbilanz und
Gesamtenthalpie -1,52 J. Es schreibt `piezo.vtu` und `piezo.json` im Bündelschema 2.
Spannung, Dehnung, E und D sind Schwerpunktwerte; die drei Energiebeiträge
sind Zellintegrale. Die elektrische Enthalpie ist keine positive Speicherenergie.

Stationäre Piezoelektrizität ist unterstützt. Tet10-Thermopiezo und
Tet10-Piezo-Dynamik sind weiterhin Folgearbeit. Der Fall ist ein analytischer
Code-Nachweis mit synthetischen Daten, keine experimentelle Validierung.
