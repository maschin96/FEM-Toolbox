# Statische Tet10-Elastizität aus Gmsh und TOML

Nach Installation aus dem Repository-Root:

```sh
python examples/tet10_elasticity/run.py --output tet10_elastic_results
```

Das Zielverzeichnis muss neu sein. `body.msh` enthält sechs geradkantige Tet10,
zwölf Triangle6 und 27 Knoten auf dem Einheitswürfel. Die ASCII-Gmsh-2.2-Datei
verwendet die Freudenthal-Zerlegung mit geteilten Kantenmittelpunkten wie die
Tet10-Wärme-/Elektrostatikbeispiele; ein Gmsh-Prozess ist nicht erforderlich.

Synthetisches Material: E=100 Pa, ν=0. Lasten: b=(-2,0,0) N/m³ und rechts
t=(2,0,0) N/m², links u=0, übrige Flächen lastfrei. Exakt gelten
u=(0,01 x²,0,0) m, εxx=0,02 x und σxx=2x Pa, alle anderen Komponenten null.
Die gespeicherte Energie beträgt 0,02/3 J. `run.py` prüft Verschiebungen,
Spannungen, Energie sowie Kraft- und Momentenbilanz gegen diese Formeln.

`solid.vtu` und `solid.json` enthalten nodale Verschiebungen/Reaktionen,
**rohe Schwerpunktwerte** für Dehnung und Spannung sowie die **über jede Zelle
integrierte** Energie. Schema 2 kennzeichnet Feldorte und Integrale getrennt.
Die Voigt-Reihenfolge ist xx, yy, zz, yz, xz, xy mit technischen Schubdehnungen.

Der Nachweis ist Code-Verifikation, keine Materialvalidierung. Thermische
Tet10-Eigendehnungen und gekoppelte Tet10-Analysen bleiben weitere Arbeit.
[Vertrag und zusätzliche Biege-/Konvergenznachweise](../../Doc/47_TET10_ELASTICITY.md).
