# Stationäre Tet10-Thermomechanik

Nach Installation aus dem Repository-Root:

```sh
python examples/tet10_thermomechanics/run.py --output tet10_thermal_expansion_results
```

Das Zielverzeichnis muss neu sein. Das mitgelieferte ASCII-Gmsh-2.2-Netz
enthält sechs geradkantige Tet10, zwölf Triangle6 und 27 Knoten auf dem
Einheitswürfel (Freudenthal-Zerlegung mit geteilten Kantenmittelpunkten).
Gmsh wird zum Ausführen nicht benötigt. Der bestehende Python-Einstieg
`solve_thermomechanics` löst Wärme und überträgt alle Knotentemperaturen an
Mechanik. Reine Thermomechanik hat keinen eigenen TOML-Model-Modus.

Synthetische SI-Daten und unabhängige Referenz:

- k=2 W/(m K), T=1 K links, Q=-4 W/m³, äußerer Wärmestrom -4 W/m² rechts,
  übrige Flächen isoliert: T=1+x² K.
- E=100 Pa, ν=0, αxx=0,01/K (übrige Komponenten null), Tref=1 K.
  Links u=0; b=(-2+2x,0,0) N/m³; rechts t=(1,0,0) N/m², sonst lastfrei.
- Exakt u=(0,01x²,0,0) m, σxx=2x-x² Pa, andere Spannungen null.
  Mechanische Energie W=0,008/3 J; Wärme-, Kraft- und Momentenbilanzen null.

Die Temperatur nahe 1 K ist ein synthetischer Manufactured-Solution-Fall zur
Trennung von Kopplungsfehlern und numerischer Auslöschung bei großen Offsets;
keine Werkstoffgültigkeit bei tiefen Temperaturen wird behauptet. Separate
Tests decken Temperaturen um 300 K und Referenztemperaturverschiebungen ab.

`run.py` prüft diese Felder und Integrale vor dem Export nach
`coupled.vtu`/`coupled.json`. Schema 2 unterscheidet rohe Schwerpunktwerte von
integrierter mechanischer Energie. Thermische Last und Energie verwenden
quadratische Temperaturinterpolation an den Integrationspunkten.

[Vertrag, Toleranzen und Vier-Netze-Nachweis](../../Doc/48_TET10_THERMOMECHANICS.md).
