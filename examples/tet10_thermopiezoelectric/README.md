# Tet10: thermisch induzierte Leerlaufspannung

```sh
python examples/tet10_thermopiezoelectric/run.py --output /tmp/tet10-thermopiezo-demo
```

Das Zielverzeichnis muss neu sein. Das Gmsh-2.2-Netz enthält sechs geradkantige
Tet10, zwölf Triangle6 und 27 Knoten im Einheitswürfel. Gmsh ist zum Ausführen
nicht erforderlich. `Model.from_file` liest alle Daten aus `model.toml`.

T=1 K links, T=2 K rechts und isolierte Seiten ergeben T=1+x. Alle mechanischen
Knoten sind festgehalten. Die linke Elektrode ist geerdet, die rechte schwebend
mit Gesamtladung null. Mit E=100 Pa, ν=0, κ=3 I F/m, eₓₓₓ=2 C/m²,
αxx=0,01/K, Tref=1 K und `electric_thermal_strain=true` folgt:

- φ=-0,01x²/3 V und D=0,
- σxx=-76x/75 Pa,
- Gesamtenthalpie 0,005/3+0,0002/9 J.

Das Skript prüft diese Felder, die Enthalpie und alle Bilanzen und schreibt
`coupled.vtu`/`coupled.json` mit Schema 2. Alle Elektroden-Mittelknoten werden
erfasst. Schaltet man den elektrischen thermischen Term aus, bleibt φ=0;
der Prüfcode bezieht sich auf die mitgelieferte aktivierte Konvention.

Synthetische SI-Parameter, keine experimentelle oder Tieftemperaturvalidierung.
Keine Wärmerückwirkung, Pyroelektrizität oder Dynamik. Details und weitere
Nachweise: [Doc/50](../../Doc/50_TET10_THERMOPIEZOELECTRIC.md).
