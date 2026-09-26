# Quadratisches elektrisches Potential mit Tet10

Nach Installation aus dem Repository-Root:

```sh
python examples/tet10_electrostatics/run.py --output tet10_electric_results
```

Das Zielverzeichnis muss neu sein. `body.msh` enthält sechs geradkantige Tet10,
zwölf Triangle6 und 27 Knoten auf dem Einheitswürfel. Die Datei verwendet die
Freudenthal-Zerlegung mit geteilten Kantenmittelpunkten wie das Tet10-Wärmebeispiel;
ASCII-Gmsh 2.2 und explizite physische Gruppen benötigen keinen Gmsh-Prozess.

`model.toml` beschreibt ein synthetisches Dielektrikum mit κ=2e-11 F/m,
Raumladung ρ=-4e-11 C/m³, φ=0 links und D·n=-4e-11 C/m² rechts, sonst D·n=0.
Exakt gelten φ=x² V, E=(-2x,0,0) V/m und D=(-4e-11 x,0,0) C/m².
Die positive gespeicherte Energie ist 4e-11/3 J; die Gesamtladungsbilanz ist null.
Das Beispiel prüft diese unabhängigen Formeln vor dem Export.

`electric.vtu` und `electric.json` enthalten nodale Potentiale/Reaktionen,
**rohe Schwerpunktwerte** für E und D sowie die **über jede Zelle integrierte**
Energie. Schema 2 kennzeichnet Feldorte und Integrale getrennt. Die Energie wird
nicht aus dem Schwerpunktwert hochgerechnet. Oberflächenladungen sind äußere
D-Flüsse; die Ladung des angrenzenden Leiters hat das umgekehrte Vorzeichen.

Dies ist Code-Verifikation mit synthetischen Daten, keine Materialvalidierung.
[Vertrag, Toleranzen und Grenzen](../../Doc/46_TET10_ELECTROSTATICS.md).
