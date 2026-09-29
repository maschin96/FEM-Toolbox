# Tet10-Wärme durch Model, Solver und Export

Nach Installation aus dem Repository-Root:

```sh
python examples/tet10_heat/run.py --output tet10_heat_results
```

Das Zielverzeichnis muss neu sein. Gmsh wird nicht benötigt: `body.msh` enthält
sechs geradkantige Tet10 und zwölf Triangle6, 27 Knoten, auf dem Einheitswürfel.
Das Netz entsteht aus der sechsfachen Freudenthal-Zerlegung eines Würfels mit
geteilten Kantenmittelpunkten; Knotenordnung und physische Gruppen liegen
explizit in der ASCII-Gmsh-2.2-Datei. Die Materialien sind synthetisch und alle
Größen SI; dies sind analytische Verifikationsfälle, keine Werkstoffvalidierung.

- `stationary.toml`: k=2 W/(m K), T=300 K links, Q=-4 W/m³, äußerer
  Wärmestrom -4 W/m² rechts, sonst isoliert. Exakt T=300+x² K und
  q=(-4x,0,0) W/m². Alle sechs Knoten jeder Dirichletfläche werden erfasst.
- `transient.toml`: rho=3 kg/m³, c=4 J/(kg K), isoliert, Q=12 W/m³,
  T(0)=300 K. Konsistente Kapazität und Backward Euler ergeben T=300+t K.
  Die Sonde (0.2,0.2,0.2) verwendet die quadratischen Formfunktionen.

`run.py` prüft Temperaturen, stationäre Schwerpunktflüsse und Leistungsbilanzen
gegen diese Formeln. `stationary.vtu/.json` und `transient/heat.pvd` samt VTUs,
Sonden-CSV und Bilanz-CSV sind die Ausgaben. **Tet10-Flüsse und Gradienten sind
nur rohe Schwerpunktwerte** (`*_at_centroid`); das Feld ist innerhalb der Zelle
ortsabhängig. Keine Knotenprojektion oder vollständige Quadraturpunktwolke.
Die Manifeste verwenden dafür Ergebnis-/Serienschema 2 und nennen die Punkte.

[Vertrag, Toleranzen und verbleibende Integrationsschritte](../../Doc/45_TET10_HEAT_INTEGRATION.md).
