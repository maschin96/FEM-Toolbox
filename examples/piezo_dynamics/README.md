# Piezo-Dynamik: Model, CSV-Anregung und ParaView

Aus dem Repository-Root:

```sh
PYTHONPATH=src .venv/bin/python examples/piezo_dynamics/run.py --output results/piezo_dynamics
```

Das Beispiel enthält **synthetische Verifikationsparameter**, keine realen
Materialdaten oder experimentellen Messungen. Der Einheits-Tet4 hat einen freien
mechanischen und einen algebraischen elektrischen Freiheitsgrad. Die Tabelle
`drive.csv` stellt einen demonstrativen, linear interpolierten Anregungsverlauf
dar und skaliert die Volumenkraft aus `model.toml`. Eigene CSV-Kurven können
über dieselbe Schnittstelle eingelesen werden; Einheiten und Spaltennamen sind
explizit. Das Beispiel berechnet außerdem einen komplexen Frequenzsweep und
einen massennormierten Eigenmodus.

Unter `transient/`, `harmonic/` und `modes/` entstehen jeweils `piezo.pvd`,
VTU-Felder, `probes.csv` und `manifest.json`. Die PVD-Datei in ParaView öffnen.
Bei der transienten Rechnung ist die Abszisse Zeit in s, beim Frequenzsweep
Frequenz in Hz und bei Moden der Modenindex. Die genaue Achsenbedeutung und
Kreisfrequenz stehen im Manifest; ParaView nennt den gemeinsamen Slider
technisch „Time“. Harmonische Felder liegen als Real-/Imaginärteil und
komponentenweise Spitzenamplitude/Phase in rad vor. Phase bei null Amplitude
ist ungültig und mit `phase_valid=0` markiert.

`tip` ist Knoten 1, `center` der Punkt `[0.25,0.25,0.25] m`. Die CSV interpoliert
die Feldwerte am Punkt. Zellfelder bleiben konstante Tet4-Werte und werden
nicht geglättet. Bei einer Änderung am Modell einen neuen Ausgabeordner wählen;
vorhandene Ordner werden nicht überschrieben.

Grenzen: keine Dämpfung, keine thermische Rückwirkung, maximal 2048 aktive
Gesamt-DOFs. Ein Frequenzraster erkennt keine Resonanzen zwischen Stützstellen.
Physikalische Gleichungen, Toleranzen und Verifikation: [Dokument 34](../../Doc/34_PIEZO_DYNAMICS.md).
