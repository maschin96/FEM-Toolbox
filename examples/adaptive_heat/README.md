# Adaptive stationäre Wärme / M10

Installiertes `fem-toolbox` (M10-Kandidat), Python 3.11–3.13; kein optionales Extra.

```sh
python -I examples/adaptive_heat/run_adaptive.py --output results/adaptive_heat
```

Sieben Stufen auf einem Einheitswürfel: glattes quadratisches Temperaturfeld und
lokalisierte Gaußschicht. Wärmeleitfähigkeit 1 W/(m K), alle Randwerte analytisch.
Adaptive Dörfler-Markierung 50 % und uniforme Verfeinerung aller Zellen werden
vom selben Eingabenetz ausgeführt. Der tatsächliche Energiefehler wird unabhängig
mit Gauss-Duffy-Regeln der Ordnung 6 und 7 integriert und auf Quadraturparität
geprüft. `comparison.json` enthält Fehler, Indikator, Effektivität, DOFs,
Netz-Hash, Bilanz und gesamte Solve-/Refine-Laufzeit. Vier weitere JSON-Dateien
enthalten alle Netze und Rohfelder für die nachvollziehbare Stufenrekonstruktion.

Für den Installations-Smoke-Test: `--steps 4`. Der volle Nachweis verwendet den
Standard von sieben Stufen; die lokalisierte Fehlergrenze 1.5 sqrt(W K) wurde
nach Pilotläufen festgelegt und im finalen Wiederholungs-/Regressionstest fixiert.
514 adaptive gegenüber 729 uniformen DOFs erreichen diese Grenze im Nachweis.
Das ist ein fallspezifischer Vorteil, keine allgemeine Beschleunigungsgarantie.

Neue Ausgabeziele verwenden, um frühere Evidenz zu bewahren. `error` oder
unerwarteter Abbruchstatus führt zu einem fehlschlagenden Beispiel. Callables
werden in diesem Skript definiert; Hashes allein erlauben keine Neuausführung.
Details: [Vertrag](../../Doc/69_ADAPTIVE_REMESHING_CONTRACT.md) und
[Workflow](../../Doc/70_ADAPTIVE_HEAT_WORKFLOW.md).
