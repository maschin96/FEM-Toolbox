# Stationärer laminarer Quad-FV-Kanal (#230)

Die vorbereiteten SI-Dateien beschreiben einen 2 m × 1 m Kanal mit 0,1 m
Extrusionsbreite, rho=1 kg/m³ und mu=1 Pa·s. Der Einlass ist 4y(1−y) m/s,
der Auslassdruck 0 Pa; oben/unten gilt no-slip. Vier wandnormale Zellen
sind ein schneller Installationsnachweis, kein Vier-Netz-Abnahmelauf.
Die Midpoint-Einlassquadratur ergibt 0,06875 kg/s.

Nach Installation des Pakets, auch aus Wheel oder sdist:

```sh
python examples/quad_flow/run.py --output quad_flow_results
```

`--input` akzeptiert eine andere FluidDomain-JSON mit relativem Netzpfad.
`report.json` enthält SI-Bilanzen, diskrete und gradientenbasierte Dissipation
mit ausgewiesener Differenz sowie Provenienz und vollständige Rohresiduen.
`fields.json` enthält Zell-U/p, Face-Massenfluss, beide Impulsresiduen,
Kontinuität und Wandtraktionen. Die Face-Reihenfolge ist die importierte
Netzreihenfolge, Wandtraktionen folgen den ausgewählten line2-Wandzeilen.
`residuals.csv` enthält alle äußeren Iterationen. Ein vorhandenes
Ausgabeverzeichnis wird vor dem Solve abgelehnt.

SIMPLE wird mit GMRES beschleunigt; die physikalischen Restgates aus Doc/71
werden auf frisch rekonstruierten Feldern geprüft. Diese API führt keine
SA-Schließung, y+-Abnahme oder FSI-Kopplung aus.
[Algorithmus und Energieinterpretation](../../Doc/74_QUAD_FV_SIMPLE.md).
