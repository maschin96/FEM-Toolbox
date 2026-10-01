# Stationärer inkompressibler Kanal mit Einlass und Druckauslass

```sh
PYTHONPATH=src .venv/bin/python examples/cfd_flow/run.py --output /tmp/fem-cfd-flow-results
```

Das Ausgabeziel muss neu sein. Alle Daten stehen in [flow.toml](flow.toml):
32x16 Zellen, 2x1 m, rho=1 kg/m³, mu=1 Pa s und Einlass u=y(1-y) m/s.
Druckauslass 0 Pa, ruhende Wände. Referenz: dp/dx=-2 Pa/m, Q=1/6 m²/s.
Der Lauf prüft Fehler <0.003 m/s und globale Massenbilanz <1e-9 kg/s.

`pressure.csv` enthält Zellfelder, `velocity.csv` versetzte Flächenfelder,
`report.json` Residuen, Iterationsstatus, Bilanzen und Provenienz. Upwind und
Picard verwenden die dokumentierten Standardwerte; für glatte Referenzen kann
die Python-API explizit zentrale Konvektion mit Peclet-Gate wählen.
[Mathematik, API-Vertrag und Vier-Netz-Nachweise](../../Doc/70_STATIONARY_FLOW_CORE.md).
Keine Turbulenzschließung, keine FSI und keine Validierung behauptet.
