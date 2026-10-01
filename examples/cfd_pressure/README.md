# Experimentelle Druckkopplung: glatte Stokes-Referenz

```sh
PYTHONPATH=src .venv/bin/python examples/cfd_pressure/run.py
```

Der Lauf löst Geschwindigkeit und Druck gemeinsam auf vier Netzen und prüft
RMS-Fehler und Konvergenzordnung gegen unabhängig hergeleitete Funktionen.
Gleichungen, SI-Einheiten, Quellen, Druck-Gauge und feste Fehlerbudgets stehen
in [Doc/69](../../Doc/69_CFD_PRESSURE_COUPLING.md). Alle vier Wände sind no-slip.
Keine Druckwerte werden aus der analytischen Referenz in den Solver eingespeist;
die Quelle folgt aus deren starken Gleichungen. Dies ist Code-Verifikation,
keine RANS-Validierung oder allgemeine Navier-Stokes-Abnahme.
