# Experimenteller DC-Spannungsteiler

Auf dem Experimental-Branch aus dem Repository-Root ausführen:

```sh
PYTHONPATH=src .venv/bin/python examples/circuit_dc/run.py
```

Nach einer Entwicklungsinstallation reicht `python examples/circuit_dc/run.py`.
Das Beispiel liest [divider.cir](divider.cir), schreibt keine Dateien und prüft:

```text
U_out = 8.000000 V
I_R1 = 0.004000 A
P_V1 = -0.048000 W
```

Handrechnung: `I=12/(1000+2000)=0.004 A`; `U_out=I*2000=8 V`.
Die Quelle liefert 48 mW, die Widerstände nehmen 16 mW und 32 mW auf.
Absolute Prüftoleranzen: 1e-12 V, 1e-14 A, 1e-13 W; für KCL und Gesamtleistung
je 1e-12 A beziehungsweise W. Bei diesen Größenordnungen liegen sie deutlich
über float64-Rundungsresten und weit unter den Referenzsignalen.

Die API ist experimentell. [Vertrag und Grenzen](../../Doc/65_EXPERIMENTAL_CIRCUIT_CORE.md):
DC-RLC und unabhängige Quellen, keine Transienten/AC/FEM- oder thermischen Ports.
