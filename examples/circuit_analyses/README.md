# Experimenteller RLC-Referenzlauf

```sh
PYTHONPATH=src .venv/bin/python examples/circuit_analyses/run.py --output /tmp/fem-circuit-analyses
```

Das Ausgabeziel muss neu sein. Nach Entwicklungsinstallation reicht `python`.
TOML und alle Parameter stehen in [circuit.toml](circuit.toml). Die Portnamen
sind Metadaten; der Lauf erzeugt kein FEM-Modell und keinen Wärmeübergang.

R=1 Ohm, L=1 H, C=1 F und ein 1-V-Sprung bei t=0 liefern
`v_C=1-exp(-t/2)*(cos(sqrt(3)*t/2)+sin(sqrt(3)*t/2)/sqrt(3))` V.
Der Lauf prüft einen vorher festgelegten Maximalfehler von 0.004 V für
Backward Euler mit dt=0.01 s und eine diskrete Energiebilanz von 1e-12 J.
Die separat getestete Zeitschrittfolge verifiziert die erste Ordnung;
kleine Bilanzreste allein sind kein Genauigkeitsnachweis.

Bei omega=1 rad/s, f=1/(2*pi) Hz, beträgt die komplexe Serienimpedanz 1 Ohm.
Bei 1 V Peak gelten I_R=1 A, V_C=-i V und P_R=0.5 W. Die absoluten
phasor-/Leistungsgrenzen von 1e-12 A/V/VA berücksichtigen nur float64-Rundung.
`transient.csv` enthält Felder, physikalische Arbeit/Verluste und numerische
Dämpfung; `report.json` enthält beide Provenienzen und Referenzkennwerte.

[Analysen-/Portvertrag](../../Doc/66_CIRCUIT_ANALYSES_AND_PORTS.md) und
[Verifikationsmatrix](../../Doc/67_CIRCUIT_CORE_VERIFICATION.md).
