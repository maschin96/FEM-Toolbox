# Experimenteller CFD-Kanaloperator

```sh
PYTHONPATH=src .venv/bin/python examples/cfd_channel/run.py
```

H=1 m, mu=2 Pa s, dp/dx=-4 Pa/m, ruhende Wände, 64 uniforme Zellen.
Die unabhängige Referenz ist u=y(1-y) m/s, Q=1/6 m²/s, tau=(2,-2) Pa.
Der erste Schnitt verifiziert ausschließlich den viskosen Impulsoperator;
keine allgemeine Druckkopplung, Turbulenzschließung oder FSI.
[Vertrag und Verifikation](../../Doc/68_EXPERIMENTAL_CFD_CORE.md).
