# Vorbereiteter Kanal: Fluidnetz und Rand-Preflight

`mesh.json` enthält das explizite ebene Quad4-Netz in Metern und seine
Volumen-/Randgruppen. `domain.json` referenziert das Netz mit SI-Fluid,
Extrusionsbreite, Einlass, Druckauslass und zwei Wänden. Die obere Wand ist
als spätere FSI-Schnittstelle markiert. Dieser Lauf prüft Eingaben und
Geometrie; er berechnet keine Strömung und verformt kein Netz.

Aus dem Repository mit eingerichteter Entwicklungsumgebung:

```sh
PYTHONPATH=src .venv/bin/python examples/fluid_mesh/preflight.py --output /tmp/fluid-preflight.json
```

Ausgabeziel muss neu sein. `--input` akzeptiert alternativ einen Pfad zu
einer Domain-JSON. Der Meshpfad bleibt relativ zu deren Verzeichnis; die
Dateien werden auch bei Start aus einem anderen Arbeitsverzeichnis gefunden.

Erwartete Werte: 16 Zellen, 16 Randfaces, minimale Zellfläche 0.125 m²,
minimaler Wandabstand 0.125 m, Kantenverhältnis 2, Schichtwachstum 1,
vorgegebener Einlassmassenstrom 0.06875 kg/s. Letzterer folgt unabhängig aus
vier Mittelpunktproben von 4y(1-y), rho=1 kg/m³ und Breite 0.1 m; keine
Normierung auf einen kontinuierlichen Zielstrom. `solver_status="not_run"`
und `wall_resolution="not_checked_without_wall_shear"` verhindern einen
unbelegten Konvergenz- oder Wandauflösungsanspruch.

[API, Schema und Grenzen](../../Doc/72_FLUID_MESH_ADAPTER.md).
Quad-FV-Operator #226 und SA-Schließung #87 folgen separat; dieses JSON
ersetzt nicht das uniforme MAC-FlowProblem/TOML-Format aus #86.
