# Transiente Wärme: Aluminium, CSV-Anregung und Messpunkte

Nach Installation aus dem Repository-Root:

```sh
python examples/transient_heat/run.py
```

Das Beispiel rechnet einen kleinen EN-AW-6060-T6-Tetraeder in SI-Einheiten.
Eine CSV-Quelle fährt von 0 auf 2 MW/m³ und zurück; die Umgebungstemperatur
springt bei 5 s von 293,15 auf 298,15 K. Konvektion ist h=20 W/(m² K).
Dieser h-Wert ist eine **Beispielannahme**, kein Materialbibliothekswert.
Das grobe Einzelelement demonstriert den Workflow, nicht die räumliche
Konvergenz eines realen Kühlkörpers.

Öffne `transient_heat_results/heat.pvd` in ParaView. Temperatur und Wärmefluss
stehen für jeden gespeicherten Zeitpunkt bereit. `probes.csv` enthält den
interpolierten Mittelpunkt und Knoten 0 für alle Schritte; `balance.csv` die
Wärmebilanz. Der Ordner darf vor dem Lauf nicht existieren.

Die Quelle lässt sich in `model.toml` beispielsweise ersetzen durch:

```toml
[loads.source]
unit = "W/m^3"
value = { kind = "step", time = 2.0, before = 0.0, after = 1000000.0 }
```

Oder als Rampe:

```toml
[loads.source]
unit = "W/m^3"
value = { kind = "ramp", start = 2.0, end = 4.0, start_value = 0.0, end_value = 1000000.0 }
```

CSV-Spalten und Einheit werden explizit angegeben. Messzeiten müssen streng
steigen; außerhalb des Messbereichs wird standardmäßig abgebrochen. Alle
Profilknoten ergänzen das Zeitraster automatisch.

[Materialquellen](../../Doc/32_MATERIAL_LIBRARY.md) ·
[API, Schema, Bilanzen und Verifikation](../../Doc/33_TRANSIENT_HEAT.md)
