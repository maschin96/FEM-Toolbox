# Model-Ergebnisse als JPEG

```sh
python -m pip install '.[visualization]'
python examples/field_visualization/run.py --output field_images
```

Aus dem Repository-Root ausführen. Für Wiederholungen ein neues Ausgabeverzeichnis
verwenden; vorhandene Bilder werden bewusst nicht überschrieben. `run.py`
zeigt `Model.from_file(...).solve().field_result` bis `plot_field`/`save_jpeg`.
Die sechs JPEGs sind bei 8×6 Zoll und 150 DPI jeweils 1200×900 Pixel groß.
Der Pfad benötigt keine grafische Sitzung. Jeder erzeugte Plot wird geschlossen.

Das eigenständige Model mit mitgeliefertem Netz stammt aus dem verifizierten
[stationären Thermopiezo-Beispiel](../thermopiezoelectric/README.md): ein Tet4,
synthetische SI-Materialdaten, vorgegebene Knotentemperaturen und einseitige
thermische Kopplung. Keine experimentelle Materialvalidierung und kein
Netzkonvergenzanspruch. Materialachsen, Lasten und Lagerung stehen in `case.toml`.

Der erste Lauf hat T=300/310/300/300 K und qx=-20 W/m². Es entstehen Bilder
für T, |u|, qx, phi und sigma_xx. Der zweite Lauf verwendet 305 K am heißen
Knoten; beide Temperaturbilder nutzen die Grenzen 300–310 K, `viridis` und
dieselbe Ansicht. Die Flächenfarbe ist der Mittelwert der Knotenwerte,
weshalb das Temperaturmaximum der Farbfläche kleiner als 310 K ist.
Beim Verschiebungsbetrag wird zuerst je Knoten die Norm gebildet;
qx und sigma_xx bleiben ungemittelte Eigentümerzellwerte.

Der [API-Vertrag](../../Doc/62_FIELD_VISUALIZATION.md) erklärt Grenzen,
Komponenten, Einheiten, aktive Knoten, Export und weitere Plotparameter.
