# Lokale Wärme-Parameterstudie

Mit installiertem Paket aus dem Repository-Root:

```sh
python examples/parameter_study/run.py --max-cases 1
python examples/parameter_study/run.py
python examples/parameter_study/run.py
```

Der erste Aufruf berechnet einen Fall und lässt den zweiten `pending`. Der zweite
setzt fort; der dritte übernimmt beide geprüften Bündel ohne Solve. Mit
`--output <verzeichnis>` ein eigenes Ziel wählen. `catalog.json` referenziert je
Fall Input, Parameter mit Einheiten, Status, Versuche und Ergebnisbündel. Das
Skript kann auch von außerhalb des Repositorys aufgerufen werden; Eingabepfade
sind relativ zum Beispiel, der Ausgabepfad relativ zum Arbeitsverzeichnis.

Ein gerader Einheitstetraeder mit Volumen 1/6 m³ besitzt drei auf 300 K fixierte
Knoten und einen freien Knoten bei x=1 m. Die Quelle ist Q=8 W/m³; Leitfähigkeiten
k=2 und 4 W/(m K) ergeben die **diskrete** Referenz
`K11=k/6`, `f1=Q/24`, also `T1=300+Q/(4k)` = 301 bzw. 300,5 K.
Natürliche Randbedingungen an übrigen Flächen sind homogen. Die Leistungsbilanz
muss innerhalb 1e-11 W schließen; Temperaturtoleranz 1e-11 K berücksichtigt
Rundung bei Subtraktion der 300-K-Grundtemperatur. Keine exakte Kontinuumslösung
für dieses Quellenproblem behaupten. Modell/Material sind synthetische Testdaten.

Für einen sichtbaren Fehlfall eine Kopie des Beispiels anlegen und in `study.toml`
`values = [2.0, 4.0, -1.0]` verwenden: der dritte Fall meldet unzulässiges Material,
die ersten bleiben gültig. Änderung der Quelle im Basismodell invalidiert die
bestehenden Ergebnisse. Ctrl-C während eines längeren Falls erzeugt
`interrupted`; der automatisierte Export-Abbruchtest in
`tests/unit/test_studies.py` erzeugt deterministisch ein halbes Bündel im zweiten
Fall und belegt, dass der erste erhalten bleibt und der zweite neu gerechnet wird.

[Vertrag, Fingerprints und Recovery](../../Doc/59_PARAMETER_STUDIES.md).
Das Beispiel enthält zwei einzelne Tet4-Läufe, keine großen Daten oder Extras.

Mit `python examples/parameter_study/report.py` entsteht anschließend eine neue
CSV-/JSON-Kennwerttabelle samt Markdown-Bericht, ohne erneuten Solve.
[Vertrag](../../Doc/60_STUDY_REPORTS.md).
