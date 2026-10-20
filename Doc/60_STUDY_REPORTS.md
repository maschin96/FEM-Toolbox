# Kennwerttabellen und Studienberichte (#60)

`fem_toolbox.study_reports` liest gespeicherte stationäre #59-Kataloge und
bestehende Ergebnisbündel. Es führt **keinen Solve** aus. Für jeden aktiven
Katalogfall entsteht genau eine Zeile; archivierte `retired_cases` sind nicht
Teil der aktuellen Studie. Fehlende/unterbrochene/ungültige/fehlgeschlagene Fälle
bleiben sichtbar. Vor der Auswertung erfolgreicher Fälle werden beide
Bündelhashes und die Provenienz gegen den Katalog geprüft.
Dies umfasst die Fall-/Versuchspfade und die gespeicherte Studienzuordnung
gegen Parameter, aufgelöstes Dokument, Fingerprints und Versionen (#219).
Fremde Bündel und ältere Bündel ohne Studienzuordnung erhalten `invalid`
mit Fehlergrund und fehlenden Kennwerten; Fortsetzung verwendet dieselbe Prüfung.

## API

```python
from fem_toolbox.study_reports import Quantity, export_study_table, read_study_table

quantities = [
    Quantity(
        "temperature",
        "point",
        "temperature_K",
        "K",
        "maximum active-node temperature",
        reduction="max",
    ),
    Quantity(
        "balance",
        "metric",
        "heat.power_balance_W",
        "W",
        "signed global power imbalance",
        absolute_tolerance=1e-11,
    ),
]
export_study_table("runs/catalog.json", "table", quantities)
table = read_study_table("table")
```

`Quantity(name, source, field, unit, definition, reduction="scalar",
component="scalar", indices=None, absolute_tolerance=None)` ist unveränderlich.
Namen sind eindeutige ASCII-Bezeichner; `case_id` und `status` sind reserviert.
Einheit und fachliche Definition sind Pflicht. Die Einheit wird gegen Feldmetadaten
bzw. den SI-Vertrag des Metriknamens geprüft, nicht automatisch konvertiert.

- `source="metric"`: Punktpfad im JSON-Bündel, z. B.
  `piezoelectric.charge_balance_C`, `piezoelectric.force_balance_N` oder
  `piezoelectric.electrodes.0.charge_C`. Listenindizes folgen der Bündelreihenfolge.
- `source="point"`: exportierter Feldname, z. B. `displacement_m`.
  Nur `active_node_ids` werden ausgewertet; inaktive Nullfüllung ist kein Messwert.
- `source="cell"`: exportierter Volumenfeldname, z. B. `mechanical_energy_J`.
  Tet4/Tet10 werden unterstützt; Randzellen-/Triangle-Padding wird ausgeschlossen.

`component` wählt `x/y/z`, eine vorhandene Voigt-Komponente oder `norm`.
Bei skalaren Feldern bleibt es `scalar`. `indices` wählt explizite globale
Knotenindizes bzw. Volumenzellzeilen (nullbasiert), keine externen Entitäts-IDs.
Eine Region wird über ihre explizit dokumentierten Indizes ausgewertet. Ohne
Auswahl werden alle aktiven Knoten/Volumenzellen benutzt.
`scalar` verlangt einen Wert; `min`, `max`, `mean`, `sum`, `max_abs` reduzieren
nach Komponentenauswahl. `mean` ist ungewichtet. Eine Summe von Feldproben ist
kein Volumenintegral; Summen von Zellenergien sind dagegen die Summe bereits
integrierter Energiebeiträge. Keine Glättung, Projektion oder implizite Vermischung
von mechanischer Energie, Kopplung und dielektrischer Enthalpie.
`origins` übernimmt Rohfeldart (`kind`) und Auswertelage einschließlich Tet10-
Schwerpunktkennzeichnung; eine Reduktion wird nicht als kontinuierliches Feld
bezeichnet.

## Dateivertrag Schema 1

Der Export erzeugt ein **neues** Verzeichnis mit `values.csv`, `metadata.json`
und `report.md`; bestehende Verzeichnisse werden abgelehnt. Ein Schreibfehler
kann ein unvollständiges Verzeichnis hinterlassen. Metadaten werden zuletzt
veröffentlicht; erst ein erfolgreicher Leser bestätigt CSV-/Metadatenkonsistenz.

CSV: UTF-8, Komma, LF, Kopf `case_id,status,<Spalten>`. Parameter heißen
`parameter:<name>` und stehen alphabetisch vor den Kennwerten in Aufrufreihenfolge.
Numerische Werte verwenden JSON-Zahlendarstellung; Parametervektoren/-strings
werden JSON-kodiert und nach CSV-Regeln gequotet. Fehlwerte sind leere Zellen,
niemals erfundene Nullen. Fallreihenfolge entspricht dem Katalog.

JSON: `table_schema_version=1`, `columns`, `rows`, `catalog_sha256`,
`csv_sha256`, `limitations`. Jede Spalte enthält Rolle (`parameter`/`quantity`),
Einheit und Definition; Kennwertspalten zusätzlich vollständige Auswahl-/Reduktions-
und Toleranzangaben. Jede Zeile enthält Fall-ID, Status, `values` (Fehlwert `null`),
`missing` (Grund je Fehlwert), `origins`, Eingabe-Fingerprints, Versionen,
Ergebnisreferenzen und bei intakten Ergebnissen die ursprüngliche Provenienz.
Die relativen Ergebnisreferenzen bleiben relativ zum ursprünglichen Studienziel;
`catalog_sha256` und `result` ermöglichen die Zuordnung. Beim Weitergeben der
Tabelle Studienverzeichnis und Quellenherkunft mit angeben. Alle Artefakte sind
bei identischen Eingaben in derselben Umgebung deterministisch.

`solver_status` übernimmt nur die Solvermeldung. `balance_checks` bewertet
**ausschließlich explizit gewählte** Kennwerte mit `abs(value) <= absolute_tolerance`
in ihrer Einheit. `balance_status` ist `passed`, `failed`, `not_available` oder
`not_checked`. Kein gewählter Test bedeutet `not_checked`, nicht bestanden.
`mesh_convergence` bleibt `not_checked`; kein Solve/Export liefert automatisch
Netzkonvergenz oder experimentelle Validierung. Definition und Schwelle stehen
in `columns`; Nutzer müssen fachlich passende Bilanzgrößen wählen.

Nicht verfügbare Größen, falsche Einheiten, unbekannte Komponenten oder unzulässige
Auswahlen ergeben einen Fehlwert mit konkretem Grund. Defekte Bündel erhalten
Zeilenstatus `invalid`. Ein beschädigter Katalog oder ungültige globale Definition
bricht den Export ab. `read_study_table(directory)` prüft Version, CSV-Hash,
Spalten/Zeilen, Fall-ID-Eindeutigkeit und jeden Wert gegen die JSON-Metadaten.
Es benötigt die ursprünglichen FEM-Dateien nicht. Die Hashprüfung ist keine
Authentifizierung gegen absichtliche Manipulation.

## Beispiel und Nachweise

```sh
python examples/parameter_study/run.py
python examples/parameter_study/report.py
```

Temperaturmaxima: 301 und 300,5 K, Leistungsbilanz kleiner 1e-11 W (Herleitung
im [Studienbeispiel](../examples/parameter_study/README.md)). Nach einem mit
`--max-cases 1` begrenzten Studienlauf bleibt der zweite Tabellenwert leer.

`tests/unit/test_study_reports.py` prüft echte #59-Läufe ohne erneuten Solve,
Status/Fehlwerte, deterministische Dateien, Hashschäden, Einheiten, Indizes,
Bilanzverletzung und gekoppelte Verschiebung/Ladung/getrennte Energien gegen
Originalresultate. Werte werden direkt serialisiert; analytische Verschiebung
wird mit `atol=1e-12` für den kleinen SI-Referenzfall verglichen. Bestehende
Numerik bleibt unverändert; keine neue Konvergenzaussage. #205 konsumiert explizit
gewählte numerische Parameter-/Zielspalten, nicht den Markdown-Bericht.
