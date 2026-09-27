# Stationäres Ergebnisbündel, Schema 1 (#140)

Status: Implementiert; Roundtrip- und Fehlerfälle sind in
`tests/unit/test_result_bundle.py` geprüft. Die Ausgabe ist ein Maschinenvertrag
für Rohfelder und zentrale Kennwerte, keine wissenschaftliche Validierung.

## Schreiben und Lesen

```python
from fem_toolbox.io import export_result_bundle
from fem_toolbox.model import Model

model = Model.from_file("examples/thermopiezoelectric/case.toml")
fields, manifest = export_result_bundle("result", model.mesh, model.solve())
```

Dies schreibt `result.vtu` und `result.json` in ein **vorhandenes** Verzeichnis.
Ein Dateipräfix ohne Endung wird verlangt. Bestehende Ausgabedateien werden
nicht überschrieben. Beide Dateien werden zuerst im selben Verzeichnis
vorbereitet; bei einem Schreibfehler vor Veröffentlichung bleibt kein halbes
Bündel, und bei einem Fehler zwischen den beiden exklusiven Dateiverknüpfungen
wird die neu geschriebene VTU-Datei entfernt. Das JSON-Manifest wird zuletzt
veröffentlicht. Die exklusiven Verknüpfungen verhindern auch ein Überschreiben
bei gleichzeitigem Export auf dasselbe Präfix.

Ohne Import eines FEM-Solvers lassen sich Metadaten, Bilanzen und der
SHA-256-Prüfwert der VTU-Datei lesen:

```sh
python -I scripts/read_result_bundle.py result.json
```

Der Leser benötigt nur die Python-Standardbibliothek. Numerische Felder in
`result.vtu` können beispielsweise mit `meshio.read` geprüft werden.

## Format

`bundle_schema_version = 1` bezeichnet diesen Exportvertrag; die
Modell-Schemaversion wird getrennt unter `model_schema_version` geführt.
`toolbox_version` und `provenance_schema_version` sind direkt im Manifest
ablesbar, auch bei direkten Solverläufen ohne `Model`.
`provenance` enthält den Modelleingabehash samt Toolbox-Version und
Solverparametern für `ModelResult`. Bei direkten Solverläufen enthält
`solver_provenance` die jeweiligen Solverdatensätze; der Bündelhash ist ein
SHA-256 über deren sortierte Eingabehashes. `files.fields_sha256` prüft die
Byte-Identität der zugehörigen VTU-Datei. Zeitstempel und lokale Pfade gehen
nicht in das Manifest ein.

`fields.point` und `fields.tetra4` beschreiben Namen, SI-Einheiten, Formen,
Komponentenreihenfolge und
`kind = "raw"`. Knotenfelder verwenden die in `active_node_ids` angegebene
Reihenfolge der Solver-Ergebnisse. Ungenutzte Netzknoten erhalten im VTU eine
Nullfüllung und sind nur über die IDs als inaktiv erkennbar. Zellfelder folgen
der Eingabereihenfolge der Tet4-Zellen; `mesh.tetra4_ids` nennt externe IDs,
falls vorhanden. Für Triangle3-Flächen wird bei zellweisen Volumenfeldern eine
Nullfüllung geschrieben und als `triangle3_padding` gekennzeichnet. Benannte
Auswahlsätze bleiben als `fem_toolbox:set:<name>`-Masken in der VTU-Datei.

| Physik | Primäre Knotenfelder | Rohfelder pro Tet4 |
|---|---|---|
| Wärme | `temperature_K` | Gradient `K/m`, Wärmestrom `W/m^2` |
| Elastizität | `displacement_m` | Engineering-Dehnung `1`, Spannung `Pa`, Energie `J` |
| Elektrostatik | `potential_V` | E-Feld `V/m`, D-Feld `C/m^2`, Energie `J` |
| Piezo | `displacement_m`, `potential_V` | Dehnung, Spannung, E, D und getrennte Energie-/Enthalpiebeiträge |

Knotenreaktionen, Lasten und Residuen werden mit feldspezifischen Namen und
Einheiten exportiert. `metrics` enthält Leistungen, Kraft-, Moment- und
Ladungsbilanzen, feldweise relative Residuen, Skalierungsgrößen und summierte
Energie. Bei Piezo-Elektroden stehen Potential und **vorzeichenbehaftete**
Gesamtladung je Name in `metrics.piezoelectric.electrodes`. Die
mechanische Energie, Kopplungsenergie und dielektrische Enthalpie bleiben
separate rohe Zellfelder; `total_enthalpy_J` ist ihre vorzeichenbehaftete
Summe. Ein erfolgreich konstruiertes Solver-Ergebnis erhält
`solver_status = "converged"`.

`limitations` benennt die grundsätzlichen Grenzen: stationär, linear, Tet4,
keine aus dem Export abgeleitete Netzkonvergenz oder experimentelle
Validierung. Im kombinierten Modell ist die Wärmekopplung einseitig.

## Numerische Prüfung

Der Roundtrip vergleicht die VTU-Rohwerte direkt mit den Ergebnisobjekten und
liest zentrale Kennwerte per isoliertem Standardbibliotheks-Skript. Die
gekoppelte Beispielrechnung aus [Doc/22](22_THERMOPIEZOELECTRIC.md) nutzt
analytisch bekannte Felder und Bilanzgrenzen von `1e-12` für Rohfelder und
`1e-10` für globale Bilanzen. Die Exportprüfungen verwenden keine zusätzliche
Numeriktoleranz: JSON-Kennwerte sind dieselben Python-Fließkommazahlen, und
VTU-Felder werden mit `assert_allclose` gegen ihre Solverarrays verglichen.

## M7-Erweiterung: Tet10-Wärme

Stationäre Tet10-Wärme verwendet Bündelschema 2, Zelltyp `tetra10` und explizit
als Schwerpunktwerte benannte Gradienten/Flüsse samt Auswertekoordinaten.
Tet4 behält Schema 1. Der eigenständige Manifestleser unterstützt beide
Versionen. [Vertrag und Nachweise](45_TET10_HEAT_INTEGRATION.md).

Tet10-Elektrostatik erweitert Schema 2 um rohe E-/D-Schwerpunktfelder und
`electric_energy_J` als `cell_integral` am `cell_volume`. Energie erhält keinen
Schwerpunktsuffix. [Details](46_TET10_ELECTROSTATICS.md).

Statische Tet10-Elastizität verwendet Schema 2 mit
`strain_1_at_centroid`/`stress_Pa_at_centroid` und den sechs Voigt-Komponenten.
`mechanical_energy_J` ist ein `cell_integral` am `cell_volume`, ohne
Schwerpunktsuffix. [Details](47_TET10_ELASTICITY.md).

Gekoppelte Tet10-Thermomechanik exportiert Wärme und Mechanik gemeinsam in
Schema 2. Beide Teilresultate müssen passende Schwerpunktkoordinaten liefern.
[Vertrag und Nachweise](48_TET10_THERMOMECHANICS.md).
