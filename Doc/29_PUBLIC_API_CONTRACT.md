# Öffentliche stationäre Python-API und Kompatibilität (#57)

Status: **Vertragsentwurf für v1.0.0** auf Basis des aktuellen Quellstands.
Die veröffentlichte Version ist weiterhin v0.5.0. Die hier genannten
Einstiegspunkte sind die für v1.0 vorgesehenen stabilen Oberflächen; die
Release-Freigabe erfolgt erst über #28. Der Vertrag umfasst lineare,
stationäre Tet4-Modelle in SI-Einheiten. Numerische Gültigkeitsgrenzen stehen
in [Doc/26](26_VALIDITY_AND_LIMITS.md).

## Umfang der stabilen Oberfläche

Der Importpfad gehört zum Vertrag. `fem_toolbox` selbst exportiert nur
`__version__`; direkte Solver werden aus ihren Modulen importiert. Weitere
erreichbare Klassen und Funktionen in `assembly`, `elements`, `fields`,
`physics` und nicht aufgeführten Untermodulen sind Forschungsbausteine und
hier nicht als stabile Anwender-API ausgewiesen.

| Importpfad | Einstiegspunkte | Eingabe und Rückgabe |
|---|---|---|
| `fem_toolbox.mesh` | `Mesh`, `inspect_boundary_topology` | Netz aus Punkten `(n,3)` in m, Tet4 `(e,4)`, optional Triangle3 `(f,3)` und benannten Indexsätzen; unveränderliche Arrays. Topologiebericht für äußere Flächen. |
| `fem_toolbox.materials` | `ThermalMaterial`, `ElasticMaterial`, `DielectricMaterial`, `PiezoelectricMaterial`, `ThermalExpansion`, `assign_region_materials` | Validierte SI-Tensoren; regionsweise Zuordnung ergibt `MaterialAssignment` in Tet4-Reihenfolge. Isotrope und, wo angeboten, Achsen-/Datenblatt-Konstruktoren sind Teil der jeweiligen Klassen. |
| `fem_toolbox.boundary` | `Dirichlet`, `Electrode` | Globale DOF-Indizes und SI-Werte; Elektrode über Dreiecksflächen mit Potential in V oder Gesamtladung in C. |
| `fem_toolbox.model` | `Model`, `ModelResult` | `Model(mesh=..., document=...)` oder `Model.from_file(path)`; `solve()` liefert `ModelResult(physics, field_result, provenance)`. `to_document()` liefert eine unabhängige Kopie. |
| `fem_toolbox.solvers.heat` | `thermal_dofs`, `solve_heat`, `HeatResult` | Ein thermisches Material oder eines pro Tet4, Dirichlet und optionale Quellen/Flächenlasten; `HeatResult`. |
| `fem_toolbox.solvers.elasticity` | `mechanical_dofs`, `solve_elasticity`, `ElasticResult` | Steifigkeit, Lagerung und optionale Kräfte, thermische Expansion/Temperatur; `ElasticResult`. |
| `fem_toolbox.solvers.electrostatics` | `electric_dofs`, `solve_electrostatics`, `ElectrostaticResult` | Permittivität, Potentialreferenz und optionale Ladungen; `ElectrostaticResult`. |
| `fem_toolbox.solvers.piezoelectric` | `piezoelectric_dofs`, `solve_piezoelectric`, `PiezoelectricResult` | Stress-Charge-Material, mechanische und elektrische Randbedingungen, optionale Elektroden/Lasten; `PiezoelectricResult`. |
| `fem_toolbox.solvers.thermopiezoelectric` | `solve_thermopiezoelectric`, `ThermopiezoelectricResult` | Wärme zuerst, dann `(u, φ)` auf demselben Netz; Paar aus Wärme- und Piezo-Ergebnis. |
| `fem_toolbox.io` | `read_gmsh`, `write_vtu`, `export_result_bundle` | Gmsh-Import nach `Mesh`, VTU-Export und `(vtu_path, json_path)` eines versionierten Ergebnisbündels. |
| `fem_toolbox.provenance` | `RunProvenance` | Deterministischer Eingabehash, Solvername, Parameter und Toolbox-Version; `to_dict()` ist JSON-kompatibel. |
| `fem_toolbox.errors` | `FEMToolboxError`, `ValidationError` und spezialisierte Unterklassen | Typisierte Fehler mit `field`; `ModelSchemaError` nennt zusätzlich den Eingabepfad. |

`solve_thermomechanics` bleibt als direkter v0.4-Einstiegspunkt nutzbar,
einschließlich `ThermomechanicalResult`; er besitzt keinen eigenen
`Model.physics`-Modus. Der Weg über `solve_elasticity(..., expansion=...,
temperature=...)` und `solve_heat(...)` ist ebenfalls möglich. Dieser
Sonderfall wird vor v1.0 nicht stillschweigend entfernt.

## Formen, Reihenfolge und Einheiten

`Mesh` verwendet nullbasierte Knoten-/Zellindizes. Aktive Knoten stehen in
aufsteigender Netzindexfolge in `result.node_ids`; inaktive Punkte haben
keine Solver-DOFs. Zellergebnisse folgen **immer** der Eingabereihenfolge von
`mesh.tetra4`. `thermal_dofs`, `mechanical_dofs`, `electric_dofs` und
`piezoelectric_dofs` erzeugen feldweise Blöcke; innerhalb eines Felds folgt
auf jeden aktiven Knoten die Komponentenreihenfolge `x,y,z` für `u`.
Direkte `Dirichlet`-DOFs müssen aus der passenden Map stammen.

| Ergebnis | Primärfeld | Rohfelder je Tet4 | Wichtige Bilanzen |
|---|---|---|---|
| `HeatResult` | `temperature (n,)` K | `gradient (e,3)` K/m, `heat_flux (e,3)` W/m² | `reactions (n,)`, `power_balance` W |
| `ElasticResult` | `displacement (n,3)` m | `strain (e,6)` 1, `stress (e,6)` Pa, `element_energy (e,)` J | `reactions (n,3)` N, `force_balance (3,)` N, `moment_balance (3,)` N m |
| `ElectrostaticResult` | `potential (n,)` V | `electric_field (e,3)` V/m, `electric_displacement (e,3)` C/m², `element_energy (e,)` J | `reactions (n,)`, `charge_balance` C |
| `PiezoelectricResult` | `displacement (n,3)` m, `potential (n,)` V | `strain (e,6)`, `stress (e,6)` Pa, `electric_field (e,3)` V/m, `electric_displacement (e,3)` C/m² | `reactions_u (n,3)` N, `reactions_phi (n,)` C, Kraft-, Moment- und Ladungsbilanz |
| `ThermopiezoelectricResult` | `heat` und `piezoelectric` wie oben | Je Teilresultat unverändert | Je Teilresultat getrennt |

Für die einzelnen Felder ist `n = len(result.node_ids)` und
`e = mesh.n_tetra4`. Die
mechanische Voigt-Reihenfolge lautet `xx, yy, zz, yz, xz, xy`; Dehnungen
verwenden technische Schubkomponenten. `E = -grad(phi)` und Wärmefluss
`q = -k grad(T)`. Äußerer Wärmestrom ist positiv nach außen;
elektrische Flächenladung meint `D · n_out`. Die Ladung einer angrenzenden
idealen Elektrode hat das umgekehrte Vorzeichen. Weitere Vorzeichen und
Gleichungen stehen in [Doc/02](02_MATHEMATICAL_MODEL.md).

Die direkte Solver-API erwartet SI-Werte ohne implizite Umrechnung. Das
versionierte [TOML-/Python-Dokument](21_MODEL_SCHEMA.md) erlaubt die dort
explizit aufgezählten Eingabeeinheiten und konvertiert sie vor dem Solve in
SI. Die `Model`-API prüft Materialregionen, Randbedingungskonflikte und
Nullräume vor der Assembly. Direkte Solver prüfen ihre jeweiligen
Voraussetzungen selbst. `ModelResult.field_result` ist derselbe Ergebnistyp
wie beim passenden direkten Solver; nur die modellweite Provenienz kommt
hinzu. Der kombinierte Modus gibt ein Paar aus `heat` und `piezoelectric`
zurück. Gemeinsame Fälle müssen dieselben Rohfelder und physikalischen
Bilanzen liefern; [der Vertragstest](../tests/unit/test_public_api_contract.py)
vergleicht alle fünf `Model.physics`-Modi mit den direkten Solvern.

## Fehler- und Ausgabe-Vertrag

Ungültige Netz-/Material-/Randdaten lösen `ValidationError` oder eine
spezifische Unterklasse aus; Aufrufer können `error.field` auswerten.
`Model` übersetzt Schema- und Zuordnungsfehler in `ModelSchemaError` mit
`path`, während physikalische Solverfehler als ihre diagnostischen
Unterklassen durchgereicht werden. Fehlende Temperaturreferenz führt zu
`ThermalNullspaceError`, fehlendes elektrisches Bezugspotential zu
`ElectricGaugeError`; unzureichende mechanische Lagerung führt zu
`ValidationError(field="supports")`. Das sind Fehler, keine automatische
Gauge-Wahl. Für Dateien gelten zusätzlich übliche I/O-Fehler;
`export_result_bundle` überschreibt keine vorhandenen Ziele.

Das [Ergebnisbündel](24_RESULT_BUNDLE.md) hält Rohfelder, Einheiten,
Kennwerte und Provenienz in `result.vtu` und `result.json` fest. Seine
`bundle_schema_version`, die `schema_version` des Model-Dokuments und die Paketversion
sind getrennte Verträge. Ein Hash belegt identische Eingaben, ersetzt aber
keinen Bilanz- oder Konvergenznachweis.

## Migration ab v0.3 und spätere Änderungen

| Anwendung ab v0.3.0 | Weg zum vorgesehenen v1.0-Vertrag |
|---|---|
| `Mesh`, Materialobjekte, `Dirichlet` und `solve_heat`/`solve_elasticity`/`solve_electrostatics` | Direkte Aufrufe bleiben möglich. Eingaben bleiben SI; DOF-Indizes weiter mit der feldspezifischen `*_dofs`-Funktion bilden. |
| Die drei veröffentlichten Referenzbeispiele | Weiter als direkte Solver-Beispiele nutzbar; `Model` ist ein zusätzlicher deklarativer Weg. |
| Wiederholbare Konfiguration als Python-Skript | Auf `Model(mesh=..., document=...)` oder TOML mit `schema_version = 1` umstellen. Materialkatalog, `regions`, Einheiten und Randbedingungen explizit angeben; keine automatische Konversion alter Skripte behaupten. |
| Eigenes VTU-Postprocessing | Für maschinenlesbare Kennwerte und Provenienz `export_result_bundle` nutzen; rohe Ergebnisarrays bleiben für eigenes Postprocessing zugänglich. |

Bis v1.0.0 können sich noch nicht freigegebene Oberflächen ändern; jede
Änderung an bereits veröffentlichtem Verhalten erhält Release Notes und
Kompatibilitätstests. Ab v1.0.0 gilt die
[SemVer-Regel](09_RELEASE_PROCESS.md): eine inkompatible Änderung der
stabilen Oberfläche benötigt eine neue Major-Version. Eine geplante
Deprecation wird mit Warnung, dokumentiertem Ersatz, Migrationsbeispiel und
mindestens einem Minor-Release Übergangszeit eingeführt; Entfernung folgt
erst in einer Major-Version. Schemadateien werden nie stillschweigend
migriert: ein inkompatibles Format erhält eine neue `schema_version` samt
expliziter getesteter Migration. Diese Regeln gelten auch für
Einheiten-/Vorzeichenänderungen und die Struktur versionierter Ergebnisbündel.
