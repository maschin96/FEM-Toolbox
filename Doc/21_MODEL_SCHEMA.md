# Versioniertes Modellformat und öffentliche Model-API (#25)

Status: Implementiert und durch Adaptertests geprüft. Die physikalischen Solver
haben ihre eigenständigen Verifikationsnachweise in Doc/10–20. Das Format
beschreibt jeweils **einen** stationären Modus: `heat`, `elasticity`,
`electrostatics` oder `piezoelectric`. Der Ablauf `T → (u, phi)` ist #26.

## Vollständiger Wärmefall

Die folgende Datei `case.toml` liegt neben einem Gmsh-Netz `mesh.msh` mit dem
Volumensatz `ceramic` und fünf Netzknoten (beispielsweise
`tests/fixtures/overlapping_groups.msh`). Aus dem Repository-Root:

```toml
schema_version = 1
units = "SI"
physics = "heat"

[mesh]
path = "mesh.msh"
length_scale_to_m = 1.0

[materials.T.ceramic]
conductivity = 2.0
unit = "W/(m K)"

[regions.T]
ceramic = "ceramic"

[[boundary.dirichlet]]
field = "T"
nodes = [0, 1, 2, 3, 4]
value = 300.0
unit = "K"
```

```python
from fem_toolbox.model import Model

model = Model.from_file("case.toml")
run = model.solve()
print(run.field_result.temperature)
print(run.provenance.input_sha256)
```

`from_file` importiert nur `.toml`. Der Gmsh-Pfad muss relativ sein und innerhalb
des Modellverzeichnisses liegen. `length_scale_to_m` ist ein positiver Faktor;
alle anderen Eingaben sind explizite SI-Werte. `Model(mesh=mesh,
document=mapping)` erlaubt denselben Schema-Vertrag mit einem bereits
eingelesenen `Mesh`. `to_document()` liefert eine unabhängige Kopie für
Bearbeitung oder Roundtrip. Die Provenienz enthält Toolbox-Version,
Schema-Version, Solverkennung und einen Hash über Dokument, Netzgeometrie,
Materialzuordnung und Solver-Eingaben.

## Tabellen und Einheiten in Version 1

`materials.<field>.<name>` definiert einen Materialkatalog;
`regions.<field>.<volume_set> = "<name>"` ordnet jedem aktiven Tet4 genau ein
Material zu. Die Feldnamen sind `T`, `u`, `phi` und `piezo` für die vier Modi.
Überlappende Volumensätze müssen dieselbe Zuordnung ergeben. Lücken,
widersprüchliche Zuordnungen und ungenutzte Materialien sind Fehler.

| Material | Werte | Einheit |
|---|---|---|
| `T` | `conductivity`: Skalar oder 3×3 | `W/(m K)` |
| `u` | `young` und `poisson` oder `stiffness`: 6×6 | `Pa` |
| `phi` | `permittivity`: Skalar oder 3×3 | `F/m` |
| `piezo` | `stiffness`: 6×6, `stress_coupling`: 3×6, `permittivity`: 3×3; optional `axes`: 3×3 | `[materials.piezo.<name>.units]` mit `stiffness = "Pa"`, `stress_coupling = "C/m^2"`, `permittivity = "F/m"` |

Randbedingungen stehen als Arrays von Tabellen unter `boundary`. Ein
Dirichlet-Eintrag benötigt `field`, `value`, `unit` und entweder `nodes`
(nullbasierte Netzknoten) oder `set` (benannter Triangle3-Satz). Mechanische
Teilkomponenten können als `components = ["x", "y"]` gewählt werden. Der Wert
ist ein Skalar, ein Wert pro Komponente oder ein Wert pro ausgewähltem DOF.

| Randbedingung | Modi | Werte und Einheiten |
|---|---|---|
| `dirichlet` | alle | `T` in `K`, `u` in `m`, `phi` in `V` |
| `flux` | Wärme | `set`, `value` in `W/m^2`; positiv nach außen |
| `convection` | Wärme | `set`, `coefficient` in `W/(m^2 K)`, `ambient_temperature` in `K` und `ambient_unit = "K"` |
| `traction` | Elastizität/Piezo | `set`, Vektor `value` in `N/m^2` |
| `surface_charge` | Elektrostatik/Piezo | `set`, `value` in `C/m^2` |
| `electrodes` | Piezo | `name`, `set`, entweder `potential` mit `potential_unit = "V"` oder `charge` mit `charge_unit = "C"` |

`[loads.source]` (Wärme) verwendet `value` in `W/m^3`,
`[loads.body_force]` (Elastizität/Piezo) einen 3-Vektor in `N/m^3` und
`[loads.volume_charge]` (Elektrostatik/Piezo) `value` in `C/m^3`.
`[solver]` ist nur bei Piezo zulässig: `method = "direct"` und
`block_scaling = true/false`; sonst gelten die dokumentierten Solver-Defaults.
Unbekannte Schlüssel und falsche Einheiten führen mit Feldpfad zu
`ModelSchemaError`. Physikalisch unzulässige oder singuläre Fälle liefern die
diagnostischen Fehler der zugrunde liegenden Solver.

## Versions- und Migrationsregel

Die ganzzahlige `schema_version = 1` ist verpflichtend. Die Toolbox akzeptiert
keine andere Schemaversion und migriert keine Eingaben stillschweigend. Eine
inkompatible Formatänderung erhöht die Schemaversion; dafür wird eine explizite,
getestete Migration mit dokumentierter Feldabbildung und Einheitenkonversion
bereitgestellt. Bestehende v1-Dateien bleiben unverändert interpretierbar.
Die Paketversion in der Laufprovenienz ist unabhängig von der Schemaversion.
