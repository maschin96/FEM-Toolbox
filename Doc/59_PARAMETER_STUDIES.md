# Stationäre Parameterstudien (#59)

Die Python-API `fem_toolbox.studies.Study` führt lokale, sequenzielle Studien
über die vorhandene `Model`-API aus. Schema 1 unterstützt `heat`, `elasticity`,
`electrostatics`, `piezoelectric` und `thermopiezoelectric` auf den bereits
unterstützten Tet4-/geradkantigen Tet10-Netzen. Sie ändert keine Feldoperatoren.
Model/TOML verwendet weiterhin SciPy; PETSc ist kein Studienparameter.

## Einstieg und Manifest

```python
from fem_toolbox.studies import Study

study = Study.from_file("examples/parameter_study/study.toml")
cases = study.resolve()  # auflösen, noch nicht rechnen
catalog = study.run("study_results")  # rechnen oder geprüfte Ergebnisse übernehmen
```

Alternativ: `Study(manifest_mapping, root=directory)`. Die Konstruktion kopiert
das Manifest. `from_file` liest UTF-8-TOML; `root` ist das Manifestverzeichnis.
Bei jeder Auflösung werden Basismodell und Netze erneut gelesen. Für Änderungen
am Studienmanifest selbst eine neue `Study` konstruieren. Rückgabewerte sind
unabhängige JSON-kompatible Kopien, keine mutable interne Konfiguration.

```toml
study_schema_version = 1
model = "model.toml"

[[parameters]]
name = "conductivity"
target = "materials.T.thermal.conductivity"
unit = "W/(m K)"
values = [2.0, 4.0]
```

Genau diese drei Wurzelschlüssel sind erlaubt. `parameters = []` bedeutet einen
unveränderten Einzelfall. Jede Parametertabelle enthält genau `name`, `target`,
`unit`, `values`. Namen sind eindeutige ASCII-Bezeichner; Zielpfade sind eindeutig.
Wertelisten sind endlich, nicht leer und ohne identische Einträge. Das kartesische
Produkt folgt der Tabellenreihenfolge; der letzte Parameter variiert am schnellsten.
Es gibt keine implizite Zufallsziehung, Zip-Semantik oder Skriptauswertung.

Zielpfade adressieren **vorhandene ganze Eigenschaften**, beispielsweise
`boundary.dirichlet.0.value`. Listenindizes sind nullbasiert, ohne Vorzeichen
oder führende Nullen. Vektoren/Tensoren werden als Ganzes mit identischer Form
ersetzt, nicht komponentenweise. Namen im Modell mit Punkten sind über diese
Pfadsyntax nicht adressierbar. Parameter können Materialwerte und Achsen,
Last-/Randwerte, Elektrodenladung/-potential, `mesh.path`, `mesh.length_scale_to_m`
und vorhandene `solver.method`, `solver.block_scaling`,
`solver.electric_thermal_strain` ersetzen. Neue Kataloge, Auswahlmengen,
Feldtypen, Topologien und Schema-/Einheitenschlüssel werden nicht parametrisiert.
Solveroptionen bei Bedarf zunächst ausdrücklich ins Basismodell schreiben.

`unit` muss zur SI-Einheit der vorhandenen Eigenschaft passen. Dimensionslose
Größen, Solveroptionen und Netzpfade verwenden `"1"`. Nicht-SI-Parameterziele
wie GPa/degC/relative Permittivität müssen vor der Studie im Basismodell in SI
umgerechnet werden; unveränderte Modelleigenschaften behalten die bestehenden
Model-Konversionen. Skalare, Boolesche Werte, Strings und Arrays müssen den
vorhandenen Typ/die Form einhalten. Physikalische Zulässigkeit prüft `Model`
pro Fall, sodass beispielsweise negative Leitfähigkeit einen sichtbaren Fehlfall
erzeugt und andere Fälle weiterlaufen.

Basismodellpfade bleiben innerhalb des Manifestverzeichnisses; Gmsh-Pfade bleiben
innerhalb des Basismodellverzeichnisses. Absolute Pfade, Ausbrüche und entsprechende
Symlinks werden abgelehnt. Netzvarianten sind vorbereitete `.msh`-Dateien.
Eingaben während eines Laufs nicht verändern: Das aufgelöste Modelldokument ist
ein Snapshot; vor dem Import werden Netzbytes erneut geprüft und aus einem
privaten Snapshot eingelesen.

## Identität und Wiederverwendung

`case_id` ist der vollständige SHA-256 der kanonischen benannten Parameterwerte,
Ziele und Einheiten. JSON hat sortierte Schlüssel, kompakte Separatoren und keine
NaN/Inf. Parameterreihenfolge beeinflusst die Ausführungsreihenfolge, nicht die
Fallidentität. Keine Parameter bedeuten genau eine konstante Fall-ID.

`fingerprints` enthält:

| Schlüssel | Inhalt |
|---|---|
| `model` | Vollständiges aufgelöstes Modelldokument ohne Netz-Dateipfad, einschließlich Netzskala, Materialien, Regionen, Lasten und Randbedingungen |
| `mesh` | SHA-256 der tatsächlichen Gmsh-Dateibytes |
| `solver` | Physik, explizite Solveroptionen und Versionsvertrag |
| `input` | Zusammenfassung der drei Hashes und aller Versionsangaben |

`versions` nennt Toolbox, NumPy, SciPy, meshio, Provenienzschema und
Studienrevision. Defaults werden durch Toolboxversion/Studienrevision gebunden;
Änderungen an ihrer Interpretation benötigen eine Revisionsänderung. Der
Versionsvertrag ist konservativ: auch ein Dependency-Update oder eine reine
Gmsh-Dateiänderung kann Neuberechnung auslösen. Er ersetzt keinen vollständigen
Umgebungs-Lock und garantiert keine bitidentische Numerik zwischen Plattformen.
Unveröffentlichte Solveränderungen erfordern neue Version/Studienrevision oder
ein neues Ausgabeziel. Absolute Ablagepfade und Uhrzeiten gehen nicht in Hashes
ein. Ein explizit parametrisierter relativer Netzpfad ist Teil des Parameterwerts
und damit der Fall-ID; sein Speicherort außerhalb des Projektbaums nicht.

Geänderte Basiseingaben behalten die Fall-ID, invalidieren aber den Inputhash.
Damit ist dieselbe Parameterkombination vor/nach Änderung nachvollziehbar.
Erfolgreiche Ergebnisse werden nur übernommen, wenn Fingerprints **und**
SHA-256 beider Bündeldateien stimmen, das Bündelschema unterstützt ist und
Bündelreferenz, VTU-Hash sowie ursprüngliche Model-Provenienz konsistent sind.
Zusätzlich müssen beide Dateien im erwarteten Fall-/Versuchsverzeichnis liegen.
Die Bündelprovenienz enthält unter `solver.parameters.study_case` die unabhängig
aus dem Fall geprüfte Parameteridentität, den Hash des aufgelösten Dokuments,
Fingerprints und Versionen. So werden auch fremde, intern konsistente Bündel
erkannt, die in den eigenen Versuchspfad kopiert wurden (#219). Die Prüfung
benötigt keinen Solve und keine ursprünglichen Eingabedateien; vollständig
kopierte Studienverzeichnisse bleiben wiederverwendbar. Frühere Bündel ohne
diese Zuordnung werden konservativ als ungültig behandelt und beim nächsten
Lauf mit erhaltener Historie neu berechnet; Berichte zeigen fehlende Kennwerte.
Dies ist Integritätsprüfung gegen unbeabsichtigte Änderungen, keine kryptografische
Authentifizierung eines fremd manipulierten Katalogs.

## Persistenter Katalog und Fortsetzung

`run(output, max_cases=None)` legt das Ausgabeziel an und liefert den gespeicherten
`catalog.json` zurück. `max_cases` begrenzt neue Versuche, nicht bereits gültige
Ergebnisse; `0` schreibt nur den aufgelösten Katalog. Wiederholter unveränderter
Lauf berechnet keinen erfolgreichen Fall neu und liefert denselben Katalog.

Der Katalog enthält `study_schema_version = 1`, geordnete `cases` und
`retired_cases` für aus dem aktuellen Manifest entfernte Fälle. Jeder Fall enthält
`case_id`, `parameters` (Wert/Einheit/Ziel), `document`, `fingerprints`, `versions`,
`status`, `attempts`, `error`, `result`, `attempt_directory`, `history`.
`history` hält Invalidierungsgrund und vorherige Fingerprints/Ergebnisreferenz;
es ist kein vollständiges Ereignisprotokoll aller Zustandsübergänge.
`retired_cases` ist ein unveränderter Archiv-Snapshot, keine Aussage aktueller
Wiederverwendbarkeit. Erneut aufgenommene Fälle werden regulär neu geprüft/gerechnet.

| Status | Bedeutung und Verhalten beim nächsten Lauf |
|---|---|
| `pending` | Noch kein Versuch; ausführen |
| `running` | Versuch begonnen; nach Prozessabbruch als unterbrochen behandeln |
| `interrupted` | KeyboardInterrupt/SystemExit oder verwaister laufender Versuch; neu versuchen |
| `failed` | Konkrete Exception mit Typ/Ursache; erneut versuchen |
| `invalid` | Eingabe geändert oder Bündel beschädigt/unvollständig; neu versuchen |
| `succeeded` | Vollständiges Bündel gespeichert; vor Wiederverwendung prüfen |

`result` ist bei nicht erfolgreichen Fällen `null`; andernfalls nennt es relative
`manifest`-/`fields`-Pfade, beide SHA-256 und `model_input_sha256`. Der unveränderte
[Ergebnisbündelvertrag](24_RESULT_BUNDLE.md) liefert Felder, Kennwerte, Einheiten
und Bilanzen an #60. Ein erfolgreicher Lauf ist kein Konvergenz-/Validierungsbeleg.

Jeder Versuch schreibt nach `cases/<case_id>/attempt-<nummer>/result.{json,vtu}`.
Existierende Versuche werden nie überschrieben. Eine unvollständige Ausgabe ohne
erfolgreichen Katalogeintrag wird nie wiederverwendet. Katalogupdates schreiben
zuerst eine temporäre Datei, synchronisieren deren Inhalt und ersetzen anschließend
den Katalog atomar. Ein beschädigter/fremdversionierter Katalog führt zum Fehler,
nicht zum stillen Neubeginn. Alte Versuche bleiben zur Diagnose erhalten; kein
automatischer Cleanup. Stromausfallgarantien jenseits der Dateisystemsemantik
sind nicht zugesagt.

Resume und Bericht validieren den vollständigen Schema-1-Fallrecord an derselben
Lesegrenze (#220): Pflichtfelder, Parameterwerte/-metadaten, endliche
Modelldokumentdaten, SHA-256-Fingerprints, Versionsangaben, relative Ergebnis-
und Versuchspfade sowie Status-/Ergebniskombinationen. Erfolgreiche Fälle
benötigen Ergebnis und Versuch; andere Zustände besitzen kein Ergebnis.
`pending` hat noch keinen Versuch, `running` keinen Fehler; `failed`,
`interrupted` und `invalid` nennen einen Fehlergrund. Historieneinträge sind
Invalidierungssnapshots mit Grund, früheren Fingerprints, optionalem Ergebnis
und Versuchszähler. Retired Records werden strukturell nach demselben Vertrag
geprüft, ohne aktuelle Bündelgültigkeit zu behaupten; wiederholte archivierte
Fall-IDs bleiben zulässig. Fachlich ungültige, aber strukturell lesbare
Modelleingaben bleiben als Fehlfälle dokumentierbar.

Vertragsverletzungen liefern `ModelSchemaError` mit konkretem Katalogfeldpfad.
Dabei bleiben Katalog und Ergebnisdateien unverändert; Resume gibt die Lockdatei
auch bei fehlgeschlagener Prüfung frei. Es gibt keine automatische Reparatur.

`.study.lock` verhindert mehrere gleichzeitige Writer. Bei normalem Abschluss,
Exceptions und Ctrl-C wird sie entfernt. Nach hartem Kill **erst prüfen, dass
kein Writer mehr läuft**, dann ausschließlich diese Lockdatei entfernen und
fortsetzen. Der Katalog erkennt verbliebenes `running`. Automatisches Lock-Stehlen
oder parallele/verteilte Ausführung sind nicht unterstützt. Katalog-I/O-Fehler
werden propagiert, fallbezogene Import-/Model-/Solve-/Exportfehler gespeichert;
KeyboardInterrupt/SystemExit werden nach Speicherung erneut ausgelöst.

## Nachweise

`tests/unit/test_studies.py` prüft Identität, Einheiten/Fehler, kartesisches Produkt,
Ortsunabhängigkeit, Material-/Rand-/Netz-/Solver-/Versionsänderungen, echte
Einzelfeld-/Kopplungsläufe, Tet10-Bündel, intakte Wiederverwendung, beschädigte JSON
und VTU, fehlende Dateien, kontrollierten Abbruch während des Exports, fehlerhafte
Materialdaten, lock-/katalogbedingte Fehler und atomare Katalogschreibfehler.
Direkte Model-Läufe werden unabhängig vom Studienexpander aufgebaut; Felder,
Kennwerte und Bilanzen werden verglichen. Rohwerte verwenden `atol=rtol=1e-12`
für die kleinen SI-Referenzfälle; Identitäten/Status und direkt serialisierte
Kennwerte werden exakt geprüft.

[Ausführbares Beispiel](../examples/parameter_study/README.md): eine freie
Temperatur am Einheitstetraeder wird gegen eine unabhängig hergeleitete diskrete
Lösung geprüft. Das ist ein Orchestrierungs-/Diskretisierungsvergleich und keine
experimentelle Validierung. Die bestehende V&V der Feldoperatoren bleibt gültig;
neue Konvergenzraten werden mit diesem Adapter nicht beansprucht.
