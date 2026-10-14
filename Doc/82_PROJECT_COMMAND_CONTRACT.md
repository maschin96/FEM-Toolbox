# Projekt- und Befehlsvertrag für M16 (#255)

Status: entschiedener Entwurf vom 14.10.2026, auf dem stabilen v1.4.0-Umfang.
Die Befehle dieses Dokuments werden erst in #256 implementiert. Aktuell gelten
die Python-APIs und Beispielskripte. #257 ergänzt Ergebnisdiagnose, #258 die
unabhängige Nutzerabnahme und #259 die Veröffentlichung. Eine vorhandene CLI,
unabhängiges Fachreview oder experimentelle Validierung wird nicht behauptet.

## Projektstruktur und Eingaben

Ein Projekt ist ein gewöhnliches Verzeichnis mit den vorhandenen
[Model-TOMLs](21_MODEL_SCHEMA.md), Gmsh-Netzen und optional einem
[Study-TOML](59_PARAMETER_STUDIES.md). Ein zusätzliches Projekt-/Modellformat
ist nicht erforderlich. Die folgenden Namen sind Beispiele, keine Pflichtnamen:

```text
project/
  model.toml
  mesh.msh
  study.toml
  sampling.json
  quantities.json
  numerical-error-bounds.json
  results/
```

Model-Schema 1/2/3, Study-Schema 1, Ergebnisbündel 1/2, Tabellen-Schema 1 und
Sampling-/Unsicherheits-Schema 1 bleiben maßgeblich. Keine automatische Migration
unbekannter Versionen. Modellphysik, SI-Einheiten und Solverparameter bleiben in
den bestehenden Dokumenten; die CLI ist ausschließlich ein Adapter.

CLI-Dateiargumente und `--output` werden relativ zum aktuellen Arbeitsverzeichnis
aufgelöst. `mesh.path` und transiente CSV-Profile bleiben relativ zum
Model-TOML; das Studienmodell relativ zum Study-TOML. Model-/Studienpfade dürfen
ihren jeweiligen Eingabebaum auch über Symlinks nicht verlassen. Eine außerhalb
des Arbeitsverzeichnisses liegende Eingabedatei ist zulässig; das ändert ihren
eigenen Auflösungsbaum nicht. Ausgaben dürfen in einem anderen Verzeichnis
liegen, aber niemals einen Eingabepfad ersetzen. Reproduktion aus einem anderen
Arbeitsverzeichnis erfordert nur angepasste CLI-Dateiargumente.

## Geplante Aufrufe und Fähigkeiten

Der Einstieg heißt vor der M14-Migration `fem-toolbox`; `python -m fem_toolbox`
führt denselben Parser aus. M14 ergänzt `femantis`/`python -m femantis` gemäß
[ADR-053](77_FEMANTIS_MIGRATION_CONTRACT.md), ohne zweite Implementierung.
`--help` und `<command> --help` zeigen Zweck, Pflichtargumente, unterstützte
Physik, Pfade, Schreibverhalten und Exitstatus. Ein optionales Backend-Argument
wird nicht angeboten: Model/TOML und Studien verwenden heute SciPy.

| Aufruf | Eingabe und vorhandener Adapter | M16-Erstumfang |
|---|---|---|
| `fem-toolbox check model model.toml` | `Model.from_file` | Die unten aufgeführten stationären Model-Physiken; Import, Material-/Schema-/Randprüfung ohne Solve |
| `fem-toolbox check study study.toml` | `Study.from_file(...).resolve()` und `Model` je aufgelöstem Fall | Stationäre kartesische Studie; jede Fall-ID und fallbezogene Preflightfehler sichtbar |
| `fem-toolbox run model.toml --output results/run-001` | `Model.from_file(...).solve()`, `export_result_bundle` | Stationärer Einzelfall mit vorhandenen SI-Feldern, Kennwerten und Bilanzen |
| `fem-toolbox study study.toml --output results/study-001 [--resume] [--max-cases N]` | `Study.run` | Stationäre Studie mit bestehenden IDs, Versuchen und Fingerprints |
| `fem-toolbox uq --model model.toml --sampling sampling.json --quantities quantities.json --numerical-error-bounds numerical-error-bounds.json --output results/uq-001 [--resume]` | `SamplingPlan`, `Quantity`, `UncertaintyStudy.run` | Gemeinsame skalare stationäre Stichproben, ausdrückliche Gesetze und Fehlerbudgets |
| `fem-toolbox report results/study-001/catalog.json --quantities quantities.json --output results/table-001` | `export_study_table`, anschließend `read_study_table` | Bestehender Studienbericht ohne Solve; keine allgemeine Vergleichsansicht |

Stationäre Model-Physiken: `heat`, `elasticity`, `electrostatics`,
`piezoelectric`, `thermopiezoelectric`; Tet4 oder geradkantige Tet10 gemäß
Supportmatrix. `thermomechanics` ist eine Python-Solverfähigkeit und kein
zusätzliches Model-TOML-Physiklabel. Transiente Wärme, Modal-/Newmark-/Harmonik
sind vorhandene Model-Fähigkeiten, benötigen jedoch separate Serienexporte und
Preflightpfade. Sie werden in der ersten CLI-Version vor dem Solve mit konkreter
Umfangsmeldung abgelehnt; vorhandene Python-Rezepte bleiben nutzbar. Ebenso
ausgeschlossen: PETSc-TOML, Adaptivität aus unveröffentlichtem M10, CFD,
Topologieoptimierung, parallele Studien und experimentelle Validierungsbehauptungen.

UQ übernimmt exakt die Parameter-/Verteilungsfelder von
`SamplingPlan.to_document()` (Schema 1), einschließlich Begründung, Seed,
Größe, Reihenfolge und Korrelation. Die dort gespeicherten NumPy-/SciPy-Versionen
beschreiben die Quelle der Spezifikation; die tatsächlich ausführende Umgebung
wird neu im erzeugten `sampling.json` protokolliert. Unbekannte Schema-/Methoden-
und Generatorwerte werden abgelehnt; keine stille Interpretation als aktuelle
Methode. #256 muss den expliziten JSON-Reader gegen diese API verifizieren.

`quantities.json` ist ein JSON-Array von benannten Argumentobjekten des
bestehenden `Quantity`-Konstruktors; erforderliche Felder sind `name`, `source`,
`field`, `unit`, `definition`, optionale Felder entsprechen dessen Signatur.
Unbekannte Schlüssel werden abgelehnt. Beispiel:

```json
[{"name":"Tmax","source":"point","field":"temperature_K","unit":"K",
  "definition":"maximum active-node temperature","reduction":"max"}]
```

`numerical-error-bounds.json` enthält genau einen Eintrag je Quantity, zum
Beispiel `{"Tmax": null}`. `null` bedeutet unbekannten numerischen Fehler,
niemals null Fehler. Endliche nichtnegative Zahlen gelten in der Quantity-Einheit.
Die Dateien sind Adapterargumente für vorhandene APIs und kein neues FEM-Modell.
Sensitivitätsoptionen bleiben zunächst in der Python-API: Die erste CLI muss
keine zusätzlichen Zentrum-/Schrittverträge erfinden. Gesetz- und
Quantity-Validierung einschließlich aller Stichproben-Zielpfade erfolgt vor
dem Erstellen des Ausgabeziels. Kein ausführbares Python in Eingabedateien.

## Preflight, Fehler und Ausgabe

`check model` konstruiert das Model und löst kein Gleichungssystem. Die bestehende
stationäre Rand-/Nullraumprüfung wird im Konstruktor ausgeführt. Ein erfolgreicher
Preflight garantiert weder Konditionierung noch Bilanz, Netzkonvergenz oder
physikalische Eignung. `check study` muss zusätzlich zur Strukturauflösung jeden
Fall als Model konstruieren: `Study.resolve()` allein prüft nicht sämtliche
physikalischen Fallfehler. Fehlfälle werden mit ID und konkretem Fehler gezeigt;
das Originalmanifest und die Quellen bleiben unverändert. Für UQ führt der
Ausführungsadapter dieselbe Vorprüfung wie `UncertaintyStudy.run` aus; eine
separate UQ-Check-Option ist kein Erstumfang.

Alle schreibenden Aufrufe erfordern ein ausdrückliches `--output`. `run` und
`report` erzeugen ausschließlich neue Verzeichnisse, auch ein vorhandenes leeres
Ziel wird abgelehnt. Einzelfallausgabe: `result.json` und `result.vtu` gemäß
[Bündelvertrag](24_RESULT_BUNDLE.md). Bericht: `values.csv`, `metadata.json`,
`report.md` gemäß [Tabellenvertrag](60_STUDY_REPORTS.md). Kein `--force`, kein
automatisches Löschen, kein atomarer Gesamtverzeichnisanspruch; unvollständige
Ausgabe nach I/O-Abbruch bleibt zur Diagnose erhalten und ist kein Erfolg.

`study` und `uq` verlangen standardmäßig ebenfalls ein neues Ziel. Nur
`--resume` erlaubt das vorhandene, zum selben Workflow gehörende Verzeichnis.
Fehlt es bei `--resume`, ist dies ein Bedienfehler. Die CLI prüft vor Mutation
die Workflowkennung (`catalog.json` bzw. `sampling.json`) und reservierte
Dateien/Unterverzeichnisse; ein fremdes oder unlesbares Ziel wird abgelehnt.
Die bestehende API validiert danach Schema, Fingerprints, Bündelhashes und
Provenienz. Eingabeänderungen invalidieren Studienfälle mit erhaltener Historie;
geänderte Sampling-Spezifikationen werden im UQ-Ziel abgelehnt. Neue NumPy-/SciPy-
Versionen können deswegen ein neues UQ-Ziel erfordern. Ein fremder Ordner wird
nicht durch ein frisches Manifest übernommen.

`.study.lock` und `.uncertainty.lock` werden niemals automatisch gestohlen.
Nach hartem Abbruch erst laufende Writer ausschließen, dann die betroffene
Lockdatei manuell entfernen; Details in den Fachverträgen. `--max-cases N`
bedeutet nichtnegatives Integer-Budget für neue Studienversuche, nicht für
wiederverwendete Erfolge. `0` schreibt nur den Katalog. Ein nicht vollständiger
Lauf ist sichtbar und bekommt Status 4. Für UQ existiert kein solches Budget.

Kein neues Provenienzschema: Model-`RunProvenance`, Bündel-Schemas, Studien-
Fingerprints/-Versionen und UQ-/Tabellen-Schemas bleiben die Nachweise. Uhrzeit
und absolute Ausgabepfade ändern deren deterministischen Kern nicht. Jede
Erfolgsmeldung nennt Artefaktpfad, vorhandene Schema-/Toolbox-Version und Umfang.
Die CLI darf fehlende Einheiten/Bilanzen oder unbekannte Fehlerbudgets nicht
durch Erfolgsvorgaben ersetzen. Ein neuer CLI-eigener persistierter Receipt
würde einen gesonderten Schemavertrag benötigen und ist hier nicht vorgesehen.

| Status | Auslöser und Ausgabe auf stderr |
|---|---|
| 0 | Preflight ohne Befund; Einzellauf samt Export; vollständige Studie; vollständige UQ ohne fehlgeschlagene angeforderte Sensitivität; intakter Bericht ohne fehlende/ungültige Werte und ohne fehlgeschlagenen gewählten Bilanztest |
| 2 | Parserfehler; `ModelSchemaError`/`ValidationError`; `ValueError` aus Gesetzen, Quantity, Schema-/Hashprüfungen oder Umfangsprüfung. Konkreter Feldpfad/Fehlertext, keine stille Reparatur |
| 3 | Dateisystem-/Mesh-I/O, `OSError`, `MeshIOError`, belegtes Ziel oder aktive Lockdatei (`FileExistsError`); beim Bündel vorhandener Pfad als konkrete Ausgabevorprüfung |
| 4 | Fertig geschriebener, aber unvollständiger Studien-/UQ-/Berichtsstand; alle betroffenen IDs und Gründe bleiben sichtbar. Auch fehlgeschlagener ausdrücklich gewählter Bilanztest |
| 5 | Nicht als Eingabefehler klassifizierter Solvefehler, etwa `RuntimeError`; keine fertige Einzelfallausgabe. In Studien gespeicherte Fallfehler ergeben dagegen 4 |
| 130 | `KeyboardInterrupt`; soweit die API es zusagt, Unterbrechungsstatus persistieren und Lock freigeben |

Zuerst spezifische Exceptions, dann allgemeinere Basisklassen behandeln:
`ModelSchemaError` ist auch `ValueError`; `FileExistsError` auch `OSError`.
Keine pauschale `except Exception: return 0`. Unerwartete interne Exceptions
bleiben diagnostizierbare Fehler mit erfolglosem Prozessstatus; sie werden
nicht als Bedienfehler etikettiert. `--help` endet mit 0 und schreibt nichts.

## Vorprüfung an drei vorhandenen Rezepten

Der Nachweis für #255 benutzt die bestehenden APIs/Skripte, weil die geplante
CLI noch fehlt. Aufrufe aus dem Repository-Root mit installiertem Paket:

```sh
python examples/parameter_study/run.py --output /tmp/m16-study --max-cases 1
python examples/parameter_study/run.py --output /tmp/m16-study
python examples/parameter_study/report.py --study /tmp/m16-study --output /tmp/m16-table
python examples/uncertainty_heat/run.py --output /tmp/m16-uq
```

Der Einzelfall wird aus `examples/parameter_study/model.toml` mit
`Model.from_file(...).solve()` und `export_result_bundle` gerechnet. Seine
freie Temperatur beträgt 301 K: am Einheitstetraeder gilt unabhängig
`K11=k/6`, `f1=Q/24`, `T1=300+Q/(4k)`, hier Q=8 W/m³, k=2 W/(m K).
Toleranzen vor Lauf: 1e-11 K für diese diskrete Temperatur, 1e-11 W für die
Leistungsbilanz; keine Aussage zum Kontinuumsfehler.

Studie: 301 und 300,5 K für k=2/4; nach Budget 1 bleibt ein Fall pending,
Fortsetzung übernimmt den ersten Versuch. Intakter Tabellenleser bestätigt
CSV-/Metadatenkonsistenz. UQ: 128 gemeinsame Stichproben, Seed 64,
`Tmax=T0+Q/8`; 1e-11 K diskreter Fehler, Sensitivitäten [1/8,1] mit 2e-10
in den jeweiligen SI-Verhältnissen. Die synthetischen Gesetze und vorher
begründete statistische Grenze stehen im bestehenden UQ-Rezept; das
Kontinuumsfehlerbudget bleibt `null`.

Zusätzlich vor Integration prüfen: dasselbe Model aus fremdem Arbeitsverzeichnis,
Model-Netzpfad außerhalb des Eingabebaums, erneuter Bündel-/Tabellenexport ins
bestehende Ziel, unbekannte Model-PETSc-Option und transiente Studie. Erwartet
sind bestehende `ModelSchemaError`-/`ValueError`-Fehler, keine Umschreibung von
Eingaben oder Ausgaben. Die neuen CLI-Ziel-/Resume-/Exitstatusregeln sind
Abnahmekriterien für #256 und werden hier nicht als implementiert bezeichnet.

Ausführungsergebnisse und CI-Commit werden im Issue-PR protokolliert. Das ist
Vertragsprüfung durch den Implementierer; unabhängige Nutzerabnahme bleibt #258.
