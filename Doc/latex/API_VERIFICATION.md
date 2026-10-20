# API-/Entwicklerband: Quellen und Kandidatenprüfung (#209)

Stand: integrierter M8-Umfang und M9-Unsicherheitsworkflow #64, Basis `267adb1`,
Paketversion 1.3.0; Dokumentationskandidat, keine M9-Releasefreigabe.

## Inventar und Abgrenzung

`api-inventory.json` inventarisiert 68 Python-Quellmodule und 197 lokal definierte
öffentliche Klassen/Funktionen, außerdem Datentypen/Konstanten und ausdrücklich
reexportierte Namen. Methoden/Konstruktoren/Properties, Dataclass-Felder,
Typannotationen, Defaults, Quell-Docstrings und direkte Raise-Stellen sind
extrahiert. Dataclass-Konstruktoren ergeben sich aus den aufgeführten Feldern;
vererbte Methoden folgen den aufgeführten Basistypen. Der Klassen-/Modulimport
selbst wird nicht ausgeführt: keine NumPy-/PETSc-/Matplotlib-/Dateidienstabhängigkeit.

`api-map.json` ordnet alle 68 Module vollständig ein und nennt je Modul
manuellen Fachvertrag, Beispiele und rationale Stabilitätsabgrenzung.
Stabile Einstiegspunkte folgen Doc/29 und den additiven Verträgen Doc/33/34,
Doc/59–62 und Doc/67. Forschungs-/interne Einträge erhalten durch ihre Aufnahme
keine neue Stabilitätszusage. Incidental Imports aus fremden Bibliotheken und
Unterstrich-Helfer sind ausgeschlossen, mit begründeter allgemeiner Regel.
Die alleinige Erreichbarkeit einer Funktion macht sie nicht zur stabilen API.

Die Referenz verbindet dieses Inventar mit manuellen Fachverträgen für
SI/Shapes/Reihenfolge, Preconditions, Mutabilität, Fehler/Nullräume, Rückgaben,
Model-/I/O-/Zeitserien-/Studien-/Surrogat-/UQ-Formate, Migration und Pflege.
Reexports verweisen auf dieselbe quellenbasierte Definition. Die eigentliche
Stabilitätszusage gilt am kanonischen Importpfad; tiefe Implementationspfade
werden nicht dadurch zusätzlich stabil. Direkte Raise-Stellen sind keine
vollständige Liste transitiver Adapter-/I/O-/Callbackfehler.

## Generierung und Driftprüfung

```sh
python scripts/api_contracts.py --all
python scripts/api_contracts.py --all --check
PYTHONPATH=src python -m pytest tests/unit/test_api_contracts.py
```

Generieren erst nach fachlichem Review von Quell-/Contractänderung; es aktualisiert
JSON und `chapters/api/inventory.tex`, keine Produktions-API. `--check` vergleicht
aktuellen AST mit JSON, prüft vollständige Modulzuordnung und existierende
Vertrags-/Beispieldateien und vergleicht das exakt erzeugte TeX. Der Handbook-
Build ruft dieselbe Prüfung vor PDF-Erzeugung auf. Das Quellenarchiv enthält die
prüfbaren Snapshots und generierte Referenz, benötigt zum erneuten PDF-Build
keinen Codeimport. Neue Quellmodule oder Exportänderungen bleiben bis zu
Inventar-/Zuordnungs-/Fachvertragsnachführung sichtbar rot.

TDD-Red/Green: zuerst fehlende Vollinventar-/Zuordnungsfunktionen, anschließend
fehlendes generiertes Kapitel. Tests prüfen echte Signatur-/Dataclass-Feldänderung,
fehlendes Reexportmodul, fehlende Modulabdeckung und die exakte TeX-Synchronisation.
Keine flächendeckenden neuen Solvertests für unveränderte Operatoren.

## Ausführbarer Entwicklungsbeitrag

`examples/developer_extension.py` im Handbucharchiv implementiert einen lokalen
linearen Wärmequellen-Callback Q=2+3x mit SI-Formvertrag und ungültigem Gradiententest.
Er verwendet ausdrücklich die Forschungs-Lastgrenze, erweitert keine produktive
API. Einheitstetraeder: unabhängige baryzentrische Integrale V/20 (offdiagonal)
und V/10 (diagonal) ergeben Lasten [13/120,2/15,13/120,13/120] W, Summe 11/24 W.
Absolute Toleranz 2e-15 W aus Größenordnung/wenigen Quadraturoperationen.

```sh
python -I Doc/latex/examples/developer_extension.py
PYTHONPATH=src python -m pytest tests/unit/test_developer_handbook_example.py
```

Der Regressionstest wurde zuerst am fehlenden Callback rot ausgeführt; anschließend
bestanden Integralreferenz und ungültiger Gradient. Die im Benutzerhandbuch
archivierten 23 vollständigen Workflows liefern weitere Minimalbeispiele;
Reader-/Schema-/Roundtriptests verwenden Produktionsparser. Signaturprüfung
ersetzt diese wissenschaftlichen Nachweise nicht.

## Lokale Abnahme am 10.10.2026

21 Inventar-/Entwicklungs-/Build-/Rezepttests bestanden; zusätzlich 192 bestehende
API-, Model-, Schema- und Readerprüfungen. Ruff, Formatprüfung (344 Dateien),
Mypy (68 Quellmodule), Inventar-Driftprüfung und Diffcheck bestanden.
Der Docker-Build erzeugte fünf PDFs und bestand alle drei absichtlichen
Negativfälle sowie den Nachbau aus dem Quellenarchiv. Der API-Band umfasst nach
der Methodenlayoutkorrektur 113 Seiten, der Gesamtband 212; kein Overfull-Befund.
Alle API-Seiten wurden als gerenderte Übersichten kontrolliert, die manuellen
Fachverträge zusätzlich in größeren Seitenbildern. Dabei gefundene Rahmen-
Überstände an Methoden wurden korrigiert; der UQ-Kurzvertrag nutzt dieselbe
getrennte Überschrift-/Signaturdarstellung. Die leeren Abbildungs-/Tabellenlisten
der Bandvorlage werden in der gemeinsamen Schlussredaktion #211 behandelt.

Das tatsächliche Quellenarchiv wurde separat entpackt. Sein Entwicklungsbeispiel
bestand mit `python -I` außerhalb des Checkouts gegen das zuvor frisch installierte
M9-Wheel, ohne `PYTHONPATH`; unabhängiges Lastintegral und Summenreferenz bestanden.
Native PETSc und die vollständige Plattform-/Paketmatrix werden vor Integration
am aktuellen PR-Commit geprüft; lokale Checks behaupten keinen solchen CI-Lauf.
Kein unabhängiges Fachreview liegt vor. V&V-Band #210, Gesamtabnahme #211 und
Releasegate #66 bleiben offen.
