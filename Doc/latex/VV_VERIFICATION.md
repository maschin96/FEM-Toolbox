# Verifikationsband #210 – Prüfprotokoll

Stand: 10.10.2026, M9-Integrationsbasis `11b64c8`, Issuebranch
`codex/210-verification-handbook`. Keine experimentelle Validierungsfreigabe.

## Quellen und wissenschaftlicher Umfang

Sechs Fachkapitel unterscheiden Code-/Lösungsverifikation und experimentelle
Validierung, ordnen 17 Nachweisgruppen zu, erläutern Normen/Bilanzen/Ordnungen,
M7-Genauigkeit/Kosten, M8-Surrogatfehler und den M9-Unsicherheitsworkflow.
Alle Angaben beziehen sich auf konkrete Verträge/Tests; historische Ergebnisse
werden ausdrücklich als solche bezeichnet. #63/#65 bleiben gemäß genehmigter
Neuplanung offen außerhalb M9. Keine Messantwort wurde durch Simulation ersetzt.

`verification-map.json` verlangt jede unterstützte Physik und die Querschnitts-
Nachweise. `verification-data.json` bindet Verträge, Tests, Auswertung und eigene
Figuren durch SHA-256. Das originale M7-JSON wird bytegleich unter `evidence/`
ins Handbucharchiv aufgenommen. Eigene Auswertungen/Plots stehen unter MIT;
fremde Messdateien werden nur mit Primärquellen und Attribution im Suchbericht
referenziert, nicht im Handbuch als eigene Daten verarbeitet.

## Red–Green–Refactor und unabhängiger Datenabgleich

Vier neue Tests wurden zunächst am fehlenden Auswertungsskript rot ausgeführt.
Sie prüfen den tatsächlichen 72-Lauf-Bestand/24 Gruppen/vier unerreichten Zielen,
fehlende/duplizierte/physikalisch fehlgeschlagene Läufe, manipulierte Medianwerte,
unvollständige Physikzuordnung/fehlende Testpfade sowie Quellen-/Tabellen-/
Figurdrift. Der Audit verwendet Standardbibliothek und unabhängig neu berechnete
Statistiken/Budgetauswahlen; keine Produktionsbenchmark-Zusammenfassung wird als
ungeprüfte Referenz übernommen. Default-Build benötigt kein Matplotlib.

Der bestehende analytische korrelierte UQ-Fall wurde mit 100000 Zeilen, Seed 641
reproduziert. Unabhängig erwartete Momente: 8 K und 13 K², Gradient [2,-1].
Tatsächlich: Mittelwert 8.007392184773114 K, Varianz 13.005819137120481 K².
Vorab bestehende 6-Standardfehler-Testkriterien und 1e-12-Gradiententoleranz
bleiben unverändert. Die illustrative deterministische Angabe 0.01 K ist keine
FEM-Diskretisierungsschranke. Seed, Versionen, vollständige Präfixstatistik und
Samplebytehash stehen in `evidence/uq-linear.json`. Echte FEM-Aufrufe werden
zusätzlich durch das integrierte 128-Zeilen-Wärmerezept geprüft.

## Abnahme

Lokale Gesamtsuite: **822 bestanden, 83 optionale PETSc-Tests übersprungen**,
113.92 s. Gezielte UQ-/Audit-/Benchmarktests: 49 bestanden; nach den finalen
Quellenänderungen 16 Audit-/Buildtests bestanden. Ruff, Formatprüfung
(347 Dateien), mypy (68 Quellmodule), Quellenprüfung und Diffcheck bestanden.

Alle fünf PDFs, drei absichtlich defekte Dokumentfälle und der vollständige
Archivnachbau bestanden. Nachweisband: 27 Seiten, Gesamtband: 232 Seiten.
Keine Overfull-Warnungen. Alle 27 Nachweisseiten wurden gerendert und visuell
geprüft; wissenschaftliche Figuren sowie Formel-/Tabellenseiten zusätzlich
in größeren Ansichten. Überlange Matrixzeilen und Kandidatenbeschriftungen
wurden korrigiert. Gemeinsame Titel-/Statusseiten kennzeichnen vier
ausgearbeitete Fachbände und die noch offene Gesamtabnahme. Native PETSc- und Plattform-/Paketprüfungen
werden am aktuellen PR-Commit vor Integration ausgewertet. Lokale PETSc-Skips
zählen nicht als Erfolg. Kein unabhängiges Fachreview liegt vor; Quellenabgleich
und eigene PDF-Prüfung werden nicht als unabhängiges Review bezeichnet.

Gesamtabnahme #211 und Release #66 bleiben gesonderte Gates.
