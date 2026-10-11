# Gesamtabnahme der deutschsprachigen M9-Dokumentation (#211)

Stand: 11.10.2026; Integrationsbasis `fe045d6`, Kandidat vor Veröffentlichung.
Paketversion wird erst in #66 geändert. Keine experimentelle Validierungsfreigabe.

## Integrierte Beiträge und Nachweisgrenzen

| Beitrag | Integration | Fach-/Ausführungsnachweis |
| --- | --- | --- |
| #206 Struktur/Build | #239/#240 | VERIFICATION.md; festgelegte Toolchain, drei negative Fälle, Archivnachbau |
| #207 Theorie | #241 | THEORY_VERIFICATION.md; zehn Kapitel, unabhängige Handrechnungen und vollständige Seitenprüfung |
| #208 Anwendung | #243 | USER_VERIFICATION.md; 23 Workflows, vollständige Eingaben, Installation/Formate/Fehlerhilfe |
| #64 Unsicherheit | #244 | Doc/67_UNCERTAINTY_WORKFLOW.md; analytische Momente, Korrelation, Sensitivität, Fehler/Resume |
| #209 API/Entwicklung | #245 | API_VERIFICATION.md; 68 Module, 197 lokale öffentliche Klassen/Funktionen, Reexports/Felder/Methoden |
| #210 Nachweise | #246 | VV_VERIFICATION.md; 17 Gruppen, 72 historische Messungen, negative Befunde und Quellenhashes |

Alle genannten PRs wurden nach erfolgreicher anwendbarer CI am jeweiligen Head
integriert. #63/#65 bleiben gemäß ADR-054 offen außerhalb M9. Ihre Datengates
sind weder erfüllt noch durch synthetische Beispiele ersetzt. Eigenprüfung und
CI sind kein unabhängiges Fachreview; ein solches liegt weiterhin nicht vor.

## Redaktioneller Vertrag

Gesamtband und vier eigenständige Bände haben vollständige Pflichtkapitel.
Einheitliche SI-/Voigt-/Fluss-/Phasorkonventionen und kanonische Importpfade
wurden mit den Quellen abgeglichen. Fachliche Herleitungen stehen im Theorieband,
Workflows im Benutzerband, Referenzsignaturen und Erweiterung im API-Band,
Nachweiskriterien und Gültigkeitsgrenzen im V&V-Band. Wiederholte Grundlagen
werden als Orientierung beibehalten; Fachverträge bleiben kanonisch.

Gemeinsame Titel nennen den tatsächlichen Commit, die Paketversion und den
Kandidatenstatus. Leere Abbildungs-/Tabellenlisten im API-Band entfallen.
Literatur, Sachregister, durchsuchbare Umlaute/Mathematik und PDF-Lesezeichen
werden automatisch geprüft. Es wird keine PDF/A-/PDF/UA-Konformität behauptet.

`coverage.json` erhält neben den unveränderten initialen Ownership-Statuswerten
für jede Quelle eine abschließende Dokumentationszuordnung. `open` im historischen
Initialstatus bezeichnet keine offene Pflichtkapitelliste. Spätere M10–M14-Pläne,
Marken-/Migrationsevidenz und vergangene Releaseberichte sind ausdrücklich
historisch bzw. zukünftig; ihre Nennung macht sie nicht zu unterstützter Physik.
Der Umfang bis M9 umfasst den generischen Unsicherheitsworkflow, jedoch keine
experimentelle Fallfreigabe. Die vollständige Quellenzuordnung steht im Anhang
und in der maschinenprüfbaren Matrix.

## Artefakt- und Rezeptprüfung

Lokale gezielte Prüfung: 24 Build-/Rezept-/API-/V&V-Tests und 5 öffentliche
API-Vertragstests bestanden. Ruff, Formatprüfung (348 Dateien), mypy
(68 Quellmodule), API-/V&V-Driftprüfung und Diffcheck bestanden. 176 lokale
Dokumentationslinkziele wurden geprüft; die persönliche lokale Löschung von
AGENTS.md bleibt erhalten und ist keine Repositoryänderung dieses Beitrags.

Alle 23 Workflows (26 Befehle einschließlich Fortsetzungen) und das
Entwicklungsbeispiel bestanden am tatsächlichen entpackten Archiv.
`recipe-verification.json` hält Befehle, Wheelhash, Hash der entpackten
Rezeptmatrix und tatsächliche Versionen fest: Python 3.13.5, Paket 1.3.0
(Kandidatenversion), NumPy 2.5.4, SciPy 1.18.1, meshio 5.3.5, Matplotlib 3.11.2.
Abhängigkeiten kompatibel. Archivstand vor abschließender reiner Anhang-
Zeilenumbruchkorrektur; Rezeptinputs unverändert und bytegeprüft.
Der erste Versuch ohne Visualisierungsextra erkannte dessen dokumentiertes
Fehlen; die neue vollständige Umgebung erfüllt die Rezeptvoraussetzungen. Die Paketprüfung verwendet
frisch gebautes Wheel mit Visualisierungsextra in einer neuen Umgebung außerhalb
des Checkouts, `python -I`, kein PYTHONPATH. Alle 23 archivierten Befehle,
Fortsetzungsbefehle, Artefaktziele und das Entwicklungsbeispiel werden ausgeführt.
Quellenarchiv wird tatsächlich entpackt; ein Kopieren aus dem Checkout ersetzt
den Archivnachweis nicht.

`python scripts/build_handbook.py --self-test` liefert fünf PDFs, extrahierten
Text, Diagnoseprotokolle, festgelegte Toolchainversionen, Eingabe-/Artefakthashes,
`handbook-source.tar.gz`, `manifest.json` und `SHA256SUMS`. Die entpackten Quellen
bauen ohne weitere Projektdateien und ohne Matplotlib alle fünf PDFs erneut;
Textinhalte werden verglichen. Keine Byteidentität zwischen Toolchains behauptet.

## Sichtprüfung der Kandidaten-PDFs

Alle 491 gerenderten Seiten wurden visuell geprüft: Gesamtband 242,
Theorie 40, Anwendung 71, API/Entwicklung 111 und Nachweise 27 Seiten.
Zusätzlich wurden die wissenschaftlichen Tabellen, Diagramme, Formeln und
der neue Quellenanhang in Einzelansichten kontrolliert. Keine abgeschnittenen
Inhalte, überlappenden Elemente oder überlaufenden Textzeilen festgestellt;
alle fünf abschließenden TeX-Protokolle sind frei von Overfull-Meldungen.
Der Build einschließlich drei negativer Prüfungen und Archivnachbau bestand.
Die Kandidatenartefakte weisen den tatsächlichen Basiscommit und den noch
geänderten Quellenzustand ausdrücklich aus. Der finale saubere Releasebuild
bleibt Aufgabe von #66.

## Veröffentlichung und verbleibende Freigabe

#211 schließt die Dokumentationsabnahme mit Kandidatenartefakten ab. #66 prüft
SemVer/Version/Release Notes, integriert aktuelles main, prüft die CI des finalen
main-Zielcommits und baut die finalen Artefakte aus einem sauberen Stand. Erst
anschließend folgen annotierter unveränderlicher Tag und GitHub Release mit
Paket-, PDF-, Quellen-, Manifest- und Prüfsummenartefakten. Keine PyPI-Publikation
ist Teil dieses Auftrags.
