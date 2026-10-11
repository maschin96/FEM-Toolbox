# M9-Handbuch und PDF-Build (#206–#211)

Status: Struktur und alle vier Fachbände sind vollständig abgenommen.
Standardbuilds sind Kandidaten; der ausdrückliche Releasebuild gehört zu #66. [Gesamtabnahme #211](ACCEPTANCE.md)
ordnet die abgeschlossenen Einzelbeiträge und Grenzen zu. Der Theorieband enthält zehn Kapitel,
SI-/Voigt-Konventionen, starke/schwache Formen, Diskretisierung, Dynamik,
Solver, Indikatoren und M8-Surrogate. Herkunft und fachliche Eigenprüfung:
[Theorie-Nachweis](THEORY_VERIFICATION.md), [Themenmatrix](theory-map.json).
`main.tex` integriert die vier separat baubaren Bände `theory.tex`,
`application.tex`, `api.tex`, `verification.tex`. Deutsche Sprache, API-Namen
unverändert. In diesem stabilen Quellstand heißt das Python-Paket weiterhin
`fem_toolbox`; die Femantis-Migration läuft getrennt in M14.

## Ein-Befehl-Build aus einem frischen Checkout

Voraussetzungen: Python 3.11–3.13, Git, laufendes Docker mit Linux-Containern,
Netzzugriff nur zum erstmaligen Aufbau des Images. Keine Python-Paketinstallation
und kein TeX auf dem Host nötig. Aus dem Repository-Root:

```sh
python3 scripts/build_handbook.py --self-test
```

Das Werkzeug baut die festgelegte Toolchain, prüft die Quellenmatrix und erzeugt
`build/handbook/`: fünf PDFs, extrahierten PDF-Text, TeX-/Fixture-Logs,
`manifest.json`, `handbook-source.tar.gz` und `SHA256SUMS`.
Standardbuilds erzeugen Kandidatenartefakte; finale Releaseartefakte werden
mit `--release` aus dem geprüften sauberen main-Stand gebaut (#66).
Die Python-Kerninstallation bleibt unverändert und benötigt kein TeX.
Bei Fehlern endet der Befehl ungleich null; ein älteres Erfolgsmanifest wird
vor dem Neuaufbau entfernt. Im Output verbleibende Dateien sind nur zusammen
mit einem aktuellen erfolgreichen Manifest als gültiger Build anzusehen.

Nur die Vollständigkeitsprüfung (ohne Docker):

```sh
python3 scripts/build_handbook.py --check-matrix
```

## Festgelegte Toolchain und Reproduzierbarkeit

`Dockerfile` fixiert das Debian-bookworm-Basisimage auf einen vollständigen
Multiarch-Digest und die Paketquelle auf den Debian-Snapshot
`20261001T000000Z`. Das Archiv liefert inklusive transitiver Abhängigkeiten
einen festen Paketbestand; `manifest.json` erfasst alle installierten
Paketversionen und Engine-Versionen. Engine: **pdfLaTeX / TeX Live 2022**,
Literatur: **BibTeX / plain**, Register: **makeindex**, Schrift: **Latin Modern**
(T1, eingebettet, ToUnicode via glyphtounicode, explizite Integral-Glyphzuordnung), PDF-Prüfung: **Poppler**.
Keine Laufzeitdownloads, kein Shell-Escape. Der eigentliche Dokumentenbuild
läuft ohne Netzwerk, als Hostbenutzer und mit schreibgeschütztem Repository.
Der Container ist ausschließlich ein Dokumentationswerkzeug.

Version stammt aus `pyproject.toml`, Commit und Quellzustand aus Git.
`SOURCE_DATE_EPOCH` stammt aus dem Commit; `FORCE_SOURCE_DATE=1` verhindert
laufzeitabhängige TeX-Datumsangaben. PDF, Buildmanifest und das im Quellenarchiv
enthaltene `metadata.tex` nennen Version, Commit und Veröffentlichungsstatus.
Ein schmutziger Checkout ist sichtbar als `dirty` gekennzeichnet; die
SHA-256-Liste der Eingabedateien im Manifest identifiziert den tatsächlichen
LaTeX-Stand. Inhaltliche Reproduzierbarkeit ist erforderlich. Byteidentität
zwischen Architekturen/Toolchains wird nicht beansprucht; Archivmetadaten sind
nicht auf Byteidentität normiert. Toolchain-Updates sind eigene überprüfte
Änderungen von Digest und Snapshot, mit erneutem PDF-/Fixture-Nachweis.

## Quellenarchiv

Das Archiv enthält alle lokalen TeX-Quellen, Bibliografie, Makefile,
Pflegeanleitung, Quellenmatrix, erzeugte Metadaten und MIT-Lizenz, ohne
Buildzwischenstände. Alle Benutzerrezeptinputs sind byteidentisch unter
`examples/` enthalten; `recipes.json` beschreibt Befehle und Sollwerte.
[Benutzer-Nachweis](USER_VERIFICATION.md) dokumentiert die isolierte Prüfung.
Eigene wissenschaftliche Abbildungen, historische Benchmarkrohwerte und der
synthetische UQ-Bericht liegen im Quellenarchiv. Die Schriften sind in den PDFs
eingebettet. Fremde Messdaten werden entsprechend dem Suchvertrag referenziert.
Nach Entpacken in derselben Toolchain: `make -C handbook` baut alle fünf PDFs.
Die Makefile dient dem Quellenarchiv; die strikten Diagnose-, PDF- und
Negativprüfungen übernimmt der oben dokumentierte Python-Befehl im Repository.
Alle Markdown-Verweise in der Matrix sind Herkunftsverweise, keine zum
TeX-Kompilieren nötigen Eingabedateien.

## Fehlerdiagnose und Verifikation

Drei absichtlich defekte Kopien des Musterbandes prüfen fehlende Eingabedatei,
undefinierten Querverweis und fehlenden Bibliografieeintrag. Jeder muss
scheitern. Positive Builds laufen pdfLaTeX, BibTeX, makeindex und zwei weitere
pdfLaTeX-Pässe. Undefinierte Zitate/Referenzen, doppelte Labels, LaTeX-Fehler
und erforderliche Referenz-Neuläufe sind fatal. Overfull-/Underfull-Layout-
warnungen bleiben im Manifest und vollständigen `.log` sichtbar und benötigen
visuelles Review; sie werden nicht pauschal verborgen.

PDF-Prüfung: durchsuchbarer Text einschließlich Umlauten, Version/Commit,
Literatur und Sachregister; PDF-Lesezeichen; eingebettete Schriften mit
Unicode-Abbildung; Seiten/Metadaten im Manifest. Visuell werden Titel,
Verzeichnisse, Symboltabelle, Gleichungen, Abbildung, Tabelle, Code,
Literatur und Index geprüft. Keine PDF/A-Behauptung. GitHub baut bei relevanten
PRs und Änderungen an main/M9, führt Negativprüfungen aus und liefert
Kandidatenartefakte. CI ersetzt das visuelle/fachliche Review nicht.

## Pflegevertrag und Kapitel hinzufügen

`coverage.json` ist die maschinenprüfbare initiale Abdeckungs-/Quellenmatrix.
Jeder bestehende Doc-Eintrag (einschließlich historischer Releases und lokaler
Nachweis-/Markenassets), README/INSTALL und jede öffentliche API-Gruppe besitzt
Zielband, kanonische Quelle, Pflegeverantwortung und Status. Die Theoriegliederung
liegt in `chapters/theory/`; der Benutzerband in `chapters/application/`. API/V&V werden in #209–#210 ausgearbeitet.
M8-Studien,
Berichte, Surrogate und Visualisierung sind explizit enthalten. M9-Datenvertrag,
Unsicherheit und Validierung bleiben geplante Inhalte im Nachweisband.

1. Ein Thema und seine kanonische Quelle in der Matrix auswählen. Mathematik:
   `Doc/02_MATHEMATICAL_MODEL.md` hat Vorrang; abweichende Entscheidungen in
   `Doc/07_RISKS_AND_DECISIONS.md` protokollieren. API: Code/Tests und die
   spezifischen Vertragsdokumente gemeinsam prüfen.
2. In `chapters/` ein Kapitel oder Untermodul erstellen und im Zielband
   einbinden. Präfixe für Labels: `eq:`, `fig:`, `tab:`, `sec:`, `code:`;
   eindeutige Namen, Bibliografie in `references.bib`, lokale Assets mit Herkunft.
3. Matrix aktualisieren; neue Doc-Dateien und öffentliche API-Gruppen werden
   bei fehlender Zuordnung abgelehnt. Historische Nachweise bleiben unverändert.
4. Kanonische Markdown-Quelle zuerst ändern; LaTeX dient bis #211 als
   redaktionell abgeleitete Darstellung. Keine automatische Synchronisation.
   Der jeweilige PR aktualisiert beide Darstellungen oder begründet offen,
   warum der noch nicht ausgearbeitete Band keine behauptete aktuelle Kopie ist.
   Zuständige Rollen sind fachliche Maintainer, API-Maintainer bzw. Release-
   Verantwortliche; sie sind keine erfundenen benannten Reviewer.
5. Build, Negativprüfungen, Inhalts-/Linkprüfung und visuelles Review ausführen;
   Quellen und wissenschaftliche Aussagen am dokumentierten Commit prüfen.
   Die ursprünglichen Statuswerte `open`, `sample`, `historical` bleiben als
   Herkunftsstand erhalten. Das ergänzende Feld `theory_coverage` kennzeichnet
   ausgearbeitete Theoriezuordnungen, ohne eine Fremdfreigabe zu behaupten.
   `theory-map.json` enthält die konkrete Kapitel-/Modul-/Testzuordnung.

## Lizenzen und Weitergabe

Eigene Quellen und Abbildung: Repository-MIT-Lizenz. Es werden keine externen
Fotos, Fontdateien oder TeX-Pakete in das Quellenarchiv vendort. Die bestehenden
Markenfotos sind in der Quellenmatrix inventarisiert, aber nicht in PDFs
verwendet; ihre Herkunft bleibt unter `Doc/brand/` dokumentiert.
Latin Modern wird über Debian/TeX Live bezogen; Lizenz **GUST Font License**:
[CTAN](https://ctan.org/pkg/lm),
[GUST](https://www.gust.org.pl/projects/e-foundry/licenses).
LaTeX und wesentliche Pakete stehen unter ihren eigenen Lizenzen, insbesondere
[LPPL](https://www.latex-project.org/lppl/); Paketdetails liegen im Container
unter `/usr/share/doc/*/copyright` und sind über die festgehaltenen Paketversionen
nachvollziehbar. Poppler (GPL), Python (PSF), Debian-Bestand behalten ihre
Lizenzen. Der Container wird hier nicht als Release-Artefakt weitergegeben;
PDFs enthalten eingebettete Latin-Modern-Schriften. Bei späterem Vertrieb eines
Containers sind dessen vollständige Lizenz-/Quellpflichten gesondert zu erfüllen.

## Benutzerbeispiele aktualisieren (#208)

Bei Änderungen unter `examples/` die gleichnamigen Dateien nach
`Doc/latex/examples/` übernehmen, SHA-256 in `recipes.json` aktualisieren und
Rezepttext/Befehle/Sollwerte prüfen. `tests/unit/test_user_handbook_sources.py`
erkennt fehlende oder abweichende Kopien und nicht eingebundene Programme.
Nicht nur Hashes erneuern: die Rezepte aus dem gebauten Wheel ausführen und
`USER_VERIFICATION.md` mit Runtime und Resultaten aktualisieren. Das PDF und
der Quellenarchiv-Roundtrip sind anschließend erneut zu prüfen. Bei Änderungen an #64 den
synthetischen Unsicherheits-Lernworkflow samt Quellenprüfung ergänzen.

## Nachweisband

[Nachweisprotokoll](VV_VERIFICATION.md), [Fachmatrix](verification-map.json) und
[quellengebundene Auswertung](verification-data.json) dokumentieren M7-Kosten,
M8-Surrogatfehler, M9-Sampling/Sensitivität und offene experimentelle Validierung.
`python scripts/verification_artifacts.py` prüft Quellen, Tabellen, Figuren und
den unveränderten 72-Lauf-Bestand; der PDF-Build führt denselben Audit aus.

## Gesamtabnahme und Quellenzuordnung

[Abnahmebericht](ACCEPTANCE.md) dokumentiert den Kandidatenumfang, Artefakt-
und vollständige visuelle Prüfung. Der Gesamtband enthält den vollständigen
Quellenanhang; `coverage.json` unterscheidet dokumentierte Fachverträge,
historische Nachweise, zukünftige Horizonte und offenen experimentellen
Folgeumfang. Initiale Ownership-Statuswerte bleiben als Planungshistorie erhalten.

## Finale Releaseartefakte

`python scripts/build_handbook.py --self-test --release` baut ausdrücklich
gekennzeichnete Release-PDFs und ein Release-Manifest. Dieser Modus verlangt
einen sauberen Checkout, eine finale Version und vollständige Commit-Provenienz.
Finale Artefakte werden aus dem tatsächlich geprüften main-Zielcommit erzeugt;
[Freigabenachweis](../68_M9_RELEASE_VERIFICATION.md) und #66 verknüpfen Tag/Release.
Der Standardbuild bleibt ein Kandidatenbuild, auch wenn die Paketversion final ist.
