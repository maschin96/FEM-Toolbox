# Versions- und Release-Prozess

Status: Entwurf

## 1. Zweck

Der Versionszyklus verbindet wissenschaftliche Reife, GitHub-Milestones,
Paketversionen und unveränderliche Git-Tags. Eine Versionsnummer bezeichnet einen
reproduzierbaren Stand mit dokumentierten Fähigkeiten und Nachweisen. Ein
Milestone gilt deshalb nicht allein durch implementierten Code als abgeschlossen.

## 2. Versionsschema

Das Projekt verwendet [Semantic Versioning](https://semver.org/) in der Form
`MAJOR.MINOR.PATCH` und Tags mit vorangestelltem `v`.

- Vor `1.0.0` erhöht jeder fachliche Milestone die Minor-Version. Die öffentliche
  API darf sich noch ändern, Änderungen werden dennoch in den Release Notes
  beschrieben.
- Ab `1.0.0` erhöht eine rückwärtsinkompatible öffentliche Änderung `MAJOR`, eine
  kompatible Funktionserweiterung `MINOR` und eine kompatible Fehlerkorrektur
  `PATCH`.
- Release-Kandidaten dürfen als `vX.Y.Z-rc.N` markiert werden. Sie erfüllen das
  finale Milestone-Tag nicht und werden durch das endgültige Tag nicht bewegt.
- Ein veröffentlichter Tag ist unveränderlich. Nachträgliche Korrekturen erhalten
  eine neue Patch-Version; Tags werden niemals umgehängt oder force-gepusht.

## 3. Geplanter Zyklus

| Zustand/Milestone | Version | Freigabebedeutung |
|---|---:|---|
| Planungsbaseline | `v0.0.1` | dokumentierte Architektur, Backlog und Release-Prozess |
| M0 — Projektfundament | `v0.1.0` | installierbares Paket und aktive Qualitätsinfrastruktur |
| M1 — Netz und Tet4-Geometrie | `v0.2.0` | verifizierte Geometrie- und Netzgrundlage |
| M2 — Ungekoppelte Einzelfelder | `v0.3.0` | verifizierte Mechanik, Wärme und Elektrostatik |
| M3 — Thermische Ausdehnung | `v0.4.0` | verifizierte sequenzielle Thermomechanik |
| M4 — Piezoelektrizität | `v0.5.0` | verifizierte quasistatische Piezoelektrizität |
| M5 — Vollständiger stationärer Workflow | `v1.0.0` | erste stabile stationäre Forschungs-API |
| M6 — Transiente und dynamische Erweiterung | `v1.1.0` | transiente Wärme und lineare gekoppelte Dynamik |
| M7 — Skalierung und weitere Elemente | `v1.2.0` | PETSc-Option, Tet10 und Genauigkeitserweiterungen |
| M8 — Parameterstudien, Feldbilder und Surrogate | `v1.3.0`* | wiederholbare Mehrfallstudien, Tet4-Feldbilder/JPEG und geprüfte polynomiale Surrogate |
| M9 — Unsicherheitsworkflow und vollständige Dokumentation | `v1.4.0`* | verifizierte Unsicherheitsfortpflanzung und vier vollständig geprüfte Fachbände; keine experimentelle Fallfreigabe |
| M10 — Kontrollierte adaptive Netzstudien | `v1.5.0`* | adaptiver stationärer Wärme-Workflow |

*Vorläufige Zielversionen der [Planungswelle 2](13_PRODUCT_VISION.md).
Keine Termine oder Freigabezusagen: Der experimentelle Folgeumfang #63/#65
setzt weiterhin geeignete Messdaten voraus; M9 wurde am 10.10.2026 explizit
auf Software und Dokumentation neu geplant (ADR-054). M10 setzt ein
positives Machbarkeitsgate voraus. Bei veränderter Releasefolge wird die
Versionszuordnung vor Veröffentlichung explizit neu entschieden.

Die Zuordnung ist die Ausgangsplanung. Falls der tatsächliche Umfang eine
inkompatible API-Änderung verlangt, wird die SemVer-Nummer vor Freigabe in einem
Issue und in diesem Dokument korrigiert. Eine niedrigere Version als nach SemVer
erforderlich ist nicht zulässig.

M11 (Topologieoptimierung) und M12 (RANS-CFD mit Fluid-Struktur-Kopplung)
sind zunächst reine Milestone-Beschreibungen ohne Zielversion oder Einzelissues.
Für diese beiden Horizonte wird das Release-Issue erst bei der Detaillierung
vor Umsetzung angelegt; vor einer Freigabe gelten alle folgenden Gates.

## 4. Release-Gate

Die Entwicklung folgt dem Milestone-Branch-Workflow aus [AGENTS.md](../AGENTS.md):
Ein `codex/milestone-m<N>` wird vom stabilen `main` abgeleitet. Jeder Issue-Branch
entsteht aus diesem Milestone-Branch und wird per geprüftem PR dorthin integriert.
`main` behält während der Entwicklung den freigegebenen Funktionsumfang.

Kurz vor Abschluss führt ein Release-PR den freigabefähigen Milestone nach `main`.
Vor diesem Merge müssen Umfang, Verifikation, Dokumentation, Installation des
Release-Kandidaten und die anwendbaren CI-/Reviewprüfungen erfüllt sein. Aktuelle
Änderungen aus `main` werden zuvor im Milestone-Branch integriert und geprüft.
Offen bleiben dürfen nur abschließende Veröffentlichungsarbeiten, keine unfertigen
Funktionen. Nach dem Merge werden die CI-Ergebnisse am tatsächlichen Zielcommit
geprüft; anschließend folgen Tag und Veröffentlichung.

Jeder zur Umsetzung detaillierte GitHub-Milestone enthält ein `release:`-Issue. Dieses Issue wird zuletzt
geschlossen. Vor dem Tag müssen alle folgenden Bedingungen erfüllt sein:

1. Alle verpflichtenden Issues des Milestones sind abgeschlossen und durch
   integrierte Pull Requests nachvollziehbar.
2. Die Definition of Done aus `AGENTS.md` und der Roadmap ist erfüllt.
3. Alle anwendbaren CI-Prüfungen laufen auf dem Zielcommit erfolgreich.
4. Patch-, Bilanz-, Benchmark- und Konvergenztests sind entsprechend dem
   Milestone dokumentiert; nicht anwendbare Prüfungen sind begründet.
5. Paketversion, Dokumentation, Beispiele und `CHANGELOG.md` nennen dieselbe
   Version und denselben Funktionsumfang.
6. Ein Release-Kandidat wurde aus einer frischen Umgebung installiert und die
   vorgesehenen Referenzbeispiele wurden daraus ausgeführt.
7. Der freizugebende Commit befindet sich auf `main` und enthält keine
   unbeabsichtigten oder geheimen Dateien.
8. Lizenzmetadaten und Lizenzdatei stimmen in Wheel, sdist, README und
   Release Notes überein; Hinweise für tatsächlich mitgelieferte
   Drittbestandteile sind geprüft.

Fehlende unabhängige Reviewer werden im Release-Issue genannt. Sie werden nicht
als erteilte Freigabe behandelt.

Für v1.0 ist zusätzlich der [Artefaktnachweis aus Issue #58](29_RELEASE_ARTIFACTS_AND_PROTECTION.md)
zu erfüllen: isolierte Wheel-/sdist-Installation und erfolgreiche CI-Prüfungen.
Branchschutz/Rulesets und deren administrativer Negativtest sind gemäß
Eigentümerentscheidung vom 23.09.2026 (ADR-020) keine Freigabekriterien.
Die Integrationsverantwortlichen prüfen die CI-Ergebnisse am Zielcommit selbst.

## 5. Tagging-Ablauf

Nach bestandenem Gate:

```sh
git switch main
git pull --ff-only
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin vX.Y.Z
gh release create vX.Y.Z --verify-tag --title "vX.Y.Z" --notes-file RELEASE_NOTES.md
```

Der konkrete Release-PR erzeugt die endgültigen Release Notes. Die temporäre
Datei `RELEASE_NOTES.md` muss nicht Teil des Repositories sein. Vor `git push`
wird geprüft, dass `git rev-parse vX.Y.Z^{commit}` dem freigegebenen
`origin/main`-Commit entspricht.

## 6. Planungsbaseline `v0.0.1`

Der Tag `v0.0.1` wird nach Integration von
[Issue #2](https://github.com/maschin96/FEM-Toolbox/issues/2) auf dem zugehörigen
`main`-Commit gesetzt. Er bestätigt nur die dokumentierte Planung, das GitHub-
Backlog und diesen Release-Prozess. Er bezeichnet keinen implementierten oder
verifizierten FEM-Solver.

## 7. Patch- und Notfallversionen

Ein Fehler in einer veröffentlichten Version erhält ein eigenes Bug-Issue mit
Regressionstest. Der Fix wird vom aktuellen unterstützten Stand abgeleitet,
regulär geprüft und als nächste Patch-Version veröffentlicht. Sicherheitskritische
Details und Zugangsdaten gehören nicht in öffentliche Issues. Eine Rückportierung
auf ältere Minor-Versionen wird je Issue ausdrücklich entschieden.

## 8. Versionsquelle im Python-Paket

Ab M0 besitzt das Paket genau eine kanonische Versionsquelle. Build-Metadaten und
`fem_toolbox.__version__` werden daraus erzeugt oder gelesen. Tests prüfen die
Übereinstimmung von Paketversion, erwartetem Release-Tag und Release-Metadaten.
Entwicklungsstände zwischen Releases werden als `X.Y.Z.devN` oder über lokale
Build-Metadaten kenntlich gemacht und nicht als finale Version veröffentlicht.

## Femantis-Namensmigration in M14 / v2.0.0 (#235)

Die bereits geplante M14-Freigabe #113 erhält die Migration zu `femantis` als
zusätzliches verpflichtendes Gate gemäß [ADR-053 und Vertrag](77_FEMANTIS_MIGRATION_CONTRACT.md).
Es entsteht kein konkurrierendes v2.0-Release. Alle v2.x-Versionen behalten den
Legacy-Importadapter; früheste Entfernung ist v3.0.0 nach separater Entscheidung.
Namens-, Registry-, Rechte-, Altimport-, Datenkompatibilitäts- und numerische
Paritätsnachweise gehören vor Veröffentlichung zum vollständigen Gate. Bereits veröffentlichte Releases und historische Artefakte bleiben unverändert.

## Dokumentationsartefakte ab M9

Vor dem Release-Merge müssen vier Fachbände und Gesamthandbuch gemäß
[Gesamtabnahme](latex/ACCEPTANCE.md) inhaltlich und vollständig visuell geprüft
sein. `python scripts/build_handbook.py --self-test` baut fünf PDFs, prüft
API-/V&V-Drift, Abdeckung, Referenzen, Schriften, drei negative Dokumentfälle
und den Nachbau aus dem tatsächlichen Quellenarchiv. Der Checkoutzustand
muss für finale Artefakte sauber sein. Der ausdrückliche Modus
`--release` lehnt schmutzige oder unvollständig zugeordnete Quellen ab und
kennzeichnet PDFs und Manifest als Release. Commit und Version werden aus dem
tatsächlichen geprüften main-Ziel gelesen.

Das Release enthält `main.pdf`, `theory.pdf`, `application.pdf`, `api.pdf`,
`verification.pdf`, `handbook-source.tar.gz`, `manifest.json` und
`SHA256SUMS` zusätzlich zu Wheel und sdist. Das Manifest führt Eingabehashes,
Toolchain, Quellzustand, Referenz-/Archivprüfungen und Layoutwarnungen auf.
Eine zusätzliche Release-Prüfsummenliste umfasst auch beide Paketartefakte;
Handbuch- und Paketprüfsummen werden nicht durch Umbenennung überschrieben.
Titel/Metadaten sowie alle nach dem Kandidaten geänderten Seiten erneut prüfen.
Quellenarchiv und PDFs werden aus demselben endgültigen Stand erzeugt.
Keine Byteidentität über unterschiedliche Toolchains und keine ungeprüfte
PDF/A- oder PDF/UA-Konformität behaupten. Experimentelle Validierung bleibt
gemäß ADR-054 außerhalb der v1.4.0-Freigabe offen.
