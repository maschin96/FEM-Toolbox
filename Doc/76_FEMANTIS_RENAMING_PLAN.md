# Femantis — vollständiger Umbenennungsplan

05.10.2026 · Planung, noch keine Code-Migration · [Issue #233](https://github.com/maschin96/FEM-Toolbox/issues/233)

[Corporate Identity](75_FEMANTIS_CORPORATE_IDENTITY.md) · [Architektur](04_SOFTWARE_ARCHITECTURE.md) ·
[API-Vertrag](29_PUBLIC_API_CONTRACT.md) · [Release-Prozess](09_RELEASE_PROCESS.md)

## 1. Ausgangslage und Ziel

Planungsbasis: stabiler `origin/main`, Commit
`603e96fdd9a0b7abf28edef4b10fad6e1ee5c777`, Paketversion 1.3.0.
Die unabhängige Markendokumentation geht per Issue-PR nach main;
aktive M12-Entwicklung wird nicht mit integriert. Der bei Arbeitsbeginn
geöffnete Branch `codex/230-quad-fv-simple` bleibt unverändert.

Eine Bestandsuche mit `rg -l 'fem_toolbox|fem-toolbox|FEM.Toolbox' --glob '!uv.lock'`
findet auf dieser Basis 259 getrackte bzw. nicht ignorierte Dateien. Das ist
kein Ersatz für die spätere vollständige Suche: Groß-/Kleinschreibung,
Punkt-/Unterstrichvarianten, binäre Artefakte und externe Dienste zusätzlich prüfen.

| Ebene | Heute | Ziel |
|---|---|---|
| Produkt | FEM Toolbox | Femantis |
| Gestaltete Wortmarke | bisher keine implementierte Wortmarke | FEMantis gemäß Leitfaden |
| Python-Distribution | `fem-toolbox` | `femantis` |
| Python-Namespace | `fem_toolbox` | `femantis` |
| Quellwurzel | `src/fem_toolbox/` | `src/femantis/` |
| GitHub-Repository | `maschin96/FEM-Toolbox` | `maschin96/Femantis`, nach Verfügbarkeitsprüfung |
| Dokumentations-/Issue-URLs | bisherige Repository-URLs | neue kanonische URLs |
| Persistierte Formatkennung | u. a. `fem_toolbox:` / `fem-toolbox-heat-series` | zunächst kompatibler Altvertrag; neue Kennungen nur versioniert |
| Alte Releases/Provenienz | historischer Name und Version | unveränderte historische Evidenz |

Der Auftrag plant die vollständige Umbenennung im Code. Die Ausführung,
Repository-Umbenennung und Veröffentlichung sind eigene Folgearbeiten.
Bestehende numerische Ergebnisse, SI-Konventionen und API-Semantik dürfen
sich durch die Umbenennung nicht verändern.

## 2. Kompatibilitätsentscheidung zur Umsetzung

Empfohlene Zielstrategie: kanonischer Namespace `femantis`; befristeter
Legacy-Namespace `fem_toolbox` als expliziter Kompatibilitätsadapter für einen
vollständigen freigegebenen Übergangszyklus. Entfernen erst in einer späteren
Major-Version mit dokumentierter Vorwarnung. Da sowohl Import- als auch
Distributionsname öffentlich sind, ist ein Major-Release vorzusehen; Nummer
und Einordnung mit der vorhandenen v2.0-Planung abstimmen, nicht automatisch
belegen oder zurücksetzen.

Der Adapter muss für jeden unterstützten öffentlichen Untermodulpfad dieselben
Objekte liefern. Kein Kopieren des Solvercodes, kein doppeltes Laden derselben
Klassen unter zwei Modulidentitäten. Exceptions, isinstance, registrierte
Objekte, relative Imports und Serialisierung sind ausdrücklich zu testen.
DeprecationWarning mit eindeutigem Migrationsziel und geeignetem stacklevel;
Warnungen in Tests gezielt prüfen, nicht global unterdrücken.

Die alte Distribution kann für die Übergangsphase ein separat geprüftes,
versionsgebundenes Metapaket werden. Der Legacy-Importadapter gehört nur in
die neue Distribution, damit beide Installationen keine konkurrierenden Dateien
besitzen. Die Publikationsberechtigung beider Paketnamen muss vorab geklärt sein.
Kein leeres Paket oder Reservierungsupload als Teil der Dokumentation.

Persistierte Kennungen zunächst erhalten: `fem_toolbox:*`,
`fem-toolbox-heat-series`, bestehende Schemafelder wie `toolbox_version` und
historische `solver_name`-Werte sind Datenschnittstellen. Altdateien bleiben
lesbar. Neue Ausgaben dürfen aktuelle Produkt-/Distributionsmetadaten nur
innerhalb des dokumentierten Schemas oder nach dessen Versionserhöhung
enthalten. Zukünftige `femantis:*`-Writer benötigen Dual-Reader und explizite
Konfliktregel: widersprüchliche alte/neue Daten ablehnen, nicht still bevorzugen.
Historische Ergebnisse werden nicht nachträglich umgeschrieben.

Diese Kompatibilitätsstrategie ist ein konkreter Umsetzungsvorschlag. Die
Architekturentscheidung wird vor abhängigem Code als ADR mit genauer
Übergangsversion und unterstützten Altpfaden beschlossen.

## 3. Vollständige Änderungsflächen

| Fläche | Änderung und Prüfnachweis |
|---|---|
| `src/` | Quellwurzel, alle absoluten/lokalen Imports, Docstrings, Meldungen und package exports; Legacyadapter separat |
| `pyproject.toml` | name, package-data/py.typed, URLs, Toolpfade, Extras und Paketmetadaten; Versionskonsistenz |
| Versions-/Provenienzmodell | `__version__`, `importlib.metadata`, `RunProvenance`, solver_name, distribution records; Altresultate lesbar |
| I/O und Schemas | `io/meshio_adapter.py` reservierte VTU-Präfixe, Ergebnisbündel, Wärme-/Piezo-Serien, Model-/Study-/Surrogat-JSON/TOML; Kompatibilität versioniert |
| Tests | Imports, Fixtures, erwartete Meldungen, Metadaten, öffentliche API, Typinformationen, Artefakt-/Installationsprüfungen |
| Beispiele/Benchmarks | Pythonimporte, Installationsbefehle, Ergebnisnamen, Bildbeschriftungen und Anleitungen; alle Rezepte erneut ausführen |
| Skripte/CI | `.github/workflows/`, `ci/`, `scripts/`, mypy-Pfade, wheel/sdist-Dateinamen, pip extras, Docker-Tags und PETSc-Checks |
| Distribution | Wheel/sdist, LICENSE/readme-Metadaten, Typdatei, Abhängigkeiten, optionale visualization- und PETSc-Pfade; saubere Umgebungen |
| Dokumentation | README, INSTALL, AGENTS, Doc, Links, API-Vertrag, zukünftige LaTeX-Handbücher; historische Evidenz markieren statt ändern |
| GitHub | Repository-Name, About, Issue-/PR-Vorlagen, Homepage/Topics, Release-/Artifact-Links, Badges und Workflows |
| Externe Integrationen | PyPI/TestPyPI falls genutzt, Containerregistry, Docs-Hosting, Webhooks, Automations, lokale Remotes und Verbraucher |
| Markenassets | Freigegebene Master, Banner/Icon, Alt-Texte und Inventar gemäß CI; wissenschaftliche Feldbilder erhalten |

Auch gegenwärtige und später integrierte CFD-/FSI-Module erfassen. Kurz vor
Migration erneut auf main und allen aktiven Milestone-Branches inventarisieren.
Historische Tags, Release-Archive, Evidenzchecksummen, DOI-/Zitationsdaten und
ursprüngliche Copyrightnamen bilden dokumentierte Ausnahmen der Restnamensuche.
Kein unkontrolliertes globales Suchen/Ersetzen in Datendateien.

## 4. Unabhängig prüfbare Arbeitspakete und Reihenfolge

Vor jedem Paket passendes Issue suchen oder erstellen. Hier sind Pakete
geplant; ihre Nummern werden erst durch GitHub vergeben. Kein künstlicher
Milestone-Branch, bevor die Release-Einordnung entschieden und Arbeit begonnen ist.

| Paket | Ergebnis | Abnahme / Abhängigkeit |
|---|---|---|
| A — Inventar und ADR | Exakte Namensliste, externe Dienste, Namens-/Rechteprüfung, Major-/Übergangsvertrag | Verantwortliche entscheiden Name, Zeitraum und Release-Ziel vor B |
| B — Namespace und Adapter | `femantis` als einzige Implementierung; getestete `fem_toolbox`-Adapter | Beide Importpfade, Objektidentität, Exceptions, Warnungen und Typen; nach A |
| C — Packaging und CI | Neue Distribution, Typdaten, Extras, Build-/Installations- und optionale Backendprüfung | Wheel/sdist in sauberen CPython-3.11–3.13-Umgebungen; nach B |
| D — Persistenz und Provenienz | Dokumentierter Alt-/Neuformatvertrag, Reader/Writer und Metadaten | Historische Fixtures unverändert lesbar, Konflikte abgelehnt, Provenienz nachvollziehbar; nach A/B |
| E — Beispiele und Dokumentation | Alle aktuellen Rezepte/Links/API-Beschreibungen aktualisiert | Ausführbare Beispiele und Publikationsabdeckung; nach B/C/D |
| F — GitHub und Release | Kanonische Repository-URLs, Markenassets, finaler Freigabebericht | Alle Gates am aktuellen Commit, unabhängiges Review; nach A–E |

Codepakete gehen vom beschlossenen Integrationsbranch aus per eigenem
`codex/<issue>-<beschreibung>`-Branch zurück in diesen. Die aktuelle reine
CI-/Planungsdokumentation ist milestone-unabhängige Repositorypflege mit
PR nach main. Die komplette Migration ist keine vorgezogene M12-Funktion.
Mit parallelen Meilensteinen die Reihenfolge koordinieren, um veraltete
Importnamen nicht nach der Migration erneut einzuführen.

## 5. TDD- und Verifikationsplan für die spätere Umsetzung

1. Vor B einen Importtest für `femantis` aus einem unabhängig aufgebauten
   Beispiel ausführen: Red muss am fehlenden Namespace scheitern.
2. Minimalen kanonischen Namespace implementieren; Green für API und
   Versionsmetadaten. Tests für sämtliche dokumentierten Untermodule ergänzen.
3. Legacytest für Klassen-/Exceptionidentität, gemischte Importreihenfolge und
   Warnung zuerst rot nachweisen; dann Adapter implementieren und refaktorieren.
4. Vor C Artefakttest auf neue Distribution, py.typed und Extras schreiben.
   Wheel und sdist jeweils in frischer Umgebung installieren, ohne Import vom
   Checkout; alte/neue Distribution einzeln und gemeinsam prüfen.
5. Vor D unveränderte historische VTU-/JSON-/TOML-Fixtures laden. Bei neuen
   Kennungen gezielte Red-Tests für Roundtrip und widersprüchliche Dualkennungen.
6. Fachliche Parität mit unabhängigen Referenzen, Patch-/Bilanz-/Benchmarktests
   und den relevanten Konvergenzsuiten sichern. Keine Referenzzahlen oder
   Toleranzen wegen Namensänderungen neu erzeugen oder abschwächen.
7. Vollständige pytest-, ruff-, Format-, mypy- und Diff-Prüfung ausführen,
   einschließlich gesonderter numerischer/optionaler Backend-CI am aktuellen
   Stand. Pfade nach Migration auf `src/femantis` und Adapter abstimmen.
8. Alle freigegebenen Beispiele, Studien-/Surrogatrezepte und Dokumentations-
   Builds gemäß ihrer aktuellen Abdeckung prüfen; keine ungeprüfte Physik zusagen.

Reine Planung benötigt keinen künstlichen Red-Test. Die obigen Red-/Green-
Nachweise gehören in die späteren Umsetzungs-PRs und sind heute nicht erbracht.

## 6. Veröffentlichung und Rückfallplan

Vor Release Namensberechtigung und Paketverfügbarkeit prüfen, alle Arbeitspakete
integrieren, main in den Integrationsbranch übernehmen und erneut prüfen.
GitHub-Umbenennung zeitlich mit Links, CI und Verbrauchern koordinieren;
Weiterleitungen und reale Clone-/Issue-/Release-URLs ausdrücklich testen.
Keine Garantie, dass externe Integrationen automatisch folgen.

Release-PR erst nach vollständiger Funktion, Dokumentation, aktuellen grünen
Prüfungen und geklärtem Review nach main. Danach CI des tatsächlichen main-
Zielcommits prüfen; erst anschließend Tag und Artefakte veröffentlichen.
Migration erklären: `pip install femantis`, `from femantis.mesh import Mesh`
sowie Übergangsverhalten des Altimports. Diese Befehle sind Zielbeispiele,
heute noch keine verfügbare Installationsanleitung.

Vor Publikation kann ein geprüfter Revert-PR den dokumentierten Ausgangsstand
wiederherstellen. Nach Veröffentlichung keine Tags verschieben oder Pakete
überschreiben; fehlerhafte Version nach Richtlinien der Registry sperren bzw.
mit Patch korrigieren. Bei Registry-Problemen stabile Altversion weiter anbieten
und Migration aussetzen. Reponame nur koordiniert zurückändern, da externe
Verweise betroffen sein können. Die Leser alter Daten bleiben erhalten.

## 7. Release-Gates

- [ ] Name und externe Zuständigkeiten geprüft; ADR/Versionsziel beschlossen.
- [ ] Kanonische Imports, Altadapter, Typen und Paketmetadaten geprüft.
- [ ] Aktuelle Restnamenliste enthält nur begründete historische/kompatible Ausnahmen.
- [ ] Alte Artefakte lesbar, neue Artefakte eindeutig versioniert und reproduzierbar.
- [ ] Numerische Parität, gesamte relevante Qualitäts- und Backend-CI bestanden.
- [ ] Beispiele, Markenassets, Dokumentation und externe Links geprüft.
- [ ] Unabhängiges Review, Release-PR und main-CI am tatsächlichen Stand abgeschlossen.
- [ ] Tag, Distribution und Release-Nachweise erst danach veröffentlicht.

Keines dieser Migrationsgates wird durch die heutige Dokumentation als
abgeschlossen behauptet.
