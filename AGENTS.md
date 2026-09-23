# Entwicklungsrichtlinien für die FEM Toolbox

Diese Vorgaben gelten für das gesamte Repository und für menschliche wie
automatisierte Beiträge. Explizite Aufgabenanweisungen haben Vorrang; spezifischere
AGENTS.md-Dateien ergänzen diese Regeln in ihren Unterverzeichnissen.

## Projektkontext und Arbeitsbeginn

- Vor Änderungen [README.md](README.md), die relevanten Dokumente unter
  [Doc/](Doc/README.md), den Git-Status und vorhandene Änderungen prüfen.
- Das Repository hat M0–M6 als v1.1.0 abgeschlossen (Stand 23.09.2026).
  Python-Paket, Tests und CI sind vorhanden; M7 ist die nächste Ausbaustufe.
  Fehlende Werkzeuge oder Prüfungen offen benennen; keine Testresultate erfinden.
- Architektur und Zielstruktur aus [Doc/04_SOFTWARE_ARCHITECTURE.md](Doc/04_SOFTWARE_ARCHITECTURE.md)
  beachten. Mathematische Konventionen aus
  [Doc/02_MATHEMATICAL_MODEL.md](Doc/02_MATHEMATICAL_MODEL.md) sind maßgeblich.
- Kleine, fachlich zusammenhängende Änderungen durchführen. Fremde Änderungen
  erhalten; keine unaufgeforderten Resets, History-Rewrites oder Force-Pushes.

## GitHub-Workflow: Issue vor Umsetzung

1. Vor einer Umsetzung nach einem passenden offenen GitHub-Issue suchen. Ein
   vorhandenes Issue verwenden oder ein neues erstellen; Duplikate vermeiden.
2. Das Issue enthält Problem und Ziel, Umfang und Nichtziele, überprüfbare
   Akzeptanzkriterien sowie Test- und Verifikationsplan. Bei Fehlern zusätzlich
   minimales Reproduktionsbeispiel, erwartetes/tatsächliches Verhalten und Umgebung
   festhalten. Passende vorhandene Labels und Meilensteine verwenden.
3. Größere Aufgaben in unabhängig prüfbare Arbeitspakete zerlegen. Fachliche
   Unklarheiten vor der davon abhängigen Implementierung auflösen.
4. Pro Milestone einen eigenen Integrationsbranch `codex/milestone-m<N>` vom
   aktuellen stabilen `main` anlegen, beispielsweise `codex/milestone-m6`.
   Für jedes Issue einen eigenen Branch
   `codex/<issue-nummer>-<kurzbeschreibung>` vom zugehörigen Milestone-Branch
   ableiten. Issue-Branches per Pull Request zurück in diesen Milestone-Branch
   integrieren; das PR-Ziel ausdrücklich prüfen. Milestone-Entwicklung wird
   nicht direkt nach `main` gemergt. Ausdrücklich beauftragte direkte Commits
   ändern diese Integrationsziele nicht.
5. Commits fokussiert halten und nach dem Muster
   `test: ...`, `feat: ...`, `fix: ...`, `refactor: ...`, `docs: ...` oder
   `ci: ...` mit Issue-Referenz formulieren, beispielsweise `docs: add workflow (#1)`.
   Vor dem Commit den gestagten Diff prüfen und nur beabsichtigte Dateien aufnehmen.
6. Der Pull Request beschreibt Problem, resultierendes Verhalten, Prüfungen und
   relevante Risiken. `Closes #<nummer>` nur verwenden, wenn alle Kriterien des
   Issues erfüllt sind; Teilbeiträge mit `Refs #<nummer>` verknüpfen.
7. Vor Integration müssen die anwendbaren CI-Prüfungen erfolgreich und
   Review-Anmerkungen geklärt sein. Unabhängiges Review anstreben; fehlende
   Reviewer transparent benennen und niemals eine Freigabe behaupten.
   Ein lokaler Commit allein schließt das Issue nicht ab.

### Milestone-Abschluss und stabiler Hauptbranch

- `main` enthält stets den aktuellen stabilen, freigegebenen Funktionsumfang.
  Neue Milestone-Funktionalität bleibt bis zur Freigabereife auf dem jeweiligen
  Milestone-Branch. Weitere Milestone-Branches werden erst bei Arbeitsbeginn
  angelegt; benötigte Vorgänger müssen zuvor stabil in `main` integriert sein.
- Der Ablauf ist `main → Milestone-Branch → Issue-Branch → Milestone-Branch`.
  Jeder Issue-PR muss vor dem Merge die anwendbaren CI- und Reviewregeln erfüllen.
  Vollständig integrierte und geprüfte Issues können nach dem Merge in den
  Milestone-Branch mit Verweis auf den PR explizit geschlossen werden; nicht
  allein auf automatische Schließung durch `Closes` vertrauen. Das Release-Issue
  bleibt bis zum vollständigen Release-Abschluss offen.
- Erst kurz vor Fertigstellung wird der Milestone als Release-PR nach `main`
  integriert. Voraussetzung sind der vollständige verpflichtende Funktionsumfang,
  abgeschlossene Verifikation, aktuelle Dokumentation sowie erfolgreiche
  anwendbare CI-Prüfungen und geklärte Review-Anmerkungen am aktuellen PR-Stand.
  „Kurz vor Fertigstellung“ erlaubt offene Veröffentlichungsarbeiten wie Tag und
  Release-Artefakte, aber keine unfertige oder ungeprüfte Funktionalität auf `main`.
- Vor dem Release-Merge Änderungen aus `main` in den Milestone-Branch übernehmen,
  Konflikte dort lösen und den resultierenden Stand erneut prüfen. Nach dem Merge
  die CI am tatsächlichen `main`-Zielcommit prüfen; erst danach gemäß
  [Release-Prozess](Doc/09_RELEASE_PROCESS.md) taggen und veröffentlichen.
- Milestone-unabhängige Repository-/Workflow-Pflege und notwendige Bugfixes der
  stabilen Version verwenden eigene Issue-Branches von `main` mit PR nach `main`.
  Relevante Korrekturen anschließend in aktive Milestone-Branches übernehmen.
  Diese Ausnahme dient nicht zur vorzeitigen Integration neuer Milestone-Funktionen.

Wenn GitHub nicht erreichbar ist, Issue-Inhalt lokal vorbereiten, die Einschränkung
melden und unabhängige lokale Arbeit fortsetzen. Keine Issue-Nummern erfinden.
Ein Auftrag zum lokalen Commit ist kein automatischer Auftrag zum Push oder Merge.

## Test-Driven Development: Red – Green – Refactor

Für neue ausführbare Funktionalität und Fehlerkorrekturen ist TDD verpflichtend:

1. **Red:** Aus einem Akzeptanzkriterium einen kleinen, deterministischen Test
   ableiten und vor der Implementierung ausführen. Prüfen, dass er wegen des
   fehlenden Verhaltens scheitert, nicht wegen einer defekten Testumgebung.
   Bei Bugfixes zuerst einen Regressionstest schreiben, der den Fehler reproduziert.
2. **Green:** Die kleinste fachlich korrekte Implementierung erstellen, die den
   Test erfüllt. Den gezielten Test und betroffene vorhandene Tests ausführen.
3. **Refactor:** Struktur und Lesbarkeit bei unverändertem Verhalten verbessern;
   anschließend die betroffenen Tests erneut ausführen.
4. Für das nächste Verhalten wiederholen. Vor Abschluss die für die Änderung
   relevante Testsuite und alle konfigurierten Qualitätsprüfungen ausführen.

Tests prüfen beobachtbares Verhalten und unabhängige fachliche Erwartungen, keine
Kopie der Implementierung. Normalfälle, Randfälle und ungültige Eingaben abdecken.
Zufallsquellen mit festen Seeds kontrollieren; Netz, Uhrzeit und externe Dienste
in Unit-Tests vermeiden. Mocking auf Infrastrukturgrenzen beschränken.

Im PR die ausgeführten Befehle, den beobachteten Red-/Green-Nachweis und eventuelle
Prüflücken kurz festhalten. Absichtlich fehlschlagende Zwischenstände dürfen lokal
entstehen; der zur Integration vorgesehene Stand muss grün sein. Tests niemals
deaktivieren oder Erwartungen abschwächen, nur um CI erfolgreich zu machen.

Für reine Dokumentationsänderungen ohne ausführbares Verhalten ist kein künstlicher
Red-Test nötig: Inhalt, Links und Konsistenz prüfen. Bei verhaltensneutralem
Refactoring vorhandene Tests zuerst ausführen und relevante Abdeckungslücken vor
dem Umbau schließen.

## Python und Architektur

- Den in README und Projektkonfiguration festgelegten Python-Support beachten.
  Entwicklungswerkzeuge und Abhängigkeiten in `pyproject.toml` zentral konfigurieren,
  sobald M0 umgesetzt wird; reproduzierbare Installation dokumentieren.
- `src/fem_toolbox/` und die geplante Teststruktur verwenden. Kleine Funktionen,
  explizite Datenflüsse und klare Modulgrenzen bevorzugen; globale Zustände vermeiden.
- Öffentliche APIs typisieren und Annahmen, Array-Formen, SI-Einheiten,
  Rückgabewerte und Fehlerfälle dokumentieren. Ungültige Daten früh mit konkreten
  Exceptions ablehnen; Fehler nicht stillschweigend unterdrücken.
- NumPy/SciPy bilden den numerischen Kern. I/O, Meshing und Visualisierung bleiben
  Adapter. Neue Abhängigkeiten fachlich begründen und ihre Lizenz prüfen.
- `pytest`, `ruff` und `mypy` sind die vorgesehenen Qualitätswerkzeuge. Nach ihrer
  Einrichtung gelten folgende Basisprüfungen aus dem Repository-Root; dokumentierte
  projektspezifische Befehle und zusätzliche CI-Checks ebenfalls berücksichtigen:

  ```sh
  python -m pytest
  python -m ruff check .
  python -m ruff format --check .
  python -m mypy src/fem_toolbox
  git diff --check
  ```

- Diese Befehle setzen eine eingerichtete Entwicklungsumgebung voraus. Fehlende
  Konfiguration nicht als bestandene Prüfung behandeln. Keine unkonfigurierten
  Quality Gates als bereits aktiv darstellen.
- Erst eine verifizierbare Referenzimplementierung erstellen; Optimierungen durch
  Profiling begründen und gegen diese Referenz prüfen.

## Wissenschaftliche Verifikation

[Doc/05_VERIFICATION_AND_VALIDATION.md](Doc/05_VERIFICATION_AND_VALIDATION.md)
und die Abnahmekriterien aus [Doc/06_ROADMAP.md](Doc/06_ROADMAP.md) gelten zusätzlich:

- Numerische Erwartungen aus analytischen Lösungen, unabhängig hergeleiteten
  Referenzen oder physikalischen Invarianten ableiten.
- Absolute und relative Gleitkommatoleranzen explizit festlegen und anhand von
  Einheiten, Größenordnung und Konditionierung begründen; keine unmotivierte
  exakte Gleichheit und kein nachträgliches Aufweiten zum Kaschieren von Fehlern.
- Für betroffene Operatoren und Felder Unit-, Patch-, Bilanz- und Benchmarktests
  ergänzen. Konvergenztests bei Änderungen an Diskretisierung oder Solververhalten
  durchführen. Geometriefehler, Starrkörpermoden und Randbedingungen berücksichtigen.
- Einzelfelder vor Kopplungen verifizieren. Kleines algebraisches Residuum allein
  genügt nicht: feldweise Residuen, Energie, Reziprozität und Erhaltung prüfen,
  soweit sie für das Modell anwendbar sind.
- Referenzdaten mit Herkunft, Parametern und Toleranzen dokumentieren. Änderungen
  von Referenzwerten fachlich begründen, nicht blind neu erzeugen.
- Änderungen an mathematischen Konventionen oder Architekturentscheidungen in
  [Doc/07_RISKS_AND_DECISIONS.md](Doc/07_RISKS_AND_DECISIONS.md) festhalten.
  „Verifiziert“ und „validiert“ gemäß den Projektdefinitionen unterscheiden.

## GitHub Actions und Sicherheit

Bei Einrichtung oder Änderung der CI:

- Pull Requests und Änderungen am Hauptbranch prüfen. Unterstützte
  Python-Versionen in einer dokumentierten Testmatrix abdecken; Installation,
  Tests, Linting, Formatierung und Typprüfung automatisieren.
- Langsame Konvergenztests gesondert ausführen, aber für relevante numerische
  Änderungen vor Integration nachweisen. Ausgelassene Tests sichtbar machen.
- Branchschutz/Rulesets sind gemäß Eigentümerentscheidung vom 23.09.2026
  (ADR-020) keine Pflicht und kein Integrations- oder Release-Gate. Die
  Integrationsverantwortlichen prüfen die anwendbaren CI-Ergebnisse am aktuellen
  Commit selbst; die oben beschriebenen PR- und Reviewregeln gelten weiterhin.
- GitHub Actions auf vollständige Commit-SHAs fixieren und Updates nachvollziehbar
  pflegen. `GITHUB_TOKEN` standardmäßig auf minimale Leserechte begrenzen und
  Schreibrechte nur gezielt für erforderliche Jobs vergeben.
- Unvertrauenswürdigen PR-Code nicht mit Secrets oder privilegiertem Token
  ausführen, insbesondere nicht über `pull_request_target`. Fremde Eingaben nicht
  direkt in Shell-Code interpolieren.
- Keine Zugangsdaten, privaten Datensätze oder lokalen Umgebungen committen.
  Abhängigkeiten regelmäßig auf Sicherheitsmeldungen prüfen; Sicherheitsbefunde
  nicht mit veröffentlichten Zugangsdaten in öffentliche Issues schreiben.

## Definition of Done

- Issue und Akzeptanzkriterien sind nachvollziehbar erfüllt.
- Relevante TDD-, Regressions- und wissenschaftliche Nachweise liegen vor;
  nicht anwendbare Prüfungen sind begründet.
- Anwendbare Qualitätsprüfungen sind bestanden; nicht ausführbare Prüfungen und
  verbleibende Risiken sind ausdrücklich benannt und gelten nicht als bestanden.
- Dokumentation, Beispiele und Entscheidungsprotokoll sind bei Bedarf aktualisiert.
- Diff ist geprüft und frei von unbeabsichtigten Dateien oder Geheimnissen.
- Abschlussbericht nennt Änderung, Prüfungen, Issue und Commit beziehungsweise PR.

## Quellen für den GitHub-Workflow

Stand der Dokumentationsprüfung: 2026-09-20. Projektspezifische Regeln oben sind
bewusste Arbeitskonventionen, keine Behauptung eines universellen Standards.

- [Issues und Pull Requests verknüpfen](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/using-keywords-in-issues-and-pull-requests)
- [GitHub Actions sicher verwenden](https://docs.github.com/en/actions/reference/security/secure-use)
