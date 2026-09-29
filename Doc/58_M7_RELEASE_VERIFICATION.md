# M7 / v1.2.0 — finaler Freigabestand (#38)

Status: **v1.2.0 veröffentlicht am 29.09.2026.**
Tag: [`v1.2.0`](https://github.com/maschin96/FEM-Toolbox/releases/tag/v1.2.0).
main-Zielcommit: `ee3465f8d3297e36c0a4709078dc5744aac2980c`.
CI-Lauf des Zielcommits: [36528677242](https://github.com/maschin96/FEM-Toolbox/actions/runs/36528677242),
19/19 Prüfungen erfolgreich.
Basis ist der integrierte Kandidat aus [PR #201](https://github.com/maschin96/FEM-Toolbox/pull/201),
M7-Commit `a2c24c886a0fd38b0a5b9586b95b9bc9a9f39201`.
Dieser Schnitt setzt die finale Version und aktualisiert die nutzerseitige
Release-Dokumentation. Produktionssolver, numerische Toleranzen und Funktionsumfang
sind gegenüber 1.2.0rc1 unverändert.

## Funktionsumfang und Nachweise

- [Supportmatrix und Installation](57_M7_RELEASE_CANDIDATE.md): Tet4/Tet10,
  stationäre Analysen und M6, SciPy/PETSc direkt/FieldSplit; PETSc nur an direkten
  Python-APIs. Die Kandidaten-Testergebnisse bleiben als historische Evidenz erhalten.
- [Tet10-Abnahme](55_TET10_ACCEPTANCE.md): Topologie, Operatoren, Model, Export,
  Vier-Netz-Nachweise, Bilanzen und Dynamik. #35 ist nach Integration von PR #199
  und #201 explizit geschlossen.
- [Performance-Abnahme](56_PERFORMANCE_ACCEPTANCE.md): drei integrierte Mess-/
  Optimierungsschritte, 159 archivierte Läufe und Grenzen; #37 ist geschlossen.
- PETSc #34 und Wärmeindikatoren #36 sind abgeschlossen; die jeweiligen
  [Backend-](38_PETSC_BACKEND.md), [FieldSplit-](39_PETSC_FIELDSPLIT.md) und
  [Indikatornachweise](43_HEAT_INDICATOR.md) bleiben maßgeblich.

Der Kandidaten-PR bestand alle 19 CI-Prüfungen. Der zunächst fehlgeschlagene
Python-3.11-sdist-Job fand beim Paketindexabruf keine SciPy-Distribution; eine
Wiederholung desselben Jobs ohne Codeänderung bestand am 29.09.2026.
[CI-Lauf](https://github.com/maschin96/FEM-Toolbox/actions/runs/36411692616).
Dies ist die Evidenz für den Kandidaten, nicht die CI am künftigen main-Zielcommit.

## Finaler Versionswechsel und Artefaktprüfung

TDD: Der Pakettest wurde zuerst auf 1.2.0 gesetzt und scheiterte erwartungsgemäß
an `1.2.0rc1 != 1.2.0`. Nach Anpassung der kanonischen Version in `pyproject.toml`
und Neuinstallation wird die installierte Version geprüft. Das Artefaktgate
verlangt ebenfalls ausdrücklich 1.2.0; die sechs frischen Wheel-/sdist-Jobs
prüfen Version, Lizenz, Betrieb ohne PETSc sowie alle 17 Beispiele je Job.
Native PETSc-real/complex-Ausführung bleibt Teil derselben CI.

Lokale Prüfung der finalen Version: **605 Tests bestanden, 83 optionale
PETSc-Tests übersprungen**. Ruff, Formatprüfung, Mypy und Diff-Prüfung bestanden.
Wheel und sdist wurden isoliert gebaut; Version 1.2.0, MIT-Metadaten und
Lizenzdateien wurden in beiden Archiven geprüft. Die nativen PETSc- und frischen
Installationsnachweise folgen in der CI des finalen PR-Stands.

Die finale Versionsnummer wurde vor dem Release geprüft. Tag und GitHub Release
wurden nach bestandener CI am main-Zielcommit erzeugt; die ausgelieferten
Distributionsmetadaten entsprechen 1.2.0. Es gab keinen PyPI-Upload.

## Integration und Release-Abschluss

- Das zum M7-Zeitpunkt aktuelle `main` war Vorfahr des Milestone-Branches;
  keine zusätzlichen main-Änderungen waren zu übernehmen.
- [Release-PR #203](https://github.com/maschin96/FEM-Toolbox/pull/203) wurde
  integriert. Die CI am tatsächlichen main-Zielcommit `ee3465f` bestand mit
  19/19 Prüfungen, darunter Tests für Python 3.11–3.13, sechs frische
  Wheel-/sdist-Installationen und native PETSc-real/complex-Jobs.
- Das annotierte Tag `v1.2.0` zeigt auf genau diesen Commit.
- Das [GitHub Release](https://github.com/maschin96/FEM-Toolbox/releases/tag/v1.2.0)
  enthält `fem_toolbox-1.2.0-py3-none-any.whl`,
  `fem_toolbox-1.2.0.tar.gz` sowie `SHA256SUMS`. Die Prüfsummen stimmen mit den
  Artefakten aus dem genannten CI-Lauf überein.
- Kein unabhängiger Reviewer ist für den Release-PR verzeichnet; dies wird
  nicht als unabhängige Freigabe dargestellt.

Issue #38 kann nach Merge dieses Statusnachtrags geschlossen werden. Der Release
und das Tag sind veröffentlicht; ein PyPI-Upload fand nicht statt.
