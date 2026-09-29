# M7 / v1.2.0 — finaler Freigabestand (#38)

Stand: 29.09.2026. Paketversion **1.2.0**, noch kein Tag und keine Veröffentlichung.
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

## Finaler Versionswechsel und Prüfung

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

Die finale Versionsnummer allein ist keine Veröffentlichung. Sie ermöglicht,
genau die späteren Distributionsmetadaten vor dem Release zu testen.
Kein PyPI-Upload, kein Tag und keine Release-Assets werden durch diesen PR erzeugt.

## Integration und noch offene Gates

Das aktuelle `origin/main` wurde am 29.09.2026 als Vorfahr des M7-Standes geprüft;
zu diesem Zeitpunkt waren keine zusätzlichen main-Änderungen zu übernehmen.
Vor dem späteren Release-Merge erneut fetch/Abgleich durchführen. Falls main
weiterläuft, Änderungen zuerst in M7 integrieren und erneut prüfen.

1. Finalisierungs-PR nach `codex/milestone-m7` integrieren, nachdem dessen
   anwendbare CI grün ist und Review-Anmerkungen geklärt sind.
2. Release-PR von M7 nach main erstellen und aktuellen Stand prüfen. Keine
   unfertige Funktionalität auf main übernehmen.
3. Nach dem Merge die CI am tatsächlichen main-Zielcommit prüfen.
4. Erst danach unveränderliches annotiertes `v1.2.0`-Tag setzen, GitHub Release
   und dazu passende geprüfte Wheel-/sdist-Artefakte mit Prüfsummen veröffentlichen.
5. Veröffentlichungsstatus in der Dokumentation nachtragen, #38 und M7 abschließen.

Ein unabhängiger Reviewer ist nicht nachgewiesen; eine solche Freigabe wird
nicht behauptet. Der [Release-Prozess](09_RELEASE_PROCESS.md) bleibt verbindlich.
