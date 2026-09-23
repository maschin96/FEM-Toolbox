# Release-Artefakte und Repository-Schutz (#58)

Stand: 23. September 2026. Die technische Artefaktprüfung ist in
[CI](../.github/workflows/ci.yml) eingerichtet. **Das administrative Gate ist
offen:** Diese Konfiguration allein schützt `main` nicht und gibt v1.0 nicht frei.

## Installation ohne Quellbaum

`Build release artifacts` baut mit `python -m build` in einer frischen
Linux-Umgebung. Der Build verwendet isolierte Build-Abhängigkeiten und erzeugt
zuerst ein sdist, daraus das Wheel. Beide Archive, SHA-256-Prüfsummen und das
eigenständige Beispiel `examples/steady_heat.py` werden sieben Tage als
`release-candidate` im jeweiligen Actions-Lauf aufbewahrt. Das ist ein
CI-Artefakt, kein veröffentlichtes Release und keine Signatur/Attestation.

Sechs nachgelagerte Jobs installieren jeweils genau ein Archiv in eine neue
virtuelle Umgebung: Wheel und sdist auf CPython 3.11, 3.12 und 3.13 unter
Linux. Diese Jobs checken das Repository **nicht** aus. Quellbaum, Tests und
editierbare Installation stehen ihnen nicht zur Verfügung. Die sdist-Installation
baut selbst mit isolierten Build-Abhängigkeiten. Prüfsummen sichern den Transfer;
`pip check`, Metadaten, Importpfad innerhalb der virtuellen Umgebung sowie
`pip freeze` machen Installation und Abhängigkeiten nachvollziehbar.

Der Solver läuft mit `python -I` außerhalb des Downloadverzeichnisses. Damit
werden auch `PYTHONPATH`, Benutzerpakete und das Skriptverzeichnis als
Importquellen ausgeschlossen. Das vorhandene analytische Wärmebeispiel prüft
`T = 400 - 40*x K`, 960 W Konvektionsleistung und die Leistungsbilanz und
schreibt ein nichtleeres VTU. Die vorhandenen absoluten Toleranzen sind
3e-12 K bzw. 1e-10 W bei relativer Toleranz null, passend zum kleinen linearen
Referenzsystem und dessen Rundungsfehlern.

Es gibt keine neue Plattformzusage: Linux/Python 3.11–3.13 ist die CI-Matrix;
macOS/ARM64/Python 3.13 ist ein lokaler Nachweis. Windows und andere Kombinationen
bleiben ungeprüft, siehe [Installation](../INSTALL.md). Änderungen an numerischen
Operatoren erfolgen hier nicht; die bestehende Gesamtsuite einschließlich
Konvergenztests bleibt separat erforderlich.

## Schutzstatus und erforderliche Integrationseinstellungen

Am 23.09.2026 wurden folgende API-Abfragen durchgeführt:

```sh
gh api repos/maschin96/FEM-Toolbox --jq '{private,permissions}'
gh api repos/maschin96/FEM-Toolbox/branches/main --jq '{protected,protection}'
gh api repos/maschin96/FEM-Toolbox/rulesets
```

Ergebnis: Repository privat, `permissions.admin=true`, `protected=false`,
`protection.enabled=false`, keine erforderlichen Statuschecks. Die Ruleset-Abfrage
antwortet HTTP 403: `Upgrade to GitHub Pro or make this repository public to enable
this feature.` Fehlende Administratorrechte sind damit nicht die beobachtete
Ursache. Sichtbarkeit und Tarif wurden nicht geändert.

Nach Bereitstellung der GitHub-Funktion muss ein Administrator Schutz für `main`
aktivieren und per API erneut belegen:

- Pull Request mit mindestens einer unabhängigen Freigabe; veraltete Freigaben
  nach neuen Commits verwerfen, offene Review-Diskussionen vor Merge klären.
- Aktueller Branch und erfolgreiche erforderliche Checks: `Quality / Python 3.13`,
  `Tests / Python 3.11`, `Tests / Python 3.12`, `Tests / Python 3.13`, alle sechs
  `Compatibility / Python <Version> / <oldest|latest>` und alle sechs
  `Artifacts / Python <Version> / <wheel|sdist>` sowie `Build release artifacts`.
- Force-Push und Branchlöschung unterbinden; Administrator-Bypass nicht als
  regulären Integrationsweg zulassen.

Der Nachweis muss die tatsächlich zurückgegebenen Checknamen und aktiven Regeln
enthalten. Ein unabhängiger Reviewer ist bislang nicht nachgewiesen.
Ein absichtlich fehlschlagender Check auf einem temporären PR muss nach Aktivierung
nachweislich den Merge blockieren. Diesen administrativen Negativtest kann ein
ungeschützter Branch nicht erfüllen; er bleibt offen. Anschließend Teständerung
entfernen und alle Checks am endgültigen Commit erneut bestehen lassen.

Die CI nutzt SHA-fixierte Actions, `contents: read` und `pull_request`, ohne
Secrets oder privilegierte PR-Ausführung. Matrixwerte sind statisch; PR-Titel
oder andere fremde Texte werden nicht als Shell-Code eingesetzt.

## Lokaler Nachweis und Freigabegrenze

Wheel und sdist wurden mit `uv build --out-dir /tmp/fem-58-dist` gebaut und
jeweils in einer separaten frischen Umgebung installiert. Das kopierte
Wärmebeispiel lief außerhalb des Repositorys mit `python -I`: 320–400 K,
960 W und Bilanzfehler 1.023e-12 W. `pip check` bestand für die sdist-Umgebung.
Als Negativkontrolle scheiterte das Beispiel trotz installierter NumPy/SciPy/
meshio-Abhängigkeiten ohne Toolbox mit `ModuleNotFoundError: fem_toolbox`.
Das anschließende Installieren des Wheels ließ denselben Test bestehen.
Dies ist ein Installationstest; es wird kein neuer physikalischer TDD-Zyklus
beansprucht.

Issue #58 und damit die Freigabe in #28 bleiben bis zum aktivierten, per API
belegten Schutz und zum nachgewiesenen Review-/Negativtest offen. Erfolgreiche
Artefaktjobs allein schließen dieses Gate nicht.
