# Release-Artefakte und Integrationsregeln (#58)

Stand: 23. September 2026. Die technische Artefaktprüfung ist in
[CI](../.github/workflows/ci.yml) eingerichtet. Branchschutz ist gemäß der
Entscheidung des Repository-Eigentümers vom 23.09.2026 keine Voraussetzung für
Integration oder Freigabe (ADR-020). Die übrigen Release-Gates gelten weiterhin.

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

## Historischer Schutzstatus

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

## Entscheidung: Branchschutz entfällt als Pflicht

Am 23.09.2026 hat der Repository-Eigentümer ausdrücklich entschieden, den
Branchschutz aus den Anforderungen zu streichen. Damit entfallen die Pflicht
zur Aktivierung von Branchschutz/Rulesets, der API-Nachweis aktiver Regeln und
der administrative Negativtest eines technisch blockierten Merges. Ein
Tarifupgrade oder eine Änderung der Repository-Sichtbarkeit ist dafür nicht
vorgesehen. Der oben erfasste API-Status bleibt als historische Beobachtung
stehen und stellt keinen offenen Freigabeblocker dar.

Die Integration erfolgt weiterhin per Pull Request. Vor dem Merge werden die
anwendbaren CI-Ergebnisse am aktuellen Commit geprüft: Quality, Tests,
Compatibility sowie Build release artifacts und alle sechs Artifacts-Jobs.
Review-Anmerkungen müssen geklärt sein; unabhängiges Review wird gemäß
`AGENTS.md` angestrebt, fehlende Reviewer werden transparent benannt. Eine
technisch erzwungene Reviewfreigabe wird nicht verlangt. Ein unabhängiger
Reviewer ist für diesen Beitrag bislang nicht nachgewiesen.

Die Einhaltung liegt damit bei den integrierenden Personen; GitHub verhindert
abweichende Integrationen nicht automatisch. Die Pflicht zu erfolgreichen
Checks und die bestehenden Regeln gegen unaufgeforderte History-Rewrites
bleiben bestehen. Siehe [ADR-020](07_RISKS_AND_DECISIONS.md).

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

Die 17 CI-Jobs des [PR-Laufs](https://github.com/maschin96/FEM-Toolbox/actions/runs/35830968433)
bestanden am Commit `7e0da81`, darunter alle sechs Artefaktinstallationen.
Issue #58 kann nach Integration des PRs mit den angepassten Akzeptanzkriterien
abgeschlossen werden. Für #28 bleiben die übrigen Kriterien des
[Release-Prozesses](09_RELEASE_PROCESS.md) erforderlich; ein ungeschützter
Branch ist kein eigenständiger Blocker mehr.
