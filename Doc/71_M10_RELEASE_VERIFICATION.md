# M10 / v1.5.0 — Implementierungs- und Freigabenachweise (#70)

Stand 19.10.2026 · Paketversion **1.5.0**. Veröffentlichungsstatus und finale
Commit-/CI-/Artefaktnachweise: [#70](https://github.com/maschin96/FEM-Toolbox/issues/70).

Die folgenden lokalen Nachweise dokumentieren den Kandidaten vom 11.10.2026;
der Abschlusslauf unten ergänzt die endgültige Version.

## Erfüllter Implementierungsumfang

- #67: [ADR-055 und positives Machbarkeitsgate](69_ADAPTIVE_REMESHING_CONTRACT.md),
  konforme Kantenstern-Bisektion ohne externe Remeshing-Lizenz, Mehrgebiets-/Randgruppenerhalt.
- #68: [Kontrollierter Wärmeablauf](70_ADAPTIVE_HEAT_WORKFLOW.md), alle fünf
  Endstatus, Material-/Randvererbung, unveränderliche Stufen und deterministisches JSON.
- #69: glatte und lokalisierte analytische Referenzen mit je sieben Stufen
  adaptiv/uniform. Unabhängige Energiefehlerintegration; eingefrorene lokalisierte
  Fehlergrenze 1.5 sqrt(W K): 514 gegenüber 729 DOFs (29.5 % weniger).

[28 Messzeilen](evidence/m10/comparison.json) und
[Umgebung, numerische Quellhashes und Rohartefakthashes](evidence/m10/environment.json)
sind versioniert. Vollständige Rohrecords werden durch das öffentliche Beispiel
regeneriert; lokale Ausgabe `results/adaptive_heat/`. Sie sind keine experimentelle
Validierung. Die Fehlergrenze wurde nach Pilotexploration und vor dem finalen
Nachweis fixiert; keine blind gewählte Grenze und kein universeller Laufzeitvorteil.

## Lokale Verifikation

- Gesamte Suite am Kandidaten: **839 bestanden, 83 optional übersprungen**.
  Optionales PETSc wird damit nicht als lokal geprüft ausgegeben.
- Ruff-Prüfung und Formatcheck bestanden; striktes mypy über 70 Quelldateien bestanden.
- Exaktes API-Inventar, Modulabdeckung und generiertes Signaturkapitel synchronisiert.
- Wheel und sdist gebaut; jeweils in separater frischer CPython-3.13-Umgebung
  installiert. Isolierte `-I`-Imports stammen aus der jeweiligen Umgebung,
  Version `1.5.0.dev1`, nicht aus dem Checkout. Vier-Stufen-Beispiel je Artefakt
  ausgeführt; NumPy 2.5.4 und SciPy 1.18.1 in den Installationsumgebungen.
- Fünf Handbuch-PDFs gebaut, drei negative Dokumentfälle abgelehnt, eingebettete
  Schriften/Lesezeichen/Referenzen geprüft und Nachbau aus dem Quellenarchiv bestanden.
  Neue Kapitel in Theorie-, Anwendung-, API-, V&V-Band und Gesamthandbuch visuell
  geprüft. Die historischen Layoutwarnungen bleiben im Build-Manifest sichtbar;
  keine pauschale vollständige Neuabnahme sämtlicher historischer Seiten behauptet.

## Veröffentlichungscheckliste

#70 bleibt offen, bis folgende Veröffentlichungsarbeiten tatsächlich erfolgt sind:

- Die verpflichtenden Arbeitspakete sind integriert und die PR-/CI-Nachweise verlinkt.
- CI einschließlich CPython 3.11–3.13, sauberer Wheel-/sdist-Installation und
  optionaler Backend-Gates läuft am tatsächlichen Integrations-/main-Ziel grün.
- Unabhängiger Reviewerstatus wird dokumentiert; derzeit liegt keine unabhängige
  menschliche Freigabe dieses neuen M10-Umfangs vor.
- Paket/CI/Release Notes werden vom Entwicklungskandidaten auf final 1.5.0 umgestellt;
  finale Installations- und Dokumentationsartefakte stammen aus demselben Commit.
- Release-PR nach main, dessen CI, unveränderlicher annotierter v1.5.0-Tag und
  Veröffentlichung gemäß [Release-Prozess](09_RELEASE_PROCESS.md) sind abgeschlossen.

Die vorhandene lokale Löschung von `AGENTS.md` und `.DS_Store` wurden nicht in
M10 aufgenommen. Die Implementierung ändert weder alte Ergebnisse noch bestehende
Solver-/Model-/Persistenzverträge; die zusätzliche JSON-Kennung hat eigenes Schema 1.

## Abschlusslauf am 19.10.2026 (#70)

Der Kandidat wird auf Paketversion 1.5.0 umgestellt. Das Repository ist öffentlich;
CI und Handbuch verwenden ephemere GitHub-gehostete Ubuntu-24.04-Runner. Die
Python-Matrix, Artefaktinstallation, PETSc real/complex und Handbuch-Negativfälle
bleiben erhalten. Main-Handbücher verlangen saubere Release-Provenienz.
Die früheren selbst gehosteten Nachweise am Commit
`8820e0075be9c0a3d6db0660cfa6cfa8949ab604` waren erfolgreich; sie ersetzen
keine Prüfung des neuen Freigabecommits. Aktuelle PR-/main-Läufe, Tag und
Veröffentlichungsartefakte werden abschließend im Release-Issue #70 verknüpft.
Unabhängiges menschliches Fachreview liegt weiterhin nicht vor.

Lokaler Abschlusslauf am finalen Versionsstand: **841 bestanden, 83 optional
übersprungen**, Ruff/Lint und Formatprüfung bestanden, striktes mypy über
70 Quelldateien bestanden, Paketabhängigkeiten konsistent und
Dokumentationsabdeckung (196 Quellen) vollständig. Die ausgelassenen Tests
betreffen optionale Backends; PETSc wird separat in der CI geprüft.
