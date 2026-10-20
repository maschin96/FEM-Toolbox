# M11 / v1.6.0 — Implementierungs- und Freigabenachweis (#83)

19.10.2026 · Freigabevorbereitung; tatsächlicher Abschluss in
[Release-Issue #83](https://github.com/maschin96/FEM-Toolbox/issues/83).

## Umfang und Versionsentscheidung

#76–#82 wurden über [PR #284](https://github.com/maschin96/FEM-Toolbox/pull/284)
im M11-Integrationsbranch integriert. Der Commit 68d5e1f hat alle 20 CI-/Handbuch-
Prüfungen bestanden. Der Eigentümer hat die anschließende Integration von M11
nach main ausdrücklich beauftragt. M11 folgt damit auf v1.5.0 und erhält
v1.6.0: additive API, unveränderte existierende Persistenz und eigene neue
Schema-1-Ergebnisse begründen den kompatiblen Minor-Schritt.

Der Umfang umfasst Tet4-SIMP, den volumengewichteten Dichtefilter, vollständige
Sensitivitäten, mehrere Lastfälle und Volumenrestriktion, deklarative
Model-basierte Studien, deterministischen JSON-/VTU-Export und separate
Vollmaterial-Geometrie-Reanalyse. Alte Solver und Beispiele bleiben kompatibel.
Kein allgemeines CAD, keine Fertigungs-, Festigkeits- oder experimentelle Freigabe.

## Numerische und methodische Gates

[Vorab festgelegter Vertrag](83_M11_TOPOLOGY_CONTRACT.md),
[Messmatrix und Einschränkungen](85_M11_VERIFICATION.md) und
[API / ausführbarer Workflow](84_M11_TOPOLOGY_WORKFLOW.md) tragen die Abnahme.
Drei Netze und zwei Radien, Einzellastvariante, Mehrlastfälle, Wiederholung,
Gradientschrittweiten, Volumen, Kräfte-/Arbeits-/Energieprüfung und je Fall
mindestens eine akzeptierte Geometrievariante sind belegt. Unskalierte
Budgetstopps, getrennte reine Schwellen und zu große Geometrien bleiben als
No-go sichtbar. Ergänzte Verbindungszellen werden explizit ausgewiesen.
Die netzdefinierte Passivmenge begrenzt die Interpretation als reine h-Konvergenz.

Die historischen 1.6.0.dev1-Dateien und Hashes bleiben unverändert. Die
numerischen Solverquellen wurden für die Finalisierung nicht geändert.
Finale Paket-/Installationsnachweise erhalten eigene Dateien und überschreiben
keine Kandidatenbelege.

## Finale lokale Prüfung von Version 1.6.0

- Gesamtsuite: **868 bestanden, 83 optionale Tests übersprungen**.
- Ruff, Formatprüfung, striktes mypy über 74 Quelldateien und Dokumentationsabdeckung bestanden.
- Vollständige Drei-Netze-/Zwei-Radien-Matrix unter 1.6.0 erneut bestanden;
  alle C-Werte identisch zum Kandidaten, Wiederholung und Einzellastfall erfolgreich.
- Wheel und sdist getrennt frisch installiert; korrekter isolierter Import,
  Version 1.6.0, MIT-Metadaten/Lizenzdatei sowie der vollständige Studien- und
  Geometrieablauf geprüft. NumPy 2.5.4/SciPy 1.18.1; Kern benötigt weder
  Matplotlib noch PETSc.
- [Finale lokale Belege und unveränderte numerische Quellhashes](evidence/m11/v1.6.0-local-verification.json).
  Dies sind Vorabprüfungen, nicht die später aus main erzeugten Veröffentlichungsartefakte.

## Integrations- und Veröffentlichungsgate

Vor dem main-Merge sind die finalen lokalen Regressionen/Qualitätsprüfungen,
Wheel-/sdist-Installationen, unveränderten Beispiele, Dokumentationsabdeckung,
fünf Handbücher samt Negativfällen/Quellenarchiv-Nachbau und aktuelle Release-PR-CI
zu prüfen. Unabhängiges menschliches Fachreview liegt nicht vor; der Status
wird gemäß Release-Prozess transparent dokumentiert und nicht als Freigabe ausgegeben.
Branchschutz ist gemäß ADR-020 kein Pflichtgate; CI am Zielcommit wird tatsächlich geprüft.

Nach dem Merge sind main-CI und Release-Handbuchmanifest am exakten Zielcommit
zu prüfen. Erst danach dürfen annotierter v1.6.0-Tag und Veröffentlichung
entstehen. #83 wird zuletzt geschlossen. Tatsächliche Commit-, Lauf-, Tag-
und Artefaktlinks werden im Release-Issue sowie der veröffentlichten
RELEASE_VERIFICATION.json gebunden. Diese Vorbereitungsdatei behauptet keinen
bereits erfolgten Tag oder eine bereits abgeschlossene Veröffentlichung.
