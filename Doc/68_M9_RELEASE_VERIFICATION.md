# M9 / v1.4.0 — Freigabenachweis (#66)

Stand: 11.10.2026. Veröffentlichung wird erst nach erfolgreicher finaler
main-CI, sauberem Paket-/Dokumentationsbuild und Tag-Zielprüfung in #66
bestätigt. Dieser Quellenbericht dokumentiert Voraussetzungen und Prüfplan;
Tag und endgültiger main-Commit werden im GitHub-Release und Issue #66
unveränderlich zugeordnet, um einen selbstreferenziellen Commit zu vermeiden.

## Umfang und abgeschlossene Voraussetzungen

Eigentümerentscheidung vom 10.10.2026 / ADR-054: verifizierter generischer
Unsicherheitsworkflow und vollständige Gesamtdokumentation. Experimentelle
Validierung bleibt #63/#65 außerhalb M9 offen. Kein synthetischer Test ersetzt
Messdaten. Verpflichtende Beiträge: #64 und #206–#211, integriert über
#239–#247. [Gesamtabnahme](latex/ACCEPTANCE.md) verbindet Fachabnahmen,
Rezeptausführung, Quellen-/API-/V&V-Prüfungen und vollständige Sichtprüfung.
[Release Notes](releases/v1.4.0.md) nennen Verhalten, Migration und Grenzen.
Ein unabhängiges Fachreview liegt nicht vor; CI und Eigenprüfung sind keine
unabhängige fachliche Freigabe.

## Paket- und Dokumentationsgates

- Vollständige Suite und Qualitätsprüfungen; CPython 3.11–3.13, neueste und
  älteste unterstützte Abhängigkeiten, native PETSc-real/complex-Parität.
- Wheel und sdist frisch installiert außerhalb des Checkouts; vorhandene
  stationäre und dynamische Beispiele sowie M9-Unsicherheitsrezept.
- Alle 23 archivierten Handbuchworkflows (26 Befehle) und Entwicklerbeispiel
  mit Visualisierungsextra, isoliertem Python und überprüften Artefaktzielen.
- Fünf vollständige PDFs, drei negative Buildfälle, vollständige Quellenmatrix,
  API-/V&V-Drift und echter Nachbau aus dem Quellenarchiv.
- Gesamtabnahme: 491 Kandidatenseiten visuell geprüft; finale Metadaten und
  gegenüber dem Kandidaten geänderte Seiten werden erneut geprüft.
- Expliziter `--release`-Build aus sauberem tatsächlichem main-Zielcommit;
  Paketversion 1.4.0, PDF-Titel, Manifest, Hashes und Release Notes konsistent.

## Release-Metadaten und Reihenfolge

1. Aktuelles main vor Integration übernehmen und resultierenden Stand prüfen.
2. Release-PR mit erfolgreicher anwendbarer CI in main integrieren.
3. Erfolgreiche CI am tatsächlichen main-Zielcommit bestätigen; finale
   Artefakte aus genau diesem sauberen Stand bauen und prüfen.
4. Annotiertes `v1.4.0` auf diesem Commit veröffentlichen; keine Umhängung.
5. GitHub Release mit Wheel, sdist, fünf PDFs, vollständigem Quellenarchiv,
   Buildmanifest, Handbuchprüfsummen und separaten Releaseprüfsummen erstellen.
6. Tatsächliche Nachweise in #66 verlinken; Release-Issue zuletzt schließen.

Keine PyPI-Veröffentlichung. Artefaktherkunft, MIT-Lizenz und externe
Quellenrechte stehen in den Release Notes. Historische M7/M8-Nachweise bleiben
historisch; es wird keine allgemeine experimentelle Validierung behauptet.

## TDD für den ausdrücklichen Releasestatus

Vier gezielte Tests scheiterten zunächst an der fehlenden Veröffentlichungs-
Metadatenfunktion. Nach Umsetzung bestehen alle 16 Build-Gatetests: Kandidaten
bleiben ausdrücklich Kandidaten; Releases benötigen sauberen Quellenzustand,
finale dreiteilige Version und vollständigen Commit. Diese Prüfung ersetzt
nicht die anschließende Kontrolle des tatsächlichen main- und Tag-Ziels.
