# Roadmap

Status: M0–M6 freigegeben; M7 konkretisiert (Stand 2026-09-23).
Spätere Meilensteine: siehe historische Produktvision und aktuelles GitHub-Backlog.

Die Roadmap priorisiert wissenschaftliches Risiko vor Funktionsmenge. Jeder
Meilenstein endet mit einer ausführbaren Demonstration und einem V&V-Nachweis.
Zusätzlich endet jeder Milestone mit einem eigenen Release-Gate und einem
annotierten Versions-Tag. Die verbindlichen Regeln stehen in
[09_RELEASE_PROCESS.md](09_RELEASE_PROCESS.md).

## M0 — Projektfundament

Zielversion: `v0.1.0`

GitHub: [Milestone M0](https://github.com/maschin96/FEM-Toolbox/milestone/1) ·
[Release-Issue #6](https://github.com/maschin96/FEM-Toolbox/issues/6)

Ergebnis: installierbares Python-Paket mit reproduzierbarer Entwicklungsumgebung.

- `pyproject.toml`, `src/`-Layout und Teststruktur
- CI für unterstützte Python-Versionen
- Formatierung, Linting, Typprüfung
- Konventionen für Einheiten, Voigt-Reihenfolge und Fehlerklassen
- Versions- und Provenienzmodell

Abnahme: Paketinstallation und leerer Smoke-Test in sauberer Umgebung.

## M1 — Netz und Tet4-Geometrie

Zielversion: `v0.2.0`

GitHub: [Milestone M1](https://github.com/maschin96/FEM-Toolbox/milestone/2) ·
[Release-Issue #10](https://github.com/maschin96/FEM-Toolbox/issues/10)

Ergebnis: validiertes internes Netzmodell und reine Tet4-Geometrieoperatoren.

- Mesh-Datenmodell, Tags und Auswahlsätze
- Gmsh-Import und VTU-Export
- Volumen, Gradienten, Flächen, Normalen
- Qualitätsindikatoren und Degenerationsprüfung

Abnahme: Geometrie-Invarianten und Rundreise Gmsh -> intern -> VTU.

## M2 — Ungekoppelte Einzelfelder

Zielversion: `v0.3.0`

GitHub: [Milestone M2](https://github.com/maschin96/FEM-Toolbox/milestone/3) ·
[Release-Issue #15](https://github.com/maschin96/FEM-Toolbox/issues/15)

Ergebnis: drei getrennte stationäre Solver.

- lineare Elastizität
- Wärmeleitung
- Elektrostatik
- gemeinsame Assembly, Lasten und Dirichlet-Behandlung
- Rohfelder und globale Bilanzen

Abnahme: Patch-Tests, analytische Benchmarks und Netzkonvergenz je Physik.

## M3 — Thermische Ausdehnung

Zielversion: `v0.4.0`

GitHub: [Milestone M3](https://github.com/maschin96/FEM-Toolbox/milestone/4) ·
[Release-Issue #18](https://github.com/maschin96/FEM-Toolbox/issues/18)

Ergebnis: sequenzielle stationäre Thermomechanik.

- anisotrope thermische Expansion
- freie und behinderte Expansion
- regionsweise Referenztemperatur
- Temperaturimport aus Wärmefeld

Abnahme: Nullspannung bei freier homogener Expansion und analytische Reaktionen
bei vollständiger Behinderung.

## M4 — Piezoelektrizität

Zielversion: `v0.5.0`

GitHub: [Milestone M4](https://github.com/maschin96/FEM-Toolbox/milestone/5) ·
[Release-Issue #23](https://github.com/maschin96/FEM-Toolbox/issues/23)

Ergebnis: monolithischer quasistatischer `u-phi`-Solver.

- Stress-Charge-Materialmodell und Achsrotation
- gekoppelte Elementblöcke
- Skalierung des Sattelpunktproblems
- direkte/inverse Benchmarks
- Elektroden mit vorgegebenem Potential

Abnahme: Reziprozität, Energietest, direkte und inverse Wirkung sowie unabhängige
Referenzrechnung.

## M5 — Elektroden und vollständiger stationärer Workflow

Zielversion: `v1.0.0`

GitHub: [Milestone M5](https://github.com/maschin96/FEM-Toolbox/milestone/6) ·
[Release-Issue #28](https://github.com/maschin96/FEM-Toolbox/issues/28)

Ergebnis: praxistaugliche Randbedingungen für Piezoaufbauten.

- Equipotential-Zwang über DOF-Kondensation oder Nebenbedingungen
- vorgegebene Gesamtladung und Floating Electrode
- kombinierter Workflow `T -> (u, phi)`
- Versioniertes TOML-Modellformat und dokumentierte Beispiele
- vollständiger Ergebnis- und Provenienzexport

Abnahme: Open-Circuit-/Short-Circuit-Benchmark, reproduzierbares End-to-End-
Beispiel und der [Vier-Netze-Konvergenznachweis mit Materialgrenze](25_THERMOPIEZO_CONVERGENCE.md).

## M6 — Transiente und dynamische Erweiterung

Freigegeben als v1.1.0; [Release-Nachweis](35_M6_RELEASE_VERIFICATION.md).

Zielversion: `v1.1.0`

GitHub: [Milestone M6](https://github.com/maschin96/FEM-Toolbox/milestone/7) ·
[Release-Issue #33](https://github.com/maschin96/FEM-Toolbox/issues/33)

Ergebnis: reproduzierbare transiente Wärmeleitung, lineare Piezoelastodynamik,
Modal-/Harmonik-Analyse. Reziproke Thermoelastizität wurde am 23.09.2026
auf Nutzerwunsch zurückgestellt.

1. #29: Masse/Kapazität, Materialdaten, Zeit-/Anfangswert- und Serienvertrag.
2. #30: Backward-Euler-Wärme samt zeitabhängigen Randwerten und Speicherbilanz.
3. #31: getrennte PRs für DAE/Modalanalyse, Newmark, Harmonik und Model/Export.
4. #32: auf Nutzerwunsch zurückgestellt; kein M6-Pflichtumfang.
5. #33: kompatibles Release mit ausführbaren Artefaktbeispielen.

Abnahme: getrennte Raum-/Zeitkonvergenz, Frequenzreferenzen, Energiebilanzen,
konsistente elektrische Zwangsbedingungen und lauffähige
v1.0-Beispiele. Model, Preflight, Provenienz und Zeit-/Frequenzexport gehören
zum Pflichtumfang. Der GitHub-Zieltermin 14.10.2026 ist ohne Kapazitätsplanung
nicht bestätigt. Details und PR-Schnitte: [Entwicklungsplan](30_M6_M7_DEVELOPMENT_PLAN.md).

## M7 — Skalierung und weitere Elemente

Zielversion: `v1.2.0`

GitHub: [Milestone M7](https://github.com/maschin96/FEM-Toolbox/milestone/8) ·
[Release-Issue #38](https://github.com/maschin96/FEM-Toolbox/issues/38)

Ergebnis: nachweisbar skalierbare Solverpfade und verifizierte Tet10-Genauigkeit.

1. #37a: Laufzeit-/Speicherbaseline vor Architektur- und Optimierungsentscheidungen.
2. #34: SciPy-Backend-Vertrag, optionales PETSc, FieldSplit und echte Backend-CI.
3. #35: geradkantige Tet10/Triangle6 durch Netz, Operatoren, M6, Model und Export.
4. #36: Fehlerindikator und Markierung zunächst für stationäre Tet4-Wärme.
5. #37b/c und #38: gemessene Optimierung, Supportmatrix und Release-Abnahme.

Abnahme: Backend-Parität, Tet10-Patch-/Konvergenztests, dokumentierter
Indikatornachweis und reproduzierbare Leistungs-/Speichermessungen bei
vergleichbarem Fehler. MPI, gekrümmte Tet10 und gekoppelte Adaptivität sind
keine impliziten Zusagen. Details: [Entwicklungsplan](30_M6_M7_DEVELOPMENT_PLAN.md).

## M8–M10 — Zweite Produktplanungswelle

Die [Bestandsanalyse und Produktvision](13_PRODUCT_VISION.md) dokumentiert
Restumfang, Abhängigkeiten und detaillierte Abnahme. Die damalige Baseline
ist historisch; M5 einschließlich #57/#58 ist abgeschlossen.
Für M6/M7 gilt nun der Entwicklungsplan vom 23.09.2026.

| Milestone | Zielversion (vorläufig) | Ergebnis | Issues / Release-Gate |
|---|---|---|---|
| [M8](https://github.com/maschin96/FEM-Toolbox/milestone/9) | v1.3.0 | reproduzierbare stationäre Parameterstudien, Kennwertberichte, Tet4-Feldbilder/JPEG und polynomiale Surrogatmodelle | #59–#61, #189, #205 / #62 |
| [M9](https://github.com/maschin96/FEM-Toolbox/milestone/10) | v1.4.0 | verifizierter Unsicherheitsworkflow und vollständige deutschsprachige Dokumentation | #64, #206–#211 / #66 |
| [M10](https://github.com/maschin96/FEM-Toolbox/milestone/11) | v1.5.0 | kontrollierte adaptive stationäre Wärme auf Tet4 | #67–#69 / #70 |

M8 setzt M5 voraus; M9 setzt M8 voraus. Nach der ausdrücklichen
Eigentümerentscheidung vom 10.10.2026 bleiben #63/#65 als experimenteller
Folgeumfang offen und sind keine v1.4.0-Pflichtgates mehr. Die generische
Unsicherheitssoftware wird mit analytischen Referenzen und klar bezeichneten
synthetischen Lehrfällen verifiziert; keine experimentelle Freigabe wird
behauptet. [Entscheidung und Recherche](66_M9_DATA_SEARCH_AND_RELEASE_OPTIONS.md).
M10 setzt #36,
M8 und ein positives Remeshing-Machbarkeitsgate voraus. Die Releasefolge bleibt
vorläufig M7 → M8 → M9 → M10; fachlich unabhängige Vorbereitung ist möglich.
Es werden keine Kalendertermine oder noch nicht belegten Validierungen zugesagt.

## M11–M12 — Grobe langfristige Milestones

| Milestone | Zielbild | Planungsstand |
|---|---|---|
| [M11 — Topologieoptimierung](https://github.com/maschin96/FEM-Toolbox/milestone/12) | Optimierer auf verifizierter FEM-Basis, zunächst lineare Strukturmechanik mit Volumenrestriktion, geprüften Sensitivitäten und Kontrolle der Netzabhängigkeit | Beschreibung ohne Einzelissues, Version oder Termin |
| [M12 — RANS-CFD und Fluid-Struktur-Kopplung](https://github.com/maschin96/FEM-Toolbox/milestone/13) | Ausschließlich RANS-Strömungsmodelle mit Turbulenzschließung und bidirektionaler Strukturkopplung; zunächst stationär, inkompressibel und isotherm | Beschreibung ohne Einzelissues, Version oder Termin |

Die [Milestone-Beschreibungen](13_PRODUCT_VISION.md) halten groben Umfang,
Abhängigkeiten, Nachweise und Grenzen fest. M12 schließt DNS, LES und hybride
Verfahren aus; dynamische FSI/URANS benötigt eine spätere Umfangsentscheidung.
Die Reihenfolge M10 → M11 → M12 ist ein Produkthorizont, keine technische
Abhängigkeit der CFD-Erweiterung vom Optimierer. Detaillierung, Architektur-
entscheidungen und Release-Issues folgen erst vor der jeweiligen Umsetzung.

## Definition of Done für jedes Arbeitspaket

Die folgenden Implementierungsnachweise gelten, soweit auf das Arbeitspaket
anwendbar. Reine Planung/Dokumentation prüft Inhalt, Links und Konsistenz;
nicht anwendbare numerische Nachweise werden begründet.

- öffentliche API und Annahmen dokumentiert
- Typen, Einheiten und Fehlersituationen definiert
- gezielte Unit-Tests vorhanden
- mindestens ein physikalischer Test für numerisches Verhalten vorhanden; bei
  reiner Dokumentation stattdessen Inhalt, Links und Konsistenz geprüft
- Bilanz-/Residualausgabe implementiert
- Beispiel ohne interne Kenntnisse ausführbar
- keine Optimierung ohne vorheriges Profiling

Zusätzlich auf Milestone-Ebene: Release-Issue abgeschlossen und Zielcommit auf
`main` annotiert getaggt. Ein einzelnes Arbeitspaket benötigt keinen eigenen
Milestone-Release.

M8-Pflichtumfang ist integriert und verifiziert. [Freigabestand v1.3.0](64_M8_RELEASE_VERIFICATION.md); keine transiente Studienausführung, Optimierung oder experimentelle Validierung.
