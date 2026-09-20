# Roadmap

Status: Entwurf

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
- YAML/TOML-Modellformat und dokumentierte Beispiele
- vollständiger Ergebnis- und Provenienzexport

Abnahme: Open-Circuit-/Short-Circuit-Benchmark und reproduzierbares End-to-End-
Beispiel.

## M6 — Transiente und dynamische Erweiterung

Zielversion: `v1.1.0`

GitHub: [Milestone M6](https://github.com/maschin96/FEM-Toolbox/milestone/7) ·
[Release-Issue #33](https://github.com/maschin96/FEM-Toolbox/issues/33)

Ergebnis: zeitabhängige Wärmeleitung und lineare Piezoelastodynamik.

- konsistente/lumped Massen- und Kapazitätsmatrizen
- Newmark- oder generalized-alpha-Verfahren
- implizite Wärmeintegration
- harmonische Analyse
- vollständig reziproke lineare Thermoelastizität

Abnahme: Eigenfrequenzen, Energiebilanz und Zeitschritt-Konvergenz.

## M7 — Skalierung und weitere Elemente

Zielversion: `v1.2.0`

GitHub: [Milestone M7](https://github.com/maschin96/FEM-Toolbox/milestone/8) ·
[Release-Issue #38](https://github.com/maschin96/FEM-Toolbox/issues/38)

Ergebnis: größere Forschungsmodelle und Genauigkeitsverbesserung.

- petsc4py-Backend und Blockpräconditioner
- parallele Assembly nach Profiling
- Tet10 als nächstes Volumenelement
- Fehlerindikatoren und perspektivisch Adaptivität

Abnahme: reproduzierbarer Skalierungstest und elementübergreifende Benchmarks.

## Definition of Done für jedes Arbeitspaket

- öffentliche API und Annahmen dokumentiert
- Typen, Einheiten und Fehlersituationen definiert
- gezielte Unit-Tests vorhanden
- mindestens ein physikalischer Test vorhanden
- Bilanz-/Residualausgabe implementiert
- Beispiel ohne interne Kenntnisse ausführbar
- keine Optimierung ohne vorheriges Profiling
- Release-Issue abgeschlossen und Zielcommit auf `main` annotiert getaggt
