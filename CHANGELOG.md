# Changelog

Alle wesentlichen Änderungen dieses Projekts werden hier dokumentiert. Das
Format orientiert sich an [Keep a Changelog](https://keepachangelog.com/) und die
Versionierung an [Semantic Versioning](https://semver.org/).

## [1.5.0] - 2026-10-19

### Added

- Konforme lokale Tet4-Kantenverfeinerung mit Erhalt von Rand-/Materialtags und benannten Gruppen (#67).
- Kontrollierter adaptiver stationärer Wärmeablauf mit Toleranz-, DOF-, Schritt-, Stagnations- und Fehlerstatus sowie deterministischem Stufenexport (#68).
- Unabhängige glatte/lokalisierte Fehlervergleiche gegen uniforme Verfeinerung, analytische Materialgrenzen- und Randwertnachweise (#69).

### Limits

- Residualindikator ist heuristisch; ausschließlich affine Tet4 und stationäre Wärme.
- Kompatible Minor-Erweiterung; vorhandene Solver- und Model-Verträge bleiben bestehen.
- GitHub-gehostete Ubuntu-Runner prüfen Python 3.11–3.13, Wheel/sdist, Abhängigkeiten, PETSc und Handbücher.
- [Release Notes](Doc/releases/v1.5.0.md); Veröffentlichungsnachweise werden in #70 verknüpft.

## [1.4.0] - 2026-10-11

### Added

- Verifizierter generischer Unsicherheitsworkflow: explizite Verteilungen, gemeinsames Sampling, normale Korrelation, Statistik, lokale Sensitivität und stationäre FEM-Studien (#64).
- Vollständige deutschsprachige LaTeX-Gesamtdokumentation mit vier Fachbänden, ausführbaren Rezepten, API-/Quellendriftprüfung und reproduzierbarem Quellenarchiv (#206–#211).

### Compatibility and limits

- Kompatible Minor-Erweiterung; bestehende Model-Schemata 1–3 und Solver bleiben bestehen. Versionsgebundene Studien werden beim Update bewusst neu gerechnet.
- Experimentelle Validierung bleibt in #63/#65 außerhalb M9 offen gemäß ADR-054. Keine pauschale Validierungsfreigabe.
- [Release Notes](Doc/releases/v1.4.0.md) und [Freigabebericht](Doc/68_M9_RELEASE_VERIFICATION.md).

### Fixed

- Rampen interpolieren große endliche Zeitgrenzen und Amplituden ohne stillen
  Zwischenüberlauf; End-/Haltewerte und konstante Profile bleiben erhalten (#221).

## [1.3.0] - 2026-10-01

### Added

- Stationäre deklarative Parameterstudien mit Fingerprints, Resume und sichtbaren Fehlern (#59).
- Rückverfolgbare CSV/JSON-Kennwerttabellen und Markdown-Berichte mit getrennten Qualitätsaussagen (#60).
- Wärme- und Piezo-Studienrezepte mit analytischen Referenzen und sauberer Artefaktreproduktion (#61).
- Optionaler stationärer Tet4-Feldbildadapter mit Headless-Rendering und JPEG-Export (#189).
- Polynomiale FEM-Surrogate mit disjunkten Trainings-/Testfällen, Fehlerbericht und JSON-Artefakten (#205).

### Compatibility and limits

- Kompatible Minor-Erweiterung; Model-Schemata 1–3 und vorhandene Solver bleiben bestehen.
- Versionsgebundene Studienfingerprints invalidieren Ergebnisse beim Update bewusst.
- Keine dynamischen Studien, garantierten Surrogatfehler oder experimentelle Validierung.
- [Release Notes](Doc/releases/v1.3.0.md) und [Freigabebericht](Doc/64_M8_RELEASE_VERIFICATION.md).

## [1.2.0] - 2026-09-29

Release v1.2.0, veröffentlicht am 29.09.2026.

### Added

- Optionale PETSc-Direkt-/FieldSplit-Backends mit nativer real/complex-Parität (#34).
- Geradkantige Tet10/Triangle6 durch I/O, stationäre Felder/Kopplungen und
  transiente Wärme/Piezo-Dynamik; quadratische Sonden und Schwerpunktfelder (#35).
- Residualindikator und deterministische Markierung für stationäre Tet4-Wärme (#36).
- Reproduzierbare Profiling-Baseline und Fehler-/Kostenvergleich mit Rohdaten (#37).
- Supportmatrix und erweiterte Kandidaten-Artefaktprüfung (#38).

### Changed

- Tet4-Flächengeometrie wird bedarfsgerecht ausgewertet; gemessene Verbesserung
  und Regressionstoleranzen sind im Performancebericht belegt (#37).

### Compatibility and limits

- Tet4-Verträge, Model-Schemata 1–3, CPython 3.11–3.13 und MIT bleiben bestehen.
  PETSc bleibt optional und wird nur über direkte Python-Solver-APIs ausgewählt.
- Keine gekrümmten/hybriden Tet10-Netze, MPI-, Adaptivitäts- oder experimentelle
  Validierungszusage. Bestehende Dynamiklimits bleiben unverändert.
- [Release-Nachweis](Doc/58_M7_RELEASE_VERIFICATION.md).

## [1.1.0] - 2026-09-23

### Added

- Konsistente und gelumpte Tet4-Masse/Kapazität, quellenbelegte thermische
  Parameterbibliothek sowie Sprung-, Rampen- und CSV-Zeitprofile (#29).
- Backward-Euler-Wärmeleitung mit Model-Schema 2, Speicher-/Randbilanz,
  ParaView-Zeitserie und Messpunkt-CSV (#30).
- Piezoelektrische Modalanalyse mit mechanischer Masse und algebraischer
  Elektrodenkondensation; ungedämpftes Newmark und komplexe Harmonik (#31).
- Model-Schema 3, vollständige dynamische Feld-/Frequenzserien mit
  Real-/Imaginärteil, Amplitude/Phase und deterministischer Provenienz (#31).
- Release-Nachweis mit getrennten Raum-/Zeitfehlern, Energie- und
  Frequenzreferenzen; frische Wheel-/sdist-Prüfung aller alten und neuen Beispiele (#33).

### Compatibility

- Stationäre Python-API, Model-Schema 1 und Ergebnisbündel 1 bleiben erhalten.
  Neue Analysearten verwenden additive Schemata 2 und 3, keine Umdeutung alter Dateien.
- CPython 3.11–3.13; unveränderte NumPy/SciPy/meshio-Abhängigkeiten und MIT-Lizenz.

### Limitations

- Keine Rückwirkung von Verformung auf Temperatur (#32 zurückgestellt), keine
  Dämpfung im Piezo-Dynamikpfad und keine experimentelle Validierung.
- Dichte Piezo-Referenz bis 2048 aktive Gesamt-DOFs; gespeicherte Zeitverläufe
  begrenzt. Tet10/PETSc und allgemeine Großmodell-Skalierung folgen später.

## [1.0.0] - 2026-09-23

### Added

- Feste und schwebende Equipotential-Elektroden mit Gesamtladungsvorgabe.
- Versioniertes SI-TOML-Modell und stabile stationäre Model-/Solver-API mit
  explizitem Kompatibilitätsvertrag und Migration von direkten Solveraufrufen.
- Einseitiger stationärer Ablauf `T → (u, φ)`, regionsweise Materialien und
  perfekte Grenzflächen; explizite elektrisch-thermische Modelloption.
- Modellweiter Preflight für Topologie, Randbedingungen, Lagerung und Gauges.
- Versioniertes VTU-/JSON-Ergebnisbündel mit Rohfeldern, Einheiten, Bilanzen,
  Provenienz und Dateiprüfsumme; reproduzierbares End-to-End-Beispiel.
- Vier-Netze-Mehrmaterialnachweis, dokumentierte Nutzungsgrenzen und
  Abhängigkeitsmatrix; isolierte Wheel-/sdist-Prüfungen auf Python 3.11–3.13.
- MIT-Lizenz und Lizenzmetadaten für die aktuellen Distributionsarchive.

### Compatibility

- Direkte Solveraufrufe bleiben erhalten; Model-Schema und Ergebnisbündel stehen
  jeweils auf Version 1. Der [API-Vertrag](Doc/29_PUBLIC_API_CONTRACT.md) definiert
  die stabile Oberfläche und die Regeln für spätere inkompatible Änderungen.

### Limitations

- Kleine lineare stationäre Tet4-Modelle mit SciPy Sparse Direct; keine
  Transienten, PETSc, Tet10, reziproke Wärmekopplung oder experimentelle Validierung.


## [0.5.0] - 2026-09-22

### Added

- Piezoelektrisches Stress-Charge-Materialmodell mit Strain-Charge-Konversion,
  Achsrotation, SI-/Voigt-Validierung und erhaltener Eingangsdatenprovenienz.
- Gekoppelte symmetrisch-indefinite Tet4-Elementoperatoren, thermische
  Eigendehnung mit explizitem elektrischen Modellschalter sowie rohe Feld- und
  Enthalpieauswertung.
- Monolithischer quasistatischer `u-φ`-Referenzsolver mit feldweiser
  Kongruenzskalierung, Sparse LU, Lagerungs- und Gauge-Prüfung, getrennten
  Residuen und globalen Kraft-, Moment- und Ladungsbilanzen.
- Analytischer direkter/inverser d₃₃-Benchmark, einfacher 1D-Leerlauf,
  eigenständiger 1D-FE-Vergleich, Reziprozitäts- und Energietests sowie
  vierstufige gekoppelte Konvergenz.
- Ausführbares Piezo-Beispiel mit direkter und inverser Wirkung und VTU-Export.

### Fixed

- Die Rundungstoleranz der globalen Piezo-Ladungsbilanz berücksichtigt bei
  physikalisch verschwindender Nettoladung die ungekürzten elektrischen
  Matrixbeiträge.

### Limitations

- Nur quasistatische lineare Tet4-Modelle und ein direkter Sparse-LU-Solver;
  allgemeine leitende Floating Electrodes, Transienten, PETSc und
  experimentelle Validierung sind nicht enthalten.

## [0.4.0] - 2026-09-22

### Added

- Sequenzielle stationäre Thermomechanik: ein Wärmefeld wird auf demselben
  Tet4-Netz an den Elastizitätssolver übergeben.
- Anisotrope thermische Eigendehnung mit Referenztemperatur, regionsweisen
  Parametern und Rotation aus orthonormalen Materialhauptachsen.
- Thermischer Elementlastterm, mechanische Spannung und Energie sowie
  deterministische Provenienz des Temperatur- und Expansionsmodells.
- Ein ausführbares Thermomechanik-Beispiel und unabhängige Patch-, Bilanz-
  und vierstufige Konvergenznachweise.

### Limitations

- Die stationäre Kopplung verläuft ausschließlich `T → u`; reziproke
  Wärmequellen, Transienten und experimentelle Validierung sind nicht enthalten.

## [0.3.0] - 2026-09-21

### Added

- Feldweise Freiheitsgradzuordnung, Sparse-Assembly und symmetrieerhaltende
  Dirichlet-Elimination als gemeinsame Referenzinfrastruktur.
- Stationäre lineare Tet4-Elastizität mit validierter isotroper und anisotroper
  Steifigkeit, Volumen-, Flächen- und Punktlasten, Reaktionen, Energien und
  deterministischer Laufprovenienz.
- Stationäre Tet4-Wärmeleitung mit anisotroper Leitfähigkeit, Wärmequellen,
  nach außen positivem Wärmestrom, Konvektion, Leistungsbilanz und
  Nullraumdiagnostik.
- Tet4-Elektrostatik mit absoluter anisotroper Permittivität, Volumen- und
  Oberflächenladungen, Gauge-Diagnostik, E-/D-Feldern, Elektrodenladungen und
  Gauß-Bilanz.
- Verifizierte Zug-, Konvektionsquader- und Plattenkondensator-Beispiele sowie
  Patch-, Benchmark- und vierstufige Konvergenznachweise für alle drei Felder.

### Changed

- Der thermische Neumann-Term folgt nun explizit der Konvention eines positiv
  nach außen gerichteten Wärmestroms; die schwache Form dokumentiert den
  korrekten negativen Lastterm.

## [0.2.0] - 2026-09-21

### Added

- Gmsh-Import für Tet4/Triangle3 mit überlappenden physischen Gruppen und
  VTU-Export für Punkt-, Zell- und Auswahldaten.
- Unveränderliches, validiertes Mesh-Datenmodell für Punkte, Tet4-Zellen,
  Triangle3-Flächen, Tags, externe IDs und überlappende Auswahlsätze.
- Verifizierte Tet4-Volumina, Formfunktionsgradienten, Orientierungskorrektur,
  Flächengeometrie und dimensionslose Qualitätsindikatoren.

### Fixed

- Interne Gmsh-Metadaten werden nicht als fachliche Auswahlsätze importiert.

## [0.1.0] - 2026-09-20

### Added

- GitHub-Milestones M0–M7 mit überprüfbaren Arbeitspaketen und Release-Gates.
- Dokumentierter Versionszyklus von der Planungsbaseline `v0.0.1` bis `v1.2.0`.
- Installierbares Python-Paketgerüst mit `src`-Layout, Teststruktur und
  Entwicklungsinstallation für den M0-Zyklus.
- Lokale Qualitätsprüfungen mit pytest, Ruff und mypy sowie eine
  sicherheitsbeschränkte GitHub-Actions-Matrix für Python 3.11 bis 3.13.
- Typisierte SI-/Voigt-Konventionen, konkrete Validierungsfehler sowie ein
  deterministisches Provenienzschema mit installierter Paketversion.

[Unreleased]: https://github.com/maschin96/FEM-Toolbox/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.5.0...v1.0.0
[0.5.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.0.1...v0.1.0
