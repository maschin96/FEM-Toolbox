# Changelog

Alle wesentlichen Änderungen dieses Projekts werden hier dokumentiert. Das
Format orientiert sich an [Keep a Changelog](https://keepachangelog.com/) und die
Versionierung an [Semantic Versioning](https://semver.org/).

## Unreleased

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

[Unreleased]: https://github.com/maschin96/FEM-Toolbox/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.0.1...v0.1.0
