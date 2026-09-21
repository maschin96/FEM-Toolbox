# Changelog

Alle wesentlichen Änderungen dieses Projekts werden hier dokumentiert. Das
Format orientiert sich an [Keep a Changelog](https://keepachangelog.com/) und die
Versionierung an [Semantic Versioning](https://semver.org/).

## Unreleased

### Added

- Unveränderliches, validiertes Mesh-Datenmodell für Punkte, Tet4-Zellen,
  Triangle3-Flächen, Tags, externe IDs und überlappende Auswahlsätze.
- Verifizierte Tet4-Volumina, Formfunktionsgradienten, Orientierungskorrektur,
  Flächengeometrie und dimensionslose Qualitätsindikatoren.

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

[Unreleased]: https://github.com/maschin96/FEM-Toolbox/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/maschin96/FEM-Toolbox/compare/v0.0.1...v0.1.0
