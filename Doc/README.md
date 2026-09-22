# Planungsdokumente

Diese Dokumente definieren die wissenschaftliche und technische Grundlage der
FEM Toolbox. Bei Widersprüchen haben die expliziten Konventionen in der
mathematischen Modellbeschreibung Vorrang; Änderungen daran benötigen einen
Eintrag im Entscheidungsprotokoll.

## Lesereihenfolge

1. [01_SCOPE_AND_REQUIREMENTS.md](01_SCOPE_AND_REQUIREMENTS.md)
2. [02_MATHEMATICAL_MODEL.md](02_MATHEMATICAL_MODEL.md)
3. [03_TET4_FORMULATION.md](03_TET4_FORMULATION.md)
4. [04_SOFTWARE_ARCHITECTURE.md](04_SOFTWARE_ARCHITECTURE.md)
5. [05_VERIFICATION_AND_VALIDATION.md](05_VERIFICATION_AND_VALIDATION.md)
6. [06_ROADMAP.md](06_ROADMAP.md)
7. [07_RISKS_AND_DECISIONS.md](07_RISKS_AND_DECISIONS.md)
8. [08_REFERENCES.md](08_REFERENCES.md)
9. [09_RELEASE_PROCESS.md](09_RELEASE_PROCESS.md)

## Produktplanung

- [Bestandsanalyse, Produktvision und Planungswelle 2](13_PRODUCT_VISION.md)

## Implementierte Einzelfelder

- [Lineare statische Elastizität: API, Grenzen und Verifikation](10_LINEAR_ELASTICITY.md)
- [Stationäre Wärmeleitung: API, Bilanzen und Verifikation](11_STEADY_HEAT.md)

- [Elektrostatik: API, Ladungsbilanz und Verifikation](12_ELECTROSTATICS.md)
- [Stationäre sequenzielle Thermomechanik](14_THERMOMECHANICS.md)
- [Piezoelektrische Materialdaten, Konversion und Achsrotation](15_PIEZOELECTRIC_MATERIAL.md)

## Veröffentlichungen

- [v0.4.0 — Sequenzielle Thermomechanik](releases/v0.4.0.md)
- [v0.3.0 — Ungekoppelte Einzelfelder](releases/v0.3.0.md)
- [v0.2.0 — Netz und Tet4-Geometrie](releases/v0.2.0.md)
- [v0.1.0 — Projektfundament](releases/v0.1.0.md)

## Statuskonvention

- **Entwurf**: fachlich geplant, Implementierung noch offen
- **Implementiert**: Code vorhanden, Verifikation noch nicht vollständig
- **Verifiziert**: definierte Code-Verifikation bestanden
- **Validiert**: gegen geeignete experimentelle oder hochvertrauenswürdige
  Referenzdaten geprüft

Die Planungsdokumente 01–09 sind Entwürfe. Implementierungs- und
Verifikationsdokumente benennen ihren jeweiligen Nachweisumfang ausdrücklich.
