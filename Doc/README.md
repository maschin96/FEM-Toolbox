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

- [Aktueller Entwicklungsstand und konkrete M6/M7-Planung](30_M6_M7_DEVELOPMENT_PLAN.md)

- [Installation, Gmsh-Einrichtung und erster Lauf](../INSTALL.md)
- [Bestandsanalyse, Produktvision und Planungswelle 2](13_PRODUCT_VISION.md)
- [Weitergabe und Drittbestandteile](27_DISTRIBUTION_INVENTORY.md)
- [Abhängigkeitsverträglichkeit und Sicherheitsupdates](28_DEPENDENCY_COMPATIBILITY.md)
- [Release-Artefakte und Integrationsregeln](29_RELEASE_ARTIFACTS_AND_PROTECTION.md)

## Implementierte Einzelfelder

- [Lineare statische Elastizität: API, Grenzen und Verifikation](10_LINEAR_ELASTICITY.md)
- [Stationäre Wärmeleitung: API, Bilanzen und Verifikation](11_STEADY_HEAT.md)

- [Elektrostatik: API, Ladungsbilanz und Verifikation](12_ELECTROSTATICS.md)
- [Stationäre sequenzielle Thermomechanik](14_THERMOMECHANICS.md)
- [Piezoelektrische Materialdaten, Konversion und Achsrotation](15_PIEZOELECTRIC_MATERIAL.md)
- [Gekoppeltes piezoelektrisches Tet4-Element](16_PIEZOELECTRIC_ELEMENT.md)
- [Monolithischer piezoelektrischer Referenzsolver](17_PIEZOELECTRIC_SOLVER.md)
- [Direkter und inverser d₃₃-Nachweis](18_PIEZO_D33_VERIFICATION.md)
- [Materialregionen und perfekte Grenzflächen](19_MATERIAL_REGIONS.md)
- [Equipotential- und Floating-Electrode-Randbedingungen](20_ELECTRODES.md)
- [Versioniertes TOML-Modellformat und Model-API](21_MODEL_SCHEMA.md)
- [Stationärer Ablauf T → (u, φ)](22_THERMOPIEZOELECTRIC.md)
- [Modellweite Randbedingungs- und Nullraumprüfung](23_MODEL_PREFLIGHT.md)
- [Versioniertes stationäres Ergebnisbündel](24_RESULT_BUNDLE.md)
- [Vier-Netze-Nachweis des stationären Mehrmaterialmodells](25_THERMOPIEZO_CONVERGENCE.md)
- [Gültigkeitsbereich und numerische Nutzungsgrenzen für v1.0](26_VALIDITY_AND_LIMITS.md)
- [Öffentliche stationäre Python-API und Migration ab v0.3](29_PUBLIC_API_CONTRACT.md)

## M6-Entwicklung

- [Tet4-Masse und Wärmekapazität](31_TET4_STORAGE.md)
- [Quellenbelegte thermische Parameterbibliothek](32_MATERIAL_LIBRARY.md)
- [Transiente Wärme, Model-Schema 2 und Zeitserien](33_TRANSIENT_HEAT.md)
- [Piezoelastodynamik und Modalanalyse](34_PIEZO_DYNAMICS.md)

## Veröffentlichungen

- [v1.0.0 — Vollständiger stationärer Workflow](releases/v1.0.0.md)
- [v0.5.0 — Piezoelektrizität](releases/v0.5.0.md)
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
