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

## M6 — v1.1.0

- [Tet4-Masse und Wärmekapazität](31_TET4_STORAGE.md)
- [Quellenbelegte thermische Parameterbibliothek](32_MATERIAL_LIBRARY.md)
- [Transiente Wärme, Model-Schema 2 und Zeitserien](33_TRANSIENT_HEAT.md)
- [Piezoelastodynamik und Modalanalyse](34_PIEZO_DYNAMICS.md)

- [M6-Release-Verifikation, Migration und Grenzen](35_M6_RELEASE_VERIFICATION.md)

## M7 — v1.2.0 (in Entwicklung)

- [SciPy-/Tet4-Performance-Baseline (#37a)](36_PERFORMANCE_BASELINE.md)
- [Interner linearer Backend-Vertrag und SciPy-Referenz (#34a)](37_LINEAR_BACKEND.md)
- [Optionaler PETSc-Direktbackend und Scalar-Parität (#34b)](38_PETSC_BACKEND.md)
- [Iterativer PETSc-Backend mit FieldSplit/Schur (#34c)](39_PETSC_FIELDSPLIT.md)

- [Tet10-/Triangle6-Topologie und I/O (#35a)](40_TET10_TOPOLOGY.md)

- [Tet10-Formfunktionen, Geometrie und Quadratur (#35b)](41_TET10_GEOMETRY.md)

- [Tet10-Feldoperatoren, Speicher und Randintegrale (#35c)](42_TET10_OPERATORS.md)

- [Residualindikator und deterministische Markierung für Tet4-Wärme (#36)](43_HEAT_INDICATOR.md)

- [Gemessene bedarfsgerechte Tet4-Flächengeometrie (#37b)](44_GEOMETRY_PERFORMANCE.md)

- [Tet10-Wärme: Model, Backend und Schwerpunkt-Export (#35d, thermischer Teil)](45_TET10_HEAT_INTEGRATION.md)

- [Tet10-Elektrostatik: Ladungen, Energie und Schwerpunkt-Export (#35d)](46_TET10_ELECTROSTATICS.md)

- [Statische Tet10-Elastizität: Biegung, Energie und Schwerpunkt-Export (#35d)](47_TET10_ELASTICITY.md)

- [Stationäre Tet10-Thermomechanik: Eigenverzerrungen und gekoppelter Export (#35d)](48_TET10_THERMOMECHANICS.md)

- [Stationäre Tet10-Piezoelektrizität: Elektroden, Model und Enthalpie (#35d)](49_TET10_PIEZOELECTRIC.md)

- [Stationäre Tet10-Thermopiezoelektrizität: T → (u, φ), Model und Enthalpie (#35d)](50_TET10_THERMOPIEZOELECTRIC.md)

- [Tet10-Modalanalyse: Masse, Elektroden und Modenexport (#35d)](51_TET10_MODES.md)

- [Tet10-Newmark: Zeitintegration, Elektroden und Zeitserien (#35d)](52_TET10_NEWMARK.md)

- [Tet10-Harmonik: komplexe Frequenzantworten und Export (#35d)](53_TET10_HARMONIC.md)

- [Genauigkeit und Kosten: Tet4/Tet10 × SciPy/PETSc (#37c)](54_ACCURACY_COST_COMPARISON.md)

- [Tet10-Gesamtabnahme: Kriterien, Analyseumfang und Grenzen (#35)](55_TET10_ACCEPTANCE.md)

- [Performance-Gesamtabnahme: Rohdatenabgleich und Grenzen (#37)](56_PERFORMANCE_ACCEPTANCE.md)

- [M7-Release-Kandidat 1.2.0rc1 und Supportmatrix (#38)](57_M7_RELEASE_CANDIDATE.md)

- [M7: finaler Release-Nachweis und Veröffentlichungsgates (#38)](58_M7_RELEASE_VERIFICATION.md)

## Veröffentlichungen

- [v1.2.0 — Release Notes, Veröffentlichung vorbereitet](releases/v1.2.0.md)

- [v1.1.0 — Transiente Wärme und Piezo-Dynamik](releases/v1.1.0.md)

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

## M8 — in Entwicklung

- [Stationäre Parameterstudien, Fallkatalog und Fortsetzung (#59)](59_PARAMETER_STUDIES.md)

- [Kennwerttabellen und Studienberichte (#60)](60_STUDY_REPORTS.md)

- [Surrogatmodelle: Training, Prüfung und Modellartefakte (#205)](61_SURROGATE_MODELS.md)

- [Stationäre Tet4-Feldvisualisierung und JPEG-Export](62_FIELD_VISUALIZATION.md)
