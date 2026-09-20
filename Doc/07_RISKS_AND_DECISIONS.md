# Risiken und Entscheidungen

Status: Entwurf

## 1. Entscheidungsprotokoll

### ADR-001: Tet4 als erstes Element

Entscheidung: Der Referenzkern startet mit dem linearen Vierknoten-Tetraeder.

Begründung: einfache affine Abbildung, konstante Gradienten, exakt integrierbare
lineare Operatoren und robuste automatische 3D-Vernetzung. Die bekannte Schwäche
in Biegung und nahe Inkompressibilität wird akzeptiert, dokumentiert und durch
Konvergenztests kontrolliert. Tet10 ist die geplante nächste Elementfamilie.

### ADR-002: Stress-Charge-Form als interne Piezo-Konvention

Entscheidung: Intern werden `C^E`, `e` und `kappa^S` verwendet.

Begründung: Diese Form führt direkt zur gewählten schwachen Form. Materialdaten in
Strain-Charge-Form (`s^E`, `d`, `kappa^T`) werden in einem separaten, getesteten
Konverter transformiert. Originaldaten und transformierte Daten bleiben in der
Provenienz erhalten.

### ADR-003: Feldweise DOFs und monolithische Piezo-Lösung

Entscheidung: DOFs sind feldweise verwaltet; Mechanik und Elektrostatik werden
monolithisch gelöst.

Begründung: Der direkte und inverse Effekt sind reziprok und stark gekoppelt.
Eine monolithische Referenz vermeidet Iterationstoleranzen einer Partitionierung
und bewahrt die Blocksymmetrie.

### ADR-004: Wärme zunächst einseitig gekoppelt

Entscheidung: Version 0.1 löst Wärme unabhängig und übergibt `T` an Mechanik und
optional das explizite thermo-piezoelektrische Materialmodell.

Begründung: Stationäre thermische Ausdehnung deckt den ersten Anwendungsfall ab.
Vollständige reziproke Thermoelastizität ist eine transiente thermodynamische
Kopplung und wird gemeinsam mit den Kapazitätstermen implementiert.

### ADR-005: SI-Einheiten ohne implizite Konversion

Entscheidung: Der Rechenkern verwendet ausschließlich SI.

Begründung: Piezoelektrische Datensätze kombinieren leicht `GPa`, `pC/N`,
relative Permittivität und Millimetergeometrien. Eine versteckte Mischung kann
plausible, aber um Größenordnungen falsche Resultate liefern.

### ADR-006: SciPy-Referenzbackend vor PETSc

Entscheidung: Zunächst SciPy Sparse Direct; PETSc folgt hinter einer Solver-
Schnittstelle.

Begründung: Für Verifikation ist geringe Komplexität wertvoller als frühe
Skalierung. Die Blockstruktur und Feldindizes werden dennoch von Anfang an
erhalten.

## 2. Technisches Risikoregister

| Risiko | Wirkung | Früher Nachweis | Gegenmaßnahme |
|---|---|---|---|
| Voigt-/Schubkonvention verwechselt | falsche Kopplung, scheinbar plausible Werte | Einzelelementtest mit einzelnen Komponenten | eine interne Konvention, getestete Konverter |
| Piezo-Vorzeichen falsch | direkter/inverser Effekt widersprüchlich | Reziprozitäts- und Energietest | schwache Form aus einem Potential ableiten |
| Materialachsen falsch rotiert | falsche Polarisationsrichtung | Rotationsbenchmark 90 Grad | Tensorrotation zentral implementieren |
| elektrische Gauge fehlt | singuläres System | Nullraumerkennung | Referenzpotential oder Nullraumsolver |
| Floating Electrode falsch modelliert | falsche Open-Circuit-Spannung | 1D-Elektrodenbenchmark | Equipotential plus Gesamtladungsbedingung |
| schlechte Blockskalierung | Solverfehler trotz korrekter Matrix | Konditions-/Residualbericht je Feld | physikalische Blockskalierung |
| degenerierte Tetraeder | instabile Gradienten | Qualitätsprüfung beim Import | ablehnen mit Element-ID und Maßzahl |
| Tet4 in Biegung/nahe `nu=0.5` | zu steife Lösung | Netz-/Elementvergleich | Warnung, Verfeinerung, später Tet10/mixed |
| unklare Gmsh-Gruppen | falsche Material-/Randzuordnung | Import-Rundreisetests | eigene Auswahlsätze und Validierung |
| gemittelte Spannung als „genau“ interpretiert | falsche lokale Schlussfolgerungen | Metadaten im Export | roh/gemittelt/projiziert kennzeichnen |

## 3. Offene Forschungsentscheidungen

Vor M5 sind festzulegen:

- Nebenbedingungstechnik für Floating Electrodes: Kondensation, Lagrange-
  Multiplikator oder globale Elektroden-DOF
- skalierungsrobuste Normen für gekoppelte Residuen
- zulässige Stabilitätsprüfung für allgemeine anisotrope Piezo-Materialdaten
- Schnittstelle für räumlich veränderliche Materialorientierungen

Vor M6 sind festzulegen:

- thermodynamisches Potential und vollständige Vorzeichenkonvention der
  thermo-piezoelektrischen Transienten
- Integrationsschema und numerische Dämpfung für Elastodynamik
- Behandlung elektrischer Freiheitsgrade im dynamischen DAE-System
