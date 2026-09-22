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

Entscheidung: Der stationäre Ausbau M3–M5 löst Wärme unabhängig und übergibt `T`
an Mechanik und optional das explizite thermo-piezoelektrische Materialmodell.

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

### ADR-007: Referenzumfang der linearen Elastizität (#12)

Isotrope Tet4-Materialien mit `nu >= 0.49` werden abgelehnt; anisotrope
Steifigkeiten werden in Globalachsen geliefert und auf Symmetrie und positive
Definitheit geprüft. Die Grenze operationalisiert ADR-001 ohne eine gemischte
Inkompressibilitätsformulierung einzuführen. Lasten sind globale Vektoren;
Druck-/Normalenadapter und Materialachsenrotation bleiben spätere Erweiterungen.
Der Referenzsolver verwendet diagonale Gleichgewichtung, Starrkörperprüfung und
LU-Pivotkontrolle. Residuen werden über `|| |K| |u| + |f| ||` normiert, damit
vorgegebene starre Bewegungen keinen künstlichen relativen Fehler erzeugen.
Einheiten, Grenzen und Nachweise stehen in [Lineare Elastizität](10_LINEAR_ELASTICITY.md).


### ADR-008: Äußerer Wärmestrom und thermischer Nullraum (#13)

Der vorgegebene Wärmestrom `q_bar = -k grad T dot n` ist positiv nach außen.
Damit lautet die schwache rechte Seite `integral w Q - integral w q_bar`.
Das frühere Pluszeichen in Abschnitt 4 des mathematischen Modells widersprach
seiner starken Form und Randdefinition und wird vor der Wärmeimplementierung
korrigiert. Konvektion `q dot n = h (T - T_inf)` liefert eine positive Randmatrix
und einen positiven Umgebungstemperatur-Lastterm. Bilanztests sichern diese
Konvention unabhängig durch analytische Wärmeleistungen ab.

Jede unverbundene thermische Komponente benötigt eine Dirichlet-Temperatur oder
Konvektion mit positivem h. Reine Neumann-Probleme werden vor Faktorisierung
mit ihrer Komponenten-Leistungsbilanz als kompatibel (additive Konstante frei)
oder inkompatibel diagnostiziert. Es wird keine automatische Referenztemperatur
erfunden. Anisotrope Leitfähigkeiten liegen in globalen Achsen vor.

### ADR-009: Elektrische Reaktionen und Einzelfeldsolver (#14)

Die bestehende negative elektrische Zeile bleibt erhalten: `-K_phi phi =
surface(D.n) - volume(rho)`. Der positive Baustein K_phi wird separat angeboten;
der freie Löseblock wird nur gemeinsam mit seiner rechten Seite vorzeicheninvertiert.
Residual und Reaktionen werden am negativen Originalsystem berechnet. Eine
Dirichlet-Reaktion beschreibt den äußeren D-Fluss des Dielektrikums; die freie
Ladung der angrenzenden ideal leitenden Elektrode hat das entgegengesetzte
Vorzeichen. Positive gespeicherte Feldenergie ist von elektrischer Enthalpie
zu unterscheiden. Die Gauge wird pro unverbundener Komponente geprüft und
niemals automatisch gewählt. Details und unabhängige Nachweise stehen in
[Elektrostatik](12_ELECTROSTATICS.md).

### ADR-010: Piezoelektrische Materialkonversion und Rotation (#19)

Intern werden `C^E`, `e` und `kappa^S` in globalen Achsen und der festgelegten
Engineering-Voigt-Reihenfolge gespeichert. Strain-Charge-Eingaben werden mit
`C^E=(s^E)^-1`, `e=d C^E` und `kappa^S=kappa^T-d C^E d^T` konvertiert.
Positive Definitheit von `s^E`, `C^E` und `kappa^S` wird geprüft; die
Permittivitätsdifferenz darf nicht stillschweigend begrenzt werden.
Materialachsen werden als Spalten einer eigentlichen orthonormalen Rotation
verstanden. Unterschiedliche Spannungs- und technische
Verzerrungstransformationen verhindern einen Schubfaktorfehler. Lokale
Eingabetensoren und globale Ergebnistensoren bleiben für die spätere
Solverprovenienz verfügbar. Details und algebraische Referenzen stehen in
[Piezo-Materialdaten](15_PIEZOELECTRIC_MATERIAL.md).

### ADR-011: Gekoppelte Elemententhalpie und thermischer D-Term (#20)

Das lokale Piezoelement verwendet die symmetrisch-indefinite Matrix mit
`K_phiu=K_uphi^T` und negativem dielektrischen Block. Mechanische thermische
Eigendehnung ist stets explizit; der Beitrag `-e epsilon_th` in `D` ist
standardmäßig **aus** und kann nur mit
`electric_thermal_strain=True` aktiviert werden. Die zugehörige elektrische
rechte Seite wird zusammen mit der mechanischen aus der gleichen
Enthalpieformulierung abgeleitet. Die lokale Energie- und Patch-Verifikation
steht in [Doc/16](16_PIEZOELECTRIC_ELEMENT.md); die globale Assembly und
Gauge-Behandlung bleiben Aufgabe von #21.

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

## 4. Planungswelle 2 (2026-09-21, #56)

Die [Produktvision](13_PRODUCT_VISION.md) ergänzt M8–M10. Die alten
Versionsverweise „0.1“ auf den stationären Zielumfang werden zu M3–M5/bis v1.0
präzisiert; Gleichungen und Vorzeichen bleiben unverändert. Der ursprüngliche
M0-Release v0.1.0 bezeichnet ausschließlich das Projektfundament.

Neue Entscheidungsbedarfe: API-/Kompatibilitätsvertrag (#57), Messdaten und
Validierungsgrenzen vor Auswertung (#63), konforme adaptive Wärme mit
Remeshing-Go/no-go (#67). Diese Planung nimmt die späteren Architektur- und
Methodenentscheidungen nicht vorweg.
