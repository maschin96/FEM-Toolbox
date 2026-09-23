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
Gauge-Behandlung sind in [Doc/17](17_PIEZOELECTRIC_SOLVER.md) umgesetzt.

### ADR-012: Feldweise Kongruenzskalierung und LU für den Piezo-Referenzsolver (#21)

Der symmetrisch-indefinite Blockoperator wird mit einer positiven diagonalen
Kongruenzskalierung aus den medianen mechanischen und elektrischen
Blockdiagonalen äquilibriert und per SciPy Sparse LU gelöst. CG wird nicht
angeboten; MINRES folgt erst mit geeignetem symmetrischem Präconditioner.
Residuen und Reaktionen werden am unskalierten Originalsystem getrennt für
Kraft (N) und Ladung (C) berechnet. Starrkörpermoden und Potential-Gauge
werden pro Netzkomponente vor der Faktorisierung geprüft. Dies definiert den
verifizierbaren Referenzpfad, keine Skalierungszusage für große Modelle.

### ADR-013: Eindeutige Materialnamen auf konformen Tet4-Zellen (#137)

Benannte Volumensätze dürfen sich überlappen, wenn sie für ein aktives Feld
denselben Materialnamen zuweisen. Unterschiedliche Namen bleiben ein Fehler,
selbst bei derzeit gleichen Tensorwerten: Materialidentität und Herkunft
werden so nicht aus zufälliger Zahlenwertgleichheit abgeleitet. Die
Materialbelegung wird pro Feld vollständig vor Assembly geprüft und in
Tet4-Zeilenfolge mit Quellgruppen berichtet. Gemeinsame Netzknoten erzwingen
stetige Primärvariablen; die schwache Form liefert das Flussgleichgewicht
ohne gesonderte Grenzflächenoperatoren. Details und Grenzen in [Doc/19](19_MATERIAL_REGIONS.md).

### ADR-014: Floating Electrode durch Potentialkondensation (#24)

Alle elektrischen Knoten einer schwebenden Elektrode teilen einen reduzierten
Potentialfreiheitsgrad. Eine Projektionsmatrix `P` bildet die ursprüngliche
symmetrische Matrix auf `Pᵀ A P` ab. Wegen des negativen elektrischen Blocks
und der Leiterladungsdefinition `Q=-sum(r_phi)` lautet die reduzierte rechte
Seite auf der Elektrode `Pᵀ b-Q`. Eine feste Elektrode verwendet dieselbe
Knotenauswahl mit Dirichlet-Werten. Ein Potentialbezug bleibt pro
elektrischer Komponente nötig; keine implizite Gauge-Fixierung. Alternative
Lagrange-Multiplikatoren sind für den Referenzpfad nicht erforderlich.
Herleitung und 1D-Verifikation stehen in [Doc/20](20_ELECTRODES.md).

### ADR-015: TOML-Schema für stationäre Model-Fälle (#25)

Version 1 verwendet TOML und explizite SI-Einheiten. Python 3.11 liest TOML
ohne zusätzliche Laufzeitabhängigkeit. Unbekannte Schlüssel und andere
Schemaversionen werden abgelehnt; inkompatible Änderungen erhalten eine neue
Version mit expliziter Migration. Der Adapter ordnet benannte Materialien und
Randbedingungen zu und übergibt sie an die verifizierten Referenzsolver.
Die vier bestehenden stationären Einzelmodi sind auswählbar; die kombinierte
thermisch-piezoelektrische Orchestrierung wurde mit #26 ergänzt. Format und API stehen in
[Doc/21](21_MODEL_SCHEMA.md).

### ADR-016: Einseitiger T-Transfer in den monolithischen Piezo-Solver (#26)

Der stationäre Ablauf löst Wärme zuerst und übergibt absolute Knotenwerte an
den vorhandenen monolithischen `u-φ`-Referenzsolver. Die thermische
Eigenverzerrung wird je Zelle aus einem explizit zugeordneten
`ThermalExpansion`-Material berechnet; der elektrische thermische Beitrag
bleibt ein ausdrücklicher Schalter. Feldweise Residuen und Bilanzen werden
als getrennte Teilergebnisse erhalten. Reziproke Wärme-Rückkopplung und
Zeitabhängigkeit gehören nicht zu diesem Modell. Nachweis in
[Doc/22](22_THERMOPIEZOELECTRIC.md).

### ADR-017: Struktureller Modell-Preflight vor Assembly (#138)

Der deklarative `Model`-Adapter prüft pro Tet4-Knotenkomponente feste
Temperatur-/Potentialbezüge und den Rang mechanischer Lagerungen, bevor ein
Feldsystem assembliert wird. Identische Dirichlet-Werte sind redundant, aber
zulässig; unterschiedliche Werte bleiben ein Fehler. Eine Elektrode und ein
vorgegebener D-Fluss auf derselben Facette sind als alternative elektrische
Randgesetze unvereinbar und werden auch im direkten Solver abgelehnt.
Lastabhängige Neumann-Kompatibilität und Kondition bleiben Diagnose der
direkten Solver. Der Preflight setzt niemals eine Gauge oder Lagerung selbst.
Vertrag und Tests stehen in [Doc/23](23_MODEL_PREFLIGHT.md).

### ADR-018: Explizite Datenblatt-Einheiten im Model-Schema (#139)

Version 1 behält `units = "SI"` als Aussage über den Solverzustand. Für
Materialtabellen sind wenige benannte, explizite Datenblatt-Einheiten
zugelassen; die Umrechnung erfolgt beim Modellaufbau vor der Assembly.
Relative Permittivität verwendet den festgelegten Wert
`ε₀ = 8.8541878188e-12 F/m`, damit ein Lauf nicht durch eine spätere
Konstantenbibliothek geändert wird. Absolute Celsius-Temperaturen erhalten
den Offset `273.15`; Temperaturdifferenzen haben keinen Offset. Beide Fälle
besitzen verschiedene Einheitenzeichen. Piezo-Tensoren werden nach der
Konversion in die bestehenden globalen Stress-Charge-Tensoren überführt.
Originalwerte, Faktoren, lokale Achsen und globale SI-Tensoren bleiben in der
Laufprovenienz nachvollziehbar. Der vollständige Vertrag steht in
[Doc/21](21_MODEL_SCHEMA.md).

### ADR-019: Rohfelder und Manifest als stationäres Ergebnisbündel (#140)

Die VTU-Datei enthält unveränderte numerische Solverfelder; das zugehörige
JSON-Manifest hält Einheiten, Entitätsbezug, Bilanzen, Einschränkungen und
Provenienz fest. Der VTU-Prüfwert koppelt beide Dateien. Ungenutzte Netzknoten
und Triangle3-Zellen bekommen explizit dokumentierte Nullfüllung, statt
scheinbar interpolierte Werte. Der Export weist vorhandene Zieldateien ab und
veröffentlicht das Manifest zuletzt. Vertrag und unabhängiger Leser stehen in
[Doc/24](24_RESULT_BUNDLE.md).

### ADR-020: Branchschutz nicht als Integrations- oder Release-Gate (#58)

Am 23.09.2026 entscheidet der Repository-Eigentümer, Branchschutz aus den
Projektanforderungen zu streichen. Die GitHub-API meldete für das private
Repository fehlenden Schutz und einen tarifbedingten HTTP 403 beim Ruleset-Abruf.
Ein Tarifupgrade oder eine Änderung der Sichtbarkeit wird nicht vorausgesetzt.

Die Aktivierung von Branchschutz/Rulesets, der Nachweis aktiver Regeln und ein
administrativer Negativtest mit blockiertem Merge entfallen als Pflicht.
Erfolgreiche anwendbare CI-Prüfungen am Zielcommit, Pull Requests und die
bestehenden Reviewregeln bleiben bestehen. Integrationsverantwortliche prüfen
sie selbst; GitHub erzwingt sie nicht. Fehlende unabhängige Reviewer werden
weiterhin genannt und nicht als erteilte Freigabe dargestellt.

Diese Entscheidung ersetzt das zuvor in #58 dokumentierte administrative Gate.
Sie schließt keinen PR und gibt keine Version frei. Nachweise und Konsequenzen
stehen in [Doc/29](29_RELEASE_ARTIFACTS_AND_PROTECTION.md).

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

## 5. M6/M7-Planung nach v1.0 (2026-09-23, #167)

Der [Entwicklungsplan](30_M6_M7_DEVELOPMENT_PLAN.md) konkretisiert die nächsten
zwei Meilensteine. Planungsentscheidungen: Backward Euler und Newmark als
erste Referenzen; geradkantige Tet10/Triangle6 als erste höherwertige Familie;
Fehlerindikatoren zunächst für stationäre Tet4-Wärme. MPI, gekrümmte Geometrie
und allgemeine gekoppelte Adaptivität benötigen spätere Entscheidungen.
Dies ändert keine implementierte Gleichung oder bestehende API.

Vor Implementierung bleiben explizite ADRs erforderlich: DAE-/Elektroden- und
Phasorvertrag (#31), thermodynamisches Potential (#32), Backend-/Scalarvertrag
(#34) sowie additive Netztopologie und Auswerteorte (#35). Die alten
M5-Entscheidungsbedarfe oben sind historisch; implementierte Entscheidungen
stehen in den Material-, Elektroden- und Solver-ADRs und Dokumenten 15–24.


## ADR-021 — M6 ohne reziproke thermomechanische Rückwirkung

Datum: 23.09.2026. Status: Nutzerentscheidung.

Der Nutzer hat die Rückwirkung von Verformung auf Temperatur ausdrücklich
aus dem aktuellen Umfang genommen. #32 bleibt als zurückgestellter Backlog
offen und wird aus dem M6-Milestone entfernt. #33 verlangt deshalb keine
reziproke Thermoelastizität und keinen Thermoelastic-Damping-Nachweis mehr.
Bestehende einseitige T→u und T→(u,phi)-Pfade behalten ihre Semantik.

Transiente Wärme verwendet konstante positive rho/c-Daten und Backward Euler
als Referenz. Vorgegebene Signalereignisse werden ins Raster eingefügt;
Quellen/Flüsse/Robin-Daten werden am linken Grenzwert des rechten
Intervallendes ausgewertet, damit ein Sprung keine frühzeitige Energie
einträgt. Dirichlet-Werte sind rechtsstetig; ihr diskreter Speicherbeitrag
wird als Reservoirreaktion bilanziert. Physikalische Dirichlet-Sprünge
benötigen zeitliche Verfeinerung oder eine endliche Anstiegsrampe.
