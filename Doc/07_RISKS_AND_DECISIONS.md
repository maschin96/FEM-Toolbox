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

## ADR-022 — Piezoelektrische DAE und mechanisches Eigenproblem (#31a)

Datum: 23.09.2026. Status: implementiert im M6-Arbeitsbranch.

Elektrische Freiheitsgrade bleiben masselos und werden nach Anwendung des
bestehenden Elektroden-/Gaugevertrags kondensiert. Der positive Schur-Operator
`K+B D⁻¹ Bᵀ` und ausschließlich mechanische konsistente Masse definieren die
Moden. Homogene Störungen, explizite Lager/Gauge und ein Größenlimit von 2048
aktiven Gesamt-DOFs begrenzen die erste dichte Referenz. Anfangspotential und
-beschleunigung sind abhängige DAE-Größen. Einzelheiten und unabhängige
Referenzen stehen in [Dokument 34](34_PIEZO_DYNAMICS.md).

## ADR-023 — Komplexer Phasorpfad ohne Änderung von LinearSystem (#31c)

Datum: 23.09.2026. Status: implementiert im M6-Arbeitsbranch.

Harmonische Felder verwenden Spitzenwertphasoren `Re(x_hat*exp(+i*omega*t))`.
Ein additiver komplexwertiger Solver hält Lasten, Randwerte und Ergebnisse
komplex; der bestehende `LinearSystem`-Vertrag bleibt reell. Frequenzen werden
explizit in Hz eingegeben, rad/s separat ausgewiesen. Der ungedämpfte Pfad
weist numerische Resonanz/Nähe zur Singularität zurück statt Dämpfung zu erfinden.
Skalierung, Grenzwert und Nachweise: [Dokument 34](34_PIEZO_DYNAMICS.md).

## ADR-024 — Interner linearer Backend-/Scalarvertrag (#34a)

Datum: 25.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die Physikadapter verantworten Randbedingungen, Elektrodenkondensation,
Nullräume und Rekonstruktion. Das Backend erhält den freien, unskalierten
Operator samt vollständiger disjunkter Feldpartition in dessen lokaler
Zeilennummerierung und positiver reeller Diagonalskalierung S. Es löst
`(S A S)y = S b`, liefert `x=S y` und prüft Residuen gegen `A x-b`.
Die Feldpartition beschreibt nach elektrischer Schur-Elimination nur das
jeweils verbleibende Feld; die rekonstruierten physikalischen Residuen bleiben
zusätzlich in den bestehenden Solvern Pflicht.

Ein Faktorisierungsprotokoll erlaubt mehrere rechte Seiten und explizite
Wiederverwendung; keine globale Backend-Auswahl oder implizite Caches.
Reelle und komplexe float64/complex128-Daten sind zulässig, ohne komplexe
Anteile abzuwerfen. SciPy implementiert Sparse-LU für allgemeine Operatoren
und dichte Cholesky für die bereits dichten positiv definiten M6-Teilprobleme.
Modale Eigenwertprobleme bleiben beim separaten SciPy-eigh-Vertrag.
Backendfehler werden nicht durch Fallback verdeckt. SolveReport nennt Backend,
Verfahren, Konvergenzgrund und unskalierte feldweise Residuen. Physikalische
Bilanzen/Einheiten und die globale Reaktionsauswertung bleiben im Adapter.
PETSc, FieldSplit, iterative Nullräume und Model-Konfiguration folgen #34b/c.


## ADR-025 — Optionaler einprozessiger PETSc-Direktadapter (#34b)

Datum: 25.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die erste native Referenz verwendet PETSc PREONLY mit SeqAIJ-LU auf COMM_SELF.
Sie teilt Datenvalidierung und Originalresiduen mit SciPy, aber keine
Faktorisierung. PETSc wird erst bei expliziter Backend-Auswahl geladen.
Ein reeller Build bildet komplexe Operatoren exakt als doppelt großes reelles
Blocksystem ab; ein komplexer Build löst nativ. Beide Varianten erhalten echte
CI-Nachweise. Die kontrollierte optionale Umgebung ist Ubuntu 24.04 mit
Python 3.12 und PETSc/petsc4py 3.19.6; der Core behält Python 3.11–3.13.

Diese Wahl liefert zunächst Scalar- und Physikparität ohne iterative Parameter
oder stillen Fallback. LU-Pivotierung ist nicht mit SuperLU identisch; Fehler
werden weitergereicht. Positive Definitheit der begrenzten dichten M6-Blöcke
wird vor der nativen LU separat geprüft. FieldSplit, Schur-Präconditionierung,
Nullraumübergabe und iterative Nichtkonvergenz folgen #34c. Keine MPI- oder
Performancezusage aus diesem Schritt. Details: [Dokument 38](38_PETSC_BACKEND.md).


## ADR-026 — FGMRES/Schur und unveränderte physikalische Nullraumprüfung (#34c)

Datum: 25.09.2026. Status: Umsetzung im M7-Issue-Branch.

Der optionale iterative Adapter verwendet Rechtspräkonditionierung mit PETSc
FieldSplit: volle Blockfaktorisierung, A11 als Schur-Näherung und native LU für
beide Teilinversen. Ein einzelnes verbleibendes Feld nutzt FGMRES/LU. Da dieser
Präkonditionierer nicht symmetrisch-positiv sein muss, ist FGMRES auch beim
symmetrisch-indefiniten Piezo-Operator angemessen; CG/MINRES werden nicht gewählt.
Diese erste Referenz verspricht keinen Laufzeit- oder verteilten Skalierungsvorteil.

Physikalische Nullräume werden weiterhin vor dem Backend diagnostiziert und
müssen durch explizite Lager/Referenzen entfernt werden. Die frühere Planung
einer Nullraumübergabe aus ADR-025 wird damit präzisiert: keine implizite
Lastprojektion, kein zusätzlicher singulärer Solve-Vertrag. Algebraische Krylov-
Konvergenz allein garantiert keine Eindeutigkeit beliebiger Eingabematrizen.

Echte PETSc-Tests prüfen native Iterationsabbrüche, Originalresiduen, Scalar-
Abbildung und Physikparität. Jeder Solver isoliert intern gelesene FieldSplit-
Optionen mit einem privaten Präfix. Ein Bericht zählt äußere Iterationen je
Solve-Aufruf; Grenzen gelten je nativer RHS. Details und Einschränkungen stehen
in [Dokument 39](39_PETSC_FIELDSPLIT.md).


## ADR-027: Additive quadratische Topologie vor Operatorfreigabe (#35a)

Datum: 25.09.2026. Status: Umsetzung im M7-Issue-Branch.

Tet10/Triangle6 verwendet VTK-/meshio-Knotenordnung und eine einzige Familie je
Netz. Bestehende Tet4-Arrays bleiben unverändert; generische Zellzugriffe dienen
zunächst nur Topologie und Metadaten. Gemeinsame Kanten werden über Eck- und
Kantenknoten geprüft. Import und Randgeometrie verlangen mittige Kantenknoten,
bevor affine Eckgeometrie verwendet wird. Model lehnt Tet10 bis #35d explizit ab.
VTU-Präsenzmasken erhalten den Unterschied zwischen fehlenden und nullwertigen
Metadaten. Mathematische Feldkonventionen bleiben unverändert.
Details: [Tet10-Topologie](40_TET10_TOPOLOGY.md).


## ADR-028: Affine Tet10-Interpolation und positive Referenzquadratur (#35b)

Datum: 25.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die erste Tet10-Geometrie beschränkt sich auf mittige Kantenknoten. Quadratische
Felder verwenden affine Eckgeometrie; gekrümmte Elemente werden nicht still
linearisiert. Triangle6 liefert tangentiale Ableitungen und lehnt Auswertung
außerhalb seiner Ebene ab. Positive kollabierte Gauss-Produkte integrieren
Volumenpolynome bis Grad 5 und Flächenpolynome bis Grad 4; Gewichte enthalten
Referenz- bzw. physische Maße und sind nicht auf eins normiert. Die vergleichsweise
hohe Punktzahl priorisiert eine einfach unabhängig prüfbare Referenz. Keine
neuen Abhängigkeiten, keine Änderung mathematischer Feldvorzeichen oder der
Tet4-API. Nachweise und Grenzen: [Dokument 41](41_TET10_GEOMETRY.md).


## ADR-029: Tet10-Operatoren und explizite positive Speicherdiagonalisierung (#35c)

Datum: 25.09.2026. Status: Umsetzung im M7-Issue-Branch.

Quadratische Operatoren verwenden einen additiven Moduleinstieg und punktabhängige
Gradienten bei unveränderter affiner Geometrie und Feldkonvention. Globale
Solverdispatch und ortsabhängige Ergebnisse folgen getrennt in #35d. Die
konsistente Tet10-Masse besitzt negative Eckzeilensummen; P1-Zeilensummen-Lumping
wird deshalb nicht übernommen. `lumping="hrz"` wählt ausdrücklich positive,
gesamtspeichererhaltende Diagonalskalierung. Standard bleibt `consistent`.
Der boolesche Tet4-Vertrag bleibt unverändert. Sechsknoten-Elektroden verwenden
die gemeinsame aktive Randtopologie. [Formeln, Quellen und Nachweise](42_TET10_OPERATORS.md).


## ADR-030: Residualindikator für stationäre Tet4-Wärme (#36)

Datum: 26.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die additive Forschungs-API in `post` verwendet Zellquellen, innere
Normalflusssprünge und explizite Neumann-/Robin-Randresiduen. Gewichtung mit
Durchmessern und kleinsten Leitfähigkeitseigenwerten liefert Größen in der
thermischen Energienorm, ohne garantierte Fehlerschranken oder Kontrastrobustheit
zu behaupten. Vollständige Triangle3-Randabdeckung und explizite Dirichlet-Flächen
verhindern das stillschweigende Weglassen von Randfehlern. Materialgrenzen
verwenden Fluss- statt Gradientensprünge. Dörfler-Markierung ist deterministisch;
Remeshing, Tet10 und gekoppelte Indikatoren bleiben getrennte Folgearbeit.
Keine Änderung der stationären Solver-, Model- oder Vorzeichenverträge und keine
neue Abhängigkeit. [Formeln und Nachweise](43_HEAT_INDICATOR.md).

## ADR-031: Tet10-Wärme mit expliziten Schwerpunktfeldern (#35d)

Datum: 26.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die thermische Integration verwendet bestehende Zeitintegratoren und Backends
mit internem Elementdispatch. `Model` gibt Tet10 zunächst ausschließlich für
stationäre und transiente Wärme frei. Quadratische Gradienten sind
ortsabhängig: Ergebnisarrays enthalten rohe Schwerpunktwerte und zusätzliche
unveränderliche `field_points`. Tet10-Exporte verwenden deshalb Schema 2 und
Feldnamen mit `_at_centroid`; Tet4 behält Schema 1. Das ersetzt keine vollständige
Feldrekonstruktion. Sondenlokalisierung verwendet affine Eckkoordinaten,
Interpolation alle zehn Formfunktionen. HRZ bleibt explizit und erhält die
Gesamtkapazität, ohne konsistente Quellmomente zu verändern. Keine neue
Abhängigkeit und keine geänderten Feldvorzeichen.
[Vertrag, Verifikation und verbleibende Schritte](45_TET10_HEAT_INTEGRATION.md).

## ADR-032: Tet10-Elektrostatik trennt Feldproben und Energieintegrale (#35d)

Datum: 26.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die stationäre Elektrostatik verwendet denselben geradkantigen Elementdispatch
und die bestehende Backend-Grenze. Die negative elektrische Gleichung und
Leiterladungsvorzeichen bleiben unverändert. E/D werden an expliziten
Schwerpunkten ausgewertet; gespeicherte Energie wird separat über die
punktabhängigen Gradienten integriert. Schema 2 unterscheidet deshalb
`raw`/`cell_centroid` von `cell_integral`/`cell_volume`; der Energiefeldname
erhält keinen Schwerpunktsuffix. Tet4 behält Schema 1. Keine neuen Abhängigkeiten
und keine Freigabe weiterer Tet10-Kopplungen.
[Nachweise und Grenzen](46_TET10_ELECTROSTATICS.md).

## ADR-033: Statische Tet10-Mechanik vor Kopplungen integrieren (#35d)

Datum: 26.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die statische Elastizität verwendet die bestehenden quadratischen Operatoren,
Lagerprüfung und Backend-Grenze. Spannungen/Dehnungen sind rohe Schwerpunktwerte
mit expliziten Auswerteorten und unveränderter technischer Voigt-Konvention.
Verzerrungsenergie wird separat über das Element integriert. Schema 2 trennt
Feldproben von Energieintegralen wie bei Elektrostatik. Tet4 behält Schema 1.
Thermische Tet10-Eigendehnungen werden bis zur eigenen Kopplungsintegration
ausdrücklich abgelehnt; die Tet4-Zellmittelung wird nicht auf P2 übertragen.
Keine neuen Abhängigkeiten oder geänderten Vorzeichen.
[Verifikation und Grenzen](47_TET10_ELASTICITY.md).

## ADR-034: Quadratische Temperaturfelder in Tet10-Thermomechanik integrieren (#35d)

Datum: 27.09.2026. Status: Umsetzung im M7-Issue-Branch.

Die stationäre einseitige Kopplung verwendet alle zehn Knotentemperaturen
in `N(x)(T-Tref)`. Die vorhandene Grad-5-Quadratur integriert thermische Lasten
und mechanische Energie; Schwerpunkt-Eigendehnung und -Spannung bleiben
ausdrücklich Feldproben. Tet4-Mittelung und Schema 1 bleiben unverändert.
Der bestehende Python-Ablauf erhält getrennte optionale Backend-Hooks;
der gemeinsame Schema-2-Export prüft die Auswerteorte beider Teilresultate.
Kein neuer Model-Modus, keine thermische Rückwirkung oder geänderten Vorzeichen.
[Nachweise und Grenzen](48_TET10_THERMOMECHANICS.md).

## ADR-035 — Stationäre Tet10-Piezoelektrizität (27.09.2026)

Der vorhandene monolithische Solver nutzt die generische aktive Zellfamilie
und die verifizierten Tet10-Operatoren. Elektrodenreduktion erfasst alle
Triangle6-Knoten. Der interne Elementadapter trennt Schwerpunktfelder von
integrierten Enthalpiebeiträgen; das additive `field_points` hält Tet4-Aufrufer
kompatibel. Schema 2 kennzeichnet die drei Energiebeiträge als Zellintegrale.
Thermopiezo und Dynamik bleiben bis zu eigenen Integrationsnachweisen gesperrt.
Begründung und analytische/Konvergenznachweise: [Doc 49](49_TET10_PIEZOELECTRIC.md).

## ADR-036 — Stationäre Tet10-Thermopiezoelektrizität (28.09.2026)

Der vorhandene Ablauf T → (u, φ) wird mit vollständigem Knotentemperaturtransfer
auf geradkantige Tet10 erweitert. Eigenverzerrung und Enthalpie werden an
Integrationspunkten ausgewertet; Feldberichte bleiben Schwerpunktwerte.
Der elektrische Thermalschalter behält seine beiden bestehenden Konventionen.
Tet4-Mittelung und -Provenienz bleiben erhalten; Tet10 hasht zusätzlich alle
Temperaturen und Expansionsdaten. Getrennte optionale Backends und Schema-2-Export
folgen dem Thermomechanikvertrag. [Details und Nachweise](50_TET10_THERMOPIEZOELECTRIC.md).

## ADR-037 — Tet10-Modalanalyse (28.09.2026)

Der gemeinsame Dynamikaufbau verwendet Tet10-Operatoren und konsistente Masse;
nur der Modalpfad wird freigegeben. Die elektrische Schur-Elimination erhält
einen explizit symmetrisierten Tet10-Dielektrikblock gegen Quadratur-Rundungsreste.
Residuen bleiben am ursprünglichen System definiert. Die Grenze von 2048 DOFs
bleibt bestehen. Modenexport nutzt Schema 2 mit Schwerpunktkennzeichnung und
P2-Messpunkten. [Analytische Nachweise und Grenzen](51_TET10_MODES.md).

## ADR-038 — Tet10-Newmark (28.09.2026)

Die bestehende ungedämpfte Newmark-Integration wird für geradkantige Tet10
freigegeben. Der gemeinsame Dynamikaufbau hält den bereits für Moden
symmetrisierten dielektrischen Block vor; Newmark verwendet ihn für die
elektrische Lastrekonstruktion wieder. Physikalische Matrizen, Residuen,
Energie-, Ereignis- und Speicherverträge bleiben unverändert. Die konsistente
P2-Masse wird nicht gelumpt. [Nachweise und Grenzen](52_TET10_NEWMARK.md).

## ADR-039 — Tet10-Harmonik (28.09.2026)

Der bestehende ungedämpfte harmonische Referenzsolver wird für geradkantige
Tet10 freigegeben. Die elektrische Lastrekonstruktion verwendet denselben
symmetrisierten Dielektrikblock wie Moden und Newmark. Originalresiduen,
komplexe Spitzenwertkonvention, Resonanzschwelle und DOF-Grenze bleiben erhalten.
Model-Schema 3 und der vorhandene Schema-2-Feldexport benötigen keine neue
Schemaversion. [Unabhängige Nachweise](53_TET10_HARMONIC.md).


## ADR-040 — Lokale stationäre Parameterstudien (29.09.2026)

#59 ergänzt einen Adapter über der vorhandenen Model-/Ergebnis-API, ohne
Feldoperatoren zu ändern. Schema 1 bildet ein endliches kartesisches Produkt
expliziter SI-Parameter auf bestehende Model-Eigenschaften ab. Fallidentität
und Eingabegültigkeit sind getrennt: Parameter definieren stabile Fall-IDs;
Dokument, Netzbytes, Solveroptionen und Versionsvertrag bestimmen Wiederverwendung.
Ein atomar ersetzter Katalog veröffentlicht den Erfolg erst nach vollständigem
Ergebnisbündel. Fortsetzung prüft beide Dateihashes und Bündelprovenienz.
Jeder Versuch erhält ein neues Verzeichnis, ein exklusiver Lock verhindert
konkurrierende Writer. Nach hartem Prozessabbruch erfolgt Lockfreigabe bewusst
manuell, ohne fremde laufende Arbeit automatisch zu übernehmen. Alte Versuche
bleiben erhalten; keine verteilte Ausführung oder automatische Bereinigung.
[API, Katalogvertrag und Nachweise](59_PARAMETER_STUDIES.md).

## ADR-041 — Studienkennwerte als CSV-/JSON-Vertrag (29.09.2026)

#60 trennt numerische CSV-Werte von vollständigen JSON-Metadaten und einem
lesbaren Markdown-Bericht. Fallidentität, Herkunft, Einheit, Auswahl und Reduktion
bleiben explizit; fehlende Werte sind null mit Grund. Bilanzbewertungen benötigen
explizite Schwellen, Netzkonvergenz wird nicht aus einem Solve abgeleitet.
[Vertrag und Nachweise](60_STUDY_REPORTS.md).

## ADR-042 — Polynomiale Surrogate mit explizitem Holdout (29.09.2026)

#205 nutzt NumPy-Least-Squares auf einer train-only skalierten vollständigen
Polynombasis vom Gesamtgrad 1/2. Optionale Ridge-Regularisierung wird als
erweitertes Least-Squares-System ohne Strafe auf den konstanten Term gelöst.
Explizite Fallpartition und Gruppierung identischer Eingaben verhindern Datenleckage.
Der #60-Adapter übernimmt SI-/Quellmetadaten; Ausschlüsse sind ausdrücklich und
vollständig protokolliert. Ein geprüftes JSON-Artefakt enthält alle für Vorhersagen
benötigten Größen, niemals ausführbare Serialisierung. Boxgrenzen und Testfehler
werden nicht als Unsicherheitsintervalle oder physikalische Validierung ausgegeben.
[Vertrag und Nachweise](61_SURROGATE_MODELS.md).

## ADR-043: optionale Matplotlib-Feldbilder (#189)

Entscheidung vom 30.09.2026: Matplotlib >=3.9,<4 als separates
`visualization`-Extra (BSD-kompatible Lizenz), expliziter Agg-Canvas ohne
pyplot-Registrierung. Tet4-Außenflächen mit nodalen Flächenmitteln oder
Eigentümerzellwerten; keine implizite Materialglättung. Tet10 und komplexe
Felder werden abgelehnt. [API, Grenzen und Nachweise](62_FIELD_VISUALIZATION.md).

## ADR-044 — Experimenteller eigenständiger DC-Schaltungskern (01.10.2026)

Auf ausdrücklichen Nutzerauftrag beginnt M13 auf einem separaten
Experimental-Branch vom stabilen M13-Integrationsstand. Erster Schnitt zu
#114/#115/#116: lineare R/C/L-Elemente und unabhängige V-/I-Quellen,
begrenzter strikt geprüfter Netzlistendialekt und Sparse-MNA mit
Zeilen-/Spaltenäquilibrierung. Passive Zweigvorzeichen, Masse `0`, C offen und
L kurzgeschlossen definieren den DC-Grenzfall. Kein CG, keine künstlichen
Leckleitfähigkeiten und keine automatische Auflösung redundanter idealer
Quellenschleifen. KCL an allen Knoten und aufgenommene Zweigleistung werden
aus der Lösung rekonstruiert und getrennt vom algebraischen Residuum geprüft.

Der Kern verwendet vorhandenes NumPy/SciPy und `RunProvenance`. Module bleiben
getrennt vom FEM-Model; die experimentelle API erhält keine v1.3-Stabilitätszusage.
Diese Entscheidung legt nur den eigenständigen DC-Schnitt fest. Vollständige
elektrische/thermische Ports und monolithische/iterative FEM-Kopplung aus #114,
Transienten und AC werden vor ihrer Implementierung entschieden. M13 ist
damit weder abgeschlossen noch freigegeben; Issues bleiben offen.
[Vertrag, Referenzen und Prüfnachweis](65_EXPERIMENTAL_CIRCUIT_CORE.md).
