# Produktvision und Planungswelle 2

Status: Planungsbaseline, keine Implementierungszusage für künftige Funktionen.
Bestandsaufnahme: 2026-09-21; Abschlussprüfung: 2026-09-22. Tracking: [#56](https://github.com/maschin96/FEM-Toolbox/issues/56).

> Historische Bestandsaufnahme: Aussagen zu offenen M3–M5-Arbeiten beziehen
> sich auf September 2026 vor v1.0. Den aktuellen Stand und die verbindliche
> Konkretisierung von M6/M7 enthält [Dokument 30](30_M6_M7_DEVELOPMENT_PLAN.md).

## 1. Evidenz und aktueller Stand

Analysiert wurde `origin/main` bei `e3382a1a84a11221c34e18c64271249784ee9b87`
(v0.3.0). Der ursprünglich ausgecheckte Branch `codex/46-readme` enthielt nur
v0.2.0; dessen README wäre als alleinige Bestandsaufnahme veraltet gewesen.
Code, öffentliche APIs, Tests, Dokumentation, Releases und alle 40 damaligen
GitHub-Issues wurden gegengeprüft. GitHub meldete M0–M2 geschlossen, M3–M7 offen,
keine offenen PRs und erfolgreiche CI auf dem genannten Hauptbranch-Commit.
Die Zählung bezieht sich auf den Stand **vor** dieser Planungswelle.

| Bereich | Nachweis im Repository | Einordnung |
|---|---|---|
| M0: Paket und Qualitätsbasis | `pyproject.toml`, `.github/workflows/ci.yml`, `conventions.py`, `errors.py`, `provenance.py` | Freigegeben als v0.1.0; Python 3.11–3.13 in Linux-CI |
| M1: Netz, Geometrie und I/O | `mesh/model.py`, `elements/tet4.py`, `io/meshio_adapter.py`, Geometrie-/Rundreisetests | Freigegeben als v0.2.0; Tet4/Triangle3, Gruppen, Qualität, Gmsh/VTU |
| M2: gemeinsame numerische Basis | `fields/`, `assembly/`, `boundary/` | Feldweise aktive DOFs, COO→CSR/CSC, symmetrieerhaltende Dirichlet-Elimination |
| M2: Elastizität | `solvers/elasticity.py`, `materials/`, `physics/elasticity.py` | Isotrop/anisotrop, Volumen-/Flächen-/Punktlasten, Spannung, Energie, Reaktionen |
| M2: Wärme | `solvers/heat.py`, Material-/Element-/Lastmodule | Stationär, anisotrop, Quellen, äußerer Wärmestrom, Konvektion und Bilanz |
| M2: Elektrostatik | `solvers/electrostatics.py`, Material-/Element-/Lastmodule | Potential, E/D, Ladungen, Gauge-Diagnostik und Energie |
| Numerische Nachweise | `tests/unit/`, `tests/patch/`, drei Dateien unter `tests/convergence/` | Analytische Fälle, physikalische Bilanzen, vier Netzstufen je Einzelfeld |
| Nutzerbeispiele | `examples/elastic_tension.py`, `steady_heat.py`, `parallel_plate.py` | Drei ausführbare Einzelfälle mit VTU-Rohfeldern |

Die Quellpfade in der Tabelle beziehen sich auf `src/fem_toolbox/`, sofern kein
anderer Präfix angegeben ist. Freigegebener Umfang und Einschränkungen stehen in
[v0.3.0](releases/v0.3.0.md) und den API-Dokumenten [10](10_LINEAR_ELASTICITY.md),
[11](11_STEADY_HEAT.md), [12](12_ELECTROSTATICS.md).

**Noch nicht implementiert:** thermische Expansion, Piezo-Kopplung,
Materialachsenrotation als allgemeine API, Floating Electrodes, deklarative
`Model`-API, Transienten, PETSc, Tet10 und adaptiver Regelkreis. Anisotrope
Matrizen werden derzeit in globalen Achsen geliefert. Der reine Verschiebungs-
Tet4 lehnt isotrope Materialien mit `nu >= 0.49` ab. Rohspannungen und Flüsse
sind elementweise konstant. Ein bestandener Code-Test ist keine experimentelle
Validierung und keine Aussage zur industriellen Einsatzreife.

## 2. Was in den vorhandenen Milestones offen ist

Die folgenden bestehenden Issues bleiben die verbindlichen Arbeitspakete;
sie werden durch die neue Welle nicht dupliziert. Zahlen schließen Release-Issues ein.

| Milestone | Offene Issues vor Welle 2 | Verbleibendes Ergebnis und Gate |
|---|---:|---|
| M3 / v0.4.0 | 3: #16–#18 | Anisotrope Eigendehnung, regionsweise Referenztemperatur, T→u; freie/behinderte Expansion, Gradientenfall und Netzkonvergenz |
| M4 / v0.5.0 | 5: #19–#23 | Materialkonversion/-rotation, gekoppelte Blöcke, skalierter u–phi-Solver; Reziprozität, Energie, direkter/inverser Effekt und zweite unabhängige Referenz |
| M5 / v1.0.0 | 6: #24–#28, #46 | Floating-/Equipotential-Elektroden, Model-Schema, T→(u,phi), Provenienz und Installation; reproduzierbarer vollständiger stationärer Workflow |
| M6 / v1.1.0 | 5: #29–#33 | Masse/Kapazität, Wärmezeitintegration, Piezoelastodynamik, harmonische Analyse, reziproke Thermoelastizität; Zeit-/Energie-/Frequenznachweise |
| M7 / v1.2.0 | 5: #34–#38 | PETSc, Tet10, Fehlerindikatoren und profilierte Skalierung; Backend-Parität, höhere Konvergenz und reproduzierbare Messung |

M6 und M7 sind breit geschnitten. Vor Beginn werden #31 (Zeit-, Eigenwert- und
Frequenzanalyse) sowie #35 (Geometrie, Quadratur, I/O, Physiken) bei Bedarf in
kleinere abhängige PRs zerlegt. Die Planung erzeugt dafür jetzt keine parallelen
Duplikat-Issues. Keine belastbare Kalenderplanung ist ohne Teamkapazität und
Aufwandsabschätzung möglich; die Reihenfolge orientiert sich an Abnahmegates.

## 3. Produktvision und Zielgruppen

**FEM Toolbox soll eine nachvollziehbare Forschungsumgebung für lineare
thermische, mechanische und piezoelektrische 3D-Fragen werden: vom Modell über
den verifizierten Solver bis zur reproduzierbaren Studie mit dokumentierter
Unsicherheit und klaren Gültigkeitsgrenzen.**

Die Vision ist aus dem vorhandenen Scope abgeleitet, noch keine durch
Nutzerinterviews bestätigte Marktpositionierung. Drei primäre Arbeitsabläufe:

1. **Forschung:** eine neue Formulierung gegen reine Element-/Assemblyfunktionen
   und veröffentlichte Referenzfälle prüfen.
2. **Lehre und Einstieg:** einen kleinen Fall aus sauberer Installation lösen,
   Randbedingungen und Bilanz verstehen und Rohfelder in ParaView betrachten.
3. **Piezo-Anwendungsforschung:** Aktor-/Sensorvarianten vergleichen und für einen
   abgegrenzten linearen Fall Prognose und Messung samt Unsicherheit bewerten.

Erfolg wird pro Release an reproduzierten Referenzfällen, vollständiger
Provenienz, erfüllten Bilanz-/Fehlergrenzen und unabhängig durchgeführten
Nutzerrezepten gemessen. Es gibt noch keine Nutzungs- oder Zeitersparnisdaten;
solche Zahlen werden nicht als Produktnachweis erfunden. Laufzeit- und
Speicherziele entstehen erst aus gemessenen repräsentativen Fällen in #37.

GUI, CAD-Kern, Cloud-Betrieb, Kontakt, Plastizität, Hysterese und allgemeine
industrielle Zertifizierung bleiben außerhalb dieser Roadmap. Topologieoptimierung
und RANS-CFD/FSI ergänzen sie als grobe Horizonte M11–M12; sie gehören nicht
zum zugesagten Umfang von M0–M10. Ein neuer Physikbereich benötigt einen eigenen Anwendungsfall,
mathematischen Entwurf und unabhängigen V&V-Plan.

## 4. Lücken vor der stabilen Version

| Befund | Konsequenz und zuständiges Issue |
|---|---|
| Direkte Solver sind nutzbar; stabile öffentliche Oberfläche und Migrationsregeln fehlen | **Neu #57:** API-/Kompatibilitätsvertrag als Input für #25; kein zweites Model-Issue |
| `main` meldet `protected=false`; CI prüft Quellinstallation unter Linux | **Neu #58:** Artefaktinstallation, Supportgrenzen und nachgewiesener Schutz/Review vor v1.0; GitHub-Einstellungen wurden hier nicht geändert |
| Gmsh-Einrichtung und vollständiges Onboarding bleiben offen | **Bestehend #46** beibehalten; Anleitung früh bearbeiten, vollständige v1.0-Abnahme mit #27 |
| Unabhängige Piezo-Referenz und Materialstabilitätsentscheidung sind Voraussetzung | **Bestehend #19/#22:** Referenzbeschaffung früh starten; keine zweite Benchmark-Aufgabe |
| M4 nennt einfache Open-Circuit-Grenzen, allgemeine Floating Electrode folgt erst M5 | #22 auf analytisch geeignete einfache Fälle begrenzen; allgemeine leitende Open-Circuit-Bedingung erst mit #24 freigeben |
| Experimentelle Daten sind bisher nicht ausgewählt | **Neu #63:** Daten-/Messvertrag mit Go/no-go vor Validierungsarbeit |

Zudem waren die ursprünglichen stationären Zielkriterien in Dokument 01 und
mehrere Modelltexte noch als „Version 0.1“ bezeichnet. Sie werden dem
stationären Ausbau bis v1.0 zugeordnet; v0.1.0 bleibt das bereits freigegebene
Projektfundament. Dies verändert keine Gleichungen oder Vorzeichen.

## 5. Neue Milestones und Arbeitspakete

### M8 — Reproduzierbare Parameterstudien (vorläufig v1.3.0)

[GitHub M8](https://github.com/maschin96/FEM-Toolbox/milestone/9).
Fachliche Voraussetzung: M5. Integration ist nach M7 vorgesehen; Entwicklung
kann nach stabiler Model-/Ergebnis-API unabhängig von PETSc beginnen.

- [#59](https://github.com/maschin96/FEM-Toolbox/issues/59): deterministisches Studienmanifest, Fall-IDs, Status und Fortsetzung.
- [#60](https://github.com/maschin96/FEM-Toolbox/issues/60): Kennwerttabellen und Berichte mit Einheiten, Bilanzstatus und Provenienz.
- [#61](https://github.com/maschin96/FEM-Toolbox/issues/61): thermisches und piezoelektrisches Studienrezept samt unabhängiger Reproduktion.
- [#62](https://github.com/maschin96/FEM-Toolbox/issues/62): Release-Gate einschließlich Abbruch-/Resume-Nachweis und Kompatibilität.

Abnahme: getrennte Einzelläufe und Studie stimmen innerhalb vorab begründeter
Toleranzen überein; Fehler werden sichtbar, veränderte Eingaben invalidieren
Ergebnisse und beide Rezepte laufen aus sauberer Installation. Kein Optimierer,
keine verteilte Ausführung und kein Dashboard als Pflichtumfang.

### M9 — Experimentelle Validierung und Unsicherheit (vorläufig v1.4.0)

[GitHub M9](https://github.com/maschin96/FEM-Toolbox/milestone/10).
Voraussetzungen: M5 und M8; geeignete Messdaten mit Nutzungsrechten. Ein linearer
quasistatischer Piezo-Fall begrenzt den Anspruch; M6-Dynamik ist nicht erforderlich.

- [#63](https://github.com/maschin96/FEM-Toolbox/issues/63): Fall, Messdatenvertrag, Unsicherheiten und vorab festgelegte Akzeptanzgrenzen.
- [#64](https://github.com/maschin96/FEM-Toolbox/issues/64): reproduzierbare Sensitivitäts-/Unsicherheitsstudien auf M8.
- [#65](https://github.com/maschin96/FEM-Toolbox/issues/65): unabhängiger Messvergleich mit Netzfehler- und Unsicherheitsbudget.
- [#66](https://github.com/maschin96/FEM-Toolbox/issues/66): Release-Gate für den abgegrenzten Validierungsworkflow.

Abnahme: vorab bestimmte Messgrößen erfüllen die festgelegten Grenzen;
Kalibrierung und Validierungsdaten bleiben getrennt. Fehlende Daten oder
ungeklärte Abweichungen halten das Gate offen. Ein negativer Befund wird
veröffentlicht, aber nicht durch nachträglich verbreiterte Toleranzen zum Erfolg
erklärt. Keine pauschale Validierung der Toolbox.

### M10 — Kontrollierte adaptive Netzstudien (vorläufig v1.5.0)

[GitHub M10](https://github.com/maschin96/FEM-Toolbox/milestone/11).
Fachliche Voraussetzungen: M7-Fehlerindikatoren (#36), M8-Provenienz und positives
Machbarkeitsgate. M9 ist in der Releasefolge vorgelagert, keine mathematische
Voraussetzung; bei Verzögerung wird die Versionsfolge per Issue neu geplant.

- [#67](https://github.com/maschin96/FEM-Toolbox/issues/67): ADR und Machbarkeit konformen Remeshings mit Erhalt von Rand-/Materialgruppen.
- [#68](https://github.com/maschin96/FEM-Toolbox/issues/68): Solve–Estimate–Mark–Remesh, Budgets, Stagnation und Fehlerabbruch.
- [#69](https://github.com/maschin96/FEM-Toolbox/issues/69): unabhängiger Fehlervergleich gegen uniforme Netzfolgen.
- [#70](https://github.com/maschin96/FEM-Toolbox/issues/70): Release-Gate für kontrollierte Fehlerreduktion.

Erster Umfang ist ausschließlich stationäre Wärme auf konformen Tet4-Netzen.
Abnahme: vier Stufen, nachgewiesene Fehlerreduktion und mindestens ein geeigneter
lokalisierter Fall mit geringerem DOF-Bedarf für dieselbe Fehlergrenze als uniforme
Verfeinerung. Rand-/Materialgruppen und Bilanzen bleiben gültig; jede Stufe ist
rekonstruierbar. Keine dynamische Zustandsübertragung, hängenden Knoten oder
gekoppelte Adaptivität. #36 liefert Indikatoren; M10 ergänzt den Regelkreis.

## M11–M12 — Langfristige Erweiterungshorizonte

Ergänzung vom 2026-09-22 auf ausdrücklichen Produktwunsch. Diese beiden
Milestones enthalten zunächst nur Beschreibungen; keine neuen Einzel- oder
Release-Issues, Versionsnummern oder Termine.

### M11 — Topologieoptimierung

[GitHub-Milestone](https://github.com/maschin96/FEM-Toolbox/milestone/12).

#### Zielbild
Die Toolbox erhält einen reproduzierbaren Topologie-Optimierer auf Basis ihrer verifizierten FEM-Solver. Einstieg ist die strukturelle Optimierung linear-elastischer Bauteile, etwa minimale Nachgiebigkeit bei vorgegebenem Materialvolumen; weitere thermische und gekoppelte Zielgrößen werden erst in einer späteren Detailplanung bewertet.

#### Grober Umfang
Design- und Nichtdesigngebiete, Volumenrestriktionen, mehrere Lastfälle, Sensitivitäten sowie ein geeigneter gradientenbasierter Optimierungsablauf. Dichtebasierte Verfahren wie SIMP mit Filterung/Projektion sind ein zu bewertender Ausgangspunkt; Mindestlängenskalen, Netzabhängigkeit und der Umgang mit Zwischenmaterial werden ausdrücklich berücksichtigt. Optimierungshistorie, Abbruchkriterien, Provenienz und exportierbare Ergebnisse gehören zum Workflow.

#### Abhängigkeiten und Nachweis
Aufbau auf verifizierter Mechanik, stabiler Model-API, Parameterstudien und skalierbaren Solvern. Vor Umsetzung werden Optimierungsverfahren und Sensitivitätsstrategie entschieden. Grobe Abnahme: unabhängig geprüfte Gradienten, reproduzierbare Referenzoptimierung, eingehaltene Restriktionen und erneute FEM-Prüfung der aus dem Dichtefeld abgeleiteten Geometrie; Netz- und Filtereinfluss dokumentieren.

#### Abgrenzung und Planungsreife
Keine automatische Fertigungsfreigabe und zunächst keine CFD-/FSI-Topologieoptimierung. Langfristiger Horizont nach M10, ohne feste Version oder Termin. Noch keine Einzel- oder Release-Issues; Detaillierung und Release-Gate folgen vor Umsetzung.

### M12 — RANS-CFD und Fluid-Struktur-Kopplung

[GitHub-Milestone](https://github.com/maschin96/FEM-Toolbox/milestone/13).

#### Zielbild
Erweiterung um dreidimensionale Strömungsmechanik auf Basis ausschließlich Reynolds-gemittelter Navier–Stokes-Gleichungen (RANS), einschließlich einer von Beginn an vorgesehenen bidirektionalen Fluid-Struktur-Kopplung (FSI). Zunächst stationäre, inkompressible, isotherme Strömungen mit dokumentierter Turbulenzschließung und Wandbehandlung; zeitabhängige gemittelte Modelle (URANS) und dynamische FSI bleiben eine explizite spätere Umfangsentscheidung innerhalb der RANS-Familie.

#### Grober Umfang
Geschwindigkeit, Druck, Turbulenzgrößen und Fluidmaterialien; Ein-/Auslass- und Wandbedingungen, Druck-Geschwindigkeits-Kopplung, nichtlineare Iteration, Konvergenzdiagnostik und Massen-/Impulsbilanzen. Diskretisierung, Stabilisierung, Grenzschichtnetze und geeignete Turbulenzmodelle werden vor Umsetzung bewertet; die vorhandene Festkörper-Tet4-Formulierung wird nicht ungeprüft als CFD-Diskretisierung übernommen.

#### Strukturkopplung
Stufenweise vom verifizierten Strömungssolver über einseitige Druck-/Schubspannungsübertragung zur bidirektionalen stationären FSI: Fluidlasten verformen die Struktur, die veränderte Geometrie wirkt auf die Strömung zurück. Konservative Lastübertragung, Verschiebungsübertragung auch zwischen unterschiedlichen Netzen, Schnittstellen-Gleichgewicht und Arbeitskonsistenz sowie Netzanpassung und robuste Kopplungsiteration sind einzuplanen. Zulässige Strukturverformungen bleiben explizit begrenzt. Für mögliche spätere dynamische FSI sind bewegte Netze (z. B. ALE), Zeitsynchronisation und Kopplungsstabilität gesondert zu entscheiden.

#### Abhängigkeiten und Nachweis
Aufbau auf Netz-/I/O-Infrastruktur, Strukturmechanik, Model-API und skalierbaren Solvern; dynamische FSI zusätzlich auf M6. Keine fachliche Abhängigkeit vom Topologie-Optimierer. Grobe Abnahme: getrennte CFD-Verifikation und Turbulenzmodell-Validierung mit unabhängigen Referenzen, Netz-/Wandsensitivität, Lastübertragungstest und bidirektionaler FSI-Benchmark mit Schnittstellenbilanzen. NASA Turbulence Modeling Resource bietet einen Ausgangspunkt für die Auswahl von RANS-Referenzfällen: [Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/).

#### Abgrenzung und Planungsreife
Keine DNS, LES oder hybriden RANS/LES-Verfahren. Kompressible Strömung, Mehrphasenmodelle, Verbrennung und gekoppelte Strömungs-Topologieoptimierung gehören nicht zum ersten Umfang. Langfristiger Horizont nach M11, ohne feste Version oder Termin; die Erweiterung benötigt einen eigenen mathematischen und architektonischen Entwurf. Noch keine Einzel- oder Release-Issues; Detaillierung und Release-Gate folgen vor Umsetzung.

## 6. Reihenfolge und Entscheidungen

```text
M2 abgeschlossen → M3 → M4 → M5 → M6 → M7  (bisherige Releasefolge)
                             └→ M8 → M9  (Studien und Validierung)

M7 / #36 + M8 → #67 Go/no-go → M10      (adaptive Wärme)
```

**Nächste Umsetzungswelle:** #16 → #17 → #18. Referenzbeschaffung für #19/#22,
API-Vertrag #57 und Installationsanleitung #46 können unabhängig vorbereitet
werden. Danach M4; M5 schließt stationäre Anwendung und Stabilisierung ab.
#58 ist vor dem v1.0-Gate verpflichtend.

Vor #16 sind Referenztemperatur, anisotrope Expansion und Temperaturübergabe auf
dem gemeinsamen Netz festzulegen. Vor #19/#21 werden Materialstabilität,
Rotation und feldweise Residualnormen entschieden. Vor #24 wird die Elektroden-
Nebenbedingungstechnik festgelegt. Vor #29–#32 müssen Potential, DAE-Behandlung
und Zeitintegrationsverfahren im ADR geklärt sein. Neue Entscheidungen ändern
[07_RISKS_AND_DECISIONS.md](07_RISKS_AND_DECISIONS.md), bevor abhängiger Code entsteht.

Nach jedem Release: Ist-Stand aktualisieren, Nutzerrezepte prüfen, Aufwand und
Risiko der nächsten Welle neu schätzen, unnötige Arbeit streichen. M8–M10 sind
geplante Produkthorizonte mit vorläufigen SemVer-Zielen, keine Termine. Eine
Umpriorisierung ändert Issue-Abhängigkeiten und Releaseplan gemeinsam.

## 7. Prüfprotokoll dieser Analyse

- GitHub: M0–M2 geschlossen und Releases v0.1.0/v0.2.0/v0.3.0 vorhanden;
  erfolgreiche CI auf `e3382a1`; Hauptbranch-Schutzstatus `false`.
- Ruff und Formatprüfung bestanden; mypy: 27 Quelldateien ohne Befund.
- `PYTHONPATH=src .venv/bin/python -m pytest -q`: **151 bestanden**. Die lokale
  Umgebung löste den Editable-Import auch nach Neuinstallation nicht auf; daher
  wurde der Quellpfad für die Tests explizit gesetzt. Das ist ein Quellbaumtest,
  kein Nachweis einer sauberen Artefaktinstallation. Der anfängliche Importfehler
  ist kein numerischer Testfehler und kein Red-Nachweis.
- Lokale Markdown-Dateilinks, GitHub-Zuordnungen und `git diff --check` geprüft.
  Reine Dokumentationsänderung: kein künstlicher TDD-Red-Test erforderlich.
- Keine experimentellen Daten ausgewertet, keine neue Physik implementiert,
  keine neuen Milestones freigegeben. Ein unabhängiges Review steht noch aus.
