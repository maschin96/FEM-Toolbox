# Entwicklungsstand und nächste Meilensteine M6/M7

Stand: 23.09.2026. Tracking: [#167](https://github.com/maschin96/FEM-Toolbox/issues/167).
Diese Planung ersetzt für M6/M7 die ältere Baseline aus Dokument 13.
Sie beschreibt geplante Arbeit, keine bereits implementierten Fähigkeiten.

M6 ist mit v1.1.0 abgeschlossen; [Release-Nachweis](35_M6_RELEASE_VERIFICATION.md).
Die nachfolgende Bestandsaufnahme ist die historische v1.0-Planungsbaseline.

## Bestandsaufnahme und Bewertung

Geprüfte Basis: `origin/main` und Release v1.0.0 bei
`c723312b63c2e5301ce39807b10e7b3fc893915a`. Der anfängliche lokale
Release-Branch hatte denselben Dateiinhalt. M0–M5 sind abgeschlossen;
M6 (#29–#33) und M7 (#34–#38) hatten jeweils fünf offene Issues ohne
abgeschlossene Arbeitspakete. Der veröffentlichte Release enthält Wheel,
sdist und Prüfsummen; die main-CI war mit 17 Jobs erfolgreich.

| Vorhanden | Beleg | Konsequenz für die nächsten Schritte |
|---|---|---|
| Stationäre Einzelfelder und Kopplungen einschließlich Elektroden | `solvers/`, Patch-/Konvergenztests, v1.0-Release | Belastbare Referenz bewahren; keine Neuentwicklung des stationären Kerns |
| Unveränderliche Materialien, Regionen, DOFs und COO-Assembly | `materials/`, `fields/`, `assembly/` | Ergänzbar, aber Masse/Kapazität und Zeitstatus fehlen |
| Strikt reelle lineare Systeme und Ergebnisarrays | `assembly/LinearSystem`, `io/result_bundle.py` | Harmonische Analyse benötigt expliziten komplexen Datenvertrag |
| Direkte SciPy-LU-Aufrufe in Feldsolvern | `solvers/heat.py`, `elasticity.py`, `piezoelectric.py` | Vor PETSc echte Backend-Grenze schaffen |
| Tet4-/Triangle3-spezifische Netze und Exporte | `mesh/model.py`, `mesh/topology.py`, `io/` | Tet10 betrifft Topologie, Randintegration, Model und Postprocessing zugleich |
| Stationäres TOML-Modell, Preflight und Ergebnisbündel v1 | `model.py`, `model_preflight.py`, Dokumente 21/23/24/29 | Transienten ergänzen statt bestehende Formate umzudeuten |
| Python 3.11–3.13, Qualitäts-, Kompatibilitäts- und Artefakt-CI | `pyproject.toml`, `.github/workflows/ci.yml` | Gute Basis; echte PETSc-Tests fehlen noch |

Die Toolbox ist eine verifizierte stationäre Forschungsbasis. Experimentelle
Validierung, belastbare Großmodell-Skalierung, Transienten, PETSc und Tet10
sind nicht nachgewiesen. Tests allein belegen keine allgemeine Einsatzreife.
Die alten Aussagen „Planungsphase“ in AGENTS.md und M0–M3-Status in der
Roadmap waren überholt; Dokument 13 bleibt als historische Analyse erhalten.
GitHub enthält inzwischen auch detaillierte spätere Horizonte bis M14;
deren Neuplanung ist nicht Teil dieses Auftrags.

## Planungsprinzip und Reihenfolge

Der Fokus bleibt erhalten: **M6 Zeit/Dynamik**, danach **M7 Skalierung und
höhere Elemente**. Beide sind mehrere Entwicklungsschritte, keine einzelnen
Feature-PRs. Bestehende Issues werden als Arbeitspakete weiterverwendet;
Unterabschnitte a–d sind getrennte PR-Schnitte, keine duplizierten Backlogs.

M6: #29 → #30 und #31 → #33. #32 ist seit der Nutzerentscheidung
vom 23.09.2026 zurückgestellt und kein v1.1-Release-Gate. Der mathematische Vertrag für #32 kann
früh vorbereitet werden; die Implementierung wartet auf die verifizierten
Einzelfelder. M7: #37a zuerst, danach #34 und #35; #36 ist fachlich unabhängig
von PETSc. Optimierung #37b/c und Release #38 folgen den Referenznachweisen.

Wichtigste Ergänzungen sind Model-/Exportintegration **innerhalb** der jeweiligen
Issues, frühe DAE-/Phasorentscheidungen und eine durchgängige Supportmatrix.
#36 wird bewusst auf Tet4-Wärme begrenzt: Das erfüllt die Voraussetzung für
M10, ohne M7 zusätzlich mit allgemeiner gekoppelter Adaptivität zu überladen.
Die erste Tet10-Freigabe verwendet geradkantige Geometrie; gekrümmte und
hybride Netze benötigen eigene Nachweise.

M6 trägt auf GitHub weiterhin den Zieltermin 14.10.2026. Ohne Kapazitäts- und
Aufwandsschätzung ist dieser nicht belastbar; er wird hier weder bestätigt
noch eigenmächtig verschoben. Zielversionen bleiben v1.1.0 und v1.2.0.
Releasegates haben Vorrang vor Kalenderdaten. Fehlende Referenzen oder
unerfüllte Gates erfordern eine sichtbare Neuplanung, keine stillen Kürzungen.

## Konkrete Arbeitspakete

### #29 — Zeitbasis, Materialien und Serienvertrag

Abhängigkeit: M5 abgeschlossen; erster M6-Schritt.

- Dichte rho [kg/m³] und spezifische Wärmekapazität c [J/(kg K)] regionsweise validieren; stationäre Materialkonstruktoren bleiben verwendbar. Konsistente und explizit gewählte gelumpte Tet4-Masse/Kapazität gegen analytische Integrale prüfen.

- Zeitachse, Anfangswerte und Last-/Randwertverläufe definieren: endliche streng steigende Zeiten, Interpolation, Sprungstellen und Auswertezeitpunkte. TOML verwendet deklarative Zeitreihen, keine ausführbaren Ausdrücke.

- Vor konkreten Solvern den additiven Model-/Ergebnisvertrag festlegen: Zustand T bzw. u/v/a, Zeitindex, SI-Einheiten, Ausgabeintervall, Provenienz und Speicherbudget. Schemas getrennt versionieren; Model-v1 und stationäres Bündel bleiben lesbar.

- Abnahme: Elementintegrale erhalten Gesamtmasse/Kapazität; ungültige Material-, Zeit- und Anfangsdaten werden deterministisch abgelehnt. Serienexport mit unabhängiger Leserprüfung einschließlich unvollständiger Läufe.

### #30 — Transiente Wärmeleitung

Abhängigkeit: Nach #29; liefert thermische Referenz für #32.

- Backward Euler als verpflichtende Referenz implementieren; Crank–Nicolson bleibt optionale Folgearbeit. Anfangstemperatur, zeitabhängige Quellen, Flüsse, Konvektion und Dirichlet-Werte berücksichtigen.

- Wärmespeicherung und Randreaktionen aus dem ursprünglichen zeitdiskreten System bilanzieren. Isolierte Körper ohne Temperaturreferenz sind bei positiver Kapazität zulässig; stationären Nullraum-Preflight nicht ungeprüft übernehmen.

- Abnahme: analytische Abkühlung und Diffusion, isolierter Körper, Mehrmaterialfall und zeitabhängiger Dirichlet-Rand. Raum- und Zeitfehler getrennt verfeinern; zeitlich erste Ordnung und stationären Grenzfall nachweisen.

- Ausführbares TOML-Beispiel mit Zeitserie und schrittweiser Energiebilanz aus installierter Distribution bereitstellen.

### #31 — Modal-, Zeit- und harmonische Piezoelastodynamik

Abhängigkeit: Nach #29; Reihenfolge 31a → 31b → 31c → 31d, jeweils separat prüfbarer PR.

- 31a: DAE-/Elektrodenvertrag und modale Analyse. Nur mechanische DOFs tragen Masse; elektrische DOFs bleiben algebraisch. Open-/Short-Circuit, Gauge, Floating-Charge und konsistente Anfangswerte festlegen. Keine künstliche elektrische Masse; keine ungeprüfte Eigenwertlösung des singulären Gesamtmassenproblems.

- 31b: Newmark mit beta=1/4, gamma=1/2 als ungedämpfte Referenz; bewegte Dirichlet-Ränder samt Geschwindigkeit/Beschleunigung und deren Trägheitsbeiträge konsistent behandeln. Numerische/physikalische Dämpfung getrennt ausweisen; generalized-alpha nicht verpflichtend.

- 31c: harmonischen Pfad mit expliziter Phasorkonvention und komplexwertigen Lasten/Ergebnissen entwickeln. Bestehendes LinearSystem akzeptiert nur reelle Werte: additive komplexe Infrastruktur oder äquivalente reelle Blockform vorab entscheiden; kein stiller Verlust imaginärer Anteile. Resonanz/Singularität diagnostizieren.

- 31d: Model-Einstieg und Zeit-/Frequenzexport einschließlich Real-/Imaginärteil, Amplitude/Phase und Provenienz integrieren. Frequenz in Hz und Kreisfrequenz eindeutig unterscheiden.

- Abnahme: unabhängige Eigenfrequenzreferenz, Open-/Short-Circuit-Vergleich, ungedämpfte Energie, zweite zeitliche Ordnung bei glatten Daten, Zeit-/Frequenzvergleich abseits Resonanz und statischer Niederfrequenzgrenzfall. Raum-, Zeit- und Frequenzrasterfehler getrennt berichten.

### #32 — Reziproke lineare Thermoelastizität (zurückgestellt)

Nutzerentscheidung 23.09.2026: Verformung wirkt derzeit nicht auf Temperatur
zurück. Der folgende Abschnitt ist zukünftiger Backlog, kein M6-Pflichtumfang.

Abhängigkeit: Nach #30 und mechanischer Zeitintegration aus #31b; Potential-/Referenzklärung vor Implementierung.

- Thermodynamisches Potential, Referenztemperatur, Entropieterm, Vorzeichen und zulässige Materialparameter im ADR und mathematischen Modell herleiten. Erste Freigabe auf mechanisch-thermische Reziprozität begrenzen; allgemeine reziproke Thermopiezo-/Pyrokopplung bleibt außerhalb.

- Gekoppelten transienten u–T-Pfad samt Model-Konfiguration, Anfangs-/Randbedingungen und Energie-/Entropiebericht implementieren. Bestehendes einseitiges T→u und T→(u,phi) behält seine Semantik.

- Abnahme: entkoppelter, quasistatischer und adiabatischer Grenzfall, unabhängige Referenz sowie konvergenter Thermoelastic-Damping-Fall. Physikalische Wärmeleitungsdämpfung von Integratordämpfung und Netzfehler trennen.

- Gate bleibt offen, wenn Potential oder unabhängige Referenz ungeklärt sind; Umfang nicht stillschweigend aus v1.1 entfernen.

### #33 — M6-Release-Gate

Abhängigkeit: Nach #29–#31 einschließlich ihrer Integrationsschritte. #32 ist zurückgestellt.

- Alle M6-Pflichtkriterien einschließlich Model-/Serien-/Frequenzexport nachweisen; alte v1.0-Python- und TOML-Beispiele bleiben unverändert lauffähig.

- Thermische Zeitordnung, mechanische Zeitordnung, Eigenfrequenzen, feldweise Residuen und Energiebilanzen mit vorab begründeten Toleranzen dokumentieren; Raum- und Zeitfehler getrennt.

- Wärme- und Piezo-Dynamik-/Harmonik-Beispiele aus frischen Wheel-/sdist-Installationen prüfen. Vollständige CI am Integrationscommit, Reviewstatus, Migration und Grenzen dokumentieren.

- v1.1.0 erst nach Gate, konsistenter Paketversion und annotiertem Tag freigeben. Bestehender Zieltermin 14.10.2026 ist ungeprüft: Aufwand und verfügbare Kapazität vor einer Terminzusage erheben.

### #34 — Optionales PETSc-Backend

Abhängigkeit: Vertrag und SciPy-Referenz zuerst; M6-Analyseverträge berücksichtigen. #37 liefert frühe Baseline.

- 34a: interne Backend-Schnittstelle für Operator, RHS, Feldindexmengen nach Randbedingungs-/Elektrodenkondensation, Skalierung und SolveReport mit SciPy implementieren. Bestehende splu-Aufrufe sind noch kein austauschbares Backend.

- 34b: optionales petsc4py mit dokumentierter getesteter PETSc-/Python-/Scalar-Konfiguration anbinden; reelle und komplexe Analysepfade explizit abdecken oder komplexe Systeme äquivalent reell darstellen. Core-Import darf kein PETSc benötigen.

- 34c: FieldSplit-/Schur-Verfahren, Nullräume, physikalische Residuen, Iterationslimit und Nichtkonvergenz diagnostizieren. Kein stiller Fallback, der Backend-Fehler verdeckt. Mindestens ein echter PETSc-CI-Job; übersprungene optionale Tests allein sind kein Nachweis.

- Abnahme: SciPy/PETSc-Parität für Einzelfelder, gekoppelte feste/schwebende Elektroden und unterstützte M6-Systeme einschließlich Harmonik. Residuen am ursprünglichen unskalierten System prüfen; DOF-Indizes nach Kondensation testen.

- Einprozess-PETSc ist Mindestumfang. MPI-/verteilte Assembly erst nach Profiling und explizitem Architekturentscheid; keine verteilte Skalierung aus einem Backendwechsel ableiten.

### #35 — Tet10 durch alle Schichten

Abhängigkeit: Nach M6-Verträgen; 35a → 35b → 35c → 35d. PETSc ist keine Voraussetzung der Elementverifikation.

- 35a: additive Tet10-/Triangle6-Topologie, Knotenordnung, Gmsh-/VTU-Rundreise, Auswahlsätze, Randtopologie und Materialregionen. Tet4-Konstruktoren und bestehende Model-Dateien kompatibel halten; hybride Tet4/Tet10-Übergänge zunächst explizit ablehnen.

- 35b: quadratische Formfunktionen, Volumen-/Flächenquadratur, Ableitungen und Orientierung. Erste Freigabe auf geradkantige Tet10 mit mittigen Kantenknoten beschränken; gekrümmte Geometrie gezielt ablehnen und dokumentieren. Kein ungeprüfter Wechsel des Geometriemodells.

- 35c: mechanische, thermische, elektrische und piezoelektrische Operatoren einschließlich M6-Masse/Kapazität, thermischer Kopplung und Randintegralen. Elektroden müssen alle sechs Flächenknoten berücksichtigen.

- 35d: Preflight, Model-Auswahl, Backend-Anbindung und Ergebnisexport integrieren. Spannungen/Flüsse sind bei Tet10 nicht elementweise konstant: Auswerteorte und Quadratur-/Roh-/Projektionswerte eindeutig kennzeichnen.

- Abnahme: Partition der Eins, Ableitungssummen, unabhängige Quadraturintegrale, Patch- und Bilanztests; glatte Manufactured Solutions mit erwarteter H1-Ordnung 2/L2-Ordnung 3 auf mindestens vier Netzen. Singularitäten nicht für Ordnungszusagen verwenden. Biegungs-/Gradientenvergleich bei gleichem Fehler und Kostenbudget; Tet10 behebt Inkompressibilitätslocking nicht automatisch.

### #36 — Fehlerindikatoren als begrenzte Forschungsbasis

Abhängigkeit: Nach verifizierter stationärer Wärme; kann unabhängig von PETSc beginnen. Anschluss an M10 #67–#69.

- Ersten residual- oder recovery-basierten Indikator ausschließlich für stationäre skalare Wärme auf Tet4 liefern. Norm, Einheiten, Randterme und Materialflusssprünge definieren; bei Tet4 nicht nur den oft verschwindenden Zellresidualterm auswerten.

- Lokale Indikatoren und deterministische Markierung als Forschungs-API exportieren; Abgrenzung zwischen heuristischem Indikator und garantiertem Fehlerschätzer dokumentieren.

- Abnahme: mindestens vier Netze mit bekannter glatter Lösung, Effektivitätsindex und lokalisierter Gradientenfall, zusätzlich Materialgrenze. Vorab vereinbarte Norm-/Trendkriterien prüfen; bloße optische Korrelation reicht nicht.

- Gekoppelte Indikatoren, Tet10-Indikatoren und automatisches Remeshing sind Folgearbeit nach eigenem Gate. M10 erhält damit seine benötigte Tet4-Wärmebasis; kein allgemeiner Adaptivitätsanspruch.

### #37 — Profiling und nachgewiesene Skalierung

Abhängigkeit: 37a als M7-Einstieg; Optimierung erst nach Messung, Abschluss nach #34/#35.

- 37a: reproduzierbare Baseline auf festgelegter Hardware erfassen: Tet4-Einzelfeld, gekoppeltes Elektrodenmodell und wiederholte M6-Schritte. Geometrie, Assembly, Randbehandlung, Faktorisierung/Solve und Export getrennt messen; Peak-Speicher, DOFs, nnz, Threads und Bibliotheksversionen angeben.

- 37b: aus Profilen Engpass und vorab messbares Verbesserungsziel wählen; Matrixkopien, COO-Speicher und Wiederverwendung konstanter Operatoren/Faktorisierungen prüfen. Ungültigmachung bei geänderten Randwertstrukturen/Zeitschritten testen. Vektorisierung/Parallelisierung nur bei belegtem Nutzen.

- 37c: Tet4/Tet10 und SciPy/PETSc bei vergleichbarer Fehlergrenze sowie mehreren Problemgrößen vergleichen. Wiederholungen und Streuung berichten; keine Hardware-Zeitgrenze als instabiles Unit-Test-Gate.

- Abnahme: numerische Parität, belegte Verbesserung im gewählten Ziel, dokumentierte Regressionen und Speichergrenzen. Ohne MPI-Nachweis nur Problemgrößenskalierung behaupten; optimierte Pfade behalten eine prüfbare Referenz.

### #38 — M7-Release-Gate

Abhängigkeit: Nach #34–#37; M6-Kompatibilität bleibt Pflicht.

- Tet4/Tet10 × SciPy/PETSc × freigegebene stationäre/M6-Analysearten als Supportmatrix dokumentieren und prüfen; nicht unterstützte Kombinationen werden früh mit konkretem Fehler abgelehnt, nicht implizit freigegeben.

- Tet10-Patch-/Konvergenz-/Bilanztests, Backend-Parität und Wärmeindikatornachweise sowie reproduzierbarer Performancebericht liegen vor. Echte PETSc-Ausführung einschließlich komplexer Konfiguration dokumentieren.

- Wheel/sdist ohne PETSc und dokumentierte optionale Installation mit PETSc testen; v1.0/v1.1-Beispiele und Model-Schemas regressionssicher. Qualitätschecks und anwendbare CI am Zielcommit grün; Reviewstatus offenlegen.

- v1.2.0 erst nach vollständig erfülltem Gate und annotiertem Tag. Keine unbelegte Leistungs-, MPI-, Adaptivitäts- oder experimentelle Validierungszusage.

## Prüfprotokoll dieser Bestandsaufnahme

- Lokale Quelltests: `PYTHONPATH=src .venv/bin/python -m pytest -q`: **286 bestanden**.
- Ruff-Lint und Formatprüfung bestanden (136 Dateien); mypy bestanden (38 Quelldateien).
- Der direkte pytest-Aufruf scheiterte lokal am Paketimport, auch nach erneuter
  editierbarer Installation. Deshalb ist der lokale Lauf ausdrücklich ein
  Quelltest mit gesetztem Importpfad, kein frischer Installationsnachweis.
  Die Ursache der lokalen Umgebung wurde in diesem Planungsauftrag nicht behoben.
- Der veröffentlichte v1.0-Artefaktnachweis und die erfolgreiche main-CI sind
  vorhandene Release-Evidenz; Artefaktinstallationen wurden hier nicht erneut geprüft.
- Reine Dokumentations-/Planungsänderung: kein Red-Test oder neuer numerischer
  Nachweis erforderlich. Relative Dokumentationslinks und `git diff --check` geprüft.
- Kein unabhängiger Reviewer in diesem Auftrag; keine Reviewfreigabe behauptet.


## Umsetzung begonnen (23.09.2026)

#29: Tet4-Speicheroperatoren, optionale thermische rho/c-Daten, kleine
quellenbelegte Bibliothek und Sprung-/Rampen-/CSV-Zeitverläufe sind implementiert.
#30: Backward Euler, Model-Schema 2, Speicherbilanz, ParaView-Serie und
Messpunkt-CSV sind implementiert und lokal verifiziert; Nachweise in [Doc/33](33_TRANSIENT_HEAT.md).
Issue-Abnahme, PR-Integration und M6-Freigabe sind davon getrennte Schritte.
Der dynamische u/v/a-/Frequenzvertrag und #31 bleiben ausstehend.
#32 ist auf Nutzerwunsch zurückgestellt (ADR-021).

## M7 begonnen (25.09.2026)

#37a: reproduzierbares Benchmarkwerkzeug und 27 lokale SciPy-/Tet4-Messläufe
für stationäre Wärme, Floating-Elektrode und transiente Wärme liegen vor;
[Messvertrag, Rohdaten und Grenzen](36_PERFORMANCE_BASELINE.md).
#37b/c und Issue-Abschluss bleiben offen.

#34a: Die interne Backend-Grenze mit SciPy-Referenz ist im Issue-Branch
implementiert; [Vertrag und Nachweise](37_LINEAR_BACKEND.md). Einzelfelder,
gekoppelte Elektroden und lineare M6-Teilprobleme nutzen sie.

#34b: Der optionale PETSc-Direktadapter und echte CI-Jobs für float64/complex128
sind im Issue-Branch umgesetzt. Konfiguration, Scalar-Abbildung, Parität und
Grenzen stehen in [Dokument 38](38_PETSC_BACKEND.md).

#34c: FGMRES mit FieldSplit/Schur, Nullraum- und Nichtkonvergenzdiagnostik sowie
zusätzliche reale/komplexe Physik- und Konvergenznachweise sind im Issue-Branch
umgesetzt; [Vertrag und Grenzen](39_PETSC_FIELDSPLIT.md). Issue-Abschluss erst
nach erfolgreicher CI, Review und Integration. Backend-Auswahl im Model bleibt
eine spätere Erweiterung.
#35a: Additive Tet10-/Triangle6-Topologie, Randzuordnung, Materialregionen und
Gmsh-/VTU-I/O sind im Issue-Branch umgesetzt; [Vertrag und Grenzen](40_TET10_TOPOLOGY.md).
#35b: Quadratische Formfunktionen, Ableitungen, affine Tet10-/Triangle6-Geometrie
und positive Quadratur sind im Issue-Branch umgesetzt;
[Vertrag und Interpolationsnachweise](41_TET10_GEOMETRY.md).
#35c: Quadratische Feldoperatoren, konsistente/HRZ-Speicher, thermische Kopplung
und Randlasten sind im Issue-Branch umgesetzt; die Elektrodenreduktion erfasst
alle sechs Flächenknoten. [Vertrag und Nachweise](42_TET10_OPERATORS.md).
Nächster geplanter Funktionsschritt: #35d (Model-/Solver-/Ergebnisintegration).


#36: Der stationäre Tet4-Wärmeindikator mit Quellen-, Flusssprung- und
Randresiduen sowie deterministischer Dörfler-Markierung ist im Issue-Branch
implementiert. Vier-Netze-Effektivität, lokalisierter Fehler und Materialgrenze
sind geprüft; [Vertrag und Grenzen](43_HEAT_INDICATOR.md). Integration und
Issue-Abschluss bleiben den CI-/Reviewgates unterworfen.


#37b: Die bedarfsgerechte Tet4-Flächenauswertung ist im Issue-Branch umgesetzt.
Bei fünf Wiederholungen je Fall/Größe sinkt die Elektroden-Medianlaufzeit für
N=6/N=10 um 20,5 %/18,7 %; beide Wärmefälle und gemessene Speicherpeaks verbessern
sich ebenfalls. [Rohdaten, Parität und Grenzen](44_GEOMETRY_PERFORMANCE.md).
#37c und der Issue-Abschluss bleiben nach vollständiger Tet10-Integration offen.

## Fortschritt #35d: thermischer Integrationsschnitt

Stationäre und transiente Tet10-Wärme sind durch Model, Preflight, lineare
Backend-Grenze und versionierten Export angebunden. Der Schnitt enthält
Gmsh-/TOML-Beispiele sowie Raum- und Zeitkonvergenznachweise. #35 bleibt für
die übrigen Feldintegrationen und Vergleichsgates offen; #37c folgt danach.
[Details und Grenzen](45_TET10_HEAT_INTEGRATION.md).

## Fortschritt #35d: elektrischer Integrationsschnitt

Stationäre Tet10-Elektrostatik ist einschließlich Model/Preflight, Backend,
Ladungsbilanz und Schema-2-Export angebunden. E/D sind Schwerpunktwerte;
Energie wird über das Element integriert. Vier-Netze-Nachweis und Gmsh-/TOML-
Beispiel ergänzen die Tests. Mechanik, Piezo, Kopplungen und Vergleichsgates
bleiben offen. [Details](46_TET10_ELECTROSTATICS.md).

## Fortschritt #35d: mechanischer Integrationsschnitt

Statische Tet10-Elastizität ist durch Model/Preflight, Backend, Lasten und
Schema-2-Export angebunden. Dehnungen/Spannungen sind Schwerpunktwerte;
Verzerrungsenergie wird integriert. Biegepatch, anisotroper Patch, Materialgrenze,
Kraft-/Momentenbilanzen und vier Netze sind nachgewiesen. Thermische Tet10-
Eigendehnungen, Piezo, Dynamik und Kostenvergleiche bleiben offen.
[Details](47_TET10_ELASTICITY.md).

## Fortschritt #35d: thermomechanischer Integrationsschnitt

Räumlich veränderliche thermische Tet10-Eigendehnungen und der bestehende
Python-Ablauf Wärme → Mechanik sind integriert. Lasten und mechanische Energie
werden punktabhängig ausgewertet, beide Backends sind separat wählbar und der
gemeinsame Schema-2-Export prüft beide Auswerteort-Arrays. Direkte Patches und
gekoppelte Vier-Netze-Konvergenz sind nachgewiesen. Kein neuer TOML-Modus;
Piezo-/Thermopiezo, Dynamik und Kostenvergleiche bleiben offen.
[Details](48_TET10_THERMOMECHANICS.md).

## Fortschritt #35d: stationärer piezoelektrischer Integrationsschnitt

Tet10-Piezo ist mit Model/Preflight, anisotropen Materialregionen, festen und
schwebenden Triangle6-Elektroden, Backends und Schema-2-Export integriert.
Quadratische gekoppelte Felder, integrierte Enthalpie, Bilanzen und vier Netze
sind nachgewiesen; [Details](49_TET10_PIEZOELECTRIC.md). Thermopiezo,
Piezo-Dynamik und Gesamt-Abnahme-/Kostenvergleiche bleiben offen.

## Fortschritt #35d: thermopiezoelektrischer Integrationsschnitt

Stationäres T → (u, φ) ist für Tet10 durch Python-API, Model, getrennte Backends
und Schema-2-Export integriert. Beide elektrischen Thermalkonventionen werden
punktabhängig integriert; analytische gekoppelte Patches, schwebende Elektroden
und Vier-Netze-Konvergenz sind geprüft. [Nachweise](50_TET10_THERMOPIEZOELECTRIC.md).
Tet10-Piezo-Dynamik und Gesamt-Abnahme-/Kostenvergleiche bleiben offen.
