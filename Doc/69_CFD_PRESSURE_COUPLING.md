# Experimentelle Druck-Geschwindigkeits-Kopplung (#86)

Verbindlicher Schnittvertrag vor Implementierung: stationäre zweidimensionale
inkompressible Stokes-Gleichungen -mu Laplace(u)+grad(p)=f, div(u)=0.
Konstante mu>0 in Pa s, f in N/m³, u in m/s und p in Pa. Kein
Konvektionsterm, keine RANS-Schließung, kein allgemeiner Navier-Stokes-Solver.

## Diskretisierung und Randbedingungen

Uniformes rechteckiges MAC-Gitter: p[ny,nx] an Zellzentren,
u[ny,nx+1] an vertikalen Flächen, v[ny+1,nx] an horizontalen Flächen.
An allen vier Wänden normale Geschwindigkeit null; tangentiale konstante
Wandwerte in Reihenfolge unten/oben/links/rechts. Nichtnullwerte an den
Ecken können unstetig sein; der Lidfall ist deshalb nur ein Bilanztest,
kein glatter Konvergenzbenchmark. Offene Grenzen werden nicht simuliert.

Innenflächen verwenden zentrale zweite Differenzen. Tangentiale Wände
verwenden Dirichlet-Ghostwerte 2 U_wall-u_first, d.h. zusätzlichen diagonalen
Beitrag mu/h² und rechte Quelle 2 mu U_wall/h². Normale Randwerte sind
explizit eliminiert. Zellweise Flussdifferenzen bilden D; G=-D^T ist der
Druckgradient an freien Geschwindigkeitsflächen. Daher ist Druckarbeit bei
D u=0 null, und kein entkoppelter alternierender Druckmodus bleibt übrig.
Die geometrisch gleich gewichteten freien Geschwindigkeitsvolumina erlauben
diese diskrete Adjungiertheit auf dem uniformen Gitter.

Das monolithische Sattelpunktsystem lautet A u-D^T p=b, D u=0.
p[0,0]=pressure_reference_Pa; die entsprechende abhängige Kontinuitätszeile
wird beim Lösen entfernt und danach ausdrücklich wieder geprüft.
Die Druckreferenz ist ein Gauge, keine physikalische Druckrandbedingung.
LU mit Zeilen-/Spaltenäquilibrierung löst die Referenzimplementierung.
Kein CG auf dem indefiniten Gesamtsystem, keine Strafkompressibilität,
keine künstliche Druckglättung. Dieser kleine Referenzsolver verwendet SciPy;
PETSc und skalierbare iterative Verfahren werden später bewertet.
Siehe unabhängiges [PETSc-Stokes-Beispiel](https://petsc.org/main/src/dm/impls/stag/tutorials/ex2.c.html)
für das versetzte Gitterprinzip; hier wird kein PETSc-Code übernommen.

## Rekonstruktion und Gates

Ergebnisse enthalten getrennte x-/y-Impulsresiduen (N/m³), Zell-Divergenz
(1/s), Zell-Volumenstromrest pro Einheitsbreite (m²/s), Status und
Dissipation/Leistungsrest (W/m). Alle Arrays sind owned/immutable.
Dissipation wird aus Geschwindigkeitsdifferenzen an Innen- und Randsegmenten
rekonstruiert, nicht allein aus u^T A u. Randarbeit wird separat berechnet;
Volumenkraftarbeit=sum(f u) dx dy. Ohne Massenerhaltung gilt eine kleine
Impulsrestnorm nicht als erfolgreiche Kopplung. Relative Grenzen besitzen
explizite absolute SI-Untergrenzen; ungültige Parameter und unaufgelöste
Bilanzen enden mit ValidationError. Eine Unknown-Grenze schützt den direkten
Referenzsolver vor unbeabsichtigt riesigen Gittern.

## Unabhängige Tests vorab

Konstantes f=(1,2) N/m³ im geschlossenen Rechteck: u=0 und p=x+2y+C Pa.
Das testet den Druckgradienten und Gauge unabhängig von der Assembly.
Ein bewegter Deckel prüft nichttriviale u/v, lokale Kontinuität, Impuls und
Leistung; eine Gaugeverschiebung darf Geschwindigkeiten nicht verändern.

Glatter Einheitsquadratfall mit psi=sin²(pi x) sin²(pi y):
u=pi(1-cos(2pi x))sin(2pi y),
v=-pi sin(2pi x)(1-cos(2pi y)), p=sin(2pi x)sin(2pi y).
f wird analytisch durch -mu Laplace(u)+grad(p) hergeleitet. Alle Wände
sind no-slip. Vier Netze 8/16/32/64: L2-RMS-Geschwindigkeitsordnung >1.8,
Druckordnung >1.5, feinster Fehler <0.01 m/s beziehungsweise Pa.
Absolute Bilanzgrenzen bei diesen O(1)-Referenzen: Divergenz 1e-9 1/s,
Impuls 1e-9 N/m³ im kleinen Lidfall, Leistung 1e-10 W/m.
Referenzdruck wird um denselben ersten Zellwert verschoben.

#84–#86 bleiben offen: allgemeine Fluidnetze, offene Randbedingungen,
Konvektion, nichtlineare Iteration und RANS/FSI fehlen. Dieser Schritt
schließt den Druckkopplungsanteil für den dokumentierten Stokes-Vertrag.

## Ausgeführter Nachweis (01.10.2026)

Red: alle 15 initialen Tests scheiterten an der fehlenden solve_stokes-API.
Green: nach Implementierung bestanden sie. Zwei ergänzte Tests prüfen das
kleinste Gitter und unauflösbare Geometrieskalierung. Nach Refactoring der
Leistungsrekonstruktion bestanden alle 31 CFD-Tests (14 Kanal, 17 Stokes).

Beispiellauf: RMS-Geschwindigkeitsfehler auf 8/16/32/64 Zellen je Achse:
0.218126 / 0.0514641 / 0.0125836 / 0.00311633 m/s.
Druckfehler: 0.0136358 / 0.00323661 / 0.000804243 / 0.000200856 Pa.
Geschwindigkeitsordnungen 2.08352 / 2.03202 / 2.01363;
Druckordnungen 2.07485 / 2.00878 / 2.00147. Maximaler Divergenzrest
auf dem feinsten Gitter 6.98617e-12 1/s. Vorab gesetzte Gates unverändert.

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
PYTHONPATH=src .venv/bin/python -m pytest -q tests/unit/test_cfd_channel.py tests/unit/test_cfd_stokes.py
PYTHONPATH=src .venv/bin/python examples/cfd_pressure/run.py
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

Gesamtsuite: **727 bestanden, 83 übersprungen**, 96.20 s. Gezielte Tests
nach Refactoring, Beispiel, Ruff, Formatprüfung und mypy bestanden.
Übersprungene optionale Prüfungen gelten nicht als bestanden. GitHub-CI
und unabhängiges Review ersetzen diese lokalen Nachweise nicht; keine
Integrationsfreigabe behauptet und #86 bleibt offen.

## Erweiterung zu #86

Dieser Nachweis bleibt der historische Stokes-Schnitt. Konvektion,
Picard-Iteration und Ein-/Auslass sind inzwischen im
[stationären Referenzkern](70_STATIONARY_FLOW_CORE.md) ergänzt.
