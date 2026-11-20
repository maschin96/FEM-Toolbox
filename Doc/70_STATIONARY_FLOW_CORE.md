# Stationärer inkompressibler Strömungskern (#86)

Dieser Vertrag ergänzt vor Implementierung den MAC-Stokes-Schnitt aus Doc/69.
Der gewählte eigenständige Referenzkern ist zweidimensional, uniform kartesisch,
stationär und isotherm: rho div(u tensor u)-mu Laplace(u)+grad(p)=f,
div(u)=0. Konstant positive rho in kg/m³ und mu in Pa s. Kein eigenständiges
laminares Release, RANS-Schließung/FSI bleiben #87/#89–#91. #86 implementiert
den laminaren Grenzfall und den nichtlinearen Operator dieser Diskretisierung.
Allgemeine 3D-/unstrukturierte Adapter sind kein stillschweigender Support.

## Rand-, Netz- und Iterationsvertrag

Versetzte MAC-Freiheitsgrade wie Doc/69. Geschlossene Wände: normale
Geschwindigkeit null, vier konstante tangentiale Wandwerte; p[0,0] Gauge.
Kanal: linker Geschwindigkeitseinlass u(y)>=0, v=0; obere/untere Wände no-slip;
rechter Druckauslass p_out mit homogenem Normalgradienten der Geschwindigkeit.
Die rechte normale Geschwindigkeit ist eine Unbekannte in einem halben
Dualvolumen. Ihr Druckgradient ist 2(p_out-p_last)/dx. Tangentiale v am
Auslass verwendet zero-gradient. Dieser Druckauslass ist der Vertrag
p=p_out und du/dn=0 für den Laplace-Impulsoperator, keine beliebige FSI-Traktion.
Der Auslass legt den Druck fest; ein zusätzlicher Gauge ist unzulässig.
Negative Auslassgeschwindigkeit bedeutet nicht unterstützte Rückströmung
und endet mit expliziter Diagnose. Einlassdaten werden nie automatisch auf
einen vorgegebenen Fluss skaliert. Kanalbreite/-höhe und Profile sind in SI.

A ist der viskose MAC-Operator. W sind Dualvolumina (Auslass-u: halbes
Zellvolumen); G=-W^-1 D^T V mit Druckzellvolumen V. Der diskrete Druckterm
bleibt damit energetisch konsistent, auch am Auslass. Im Kanal werden alle
Kontinuitätszeilen gelöst, bei geschlossenen Wänden wird eine abhängige
Zeile entfernt und danach wieder überprüft.

Konvektive Dualflächenströme aus arithmetisch interpolierter normaler
MAC-Geschwindigkeit; jeder gemeinsame Fluss wird mit entgegengesetzten
Vorzeichen verwendet. Standardmäßig wählt First-order donor-cell/upwind den transportierten
Geschwindigkeitswert nach dem Vorzeichen des eingefrorenen Flusses. Diese
bewusst dissipative Referenz unterstützt beide Flussrichtungen im Inneren;
keine künstliche Druckstabilisierung. Konservative Form statt ungeprüftem
punktweisen u grad(u). Upwind hat bei glatten Lösungen erste Ordnung;
Diffusion/Stokes bleiben zweiter Ordnung. Zusätzlich ist `central` explizit
wählbar (an jeder Dualfläche |Pe|<2 erforderlich). `linear-upwind` ergänzt
die Donor-Rekonstruktion um eine Minmod-begrenzte lineare Steigung, als
verzögerte Korrektur in der Picard-Quelle. Gemeinsame Flüsse und Korrekturen
bleiben konservativ. Der Limiter ist Teil des nichtlinearen Operators und
wird für den endgültigen Residualtest neu ausgewertet. Er verhindert neue
Rekonstruktionsextrema, verspricht aber keine universelle Konvergenzordnung.
Die Verfahren werden nie stillschweigend gewechselt.
Primärreferenz für Upwind und Zentral-Péclet-Grenze:
[NIST FiPy Numerical Schemes](https://pages.nist.gov/fipy/en/3.4.5/numerical/scheme.html).

Picard: Konvektionsflüsse einfrieren, monolithischen Oseen-Schritt lösen,
Geschwindigkeit UND relativen Druck mit 0<alpha<=1 relaxieren. Der
anfängliche Zustand ist die konsistente Stokes-Lösung. Jeder Schritt wird
gegen den neu rekonstruierten NICHTLINEAREN Operator geprüft. Änderung des
Iteranden allein ist kein Abbruchkriterium. Exakte Massenbedingung bleibt
unter Relaxation zwischen konsistenten Zuständen erhalten. Maximalzahl und
Unknown-Grenze sind explizit; fehlende Konvergenz ist ein Fehler mit letzter
Residualdiagnostik. Keine zugesagte Konvergenz für beliebige Reynolds-Zahlen.

## Diagnostik, unabhängige Gates und Ressourcen

Getrennte x/y-Impulsresiduen, lokale Divergenz, globale Ein-/Auslassmassenströme,
gewichtete globale diskrete Impulsreste und Iterationshistorie. Residuen
in N/m³ und kg/s pro Einheitsbreite. Rekonstruktion aus finalen Feldwerten
und gemeinsamen konvektiven Flüssen; kein ausschließliches LU-Residuum.
Viskose Dissipation aus Differenzen und Randsegmenten, Quellen-/Randarbeit,
Druckarbeit und konvektive Leistung getrennt auswerten. Upwind-Leistung ist
kein pauschaler physikalischer Dissipationswert. Leistungsrest ist ein
zusätzlicher Gate; die Energieidentität ersetzt keine Genauigkeitsprüfung.
Absolute und relative Toleranzen sind konfigurierbar und in SI ausgewiesen.

Vorab festgelegte Referenzen:

- Voll entwickelter Kanal H=1 m, L=2 m, mu=1 Pa s, rho=1 kg/m³,
  u=y(1-y) m/s, p_out=0: dp/dx=-2 Pa/m. Inletprofil an Zellzentren.
  Vier Raster mit ny=8/16/32/64, nx=2ny: feinster Geschwindigkeitsfehler
  <0.001 m/s, Druckgradientfehler <0.02 Pa/m; Flussabweichung <1e-9 relativ.
- Nichtlinearer glatter no-slip-Fall psi=sin²(pi x)sin²(pi y), p=sin(2pi x)sin(2pi y),
  analytische Quellen einschließlich rho(u grad)u. rho=0.1, mu=1.
  Zentrale Konvektion: vier Raster 16/32/64/128; feinste RMS-Geschwindigkeit
  und Druck <0.01 in m/s bzw. Pa, Ordnung >1.8. Linear-Upwind bei rho=1:
  dieselben vier Raster, beide Fehler <0.1, Ordnung >0.7. Quellen und
  analytische Felder werden jeweils mit der zugehörigen Dichte hergeleitet.
- Nichttriviale Deckelströmung, Gaugeinvarianz, Nullfall, ungültige Netze/BCs,
  Nichtkonvergenz und Auslassrückströmung. Bilanztests an kleinen O(1)-Fällen:
  lokale Divergenz <1e-9 1/s, Impulsrest <1e-7 N/m³, Leistungsrest <1e-8 W/m.

Keine neuen Bibliotheken; NumPy/SciPy-Referenz, reproduzierbare JSON-Provenienz
und TOML-Konfiguration ohne ausführbare Serialisierung. TOML beschreibt das
gewählte uniforme Netz, Fluid und Randtyp; es ist kein Gmsh-Tet-CFD-Import.
Der vollständige RANS-/FSI-Vertrag #84 und die allgemeineren Adapter aus #85
bleiben separat; ihre für diesen Kern benötigten Festlegungen stehen hier.

## API und Konfiguration

`CartesianFluidGrid(cells, size_m)` und `IncompressibleFluid(rho, mu)` sind
validierte unveränderliche Modelle. `ClosedWalls((U_bottom,U_top,V_left,V_right))`
oder `ChannelBoundary(inlet_profile, outlet_pressure_Pa)` legen den Randtyp
fest; die Kanalprofil-Länge muss ny entsprechen. `solve_stationary_flow`
verwendet standardmäßig Donor-Cell, alpha=0.7, maximal 200 Schritte,
relative Toleranz 1e-10, absolute Toleranzen 1e-12 in N/m³, kg/s und W/m,
maximal 50000 Unbekannte. Alle Kontrollen sind explizite Python-Parameter.
`convection_scheme` unterstützt `donor-cell`, `linear-upwind`, `central`.

Schema 1: exakt `flow_schema_version`, `grid`, `fluid`, `boundary`.
Grid enthält `cells`/`size_m`; Fluid `density_kg_m3`/`viscosity_Pa_s`.
Boundary: `type="channel"`, `inlet_velocity_m_s`, `outlet_pressure_Pa`;
oder `type="closed"`, `tangential_velocities_m_s`.
Unbekannte Schlüssel, fehlende Werte und inkompatible Profile werden abgelehnt.
`FlowProblem.from_file` lädt TOML, `to_document` liefert detached Daten,
`solve` startet den Standardlauf ohne Volumenkräfte. Zusätzliche Kräfte und
Solverkontrollen gehören zur Python-API und sind in der Laufprovenienz enthalten.
Siehe [ausführbaren TOML-/CSV-/JSON-Lauf](../examples/cfd_flow/README.md).

## Akzeptanzmatrix #86

| Kriterium | Nachweis |
| --- | --- |
| Analytischer laminarer Referenzfall, Fehler/Konvergenz/Bilanzen | Poiseuille mit echtem Einlass/Druckauslass auf vier Netzen; glatte nichtlineare Manufactured Solution für zentrale und begrenzte Upwind-Konvektion |
| Druck-Gauge, Ein-/Auslass und Wände eindeutig | Zwei explizite Randtypen, Gauge nur bei ClosedWalls, outlet pressure bei ChannelBoundary; Gaugeinvarianz und inkompatibler Gauge getestet |
| Feldresiduen, Massen-/Impulsbilanzen und Iterationsstatus | Aus finalen Feldwerten rekonstruierte lokale Kontinuität, viskose Flüsse und Druckgradienten; nichtlinear neu rekonstruierte Konvektionsflüsse; getrennte Impulsfelder, globale Flüsse, Impulsreste und Historie |
| Nichtkonvergenz, Parameter/Netzfehler diagnostiziert | Maximalzahl mit letztem Residuum, Backflow, Péclet-Grenze, Profilform, positive SI-Daten, Unbekanntenlimit und Wertebereichsfehler |

#86 betrifft diesen gewählten stationären MAC-Kern. Es ist kein M12-Release,
kein allgemeiner CFD-Import aus Tet4/Tet10 und keine Turbulenzmodellvalidierung.
Referenzen bestätigen die Diskretisierung für ihren angegebenen Bereich,
nicht jede Reynolds-Zahl oder beliebige Einlass-/Auslassgeometrie.

## Ausgeführter TDD-/Verifikationsnachweis (01.10.2026)

- 11 initiale Strömungs-/Rand-/Iterations-Tests scheiterten an der fehlenden API;
  nach Implementierung bestanden sie. Zusätzliche Fälle prüfen Rückströmung,
  Gaugeinvarianz, die Stokes-Grenze und numerische Fehlerzustände.
- 8 Konfigurations-/TOML-Tests scheiterten am fehlenden FlowProblem und wurden grün.
- Der Artefaktregressionstest scheiterte am fehlenden Kanalskript; danach grün.
- Die erste Donor-Cell-MMS zeigte konkurrierende O(h)/O(h²)-Fehler und verfehlte
  das Druckbudget bei rho=1 auf dem feinsten 128er Netz (0.2266 Pa statt <0.1).
  Das Budget wurde nicht aufgeweitet. Explizite Minmod-Linear-Upwind-Rekonstruktion
  wurde als neues Verfahren implementiert: zuerst am fehlenden Verfahren rot,
  danach dieselben Fehlergrenzen und vier Netze grün. Donor-Cell bleibt ein
  ausdrücklich erster Ordnung dissipatives Verfahren und bekommt keine falsche
  Hochgenauigkeitsbehauptung. Die zentrale Variante verwendet feinere
  16/32/64/128-Netze für die asymptotische zweite Ordnung; ihre Fehlergrenze
  wurde von 0.03 auf 0.01 verschärft, nicht aufgeweitet.
- Unabhängiger affiner Innen-Patch u=x/v=-y auf 8/16/32/64 Netzen bestätigt
  analytisch Konvektion (x,y), den Donor-Cell-Leitfehler +/-h/2 und die exakte
  Rekonstruktion durch Central/Linear-Upwind im glatten Patch.

Alle 65 CFD-Tests einschließlich drei Vier-Netz-Fälle und Beispielregression
bestehen. Gesamtsuite: **761 bestanden, 83 übersprungen**, 137.84 s.
Ruff, Formatprüfung, mypy und Diff-Prüfung bestanden. Übersprungene optionale
Prüfungen gelten nicht als bestanden. Der TOML-Kanallauf konvergiert in zehn
Iterationen, maximaler Geschwindigkeitsfehler 0.000795031 m/s.

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
PYTHONPATH=src .venv/bin/python examples/cfd_flow/run.py --output /tmp/fem-cfd-86-complete
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

CI prüft den CFD-Kern in der Python-3.11–3.13-Testmatrix, zusätzlich auf ältesten
und aktuellen Bibliotheksständen. Wheel-/sdist-Jobs starten den TOML-Kanal ohne
Checkout des Quellcodes aus dem installierten Artefakt. Diese Konfiguration ist
kein vorweggenommenes bestandenes CI-Ergebnis; maßgeblich ist der aktuelle PR-Lauf.
Manuelle Diffprüfung erfolgt in dieser Sitzung; ein unabhängiger Reviewer ist
nicht verfügbar und eine unabhängige Freigabe wird nicht behauptet. Integration
geht ausschließlich in codex/milestone-m12 nach erfolgreichen CI-Prüfungen;
kein M12-Release und kein Merge nach main.
