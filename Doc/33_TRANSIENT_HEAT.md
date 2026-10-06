# Transiente Wärmeleitung und Zeitserien

Status: Bestandteil von v1.1.0; Verifikation und Release-Gate siehe
[Dokument 35](35_M6_RELEASE_VERIFICATION.md). Keine reziproke thermomechanische Rückwirkung
(ADR-021), keine nichtlinearen/temperaturabhängigen Materialien, kein Phasenwechsel.

## Modell und Bilanz

Für konstante positive Dichte rho und spezifische Wärmekapazität c:

```text
rho c dT/dt - div(k grad T) = Q
C (T[n] - T[n-1])/dt + (K + H[n]) T[n] = f[n] + r[n]
```

`C` ist die konsistente Tet4-Kapazität, alternativ ausdrücklich zeilensummiert
mit `lumped=True`. `H` enthält Konvektion. Der Referenzintegrator ist Backward
Euler (erste Ordnung); jeder Schritt wird mit diagonal skaliertem SciPy Sparse
LU gelöst. Positive Kapazität macht auch isolierte Körper ohne künstliche
Temperaturreferenz lösbar. Elektrische oder mechanische Zeitintegration ist
nicht Bestandteil dieses Wärmesolvers.

Das ursprüngliche zeitdiskrete Residuum ist
`C ΔT/dt + (K+H)T - f`. An festen Temperaturknoten ist es die Wärmeleistung
**in** den Körper aus dem Reservoir. Die Bilanz in W lautet:

```text
source_power + dirichlet_inflow - prescribed_outflow
- convection_outflow - storage_power = balance
```

`stored_energy_change = 1ᵀ C (T[n]-T[0])` in J misst die Speicheränderung.
Reaktionen bei t0 sind noch nicht definiert: Initialframe enthält Nullen,
das Manifest kennzeichnet sie als nicht ausgewertet. Ein reines stationäres
`K T-f` wäre für zeitabhängige Temperaturvorgaben falsch.

## Python-API und Zeitvertrag

```python
from fem_toolbox.materials.catalog import thermal_material
from fem_toolbox.solvers.transient_heat import solve_transient_heat
from fem_toolbox.time import Step

result = solve_transient_heat(
    mesh,
    thermal_material("aluminium_6060_t6_20c"),
    times=[0, 1, 2, 3, 4],
    initial_temperature=293.15,
    source=Step(2, before=0, after=1e6),  # W/m³
    probes={"corner": 0},
    output_every=2,
)
```

- Alle Zeitwerte sind Sekunden, Temperaturwerte Kelvin. Zeitachsen müssen
  mindestens zwei endliche, streng steigende Werte enthalten. Ein Start bei
  einem anderen endlichen Zeitpunkt als null ist zulässig.
- `Step(time,before,after)` ist rechtsstetig. `Ramp(start,end,start_value,end_value)`
  interpoliert linear und hält außerhalb des Rampenintervalls die Endwerte.
  Alle vier Rampenparameter und Abfragezeiten müssen endlich sein; `end > start`
  ist erforderlich. Auch endliche Zeitgrenzen mit überlaufender Differenz werden
  durch gemeinsame Skalierung korrekt interpoliert (#221). End-/Haltewerte
  und konstante Amplituden werden exakt erhalten; die Amplitudenrechnung
  vermeidet überlaufende Zwischendifferenzen bzw. gleichnamige Summen.
  Nichtendliche Eingaben und ungeordnete Grenzen führen zu `ValidationError`.
  Dies erweitert nicht den Solverzeitbereich: Integrationsraster benötigen
  weiterhin endliche darstellbare positive Zeitschritte.
- `TimeSeries(times,values,interpolation="linear",extrapolation="error",unit=...)`
  speichert unveränderliche Messwerte. `previous` wählt rechtsstetige
  stückweise konstante Werte; `hold` muss für Extrapolation explizit gewählt werden.
- `TimeSeries.from_csv(path,value_column=...,unit=...)` liest UTF-8-CSV,
  standardmäßig Komma und Zeitspalte `time_s`. Zeiten sind bereits in s,
  Messwerte bereits in der angegebenen SI-Einheit; keine implizite Umrechnung.
  Duplikate, nichtendliche Werte, fehlende Spalten und unbekannte Modi werden abgelehnt.
- Alle inneren Profilknoten werden automatisch ins Integrationsraster eingefügt.
  Zeitauflösung bleibt Aufgabe des Modells; Ereignispunkte allein ersetzen keinen
  Zeitkonvergenznachweis.
- Quellen, Flüsse und Konvektion verwenden am rechten Intervallende den **linken
  Grenzwert**. Ein Quellensprung bei t=2 heizt deshalb nicht schon das Intervall
  vor t=2. Vorgegebene Temperaturen sind am Endpunkt rechtsstetig.
- `TemperatureBC(nodes, value)` verwendet aktive nullbasierte **Netzknoten**;
  `value` ist ein skalarer Zeitverlauf in K. Unterschiedliche räumliche
  Randtemperaturen können durch getrennte Knotengruppen beschrieben werden.
  Anfangswerte sind ein Skalar oder ein Vektor in sortierter aktiver Knotenreihenfolge.
  Sie müssen bei t0 mit Randwerten auf 1e-10 K übereinstimmen; keine stille Projektion.
- `TransientConvection(coefficient, ambient)` unterstützt Zeitprofile für h und
  Umgebungstemperatur. Fluss- und Konvektionsmappings adressieren äußere
  Triangle3-Zeilen; positiver Fluss entzieht Wärme.
- Räumlich-zeitliche Python-Quellen/Flüsse können `f(x,t)` sein. Ihre Sprünge
  sind explizit in `times` einzutragen. Der linke Grenzwert wird für Callables
  durch den unmittelbar vorherigen darstellbaren Zeitpunkt angenähert.
  Integration wie im stationären Kern:
  Quadraturgrad zwei, affine räumliche Lasten exakt.
- Konsistente Kapazität garantiert keine diskrete Positivität/Monotonie.
  Ideale Temperatursprünge benötigen Verfeinerung bzw. eine endliche Rampe;
  mögliche Überschwinger sind kein zeitlich aufgelöster physikalischer Nachweis.

## Model-Schema 2

`Model.from_file("model.toml").solve()` unterstützt zusätzlich
`physics="transient_heat"`, `schema_version=2`, `units="SI"`.
Schema 1 und der stationäre Bündelvertrag bleiben unverändert. Schema 2 ist
hier ausschließlich der transiente Wärmevertrag; u/v/a und Frequenzen verwenden
[Schema 3](34_PIEZO_DYNAMICS.md), ohne elektrische Massen.

Das vollständige Beispiel steht unter
[examples/transient_heat](../examples/transient_heat/README.md).

- `materials.<name>.preset` wählt einen unveränderlichen Katalogeintrag.
  Eigene Materialien verwenden `conductivity` (Skalar oder 3×3), `density`,
  `specific_heat` und `units={conductivity="W/(m K)",density="kg/m^3",specific_heat="J/(kg K)"}`.
  Preset und Einzelparameter dürfen nicht gemischt werden.
- `regions` ordnet alle Volumensätze eindeutig Materialnamen zu. Vorhandene
  Regionsprüfung erkennt fehlende Zellen und widersprüchliche Überlappungen.
- `time={values=[...],unit="s"}`, `initial_temperature={value=...,unit="K"}`.
- `loads.source={value=...,unit="W/m^3"}`; `value` ist Zahl oder Profil-Tabelle.
- `boundary.temperature` ist eine Tabellenliste mit `nodes` oder `set`,
  `value` und `unit="K"`. `boundary.flux` enthält `set,value,unit="W/m^2"`.
  `boundary.convection` enthält `set,coefficient,ambient_temperature`,
  `unit="W/(m^2 K)"` und `ambient_unit="K"`.
- Profile: `kind="step"` mit `time,before,after`; `kind="ramp"` mit
  `start,end,start_value,end_value`; `kind="series"` mit `times,values` und
  optional `interpolation,extrapolation,unit`.
- `kind="csv"` benötigt `path,value_column,unit`; optional `time_column`,
  `delimiter`, `interpolation`, `extrapolation`. Nur `from_file` lädt CSVs,
  ausschließlich unterhalb des Modellverzeichnisses. Die aufgelösten Werte
  werden in `to_document()` und Provenienz eingebettet; spätere Dateiänderungen
  verändern ein bereits geladenes Modell nicht.
- `probes.<name>` enthält genau `node` oder `point=[x,y,z]` in m. Punkte werden
  mit Tet4-Formfunktionen interpoliert, nicht auf den nächsten Knoten gerundet.
  Außenliegende Punkte werden abgelehnt; auf Elementgrenzen wird das erste
  passende Element gewählt (stetiges Temperaturfeld).
- `solver.method="backward_euler"`, `solver.lumped=false`.
- `output.every` ist ein positiver Schrittindex-Abstand auf dem ergänzten Raster;
  erster und letzter Zustand werden immer gespeichert. `max_history_bytes`
  begrenzt geschätzte numerische Historienarrays (Standard 256 MiB), nicht den
  gesamten Prozessspeicher oder Sparse-LU-Speicher. Messpunkte und Bilanzen
  werden unabhängig vom Feldintervall bei jedem Integrationsschritt erfasst.

## ParaView und Messpunkte

`export_heat_series(path, mesh, result.field_result)` aus
`fem_toolbox.io.heat_series` erzeugt ein **neues** Verzeichnis:

- `heat.pvd`: in ParaView öffnen, dann Play/Zeitschieber verwenden;
- `heat_*.vtu`: Temperaturen, Reaktionen, Residuen, Temperaturgradienten und Wärmefluss;
- `probes.csv`: `time_s` und `<Messpunkt>_K` bei jedem Integrationszeitpunkt;
- `balance.csv`: Schrittweite, Leistungen, Speicheränderung und Residuen;
- `manifest.json`: separater Serien-Schema-Identifier mit Version 1,
  Einheiten, Zeit-/Dateizuordnung, SHA-256 der VTUs, Netzfingerprint und Provenienz.

Inaktive Knoten und nicht belegte Triangle3-Volumenfelder sind nullgefüllt.
Die Masken `active_temperature` und `valid_volume_field` unterscheiden diese
von echten Nullwerten. Keine geglätteten/scheinbar nodalen Wärmeflüsse.
Ergebnisse eines anderen Netzes werden abgelehnt.

Für schrittweisen Export kann `HeatSeriesWriter(path,mesh)` mit
`solve_transient_heat(...,on_frame=writer.write)` verwendet werden. Erst
`writer.finish(result)` markiert den Lauf vollständig und ergänzt CSV/Provenienz.
Ein Abbruch hinterlässt `status="incomplete"` und lesbare bereits geschriebene
Frames. Die Rückgabe hält weiterhin die ausgewählten Frames im Speicher;
dies ist kein unbegrenzter Out-of-core-Solver und kein Restart-Format.

## Nachweise

Red: fehlende Material-/Zeit-/Solver-/Exportmodule und abgewiesener schema-2-
Modelleinstieg wurden vor Implementierung beobachtet. Ein zusätzlicher
Regressionstest reproduzierte den Export auf ein falsches Netz; der neue
Fingerabdruck verhindert ihn.

- Analytische Speichermatrizen, starre Translation, Knotenumordnung und Skalierung.
- Isoliertes Heizen mit Quellensprung, Mehrmaterialspeicherung, bewegte
  Temperaturvorgaben mit freien inneren Knoten, zeitabhängige Flüsse/Robin-Daten.
- Gleichmäßige Abkühlung des regulären Tetraeders: exakte Rate 1/s,
  10/20/40/80 Schritte bis 1 s. Temperaturfehler 0,176638 / 0,090100 /
  0,045512 / 0,022873 K; Ordnungen 0,971 / 0,985 / 0,993.
- Vier Netze n=2,4,6,8 mit unabhängig vorgegebener quadratischer Raumlösung
  und zeitlich linearer Temperatur. Elliptische Anfangsprojektion vermeidet
  einen Startfehler; Backward Euler ist für die lineare Zeitabhängigkeit exakt.
  Beobachtete L2-Ordnungen 2,0 und H1-Ordnungen 1,0.
- Stationärer Abkühlgrenzfall, Speicher-/Randbilanz und freie Residuen.
- Unabhängiges Einlesen von PVD (XML), VTU (meshio) und CSV; Abbruchtest,
  Netzzuteilung und Wiederholbarkeit der TOML-/CSV-Provenienz.

Toleranzen: Zeitordnungen 0,9–1,1; Raumordnungen asymptotisch L2 1,8–2,2,
H1 0,85–1,15. Für analytisch exakte kleine SI-Testfälle 1e-9 W Bilanz- bzw.
2e-11 K Temperaturtoleranz; diese decken LU-/Summationsrundung ab. Der Solver
verwendet freie relative Residuen ≤1e-9 und eine Leistungsgrenze aus 1e-9
mal Leistungsskala plus 100 Maschinen-Epsilon mal algebraischer Rundungsskala.
Diese Nachweise sind Code-Verifikation, keine experimentelle Validierung.


### Ausgeführte lokale Prüfungen (23.09.2026)

- `PYTHONPATH=src .venv/bin/python -m pytest -q`: 336 Tests bestanden.
- Ruff Lint/Format, mypy (44 Quelldateien), `git diff --check`: bestanden.
- Wheel und sdist mit `uv build` erstellt; jeweils in frische Python-3.13-
  Umgebungen installiert (NumPy 2.5.3, SciPy 1.18.1, meshio 5.3.5).
- Transientes TOML/CSV-Beispiel und unverändertes stationäres Wärmebeispiel
  mit `python -I` außerhalb des Repositories aus beiden Installationen ausgeführt.
  Maximales transientes Bilanzdefizit: 1,473e-13 W.
- Der lokale editable-Import wird durch ein macOS-Hidden-Attribut der `.pth`
  Datei übersprungen; Quelltests verwenden deshalb explizit `PYTHONPATH=src`.
  Die isolierten Wheel-/sdist-Importe funktionieren ohne diesen Zusatz.

CI für Python 3.11–3.13 ist um das transiente Artefaktbeispiel ergänzt.
Ein lokaler Nachweis ersetzt weder diese Matrix noch unabhängiges Review.

## M7-Erweiterung: Tet10-Wärme

Der bestehende Backward-Euler-Pfad unterstützt geradkantige Tet10 und
Triangle6, einschließlich Model, quadratischer Sondeninterpolation und
optionalem expliziten HRZ (`solver.lumping = "hrz"`). Der Tet4-Vertrag bleibt
unverändert. Tet10-Zeitserien verwenden Schema 2 mit ausdrücklich benannten
Schwerpunktfeldern. [Vertrag und Nachweise](45_TET10_HEAT_INTEGRATION.md).
