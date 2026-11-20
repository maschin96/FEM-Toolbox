# Fluidnetz, Wandqualität und Rand-Preflight (#85)

Der experimentelle Adapter implementiert den ersten, ebenen Schnitt des
[M12-Vertrags](71_M12_RANS_FSI_CONTRACT.md). `FluidMesh2D` und `FluidDomain`
sind Daten-/Geometrie-APIs ohne Strömungssolver, SA-Transport oder Netzbewegung.
Quad-FV-Operator #226 und RANS #87 folgen separat. Der uniforme MAC-Kern und
seine bestehende FlowProblem-API bleiben unverändert. Hex8-/3D-Adapter bleiben
für die zweite M12-Stufe offen; dieser Schnitt ist keine 3D-Milestone-Freigabe.

## Einstieg

```python
from fem_toolbox.cfd import FluidDomain, FluidMesh2D

mesh = FluidMesh2D.from_file("examples/fluid_mesh/mesh.json")
domain = FluidDomain.from_file("examples/fluid_mesh/domain.json")
print(domain.inlet_mass_flow_kg_s)  # 0.06875 kg/s, vorgegebener Einlass
print(domain.wall_quality.wall_distance_m.min())  # 0.125 m
```

`FluidMesh2D(points, quad4, line2, boundary_sets, volume_sets=None)` ist die
entsprechende direkte API; `from_document(mapping)` akzeptiert das JSON-Schema.
`to_document()` liefert eine unabhängige JSON-kompatible Kopie des Netzeingangs.
`FluidDomain.from_document(mapping, mesh=mesh)` prüft die Domäne ohne I/O.
`from_file` liest Domain-JSON und das relativ dazu referenzierte Mesh-JSON;
absolute Pfade, `..` und aus dem Eingabebaum führende Symlinks werden abgelehnt.
Keine implizite Einheitenskalierung, automatische Reparatur oder Netzgenerierung.

## Mesh-JSON Schema 1

Pflichtfelder: `fluid_mesh_schema_version=1` (Integer, kein Bool),
`length_unit="m"`, `points`, `quad4`, `line2`, `boundary_sets`.
Optional: `volume_sets`, ansonsten leere Auswahlmenge. Unbekannte Felder
werden abgelehnt. Alle Indizes sind nullbasiert; Diagnosen nennen Zeilenindizes.

| Feld | Form / Vertrag |
| --- | --- |
| points | float64[npoint,2], endliche Meter, mindestens vier eindeutige Koordinaten; alle Punkte verwendet |
| quad4 | int64[ncell,4], vier verschiedene Knoten, gegen Uhrzeigersinn, positive konvexe Zelle |
| line2 | int64[nboundary,2], zwei verschiedene Knoten; jede Außenkante genau einmal, keine innere oder fremde Kante |
| boundary_sets | Namen -> nicht leere eindeutige line2-Zeilenindices; Reihenfolge bleibt erhalten |
| volume_sets | Namen -> nicht leere eindeutige quad4-Zeilenindices; reine Auswahlsätze, keine zweite Fluidmaterialzuweisung |

Gleitkomma-/Bool-Indizes, doppelte Zellen/Koordinaten, mehrfach besessene
Kanten, inkonsistente interne Orientierung, Überlappung und Teilkantenkontakt
(T-Knoten) werden diagnostiziert. Das Netz muss über gemeinsame vollständige
Faces zusammenhängen. Löcher mit vollständigem Rand sind möglich; getrennte
Fluidkomponenten und ungenutzte Punkte werden im ersten Adapter abgelehnt.
Benannte Auswahlen dürfen überlappen; physikalische Randzuweisungen müssen
später genau einmal je Randface erfolgen. Line2-Richtung bleibt gespeichert;
die äußere Normale wird aus der CCW-Zelle abgeleitet und nicht aus einer
möglicherweise umgekehrten Eingabezeile erraten.

Maximal 50000 Zellen im Referenzadapter. Überlappungsprüfung verwendet eine
Bounding-Box-Sweepprüfung und Separating Axes; T-Knotenprüfung einen
SciPy-cKDTree. Keine zusätzliche Bibliothek. Diese Referenz ist nicht als
profilierter Großnetzimport optimiert. Nicht darstellbare Geometrie wird
abgelehnt; Bool-Werte werden auch in gemischten numerischen Listen erkannt.

## Geometrie und Netzgates

`cell_areas_m2`, polygonale `cell_centres_m` und `cell_aspect_ratio` besitzen
einen Wert pro Quad4. Aspect Ratio ist max/min Kantenlänge; große Anisotropie
allein wird nicht als Fehler bewertet. Bilineare Jacobians sind an vier Ecken
und vier Gaußpunkten positiv. Zentroiden stammen aus unabhängiger
flächengewichteter Dreieckszerlegung des konvexen Polygons.

`faces` enthält alle inneren und äußeren Kanten. `face_owners` ist der erste
Zelleigentümer, `face_neighbours` der zweite bzw. -1 am Rand.
`face_centres_m`, `face_lengths_m`, `face_normals` verwenden diese Reihenfolge;
Normale zeigt vom Eigentümer nach außen. `boundary_face_indices` ordnet die
Eingabe-line2-Zeilen diesen Faces zu. `boundary_normals` folgt line2.
Diese Abbildung ist die Basis des späteren konservativen FV-Operators.

Geometriegates sind Nichtorthogonalität <=60 Grad und Skewness <=0.25 sowie
positive Eigentümer-Normaldistanz. Winkel: Eigentümer-Nachbar-Verbindung zu
Face-Normale; am Rand Eigentümer-zu-Face-Verbindung. Skewness: Abstand des
Facezentrums zum Schnitt derselben Verbindung mit der Facelinie, dividiert
durch Verbindungslänge. Vollständige Arrays heißen `nonorthogonality_deg` und
`skewness`; maximale Werte über die entsprechend benannten Properties.
Für O(1)-Referenzen sind Geometrie-/Distanzbudgets abs=1e-14 m bzw. m²,
Winkelbudget abs=1e-12 Grad, jeweils rtol=0. Das deckt reine float64-Rundung,
keinen Diskretisierungsfehler oder nachträglich aufgeweitete Testgrenzen.

## Wandmaße und spätere y+-Prüfung

`mesh.wall_quality(["bottom", "top"])` verlangt mindestens eine vorhandene
Wandgruppe und liefert unveränderliche `WallQuality`:

- `boundary_indices`: sortierte eindeutige line2-Wandzeilen.
- `wall_distance_m`: minimale euklidische Entfernung jedes Zellzentrums zu
  allen ausgewählten Wandsegmenten einschließlich Endpunkten. Kein Abstand
  nur zu Wandknoten, entlang Netzlinien oder ausschließlich zur Eigentümerwand.
- `first_cell_distance_m`: Normalabstand der Eigentümerzelle zur Wandface.
- `layer_growth`: größtes Verhältnis aufeinanderfolgender projizierter
  Zellschichtdicken je ausgewählter Wandface, mindestens 1.

Die Schichten werden im strukturierten Quad-Netz entlang jeweils
gegenüberliegender Faces bis zum Rand verfolgt; Dicke ist die Projektion der
Ein-/Austritts-Facemittelpunkte auf die ursprüngliche Wandnormale. Zyklen,
verschwindende Dicken und Wachstum >1.2 werden abgelehnt. Am Grenzwert gilt
lediglich ein float64-Rundungsbudget von 64 eps. Dies ist ein konservativer
Prüfvertrag für diesen Adapter, keine allgemeine Grenzschichtsuche auf
unstrukturierten/gefalteten Netzen. Bei FluidDomain übernehmen genau die
physikalisch als `wall` zugewiesenen Faces diese Rolle.

`walls.y_plus(density_kg_m3=..., viscosity_Pa_s=..., wall_shear_Pa=...)`
verlangt einen endlichen vorzeichenbehafteten skalaren Wandschubwert pro
`boundary_indices`-Eintrag. Es berechnet `d1 sqrt(abs(tau_w)/rho)/(mu/rho)`
und lehnt y+>1 mit Wand-ID und Wert ab. Keine Schubspannung wird geraten;
die Funktion führt keinen Solve aus. Ohne tatsächliche finale Schubwerte
bleibt y+ ausdrücklich ungeprüft. Netzänderungen verlangen eine neue
FluidMesh2D-/WallQuality-Instanz; dies ist keine FSI-Netzbewegung.

## Domain-JSON Schema 1 und Randbedingungen

Pflicht: `fluid_domain_schema_version=1`, `mesh` als relativer Pfad,
`width_m>0`, `fluid`, `boundaries`. Optional: `pressure_reference_cell`.
`fluid` enthält genau positive endliche `density_kg_m3` und `viscosity_Pa_s`
(kg/m³ und Pa s). Die existierende `IncompressibleFluid` wird wiederverwendet.
Width ist zur Extrusion der ebenen Massenströme in kg/s Pflicht, kein Default.

`boundaries` ordnet ausgewählten `boundary_sets`-Namen strikte Tabellen zu:

| type | Weitere Schlüssel / Werte |
| --- | --- |
| wall | optional `fsi` (Bool, Default false); stationäre no-slip-Wand, SA-Wert null |
| inlet | `velocity_m_s`: (2,) oder (nselection,2), `nu_tilde_m2_s`: Skalar oder (nselection,), endlich und nichtnegativ |
| pressure_outlet | `pressure_Pa`: endlicher Skalar; festgelegter Druckauslassvertrag aus Doc/71 |
| symmetry | keine zusätzlichen Schlüssel; Normalgeschwindigkeit null, tangentialer/SA-Normalgradient null |
| periodic | `partner`: anderer Gruppenname, `translation_m`: endlicher nicht null (2,)-Vektor, `pressure_drop_Pa`: endlicher Skalar |

Genau eine physikalische BC je line2-Face; FSI markiert ausschließlich eine
Wand. Fehlende Zuordnung, Überdeckung durch zwei BCs, unbekannte Typen oder
zusätzliche Felder werden vor dem Solve abgelehnt. Mindestens eine Wand ist
für diesen wandaufgelösten Adapter Pflicht. Einlassgeschwindigkeit darf an
keiner Face nach außen zeigen und benötigt insgesamt positiven Zufluss.
`inlet_mass_flow_kg_s=-rho width sum(U·n_face length_face)` ist der vorgegebene
Einlassfluss; das ist keine ausgerechnete Ein-/Auslassbilanz.

Genau eine Druckauslassgruppe ohne zusätzlichen Gauge oder ein gültiger
Zell-Gauge ohne Druckauslass; mehrere Faces derselben Auslassgruppe sind
zulässig. Einströmung ohne Druckauslass ist inkompatibel, da sonst nur
undurchlässige oder paarweise periodische Grenzen verbleiben.

Periodische Paarung ist reziprok, mit entgegengesetzter Translation und
Druckdifferenz: `p(partner)=p(source)-pressure_drop_Pa`. Übersetzte Segmente
müssen geometrisch kongruent sein, ihre **äußeren** Normalen entgegengesetzt;
line2-Endpunktreihenfolge darf umgekehrt sein. Geometrisches Matchbudget
1e-10 mal größte Gebietsausdehnung, Normalenbudget 1e-10. Eindeutige Paarung
wird unabhängig von Gruppenreihenfolge unter `periodic_partners` als
line2-ID -> Partner-ID bereitgestellt. Physische Periodizitäts-/Druckoperatoren
gehören zu #226, nicht zum Import.

Alle gespeicherten Arrays besitzen kopierte immutable Bytepuffer, Gruppen/BCs
unveränderliche Mapping-Container. Schema-/Geometriefehler sind konkrete
`ModelSchemaError`/`ValidationError` mit Pfad und betroffener Entität;
Dateifehler werden als `ModelSchemaError` verknüpft. Keine stillen Reparaturen.

## Beispiel und Nachweise

[Vorbereiteter Kanalfall](../examples/fluid_mesh/README.md): 16 Quad-Zellen,
16 Randfaces, L=2 m, H=1 m, b=0.1 m, rho=mu=1 in SI. Vier Einlasssamples von
4y(1-y) ergeben den diskreten vorgegebenen Fluss 0.06875 kg/s. Der analytische
kontinuierliche Fluss ist 1/15 kg/s; die Sampledifferenz wird nicht versteckt
oder auf einen Zielwert skaliert. Wandabstand mindestens 0.125 m, Schichtwachstum
1. Kein Strömungsergebnis und kein y+-Nachweis ohne tatsächliche Schubwerte.

Tests prüfen unabhängige Rechteck-/Schergeometrie, Normale/Orientierung,
Segmentendpunktabstände, Anisotropie, Schichtwachstum, y+, T-Knoten,
Überlappung, ganze/falsche Indizes, vollständige Gruppen/BCs, periodische
Paarung, Gauge, SI-Daten, immutable Kopien und Dateipfade. Beispielregression
startet den Dateilauf aus einem fremden Arbeitsverzeichnis. Keine neuen
Solveroperatoren: neue Raum-/Solverkonvergenz ist nicht anwendbar; vorhandene
CFD-Verifikation bleibt erforderlich. TDD- und Qualitätsresultate werden
am aktuellen PR-Stand dokumentiert, nicht vorweggenommen.
