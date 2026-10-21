# M12 — konservative Schnittstelle und einseitige Strukturantwort (#89)

Stand 21.10.2026, experimentelle API in `fem_toolbox.coupling`.
Dies schließt keine bidirektionale FSI- oder M12-Freigabe ein.
Der [M12-Vertrag](71_M12_RANS_FSI_CONTRACT.md) bleibt maßgeblich.

## Geometrie und Überdeckung

`build_interface_map` nimmt ebene Fluidsegmente `[nsegment,2,2]` in m,
äußere Fluid-Einheitsnormalen `[nsegment,2]`, die explizite Breite b in m,
einen Tet4-`Mesh` mit benanntem Triangle3-Facettensatz und Lref>0 entgegen.
Die Fläche entsteht durch Extrusion über z von `span_origin_m` bis
`span_origin_m+b`. `build_domain_interface` übernimmt Segmente, Normalen und
Breite aus einer `FluidDomain`; die benannte Wand muss ausdrücklich `fsi=true`
tragen. Struktur-Normalen stammen vom eindeutigen Tet4-Eigentümer und müssen
den Fluidnormalen entgegengesetzt sein. Die Reihenfolge gespeicherter
Segmentendpunkte oder Triangle3-Knoten definiert keine physikalische Normale.

Jedes strukturelle Dreieck wird in die lokale Segment-/Spannweitenebene
projiziert und mit dem Fluidrechteck geschnitten. Eigener konvexer Polygonclip
und eine Fächertriangulierung erzeugen die Überdeckung; es gibt keine neue
Laufzeitabhängigkeit. Je Teildreieck integriert die symmetrische Dreipunktregel
Polynome vom Grad zwei exakt. Damit werden affine Traktion mal P1-Formfunktion
und deren Momente exakt integriert, bis auf Gleitkommafehler.

Die Zuordnungstoleranz beträgt 1e-10 Lref, das relative Flächengate 1e-10.
Alle Fluidflächen und alle ausgewählten Strukturflächen müssen vollständig
überdeckt sein. Doppelte Fluidintervalle und überlappende Strukturpolygone
werden zusätzlich direkt geprüft: Gesamtflächen allein würden symmetrisch
verteilte Lücken und Doppelungen unter Umständen nicht erkennen.
Keine Extrapolation, keine Ergänzung von Lücken und kein nächster Nachbar.
Falsche Breite, Orientierung, fehlende Auswahl und inkompatible Flächen sind
konkrete Fehler. Standardbudget: 200000 Quadraturpunkte, jeweils höchstens
2000000 geometrische Kandidaten-/Polygonpaarprüfungen; Budgetüberschreitung
ist kein akzeptiertes Ergebnis. Tet10/Triangle6 sind hier nicht unterstützt.

## Gleicher Operator für Weg und Kraft

Die skalare dünn besetzte Matrix W enthält die baryzentrischen P1-Gewichte.
Für Vektoren gilt H=kron(W,I3). Konstanten und Koordinaten werden beim Aufbau
explizit auf Reproduktion geprüft. Öffentliche Felder sind unveränderlich;
`interpolation_matrix` liefert eine unabhängige Kopie von H.

```text
u_q = H u_s                         [m]
F_s = H^T F_q                       [N]
F_q = t_q A_q                       [N]
u_s^T F_s = u_q^T F_q               [J]
sum F_s = sum F_q                   [N]
sum x_s cross F_s = sum x_q cross F_q [N m]
```

`interpolate` und `transfer_forces` erlauben beliebige virtuelle 3D-Vektoren
zur Arbeitsprüfung. `audit` berichtet rohe Kraft-, Moment- und Arbeitsfehler
und prüft das unveränderte relative Gate 1e-10 mit absoluten Floors
1e-12 N, N m und J. Die Bezugsgröße berücksichtigt absolute Einzelbeiträge,
damit sich aufhebende Kräfte nicht durch eine Nullresultierende normalisiert
werden. Arbeitsgleichheit prüft das Mapping, keinen dynamischen Energieaustausch.

`loads_from_fluid_traction` erwartet sigma_f*n_f in Pa, `[nsegment,2]`,
stückweise konstant wie beim FV-Wandschub. Das Strukturvorzeichen ist
`t_s=-t_f-p_external*n_f`. Fluid- und externer Druckanteil bleiben separat
berichtet. Ein Druckdatumwechsel erhält Lasten nur bei gleichem Wechsel von
p_external. Der Aufrufer wählt die Fluidtraktionen in derselben Reihenfolge
wie die Segmente; die Solver-Wandreihe folgt `domain.wall_quality.boundary_indices`.
Das Ergebnis enthält Strukturknotenkräfte in N, keine Druckwerte als Punktkräfte.

## Ebene Bewegung und bestehende Mechanik

`segment_endpoint_displacements` liefert ebene absolute Wege
`[nsegment,2 endpoints,2]` in m. Es prüft z-Verschiebungen und
Spannweitenvariation gegen 1e-8 Lref sowie die Darstellung der strukturellen
Spur durch lineare Fluidsegmente. Für jeden Überdeckungsteil werden aus den
drei Gaußwerten die affinen Abweichungen an allen Ecken rekonstruiert;
so können schmale Spitzen zwischen Quadraturpunkten nicht verborgen bleiben.
Nicht aufgelöste P1-Knicke werden abgelehnt;
eine beliebig gekrümmte verformte Schnittstelle wird nicht behauptet.
`interpolate` bleibt für allgemeine virtuelle Arbeitsproben uneingeschränkt.
Harmonische innere Netzbewegung ist das getrennte Arbeitspaket #90.

`solve_one_way(mapping, model, traction_Pa, external_pressure_Pa=...)` benutzt
den vorhandenen statischen `solve_elasticity` und `Model.elasticity_inputs()`.
Vorhandene Körperlasten, Traktionen und Dirichletwerte bleiben erhalten;
übertragene Knotenkräfte werden hinzugefügt. Koordinaten und DOF-Reihenfolge
müssen zur Schnittstelle passen. Kein neues Strukturelement und keine Änderung
der stabilen Model-Schemata. Die bestehende Gleichgewichtsprüfung bleibt aktiv;
zusätzlich werden max|u|/Lref<=0.02 und max||grad u||_F<=0.02 geprüft.

## Unabhängige Prüfungen

Versionierte Fälle stehen in `tests/unit/test_interface_map.py`:

- Nicht passende Netze: Fluidsegmente mit Längen 0.7 und 2.3 m über einer
  Fläche von 3 mal 4 m, Struktur-Triangulierung aus 48 Tet4.
  Fluidtraktion (-3,-7) Pa und p_external=5 Pa ergeben die Nettokraft
  (36,24,0) N und das Moment um den Ursprung (-48,72,0) N m.
- Affine Traktion t=(x,2+z,0) Pa ergibt unabhängig integrierte Kraft
  (18,48,0) N und Moment (-112,36,54) N m. Beliebige virtuelle Wege
  mit festem Zufallsseed prüfen die adjungierte Arbeitsgleichheit.
  Affine Vektorfelder werden exakt reproduziert.
- Rotation der Geometrie und Umkehr der Segmentreihenfolge erhalten die
  entsprechend gedrehten physikalischen Kräfte. Eine herausgegebene
  Matrixkopie kann die ursprüngliche Abbildung nicht verändern.
- Der einseitige Fall benutzt auch einen echten Fluid-Solve: geschlossener
  3 mal 1 m Kanal mit Breite 4 m, Fluiddruck 7 Pa, externer Druck 5 Pa,
  Struktur darüber mit E=1000 Pa und nu_s=0. Seitlich geführte Struktur,
  obere y-Knoten fest. Erwartet werden u_y=0.002(2-y) m, Reaktion
  (0,-24,0) N sowie verschwindende Kraft- und Momentbilanz.
  Eine zweite Variante prüft die Lastabbildung ohne Fluidsolver getrennt.
- Lücken, falsche Breite, falsche Normale, z-/spannweitenabhängige Bewegung,
  symmetrische Lücken mit Doppelungen, doppelte strukturelle Körper, ungültige
  Breite und Ressourcenüberschreitung sind Negativfälle.

TDD: Zunächst fehlender Kopplungsimport; danach fehlende einseitige API und
ein reproduzierter Nichtfehler bei symmetrischen Lücken/Doppelungen.
Explizite Intervallüberlappungsprüfung behebt diesen Fehler. Alle 20 gezielten
Tests bestanden vor dem zusätzlichen Spitzentest. Dieser reproduzierte eine
zwischen Proben verborgene Spannweitenspitze; die anschließende Eckprüfung
behebt den Fehler. Alle 21 gezielten Tests bestehen.
Unabhängiges menschliches Fachreview liegt nicht vor.
Netzbewegung, bidirektionaler Fixpunkt, FSI-Netzfolgen und deklarative
Release-Workflows bleiben separate verpflichtende Arbeitspakete.
