# Geradkantige Tet10-Geometrie und Quadratur (#35b)

M7 ergänzt quadratische Lagrange-Interpolation für affine Tetraeder und ihre
planaren Dreiecksflächen. Die [Knotenordnung aus #35a](40_TET10_TOPOLOGY.md)
bleibt unverändert. Feldoperatoren, globale Tet10-Solver und Model-Integration
sind weiterhin Gegenstand von #35c/d.

## Formfunktionen und Ableitungen

Für baryzentrische Eckfunktionen L_i gilt:

```text
Ecke i:    N_i  = L_i (2 L_i - 1)
Kante ij:  N_ij = 4 L_i L_j

grad N_i  = (4 L_i - 1) grad L_i
grad N_ij = 4 (L_i grad L_j + L_j grad L_i)
```

Dies ist die quadratische nodale Lagrange-Basis, siehe
[DefElement: Degree 2 Lagrange on a tetrahedron](https://defelement.org/elements/examples/tetrahedron-lagrange-equispaced-2.html).
Die dortige lokale Nummerierung ist nicht der Toolbox-Vertrag: hier gelten
Tet10-Kanten 01,12,20,03,13,23 bzw. Triangle6-Kanten 01,12,20.

`Tet10Geometry(coordinates, entity_id=None, degeneracy_tolerance=1e-12)`
erwartet ein (10,3)-Array in Metern. `shape_values(point)` liefert zehn Werte,
`shape_gradients(point)` ein (10,3)-Array in 1/m und `gradient_matrix(point)`
dessen Transponierte. Die Auswertung bezieht sich immer auf die Eingabereihenfolge.
Im Gegensatz zu Tet4 hängen die Gradienten vom Auswertepunkt ab.

Die vier Eckpunkte definieren eine affine Geometrie. Mittige Kantenknoten werden
mit der relativen Kantentoleranz 1e-10 aus #35a geprüft; gekrümmte Elemente werden
abgelehnt. Die dimensionlose Tet4-Entartungsprüfung und Qualitätsmaße gelten auch
hier. `signed_volume` erhält das Eingabevorzeichen, `volume` ist positiv.
`oriented_node_order`, `oriented_coordinates` und `oriented_connectivity` führen
alle zehn Knoten konsistent in eine positive Eckorientierung über. Eingaben
werden kopiert, veröffentlichte Koordinaten und Ergebnisarrays sind unveränderlich.

`faces` liefert vier `Tet10Face`-Objekte mit sechs nach außen orientierten
lokalen Knoten, Gegenüber-Ecke, Fläche, Schwerpunkt und äußerer Einheitsnormaler.
`Triangle6Geometry` erwartet (6,3)-Koordinaten und verwendet die Eingabeorientierung
für die Normale. Ihre Gradienten sind tangentiale Vektoren im dreidimensionalen
Raum. QR-Zerlegung berechnet die affine Rückabbildung ohne quadrierte Kondition
über Normalengleichungen. Flächen mit doppelter Fläche / maximaler Kantenlänge²
kleiner oder gleich 1e-12 werden abgelehnt.

Auswertung außerhalb des Tetraeders bzw. innerhalb der Dreiecksebene ist
polynomiale Extrapolation. Triangle6 lehnt Punkte mit Ebenenabstand größer als
1e-10 der maximalen Kantenlänge ab; es gibt keine stille Projektion beliebiger
Punkte auf die Fläche. Nicht endliche, komplexe oder falsch geformte Eingaben
werden ausdrücklich abgelehnt.

## Quadraturvertrag

`tetrahedron_quadrature()` und `triangle_quadrature()` liefern `QuadratureRule`
mit `points`, `weights` und garantiertem Gesamtpolynomgrad `degree`. Die
Fabrikfunktionen liefern unveränderliche Arrays. Referenzpunkte verwenden
kartesische Simplexkoordinaten, keine vollständigen baryzentrischen Koordinaten.
Gewichte summieren sich auf 1/6 im Referenztetraeder und 1/2 im Referenzdreieck.

Die Referenz nutzt positive Gauss-Legendre-Produkte auf [0,1] und die Abbildungen

```text
Tetraeder: (a,b,c) -> (a, (1-a)b, (1-a)(1-b)c), J = (1-a)²(1-b)
Dreieck:   (a,b)   -> (a, (1-a)b),              J = (1-a)
```

Ein n-Punkt-Gauss-Verfahren integriert eindimensionale Polynome bis Grad 2n-1
exakt ([NIST DLMF §3.5](https://dlmf.nist.gov/3.5)). Nach Einsetzen eines
Simplexmonoms vom Gesamtgrad d hat der Tetraederintegrand höchstens die Grade
d+2, d+1 und d. Vier Punkte je Richtung ergeben daher 64 positive Punkte mit
Exaktheit bis Gesamtgrad 5. Im Dreieck genügen drei Punkte je Richtung für neun
positive Punkte und Gesamtgrad 4. Insbesondere sind Produkte zweier quadratischer
Formfunktionen erfasst. Dies ist eine nachvollziehbare Referenzregel, keine
Behauptung einer minimalen Punktzahl oder optimierten Laufzeit.

Die jeweilige Methode `geometry.quadrature()` bildet Punkte in physische
Meterkoordinaten ab und multipliziert die Gewichte mit 6V bzw. 2A. Physische
Gewichte haben somit die Einheit m³ bzw. m² und positive Summen V bzw. A,
unabhängig von der Eingabeorientierung. Eine krummlinige Jacobiauswertung oder
Integration nichtpolynomialer Materialdaten ist damit nicht freigegeben.

## Nachweise

- Kronecker-Eigenschaft an allen Knoten, Partition der Eins und Ableitungssumme.
- Reproduktion aller kartesischen Monome bis Grad 2 mit unabhängigen analytischen
  Ableitungen; Triangle6 gegen den tangentialen Projektor I - n nᵀ.
- Quadratur jedes Referenzmonoms bis zum angegebenen Grad gegen
  `Integral x^a y^b z^c = a! b! c! / (a+b+c+3)!` bzw. die zweidimensionale Form.
- Alle 24 Eckpermutationen, positive Umorientierung, Flächenspur und geschlossene
  Normalenbilanz. Einheitenprüfung mit Skalen 1e-6, 1 und 1e6.
- Vier Würfelnetze mit n=1,2,3,4 Unterteilungen: Interpolation von
  `f=x³+y³/2+z³/4`, unabhängig mit fünf Gauss-Punkten je Richtung überintegriert.

| n | L2-Fehler | H1-Seminormfehler |
|---|---:|---:|
| 1 | 0.03952847 | 0.25617377 |
| 2 | 0.00494106 | 0.06404344 |
| 3 | 0.00146402 | 0.02846375 |
| 4 | 0.00061763 | 0.01601086 |

Alle drei Verfeinerungsraten sind L2=3 und H1=2. Dieser Nachweis betrifft die
Interpolation eines glatten Felds, keinen gelösten PDE-Benchmark. Algebraische
Identitäten auf gut konditionierten Elementen verwenden 2–3e-14 absolute bzw.
relative Toleranz; Monomquadratur zusätzlich 1e-16 absolut für kleine Momente.
Der Ratenvergleich erlaubt 1e-9 für akkumulierte Auslöschung. Diese Grenzen
werden nicht als Genauigkeitsgarantie für beliebig schlecht konditionierte
Geometrien verstanden. Physikalische Patch-, Energie-, Bilanz- und Solver-
Konvergenznachweise bleiben Teil von #35c/d. Keine experimentelle Validierung.
