# Stationäre Wärmeleitung (Issue #13)

Status: Code- und Lösungsverifikation für die beschriebenen Referenzfälle;
keine experimentelle Validierung. Entwicklungsstand nach v0.2.0, noch kein
vollständiges M2-Release. Kein transienter oder thermoelastischer Solver.

## Gleichung und Vorzeichen

```text
-div(k grad T) = Q             im Volumen
q = -k grad T                 Wärmeflussvektor
q dot n = q_bar               vorgeschriebener äußerer Wärmestrom
q dot n = h (T - T_inf)        Konvektion
```

Positive Quellen führen Wärme zu; positives `q_bar` führt Wärme ab. Die schwache
Form lautet mit Testfunktion w:

```text
integral grad(w) dot k grad(T) + integral_Robin w h T
= integral w Q - integral_Neumann w q_bar + integral_Robin w h T_inf.
```

Die Randdefinition mit äußerer Normaler erzwingt das Minuszeichen vor `q_bar`.
[ADR-008](07_RISKS_AND_DECISIONS.md) korrigiert den früheren Widerspruch zwischen
starker und schwacher Form in [Dokument 02](02_MATHEMATICAL_MODEL.md).

## API und Einheiten

`ThermalMaterial(k)` aus `fem_toolbox.materials` enthält die unveränderliche
Leitfähigkeitsmatrix `(3,3)` in W/(m K), bereits in globalen Achsen.
`ThermalMaterial.isotropic(k)` akzeptiert einen positiven Skalar. Endlichkeit,
Symmetrie (relativ `1e-12` zum größten Eintrag) und positive Definitheit werden
geprüft. Kleine akzeptierte Asymmetrie wird gemittelt; Eigenwertverhältnisse
kleiner/gleich `1e-12` werden abgelehnt. Es gibt keine implizite Materialrotation
oder Temperaturabhängigkeit.

`fem_toolbox.elements.heat.conductance_matrix` liefert exakt
`V G.T k G` in W/K, Größe `(4,4)`, in der Eingabeknotenreihenfolge. Der konstante
Temperaturmodus liegt im Nullraum. Beide Elementorientierungen sind zulässig.

`thermal_dofs(mesh)` und `solve_heat(...)` liegen in `fem_toolbox.solvers.heat`.
Nur aktive Tet4-Knoten erhalten einen DOF, aufsteigend nach internem Knotenindex.
Das Feld heißt `T`; `Dirichlet(dofs.indices("T", nodes), values)` beschreibt
Temperaturvorgaben. Werte müssen endlich sein und dieselbe Temperaturreferenz
wie die Umgebung benutzen, normalerweise absolute Kelvin. Die lineare Algebra
lässt auch konsistent referenzierte Temperaturdifferenzen zu; sie prüft keine
thermodynamische Positivität. Ein Material gilt homogen; eine Materialsequenz
ordnet jeder Tet4-Zeile explizit ein Material zu.

Optionale Eingaben:

- `source`: Skalar in W/m³ oder Funktion `x -> Q(x)`;
- `outward_fluxes`: Zuordnung von äußeren Triangle3-Zeilen zu W/m²-Skalaren oder
  Funktionen `x -> q_bar(x)`;
- `convection`: Zuordnung von äußeren Triangle3-Zeilen zu
  `Convection(coefficient=h, ambient=T_inf)` aus `fem_toolbox.physics.heat`.
  h ist je Fläche konstant, endlich und nichtnegativ, in W/(m² K); die
  Umgebungstemperatur ist ein Skalar oder eine Funktion.

Funktionen erhalten Punkte `(3,)` in Metern und liefern reelle endliche Skalare.
Konstante Lasten werden exakt integriert. Für Funktionen verwenden die
Lastoperatoren Quadraturgrad 2: vier Punkte im Tetraeder, drei im Dreieck.
Damit sind affine Quellen, Flüsse und Umgebungstemperaturen exakt integriert.
Höhere Ortsabhängigkeit verlangt einen zusätzlichen Quadraturkonvergenznachweis.
Die Robin-Matrix ist die exakte konsistente Dreiecksmatrix
`h A/12 * [[2,1,1],[1,2,1],[1,1,2]]`, kein diagonales Lumping.

Vorgegebener Fluss und Konvektion auf derselben Fläche addieren sich explizit.
Nicht angegebene Außenflächen sind adiabatisch. Interne oder nicht zum
Tet4-Rand gehörende Dreiecke werden als thermische Randflächen abgelehnt.
Dreiecksumordnung ändert das Vorzeichen des vorgegebenen äußeren Flusses nicht.

## Eindeutigkeit und Ergebnis

Jede knotenzusammenhängende Komponente benötigt mindestens einen vorgegebenen
Temperaturwert oder Konvektion mit h>0. Andernfalls meldet der Solver vor der
Faktorisierung `ThermalNullspaceError` mit Komponentenindex, `net_power` (W),
`tolerance` (W) und `compatible`:

- kompatibel: die Quellen-/Flussbilanz verschwindet, die additive
  Temperaturkonstante bleibt frei;
- inkompatibel: ohne Wärmeabfuhr/-zufuhr kann kein stationäres Feld existieren.

Es wird keine versteckte Referenztemperatur gesetzt. Die erste unreferenzierte
Komponente wird gemeldet; globale Aufhebung zwischen getrennten Körpern genügt
nicht. Kompatibilitätstoleranz: `1e-10 * sum(abs(component_rhs))`; ein lastfreier
Fall wird exakt als kompatibel erkannt. h=0 beseitigt keinen Nullraum.

Die symmetrieerhaltende Dirichlet-Elimination erhält das Originalsystem. Der
freie Block wird diagonal skaliert und mit SciPy SuperLU gelöst. Pivots
kleiner/gleich `1e-12` werden abgelehnt; dies ist eine Schutzprüfung, keine
vollständige Konditionszahlanalyse. Unzulässig schwache Konvektion oder extrem
schlecht konditionierte Geometrien können daher eine Fehlermeldung erzeugen.

`HeatResult` enthält unveränderliche Arrays `node_ids`, `temperature` in K,
`gradient` in K/m, `heat_flux` in W/m² sowie `reactions` und `residual` in W.
Temperaturen und nodale Leistungen stehen in aktiver Knotenreihenfolge,
Gradienten und Flüsse roh je Tet4 in Eingabereihenfolge. Das Residuum ist
`K T - f`. Reaktionen sind an freien DOFs null; positiv bedeutet am
Dirichlet-Knoten durch den Thermostaten zugeführte Wärmeleistung.

Die globalen Leistungsgrößen sind alle in W:

```text
power_balance = source_power + dirichlet_inflow
              - prescribed_outflow - convection_outflow.
```

`convection_outflow` integriert `h(T-T_inf)` mit derselben konsistenten
Randdiskretisierung. Positive Werte führen Wärme ab. Die rohe, elementweise
Flussrekonstruktion ist auf allgemeinen P1-Netzen nicht lokal flusskontinuierlich;
die globale Bilanz basiert deshalb auf integrierten Lasten, Robin-Flüssen und
Dirichlet-Reaktionen, nicht auf geglätteten Flusswerten.

Die freie Residualnorm wird durch `|| |K| |T| + |f| ||_2` geteilt; Grenzwert
`1e-9`. Zusätzlich muss die Leistungsbilanz innerhalb `1e-9` der Summe absoluter
nodaler Quell-, Rand- und Reaktionsleistungen liegen, ergänzt um das
Rundungsbudget `100 eps * sum(|K| |T| + |f|)`. Das berücksichtigt die Auslöschung
bei großen konstanten Temperaturanteilen. Die Bilanz ist zusätzlich in den
Referenztests gegen unabhängige analytische Leistungen geprüft.

`RunProvenance` hasht Koordinaten, Konnektivität, Materialien, integrierte
Quell-/Fluss-/Robin-Beiträge und aufgelöste Dirichlet-Werte als kanonisches JSON
(`heat-discrete-v1`). Solver, Vorzeichen, Quadratur und numerische Toleranzen
werden aufgezeichnet. Callable-Quelltexte werden nicht archiviert; die
vollständige Originalmodellarchivierung folgt mit #27.

## Nachweise und Toleranzen

- Handintegrierte Tet4-Leitmatrix, Rang 3, konstanter Nullmodus, Materialfehler,
  lokale Knotenumordnung und Rotation eines anisotropen Materials samt Geometrie.
- Exakte Integrale von konstanten/affinen Volumen- und Flächenlasten sowie
  Robin-Matrix und Umgebungstemperaturlast.
- Anisotroper affiner Patch `T=300+x+2y-3z` auf einem innen verzerrten Würfel:
  konstante Gradienten/Flüsse und globale Bilanz.
- Einheitswürfel, k=2, T(0)=400, rechter äußerer Fluss 20 W/m²:
  `T=400-10x`, linke Wärmezufuhr 20 W, rechte Abfuhr 20 W.
- Konvektionsquader: `q=(T_left-T_inf)/(L/k+1/h)`;
  Einheitswürfel k=2, h=4, T_left=400, T_inf=300 liefert q=400/3 W/m².
  Quader 2×3×4 m mit denselben Parametern liefert q=80 W/m², Leistung 960 W
  und `T=400-40x`. Die nicht erwärmten Seitenflächen sind adiabatisch.
- Zwei gleich dicke Schichten im Einheitswürfel, k1=2, k2=4, Endtemperaturen
  400/300 K: q=800/3 W/m², stetige Temperatur und stückweise konstante Gradienten.
- Ausschließlich Konvektion ohne Quelle: konstante Umgebungstemperatur als
  eindeutige Lösung; reine Neumann-Fälle mit null, kompatibler und inkompatibler
  Last sowie getrennten Komponenten; h=0 als expliziter Nullraumfall.

Rundungstoleranzen bei diesen SI-Referenzparametern: bis `2e-15` für kleine
Elementmatrizen, `3e-12 K` für Temperaturen um 300–400 K, `5e-11 W/m²` für
rekonstruierte Flüsse und `1e-10 W` für Benchmarkleistungen. Rotierte Geometrien
verwenden zusätzlich `rtol=2e-13, atol=2e-14` für die Leitmatrix. Die Grenzen
berücksichtigen Differenzenbildung aus absoluten Temperaturen und Akkumulation,
ohne Vorzeichen- oder Materialfehler zu verdecken.

## Netzkonvergenz

Manufactured Solution auf dem Einheitswürfel:
`T=300+x(1-x)+y(1-y)/2+z(1-z)/4`, Leitfähigkeit
`[[2,.3,.2],[.3,3,.1],[.2,.1,4]] W/(m K)`, Q=9 W/m³, exakte Temperaturen auf
dem gesamten Rand. Die unabhängigen Fehlerintegrale verwenden eine
Duffy-transformierte 4×4×4-Gaußregel; der quadrierte quadratische
Temperaturfehler wird dabei exakt integriert.

| Unterteilungen je Achse | Tet4 | L2-Fehler (K m^(3/2)) | Energienormfehler (sqrt(W K)) |
|---|---:|---:|---:|
| 2 | 48 | 0.07597742720 | 0.50000000000 |
| 4 | 384 | 0.01899435680 | 0.25000000000 |
| 6 | 1296 | 0.00844193636 | 0.16666666667 |
| 8 | 3072 | 0.00474858920 | 0.12500000000 |

Energienorm: `sqrt(integral grad(error).T k grad(error))`; sie bezeichnet hier
eine mathematische Fehlernorm, keine gespeicherte Wärmeenergie. Alle drei
Übergänge ergeben Ordnung 2 bzw. 1. Auf den letzten beiden Übergängen gelten
die asymptotischen Prüffenster [1.8,2.2] bzw. [0.85,1.15]. Die globale Wärmeabfuhr
beträgt 9 W; wegen Summation über bis zu 3072 Elemente gilt `atol=1e-9 W`.

## Reproduktion

Nach Installation der Entwicklungsabhängigkeiten im Repository-Root:

```sh
python -m pytest tests/unit/test_heat.py tests/patch/test_heat_patch.py
python -m pytest tests/convergence/test_heat_convergence.py -s
python examples/steady_heat.py
```

Das Beispiel schreibt `steady_heat.vtu`; Temperatur und Reaktionsleistung sind
Punktfelder, Temperaturgradient und Wärmefluss rohe Volumenzellfelder.
Triangle3-Werte sind Null-Platzhalter. Bei unbenutzten Mesh-Punkten müssen aktive
Ergebniswerte vor dem Export über `node_ids` zurückabgebildet werden. Die kurze
Netzstudie läuft ohne Ausnahmen in der regulären Testsuite und CI-Matrix.
