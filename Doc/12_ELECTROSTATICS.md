# Elektrostatik auf Tet4 (Issue #14)

Status: Code- und Lösungsverifikation für die folgenden Referenzfälle;
keine experimentelle Validierung. Entwicklungsstand nach v0.2.0. Keine
piezoelektrische Kopplung, Floating Electrodes oder vollständige Model-API.

## Gleichungen und Ladungsvorzeichen

Es gelten `E=-grad(phi)`, `D=kappa E` und `div(D)=rho`. Die projektweit
festgelegte schwache Form bleibt unverändert:

```text
integral grad(psi) dot D = integral_boundary psi omega - integral_volume psi rho
omega = D dot n_out
A phi = f,   A = -K_phi,   K_phi = V G.T kappa G
f = surface_charge_load - integral N rho.
```

`dielectric_matrix` liefert den **positiven** Baustein `K_phi` in F. Der Solver
assembliert dessen negatives Vorzeichen. Nur für die Skalierung/Lösung des
freien Blocks werden beide Gleichungsseiten gemeinsam mit -1 multipliziert;
Reaktionen und Residuen stammen weiterhin aus der ursprünglichen negativen
Gleichung. Das hält die spätere symmetrische Piezo-Blockkonvention konsistent.

`surface_charges` gibt **omega=D.n_out** am Dielektrikum in C/m² vor. Das ist
nicht die freie Ladung einer angrenzenden ideal leitenden Elektrode: deren
Ladung hat bei feldfreiem Leiterinneren das umgekehrte Vorzeichen. Auf
Potentialrändern gilt deshalb `electrode_charges = -reactions`.
Der Einflächen-Kondensatortest sichert diesen Unterschied unabhängig ab.

## API und Grenzen

- `DielectricMaterial(permittivity)` aus `fem_toolbox.materials` enthält eine
  unveränderliche symmetrische positiv definite Matrix `(3,3)` in **absoluten
  F/m**, in globalen Achsen. `isotropic(epsilon)` akzeptiert einen positiven
  endlichen Skalar. Relative Permittivität muss vorher ausdrücklich in F/m
  umgerechnet werden; es gibt keine implizite Multiplikation mit epsilon_0.
- Symmetrietoleranz: `1e-12` relativ zum größten Eintrag, akzeptierte
  Rundungsasymmetrie wird gemittelt. Eigenwertverhältnisse kleiner/gleich `1e-12`
  werden abgelehnt. Keine Materialrotation oder Feldabhängigkeit wird impliziert.
- `fem_toolbox.elements.electrostatics.dielectric_matrix` liefert exakt
  integriertes `V G.T kappa G`, `(4,4)`, in ursprünglicher Knotenreihenfolge.
- `electric_dofs(mesh)` und `solve_electrostatics` liegen in
  `fem_toolbox.solvers.electrostatics`. Das Feld `phi` trägt einen DOF je aktivem
  Tet4-Knoten in aufsteigender interner Knotenreihenfolge; unbenutzte Punkte
  erzeugen keine DOFs. Vorgaben: `Dirichlet(dofs.indices("phi", nodes), values)`
  mit Werten in V.
- Ein Material gilt homogen; eine Materialsequenz weist jeder Tet4-Zeile ein
  Material zu. `volume_charge` ist rho in C/m³, `surface_charges` eine Zuordnung
  von äußeren Triangle3-Zeilen zu omega in C/m². Nicht angegebene Ränder haben
  D.n=0. Interne/nicht zum Tet4-Rand gehörende Flächen werden abgelehnt.
- Dichten sind endliche Skalare oder Funktionen `x -> scalar`, x `(3,)` in m.
  Grad-2-Quadratur (vier Punkte im Tetraeder, drei im Dreieck) integriert affine
  Ladungsverteilungen exakt. Allgemeinere Ortsabhängigkeit benötigt zusätzliche
  Quadraturverifikation. Umordnung eines Dreiecks ändert das vorgegebene äußere
  Ladungsvorzeichen nicht.

## Gauge, Solver und Ergebnis

Jede knotenzusammenhängende Komponente benötigt ein vorgegebenes Potential.
`ElectricGaugeError` meldet die erste unreferenzierte Komponente mit Index,
`net_charge` in C, `tolerance` in C und `compatible`. Für reine Ladungsränder
muss Volumenladung minus vorgeschriebenem äußeren D-Fluss verschwinden.
Kompatibilitätstoleranz: `1e-10 * sum(abs(component_rhs))`; ohne Ladungen ist
sie null. Ein kompatibles Problem hat eine freie additive Potentialkonstante,
ein inkompatibles keine Lösung. Keine automatische Erdung oder versteckte
Gauge-Wahl. Getrennte Körper werden einzeln geprüft.

Der freie Block wird diagonal skaliert und mit SciPy SuperLU gelöst; Pivots
kleiner/gleich `1e-12` führen zur Fehlermeldung. Das ist ein Schutz vor numerisch
singulären Systemen, kein vollständiger Konditionszahlnachweis.

`ElectrostaticResult` enthält unveränderliche Arrays:

- `node_ids`, `potential` `(nactive,)` in V;
- `electric_field`, `electric_displacement` `(ncells,3)` in V/m bzw. C/m²;
- `element_energy` `(ncells,)` und `energy` in J: positive gespeicherte Energie
  `0.5 integral E dot D`, nicht die negative elektrische Enthalpie;
- `residual`, `reactions`, `electrode_charges` `(nactive,)` in C. Reaktionen sind
  an freien DOFs null, das vollständige Residuum `A phi-f` bleibt verfügbar.

Die Gesamtladung folgt aus der diskreten Gauß-Bilanz:

```text
charge_balance = prescribed_surface_charge + dirichlet_outward_charge - volume_charge
```

Alle drei Größen und die Bilanz sind in C. Die freie Residualnorm wird durch
`|| |A| |phi| + |f| ||_2` geteilt, Grenzwert `1e-9`. Die Ladungsbilanz wird
zusätzlich gegen `1e-9` der absoluten nodalen Last-/Reaktionssummen plus
`100 eps * sum(|A| |phi| + |f|)` Rundungsbudget geprüft. Rohe P1-D-Felder sind
nicht allgemein lokal flusskontinuierlich; die globale Bilanz verwendet daher
integrierte Ladungen und Reaktionen.

Bei rho=omega=0 gilt `2 U = sum(phi * electrode_charges)`.
Auf einer vorgegebenen Kondensatorplatte erhält man Q durch Summieren ihrer
nodalen Elektrodenladungen, C=Q/Delta_phi. Die API nimmt keine automatische
Plattenzuordnung vor und modelliert keine Floating Electrode.

`RunProvenance` zeichnet das diskrete Modell als kanonisches JSON auf:
Koordinaten, Konnektivität, Materialien, integrierte Ladungslasten und aufgelöste
Potentialvorgaben (`electrostatics-discrete-v1`). Solver, Vorzeichen, Quadratur
und numerische Toleranzen sind gespeichert. Callable-Quelltexte werden nicht
archiviert; die umfassende Eingabearchivierung bleibt Teil von #27.

## Unabhängige Referenzen

- Handintegrierte Leitmatrix, Rang drei und konstantes Potential als Nullmodus;
  Materialfehler, Kopien/Unveränderlichkeit und Knotenumordnung.
- Konstante/affine Volumen- und Oberflächenladungsintegrale mit expliziten
  Vorzeichen, unabhängig von den Solverergebnissen.
- Anisotroper affiner Patch `phi=x+2y-3z` auf innen verzerrtem Würfelnetz,
  einschließlich konstanter E-/D-Felder, Energie und Randarbeit. Ein Offset von
  100 V lässt das Feld unverändert. Nur eine Potentialreferenz von 25 V ohne
  Ladungen ergibt überall konstantes Potential und Nullfeld.
- Plattenkondensator: Abstand 0.002 m, Fläche 0.03×0.04 m², epsilon=2e-11 F/m,
  Spannung 100 V, seitlich D.n=0. Exakt: C=12 pF, Q_high=1.2 nC,
  U=60 nJ und E=(-50000,0,0) V/m. Dies ist das ideale Parallelplattenmodell
  ohne Randfelder, keine Simulation eines unbeschränkten Außenraums.
- Einflächenfall: linker Potentialrand phi=0, rechter äußerer D-Fluss 2e-11 C/m²,
  epsilon=2e-11 F/m im Einheitswürfel. Exakt phi=-x V, D=(2e-11,0,0) C/m²;
  linke äußere Reaktion -2e-11 C, angrenzende Leiterladung +2e-11 C.
- Zwei gleich dicke Dielektrika, epsilon1=2e-11, epsilon2=4e-11 F/m,
  Einheitstiefe/-fläche, Endpotentiale 1 und 0 V: C=8e-11/3 F, stetiges D_x
  bei stückweise verschiedenem E_x.
- Reine Ladungsfälle mit verschwindender, ausgeglichener und inkompatibler
  Gesamtladung; getrennte Komponenten und widersprüchliche Potentialvorgaben.

Absolute Rundungstoleranzen sind den SI-Größen angepasst: `2e-26 F` bei
Elementmatrizen O(1e-11 F), `2e-14 V` im Einheits-Patch, `3e-12 V` bei 100 V
und schmalem Kondensatorspalt, `2e-22 C` für Q O(1e-9 C) sowie `2e-20 J` für
U O(1e-8 J). Die unterschiedlichen absoluten Toleranzen sind keine gemeinsamen
einheitenlosen Schwellwerte; sie berücksichtigen Skalierung und Rundung.

## Netzkonvergenz

Manufactured Solution `phi=x(1-x)+y(1-y)/2+z(1-z)/4`,
kappa=`1e-11 * [[2,.3,.2],[.3,3,.1],[.2,.1,4]] F/m`, rho=9e-11 C/m³,
exakte Potentiale auf dem gesamten Würfelrand. Die unabhängigen L2- und
Energienormintegrale verwenden eine Duffy-transformierte 4×4×4-Gaußregel.

| Unterteilungen je Achse | Tet4 | L2-Fehler (V m^(3/2)) | Energienormfehler (sqrt(J)) |
|---|---:|---:|---:|
| 2 | 48 | 0.07597742720 | 1.58113883008e-6 |
| 4 | 384 | 0.01899435680 | 7.90569415042e-7 |
| 6 | 1296 | 0.00844193636 | 5.27046276695e-7 |
| 8 | 3072 | 0.00474858920 | 3.95284707521e-7 |

Energienorm=`sqrt(integral grad(error).T kappa grad(error))`. Die beobachteten
Ordnungen sind 2 und 1 auf allen drei Übergängen. Auf den letzten beiden
Übergängen gelten die Prüffenster [1.8,2.2] bzw. [0.85,1.15]. Die integrierte
äußere D-Ladung ist 9e-11 C; für die gesamte Netzstudie gilt `atol=1e-20 C`.
Alle Konvergenztests laufen in der regulären CI-Matrix.

## Reproduktion

Nach Installation der Entwicklungsabhängigkeiten im Repository-Root:

```sh
python -m pytest tests/unit/test_electrostatics.py tests/patch/test_electrostatics_patch.py
python -m pytest tests/convergence/test_electrostatics_convergence.py -s
python examples/parallel_plate.py
```

Das Beispiel schreibt `parallel_plate.vtu` mit Potential und Elektrodenladung
als Punktfeldern sowie E, D und Energie als Volumenzellfeldern. Bei Netzen mit
unbenutzten Punkten müssen aktive Ergebnisse über `node_ids` zurückabgebildet
werden. Es gibt keine Glättung oder implizite Projektion der Rohfelder.
