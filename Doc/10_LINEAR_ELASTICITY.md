# Lineare statische Elastizität (Issue #12)

Status: Code- und Lösungsverifikation für die unten beschriebenen Fälle bestanden;
keine experimentelle Validierung. Diese API gehört zum Entwicklungsstand nach
v0.2.0. Der M2-Gesamtumfang ist noch nicht freigegeben.

## Modell und API

Es gelten kleine Verzerrungen, globale kartesische Achsen und SI. Die
Voigt-Reihenfolge lautet `xx, yy, zz, yz, xz, xy`, mit technischen
Schubverzerrungen. `ElasticMaterial(C)` akzeptiert eine symmetrische positiv
definite Steifigkeitsmatrix `(6,6)` in Pa. Anisotrope Daten müssen bereits in
Globalachsen vorliegen. `ElasticMaterial.isotropic(E, nu)` erzeugt die isotrope
Matrix. `E > 0` und `-1 < nu < 0.49` sind erforderlich; die obere Schranke lehnt
nahezu inkompressible, für den Verschiebungs-Tet4 ungeeignete Eingaben ab.

Symmetrie wird relativ zum größten Koeffizienten mit `1e-12` geprüft; akzeptierte
Rundungsasymmetrie wird gemittelt. Ein Eigenwertverhältnis kleiner/gleich `1e-12`
wird als numerisch nicht auflösbar abgelehnt. Diese Prüfung ersetzt keine Prüfung
physikalisch sinnvoller Materialgrößen und keine Locking-Studie für Anisotropie.

`strain_matrix(geometry)` und `stiffness_matrix(geometry, material)` im Modul
`fem_toolbox.elements.elasticity` liefern `B` und `V B.T C B` in der unveränderten
lokalen Knotenreihenfolge; beide Orientierungen sind zulässig.

`mechanical_dofs(mesh)` und `solve_elasticity(...)` liegen im Modul
`fem_toolbox.solvers.elasticity`. Die DOFs umfassen nur Tet4-Trägerknoten,
aufsteigend nach internem Knotenindex, je Knoten `[ux, uy, uz]`. Vorgaben sind
`Dirichlet(dofs.indices("u", nodes, components), values)`; Werte sind in Metern
und müssen als flacher Vektor zu den Indizes passen. Nichtnullige Vorgaben sind
zulässig. Ein Material gilt homogen; eine Materialsequenz weist jeder Tet4-Zeile
explizit ein Material zu.

Optionale Lasten:

- `body_force`: Vektor `(3,)` in N/m³ oder Funktion `x -> vector`;
- `tractions`: Zuordnung von Triangle3-Zeilen zu globalen Vektoren in N/m² oder
  Funktionen; nur äußere Tet4-Flächen sind zulässig;
- `point_loads`: Zuordnung von aktiven Mesh-Knoten zu Vektoren `(3,)` in N.

Funktionsargumente sind räumliche Punkte `(3,)` in Metern. Konstante Lasten
werden exakt integriert. Für Funktionen gilt Quadraturgrad 2 mit vier Punkten
im Tet4 und drei Punkten im Triangle3; affine Lastverteilungen sind damit exakt.
Höhere Ortsabhängigkeit benötigt eine eigene Quadraturkonvergenzprüfung.
Traktionen sind vorzeichenbehaftete **globale Vektoren**, keine Druckskalare;
Flächenumordnung ändert ihre physikalische Richtung nicht. Alle Beiträge addieren
sich, auch auf vorgeschriebenen Freiheitsgraden. Punktlasten verursachen im
Kontinuum lokale Singularitäten; der analytische Zugtest verwendet Flächenlasten.

## Solver und Ergebnis

Vor der Faktorisierung wird auf jeder knotenzusammenhängenden Komponente geprüft,
ob alle sechs Starrkörpermoden durch Vorgaben unterdrückt sind. Zentrierte und
längennormierte Koordinaten machen die Rangprüfung unabhängig von Translation
und einheitlicher Skalierung; ihre dimensionslose Rangtoleranz ist `1e-10`.
Die symmetrieerhaltende Elimination baut auf #11 auf. SciPy SuperLU löst den
freien Block nach diagonaler Gleichgewichtung. Pivots kleiner/gleich `1e-12`
werden abgelehnt, um zusätzliche Mechanismen (z.B. punktverbundene Körper) oder
schlecht konditionierte Systeme nicht stillschweigend zu lösen. Das ist eine
numerische Schutzprüfung, kein vollständiger Konditionszahlnachweis.

`ElasticResult` enthält unveränderliche Arrays:

- `node_ids`, `displacement`: aktive Knoten und Verschiebung `(nactive,3)` in m;
- `strain`, `stress`: rohe elementweise Werte `(ncells,6)` bzw. Pa;
- `element_energy`, `energy`: Energie je Element bzw. Summe in J;
- `applied_loads`, `reactions`, `residual`: `(nactive,3)` in N;
- `force_balance`, `moment_balance`: Gesamtbilanz in N bzw. N m.

Reaktionen sind auf freien DOFs null; das Residuum `K u - f` bleibt vollständig
verfügbar. Die relative freie Residualnorm verwendet als Nenner
`|| |K| |u| + |f| ||_2`. Diese Rückwärtsfehler-Skalierung bleibt für starre
Bewegungen aussagekräftig, bei denen sich innere Kräfte aufheben. Der Grenzwert
ist `1e-9`; physikalische Bilanzen werden zusätzlich berichtet und in den
Referenztests geprüft. Die Momentenbilanz bezieht sich auf den globalen Ursprung.

Jeder Lauf liefert außerdem `provenance: RunProvenance`. Der deterministische
SHA-256 bezieht sich auf das diskrete Modell (Koordinaten, Tet4-Konnektivität,
Materialmatrizen, integrierte Lasten und aufgelöste Dirichlet-Vorgaben) als
kanonisches JSON. Funktionsquelltexte werden nicht serialisiert: ihr ausgewerteter
Lastvektor ist maßgeblich. Die Parameter benennen das Eingabeschema
`elasticity-discrete-v1`, Solver, Quadratur und numerische Toleranzen. Die
vollständige Archivierung der ursprünglichen Modellbeschreibung folgt mit #27.

## Unabhängige Nachweise

Die Tests verwenden algebraische und analytische Erwartungen, keine gespeicherten
Solverresultate:

1. Isotropie, reine Schubverzerrung, anisotrope Energie und Materialfehlerfälle.
2. Sechs Starrkörpermoden, Symmetrie, Rang sechs der Elementsteifigkeit sowie
   Knotenumordnung einschließlich umgekehrter Orientierung.
3. Affiner anisotroper Patch auf einem innen verzerrten Würfelnetz: exakte
   Verzerrung/Spannung, Energie sowie Kraft-, Momenten- und Arbeitsbilanz.
4. Zugwürfel: Kantenlänge 1 m, E=1000 Pa, nu=0.25, Zugtraktion 10 Pa.
   Symmetriebedingungen `ux=0` auf x=0, `uy=0` auf y=0, `uz=0` auf z=0 erlauben
   freie Querkontraktion. Exakt gilt `u=(0.01x,-0.0025y,-0.0025z)`,
   `sigma_xx=10 Pa`, restliche Spannungen null, Energie 0.05 J und
   Gesamtlagerreaktion `(-10,0,0) N`.
5. Quadratische Manufactured Solution `u=(x²,0,0)` auf dem Einheitswürfel,
   E=100 Pa, nu=0.25: `epsilon_xx=2x`, `sigma=(240x,80x,80x,0,0,0) Pa`,
   Körperlast `(-240,0,0) N/m³`, exakte Verschiebung auf dem gesamten Rand.
   Die Fehlernormen werden unabhängig mit einer Duffy-transformierten
   4×4×4-Gaußregel integriert; der quadratische Verschiebungsfehler wird dabei
   exakt integriert. Energienorm: `sqrt(integral delta_epsilon.T C delta_epsilon)`.

| Unterteilungen je Achse | Tet4 | L2-Fehler | Energienormfehler |
|---|---:|---:|---:|
| 2 | 48 | 0.04564354646 | 3.16227766017 |
| 4 | 384 | 0.01141088661 | 1.58113883008 |
| 6 | 1296 | 0.00507150516 | 1.05409255339 |
| 8 | 3072 | 0.00285272165 | 0.79056941504 |

Die L2-Norm hat die Einheit m^(5/2), die Energienorm sqrt(J). Aus den tatsächlichen
Verfeinerungsverhältnissen ergibt sich auf allen drei Übergängen Ordnung 2 bzw.
1. Der Test verlangt auf den letzten beiden Übergängen [1.8,2.2] bzw. [0.85,1.15].
Die Fenster erlauben asymptotische Restfehler, aber keine falsche P1-Ordnung.

Absolute Testtoleranzen beziehen sich auf diese Parameter: bis `5e-14` für
Elementalgebra bei Längen O(1 m) und E=120 Pa, `3e-14` für Patchverzerrungen und
`5e-13 Pa` für Patchspannungen, `2e-12 N` bzw. `2e-12 N m` für Zugbilanzen.
Im vierstufigen Konvergenztest werden wegen Akkumulation über bis zu 3072 Elemente
`1e-9 N` bzw. `1e-9 N m` zugelassen. Diese Grenzen sind Rundungsbudgets, keine
nachträgliche Anpassung fehlerhafter Referenzwerte.

## Reproduktion

Nach Installation mit Entwicklungsabhängigkeiten im Repository-Root:

```sh
python -m pytest tests/unit/test_elasticity.py tests/unit/test_elastic_loads.py
python -m pytest tests/patch/test_elasticity_patch.py
python -m pytest tests/convergence/test_elasticity_convergence.py -s
python examples/elastic_tension.py
```

Das Beispiel schreibt `elastic_tension.vtu`. Auf Triangle3-Zellen enthalten die
Spannungs-/Verzerrungsfelder ausschließlich Null-Platzhalter; mechanische
Rohwerte liegen nur auf den Volumenzellen. Für Netze mit unbenutzten Punkten muss
beim Export über `result.node_ids` auf alle Mesh-Punkte zurückabgebildet werden.
Die kurze Konvergenzstudie läuft in der regulären Testsuite und damit in allen
CI-Python-Versionen; keine Konvergenztests werden ausgelassen.

## M7-Erweiterung: Tet10

Der direkte Solver und Model-Schema 1 unterstützen zusätzlich geradkantige
Tet10 mit Triangle6-Traktionen. Dehnung/Spannung sind explizite Schwerpunktwerte;
Energie bleibt ein Zellintegral. Thermische Tet10-Eigendehnungen sind inzwischen
[punktabhängig integriert](48_TET10_THERMOMECHANICS.md). [Vertrag, Biegepatch und Konvergenz](47_TET10_ELASTICITY.md).
