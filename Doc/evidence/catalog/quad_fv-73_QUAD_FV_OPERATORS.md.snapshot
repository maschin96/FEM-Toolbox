# Konservative Quad-FV-Grundoperatoren (#229 / #226a)

Stand: 05.10.2026. Experimentelle Operator-API auf dem
[Fluidnetzadapter](72_FLUID_MESH_ADAPTER.md), getrennt vom uniformen
MAC-Kern. #226 ist gemäß Entwicklungsrichtlinien in unabhängig prüfbare
Pakete geteilt: #229 liefert diese Grundlage; #230 ergänzt SIMPLE,
Transportschemes, physikalische Randbedingungen und sämtliche ursprünglichen
Strömungs-/Bilanzgates. #226 und die Freigabe von #87 bleiben bis dahin offen.
Ein Operator-Patch ist kein konvergierter Strömungslauf und keine Validierung.

## Öffentliche API und Daten

`fem_toolbox.cfd.QuadFVGeometry(mesh, width_m=..., periodic_pairs=None,
max_cells=50000)` besitzt die validierte Quad4-Geometrie, extrudierte SI-Volumina
und Flächen. `QuadFVGeometry.from_domain(domain)` übernimmt Breite und periodische
Paarung aus `FluidDomain`; es übernimmt keine physikalische BC-Orchestrierung.
Keine Datei-I/O, neuen Abhängigkeiten oder Änderung der FEM-Model-Schemata.

| Eingabe / Ausgabe | Form und Bedeutung |
| --- | --- |
| Zellskalar / Randwerte | `(ncell,)` / `(nline2,)` |
| Zellvektor / Randwerte | `(ncell,k)` / `(nline2,k)` |
| `gradient(...)` | `(ncell,2)` bzw. `(ncell,k,2)`, letzte Achse x/y |
| `face_values(...)` | `(nface,)` bzw. `(nface,k)` |
| `FaceGradient.gradient` | `(nface,2)` bzw. `(nface,k,2)` |
| `FaceGradient.compact_normal`, `nonorthogonal_normal` | `(nface,)` bzw. `(nface,k)` |
| `cell_balance(face_flux)` | `(ncell,)` bzw. `(ncell,k)`, integrierte Bilanz |
| `StressFlux.stress_Pa`, `traction_Pa` | `(nface,2,2)` bzw. `(nface,2)` |
| `StressFlux.face_force_N`, `cell_force_N` | `(nface,2)` bzw. `(ncell,2)` |
| `rhie_chow_flux(...)` | `(nface,)`, kg/s |

Alle Facefelder folgen `mesh.faces`, Randwerte folgen der ursprünglichen
`mesh.line2`-Reihenfolge. Zell- und Randwerte müssen dieselbe Komponentenform
haben. Alle gespeicherten öffentlichen Arrays und zurückgegebenen Feldarrays
besitzen unveränderliche eigene Bytepuffer; das Sparse-Matrix-Ergebnis ist eine
abgetrennte veränderliche CSR-Matrix. Mutieren dieser Matrix beeinflusst keinen
späteren Aufruf. Felder werden weder aus Dateien geladen noch automatisch aus
Randtypen abgeleitet. Auch Neumann-/Symmetrie-Randwerte muss der spätere
BC-Orchestrator konsistent rekonstruieren; es gibt keinen stillen Nullwert.

Beispiel eines affinen Skalarpatches mit dem vorbereiteten Netz:

```python
from fem_toolbox.cfd import FluidDomain, QuadFVGeometry

domain = FluidDomain.from_file("examples/fluid_mesh/domain.json")
fv = QuadFVGeometry.from_domain(domain)
x = domain.mesh.cell_centres_m
b = domain.mesh.face_centres_m[domain.mesh.boundary_face_indices]
gradient = fv.gradient(2 + 3 * x[:, 0] - 4 * x[:, 1], 2 + 3 * b[:, 0] - 4 * b[:, 1])
# gradient[cell] = (3, -4) an allen Zellen, bis auf float64-Rundung.
```

## Geometrie, Rekonstruktion und Erhaltung

Volumen `V=b*cell_area`, Facefläche `A=b*face_length`. Normale zeigt vom
ersten Eigentümer P nach außen. Innere Faces haben Nachbar N; periodische
Nachbarzentren werden durch die vorgegebene Paartranslation in P-Koordinaten
zurückversetzt. `periodic_pairs` am Konstruktor verwendet reziproke **line2**-IDs;
`periodic_partners` am Ergebnis verwendet **mesh.faces**-IDs. Die kongruenten
Segmente besitzen entgegengesetzte äußere Normalen. Pro Paar ist die kleinere
Face-ID aktiv, die größere trägt denselben Zustand und entgegengesetzten Fluss.
`active_faces` und `neighbours` machen diese Zuordnung nachvollziehbar.

Gewichtete Least Squares minimiert pro Zelle
`sum (phi_Q-phi_P-grad(phi)_P·delta_x)^2/|delta_x|^2`.
Q ist ein innerer/periodischer Nachbar oder ein vorgegebener Randface-Sample.
Die zwei räumlichen Normalgleichungen werden vorab assembliert; die Anwendung
erfolgt mit zwei Sparse-Matrizen. Der Fit ist linear-exakt. Bei Rangverlust
oder Konditionszahl >1e12 wird die Stencil als unauflösbar diagnostiziert;
anisotrope orthogonale Wandzellen werden durch die Richtungsgewichtung weiterhin
unterstützt. Kein quadratischer Fit oder impliziter Schemewechsel.

Mit `delta_x=x_N-x_P`, `d_n=delta_x·n>0` und
`t=(x_f-x_P)·n/d_n` interpoliert `face_values` die beiden linearen
Zentrum-zu-Face-Extrapolationen mit Gewichten `1-t,t`. Dadurch bleibt auch
auf schiefen Faces ein affines Feld exakt. Physische äußere Faces übernehmen
ihren expliziten Rand-Sample.

`face_gradient` korrigiert den interpolierten Gradienten `g_bar` entlang n:

```text
g_f = g_bar + n*(delta_phi-g_bar·delta_x)/d_n
g_f·n = delta_phi/d_n - g_bar·(delta_x-d_n*n)/d_n
```

An physischem Rand ersetzt der Randface-Sample N. Der erste normale Term
ist kompakt, der zweite nichtorthogonal und separat inspizierbar. Beide
verwendeten Distanzen/Normalen sind identisch mit der Druckkorrektur. Das
Iterieren einer nichtorthogonalen Druckkorrektur gehört zu #230.

`cell_balance` addiert den gemeinsamen integrierten Fluss in P und subtrahiert
ihn in N. Periodische Partner werden dabei genau einmal verarbeitet; ein
extern vorgegebenes Paar muss bitweise entgegengesetzte Flüsse tragen, sonst
entsteht eine Diagnose statt zweier inkonsistenter Facewerte. Eine skalare
Bilanzsumme ist der Fluss über die verbleibenden physikalischen Außenfaces.
Für volumetrische Residuen teilt der Aufrufer durch `cell_volumes_m3`.

## Voller Spannungstensor

`stress_flux(U, Ub, viscosity_Pa_s=mu, boundary_viscosity_Pa_s=mub,
pressure_Pa=p, boundary_pressure_Pa=pb)` verlangt U:(ncell,2), Ub:(nline2,2).
mu/mub sind positive Skalare oder zell-/randweise positive Werte. p/pb werden
zusammen angegeben; fehlen beide, wird ausschließlich viskose Spannung berechnet.
Viskosität und Druck werden mit derselben linear-exakten Face-Rekonstruktion
interpoliert. Nichtpositive extrapolierte Faceviskosität führt zum Fehler.

```text
sigma_f = -p_f I + mu_f*(grad(U)_f+grad(U)_f^T)
traction_f = sigma_f n_f                         [Pa]
force_f = A_f traction_f                        [N]
cell_force = sum outward force_f = int div(sigma) [N]
```

Damit bleiben alle normalen, Scher- und Kreuzterme bei variabler mu erhalten.
Das ist die positive Spannungsdivergenz; im Impulsresiduum wird sie subtrahiert.
Kraft auf eine Struktur ist `-force_f`. Der Operator enthält weder Körperkraft
noch Konvektion und meldet folglich keine gesamte Impuls-/Leistungsbilanz.

## Kompakter Druckblock und Rhie-Chow

Beide Druckmethoden verlangen `momentum_diagonal_kg_s` als positiven Skalar
oder `(ncell,)`-Feld der **unrelaxierten integrierten** Impulsdiagonale a.
`d_P=V_P/a_P`, lineare Mobilitätsinterpolation liefert d_f. Relaxationsbeiträge
werden hier nicht erraten; ihre konsistente Rückrechnung gehört zu #230.

`pressure_correction_matrix(..., density_kg_m3=rho,
pressure_boundary_indices=())` assembliert den positiven kompakten Block
`-div(rho*d*grad(pc))` mit `c_f=rho*A_f*d_f/d_n`. Innere Beiträge sind
`[[c,-c],[-c,c]]`; eine ausgewählte äußere line2-ID bedeutet pc=0 und liefert
einen positiven Eigentümer-Diagonalbeitrag. Andere äußere Korrekturflüsse sind
fest null. Periodische Faces dürfen keine zusätzliche Druckvorgabe bekommen.
Es wird **keine Gaugezeile** eliminiert. Ohne Druckauslass hat der verbundene
Graph genau den konstanten Nullraum; mit Auslass ist dieser beseitigt.
Der spätere Solver setzt gegebenenfalls die einzelne Gauge und prüft danach
auch die ersetzte Kontinuitätszeile. Nichtorthogonale Anteile sind explizit.

`rhie_chow_flux(U,Ub,p,pb,...)` benutzt dieselbe Mobilität und Geometrie:

```text
m_f = rho A_f [U_f·n + d_f*(interpolate(grad p)·n-grad_face(p)·n)]
```

Innere/periodische und explizit ausgewählte Druckfaces erhalten die Korrektur;
sonstige Außenfaces behalten vorgegebene Geschwindigkeitsflüsse. Ein affiner
Druck erzeugt keine künstliche Korrektur; ein alternierender Druck wird im
Facefluss sichtbar. Der kompakte Druckkorrekturblock besitzt nur den konstanten
Nullraum; diese Aussage ist keine Behauptung über den isolierten Rhie-Chow-
Zusatzterm, der absichtlich affine Drücke annihiliert. Periodischer Druck
bezeichnet den periodischen Anteil; der vorgegebene Druckabfall ist später eine
explizite Körperkraft. Ein Flussaufruf löst weder Kontinuität noch Impuls.

## Nachweise und Grenzen

Vorab gesetzte Gates stehen in #229 und ergänzen den
[vollständigen Zielvertrag](71_M12_RANS_FSI_CONTRACT.md), ohne dessen
Strömungsgates zu ersetzen. Affine Skalar-/Vektor-/Druckpatches auf Rechtecken,
20°-Scherung und deterministisch gestörten Quads bestehen mit atol=1e-11,
rtol=0 für O(1)-SI-Größen. Starre Translation/Rotation erzeugen Nullspannung.
Für U=A x, A=((2,3),(-1,4)), mu=2+x+2y und p=5+3x-y folgt unabhängig
`div(sigma)=(5,19)`; bei b=3 m und Gebietsfläche 1 m² ist die Gesamtkraft
(15,57) N. Zell- und Außenbilanz stimmen mit dieser Rechnung.

Eigenwerte eines 4x4-Zell-Druckblocks zeigen genau einen konstanten Nullmodus;
Auslass beseitigt ihn, periodische Paarung erhält ihn. Alternierende Drücke
werden zusätzlich im tatsächlichen Rhie-Chow-Fluss geprüft. Eine unabhängige
Zwei-Zellen-Rechnung mit V=(1,1), a=(2,4), rho=3, A=2, d_n=.5 ergibt c=4.5.

Für phi=sin(pi*x)sin(pi*y) auf 16/32/64/128 Rechtecknetzen ist der RMS-Fehler
beider Gradientenkomponenten über **alle** Zellen einschließlich Rand:

| Zellen pro Richtung | RMS-Gradientenfehler | Ordnung |
| --- | --- | --- |
| 16 | 0.0100737506659 | — |
| 32 | 0.00252208152193 | 1.997914 |
| 64 | 0.000630748317520 | 1.999479 |
| 128 | 0.000157701328510 | 1.999870 |

Feinster Fehler <0.02 und jede Ordnung >1.4 erfüllen die vorab gesetzte separate
Rekonstruktionsprüfung. Es wurde keine Lösungsmatrix für Poiseuille, nichtlineare
MMS oder variable-mu-Kanalströmung gelöst; deren Vier-Netz-Gates bleiben #230.
Maximal 50.000 Zellen gemäß Netzadapter, weitere konfigurierbare Zellgrenze,
konkrete Form-/Material-/Mobilitäts-/Paarungsdiagnosen, keine automatische
Randreparatur. TDD- und Qualitätsprüfung sind im zugehörigen PR dokumentiert.

## Ausgeführte Prüfung (05.10.2026)

- Red: fehlender öffentlicher `QuadFVGeometry`-Import; danach Green der
  analytischen Operatorpatches. Zusätzlicher Red-Test: leeres zweidimensionales
  Druckrandindex-Array wurde angenommen; nach expliziter 1D-Formprüfung Green.
- 56 neue Unit-/Vier-Netz-Tests bestanden; 111 Operator-/MAC-Tests im
  Zwischenstand bestanden. Finale Gesamtsuite: 873 bestanden, 83 übersprungen
  (140.62 s); übersprungene Tests gelten nicht als bestanden.
- `python -m ruff check .`, `python -m ruff format --check .` (332 Dateien),
  `python -m mypy src/fem_toolbox` (77 Quelldateien) und `git diff --check`
  bestanden mit der vorhandenen `.venv`.
- Der dokumentierte Datei-/API-Patch liefert exakt (3,-4) auf dem
  Kanalbeispiel; 232 lokale Links der betroffenen Dokumente geprüft.
- GitHub-Matrix und unabhängiges Review sind eigene Integrationsgates;
  diese lokale Prüfung behauptet keine Review-Freigabe.
