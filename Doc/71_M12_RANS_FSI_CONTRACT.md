# M12-Vertrag: stationäre RANS und bidirektionale FSI (#84)

Entscheidungsstand: 04.10.2026, ADR-049. Dieser Vertrag legt die Entwicklung
für #85, #226 und #87–#94 fest; er behauptet keine implementierte oder validierte
RANS-/FSI-Funktion. Der bestehende [uniforme MAC-Kern](70_STATIONARY_FLOW_CORE.md)
aus #86 bleibt eine unabhängige Referenz. Seine dokumentierten Resultate werden
nicht auf neue Netze oder Turbulenz übertragen. Integration: Issue-Branches
nach `codex/milestone-m12`; Release nach `main` erst über #94 und dessen Gates.

## 1. Umfang und gewählte Verfahren

Zunächst ebene 2D-Mittelströmung, konstante positive Dichte und molekulare
Viskosität, stationär, inkompressibel, isotherm. FSI koppelt diese Strömung mit
der vorhandenen 3D-Kleinverformungsmechanik eines explizit extrudierten Körpers
konstanter Breite b. Eine ebene Strömung genügt nur bei nachgewiesener
Spannweiteninvarianz; allgemeine 3D-Fluide und beliebige 3D-Verformung gehören
nicht zu dieser Abnahme. Breite, seitliche Strukturbindungen und externe
Druckbelastung sind Pflichtdaten, keine impliziten Einheitsannahmen.

Entscheidungen:

- Eine Schließung: SA-noft2 mit der unten vollständig benannten S-Modifikation,
  ohne Trip, Transition, Rauheit, Rotation/Krümmungskorrektur oder Wandfunktion.
- Körperangepasste, strukturierte, konvexe Quad4-Fluidzellen mit geraden Kanten;
  konservative zellzentrierte Finite Volumen, Rhie-Chow-Flussinterpolation und
  SIMPLE-Druckkorrektur. Deformation erhält Konnektivität; kein Remeshing.
- Begrenztes lineares Upwind für Transport; Donor-Cell als explizit gewählte
  Referenz. Gradienten aus linear exakter Least-Squares-Rekonstruktion;
  konsistente nichtorthogonale Diffusionskorrektur. Kein stiller Schemewechsel.
- Partitionierte stationäre Dirichlet-Neumann-FSI mit konstanter Relaxation,
  transponierter Lastabbildung und harmonischer Fluidnetzverschiebung.

Die kollokierte Erweiterung ersetzt den MAC-Operator nicht. Die neue
Druckkopplung ist für körperangepasste Netze gewählt; ihre Prüfung gehört zu
[Operator-Issue #226](https://github.com/maschin96/FEM-Toolbox/issues/226),
bevor #87 auf solchen Netzen akzeptiert wird. Begründung
und Primärimplementierung für SIMPLE/Rhie-Chow: [INL OpenPronghorn](https://openpronghorn.inl.gov/source/executioners/SIMPLE.html).
Festkörper-Tet4-/Tet10-Formfunktionen dienen ausschließlich der Struktur und
der Schnittstelleninterpolation. Kein gleichrangiges u/p-Tet4 ohne eigene
Stabilitätsanalyse, keine Druckstrafe. DNS/LES, URANS, dynamische FSI,
Kompressibilität, Mehrphasenströmung und Kontakt bleiben ausgeschlossen.

## 2. Gleichungen, Vorzeichen und schwache Formen

U bezeichnet die mittlere Fluidgeschwindigkeit (m/s), d_s die
Strukturverschiebung (m). Die Symmetrisierung ist S(U)=(grad U+grad U^T)/2.
Druck p ist der mechanische Modelldruck in Pa; SA liefert kein k und damit
keinen separat rekonstruierbaren isotropen Reynolds-Spannungsanteil.

```text
nu = mu/rho                         [m²/s]
nu_t = nu_tilde f_v1                 [m²/s]
mu_eff = mu + rho nu_t               [Pa s]
sigma_f = -p I + 2 mu_eff S(U)       [Pa]
div U = 0
rho div(U tensor U) - div sigma_f = f [N/m³]
```

Bei variabler Wirbelviskosität ist `-mu_eff Laplace(U)` kein Ersatz für
`-div(2 mu_eff S(U))`. Scher- und Kreuzterme sowie Gradienten von mu_eff
müssen im Operator und bei der Traction-Auswertung übereinstimmen.

Für Geschwindigkeitstest v mit v=0 auf Dirichlet-Grenzen und Drucktest q:

```text
-int rho grad(v):(U tensor U) + int 2 mu_eff S(v):S(U)
-int p div(v) = int v·f + int_boundary v·(sigma_f n_f - rho U (U·n_f))
int q div(U) = 0
```

Alle Volumenintegrale liegen auf dem aktuellen Fluidgebiet. n_f zeigt nach
außen aus dem Fluid. Die konservative Randform enthält den konvektiven Term;
bei natürlichen Grenzen darf er nicht verschwinden. FV integriert dieselben
Bilanzen zellweise und verwendet jeden inneren Fluss genau einmal mit
entgegengesetzten Eigentümervorzeichen. Die Struktur behält die bestehenden
[SI-/Voigt-Konventionen](02_MATHEMATICAL_MODEL.md).

Randvertrag:

| Rand | Geschwindigkeit und Druck | SA-Variable |
| --- | --- | --- |
| Einlass | vorgegebene endliche U, kompatibler positiver Nettofluss; Drucknormalgradient aus Impulsbilanz | explizites nichtnegatives nu_tilde-Profil in m²/s |
| Druckauslass | p=p_out, normaler Geschwindigkeitsgradient null; keine zusätzliche Gauge | normaler Gradient null |
| Feste oder FSI-Wand | U=0 am aktuellen stationären Rand, Drucknormalgradient aus Impulsbilanz | nu_tilde=0 |
| Symmetrie | U·n=0, tangentialer Normalgradient null, dp/dn=0 | normaler Gradient null |
| Periodisches Kanalpaar | identische U, nu_tilde und periodischer p-Anteil; vorgegebener Druckabfall als f_x | periodisch |

Der Druckauslass ist eine definierte Randapproximation, keine frei wählbare
FSI-Traktion. Rückströmung an einem solchen Auslass führt zum Fehler mit
Face-ID/Fluss; keine Ersetzung durch erfundene Einlassdaten. Unvereinbare
Geschwindigkeitsvorgaben/Nettoflüsse werden vor dem Solve abgelehnt. Ohne
Druckauslass wird exakt ein Zellwert des periodischen/geschlossenen Drucks
fixiert; alle Kontinuitätsbilanzen einschließlich einer gegebenenfalls ersetzten
Zeile werden anschließend geprüft. Mehrere Druckreferenzen sind unzulässig.

## 3. SA-noft2: Gleichung, Konstanten und Schutzmaßnahmen

Mit a=nu_tilde>=0 in m²/s und Wandabstand d>0 in m:

```text
div(U a) = P - D + div(K grad a) + (c_b2/sigma)|grad a|²
K=(nu+a)/sigma; P=c_b1 S_tilde a; D=c_w1 f_w (a/d)²
chi=a/nu; f_v1=chi³/(chi³+c_v1³)
f_v2=1-chi/(1+chi f_v1)
Omega=sqrt(2 W:W); W=(grad U-grad U^T)/2
S_bar=a f_v2/(kappa² d²)
S_tilde=Omega+S_bar                                if S_bar>=-c_2 Omega
S_tilde=Omega+Omega(c_2² Omega+c_3 S_bar)/((c_3-2c_2)Omega-S_bar) otherwise
r=min(a/(S_tilde kappa² d²),10); r=10 if S_tilde=0
g=r+c_w2(r⁶-r); f_w=g[(1+c_w3⁶)/(g⁶+c_w3⁶)]^(1/6)
c_b1=0.1355; c_b2=0.622; sigma=2/3; kappa=0.41
c_v1=7.1; c_w2=0.3; c_w3=2; c_2=0.7; c_3=0.9
c_w1=c_b1/kappa²+(1+c_b2)/sigma
```

Modellquelle: [TMR SA-noft2 und Note 1(c)](https://tmbwg.github.io/turbmodels/spalart.html),
S-Modifikation: [Allmaras, Johnson, Spalart (2012)](https://www.iccfd.org/iccfd7/assets/pdf/papers/ICCFD7-1902_paper.pdf).
Diese Form hat keinen f_t2- oder Trip-Term. Der dokumentierte r-Grenzwert und
die S-Modifikation sind Teil des Modells. Alle SA-Gleichungsterme haben
m²/s²; nu_tilde ist weder turbulente Energie k noch dynamische Viskosität.
Für Skalar-Testfunktion w mit homogener Dirichlet-Testspur folgt:

```text
-int a U·grad(w) + int K grad(w)·grad(a)
 = int w(P-D+c_b2/sigma |grad a|²)
   + int_boundary w(K grad(a)·n_f - a U·n_f)
```

An Wand-Dirichlet-Punkten wird a=0 gesetzt; volumenbezogene Quellen werden
nur an Zellzentren mit positivem d ausgewertet. Keine Auswertung von a/d bei
d=0. Wandabstand ist die minimale euklidische Entfernung zu den aktuellen
Wandsegmenten, einschließlich Endpunkten; nächster Knoten oder Entfernung
entlang einer Netzlinie genügt nicht. Nach jeder Netzbewegung neu berechnen.

Implementierung #87 nach #226: Destruktion bei eingefrorenem f_w als nichtnegativen
impliziten Sink behandeln; Diffusionskoeffizient einfrieren. Produktion und
Gradientquadratterm im finalen nichtlinearen Residuum neu auswerten.
Line-Search verkleinert den gemeinsamen Aktualisierungsschritt bis a>=0 und
alle Größen endlich sind (Start alpha=0.5, Halbierung, höchstens 20 Versuche).
Sonst konkrete Nichtkonvergenz mit Minimum, Zell-ID und letztem Residuum.
Kein Clipping negativer a-Werte, keine harte Obergrenze für nu_t/nu und kein
Wechsel auf SA-neg. Schutz gegen Division/Überlauf ist diagnosepflichtig.
Der laminare Modus wird ausdrücklich gewählt, nicht aus kleinem a erraten.

Alle Inletwerte sind Pflicht; bei der unten gewählten Kanalreferenz entfällt
ein Einlass durch Periodizität. Geeignet sind zunächst glatte, wandgebundene,
voll turbulente Strömungen im geprüften Bereich. Transition, starke Rotation,
Rauheit und allgemeine stark abgelöste Strömungen werden nicht validiert.

## 4. Netz, Diskretisierung und nichtlinearer Abbruch

Fluidnetz-Adapter #85: UTF-8-JSON `fluid_mesh_schema_version=1`,
`length_unit="m"`, points[n,2], quad4[ncell,4] gegen den Uhrzeigersinn,
explizite line2[nface,2] und benannte Randgruppen als Face-Indizes.
Alle äußeren Kanten genau einmal, keine innere Kante als Rand, keine Lücken,
mehrfachen Eigentümer oder sich schneidenden Zellen. Genau ein physikalischer
Randtyp je Face; FSI ist zusätzlich eine Wandmarkierung. Periodische Partner
sind explizit und kongruent nach Translation; ihre äußeren Normalen sind
entgegengesetzt. Gespeicherte Segmentendpunktreihenfolge darf abweichen
(Präzisierung ADR-050).
JSON und TOML laden Daten, keinen ausführbaren Code. Vorbereitete Netze und
Import ohne interne Manipulation sind Pflicht; kein automatischer Generator.

Geometriegates vor jedem Solve: positive Fläche, positiver Quad-Jacobian an
allen Ecken und 2x2-Gaußpunkten, konvexe Zellen, positive Eigentümerdistanz
zur Face, maximale Nichtorthogonalität 60 Grad, Skewness <=0.25.
Nichtorthogonalität ist der Winkel zwischen Eigentümer-Nachbar-Verbindung
und Flächennormale; Skewness die Distanz des Face-Mittelpunkts zum Schnitt
dieser Verbindung mit der Facelinie, geteilt durch Verbindungslänge.
An Randfaces wird die Zellzentrum-zu-Face-Verbindung verwendet.
Anisotrope Wandzellen dürfen großes Seitenverhältnis haben; ausschließlich
hieraus wird keine Entartung abgeleitet. Wandnormaler Abstand der ersten
Zellzentren, minimale Fläche, Winkel, Skewness und Streckung werden berichtet.
Nachrechnung verlangt max y+=d_1 sqrt(|tau_w|/rho)/nu<=1 an allen no-slip-Wänden;
benachbarte wandnormale Schichtdicken dürfen höchstens Faktor 1.2 wachsen.
Ein Vorabschätzwert ersetzt das finale y+-Gate nicht. Verletzung ergibt einen
nicht akzeptierten Lauf mit Face-/Zell-IDs, keine automatische Wandfunktion.
Diese Zahlen sind konservative Projektgrenzen, keine universellen Modellgarantien.

Rhie-Chow: Massendurchsatz der Face aus interpoliertem Impulsprädiktor und
kompaktem Druckunterschied, einschließlich derselben Impulsdiagonale und
nichtorthogonalen Korrektur wie im SIMPLE-Schritt. Auf orthogonalen Faces:
`U_f=interpolate(U)+d_f(interpolate(grad p)-grad_face p)`,
`d_P=V_P/a_P` aus der integrierten unrelaxierten Impulsdiagonale;
Relaxationsbeiträge müssen konsistent rückgerechnet werden. Keine unabhängige
Interpolation von Massenstrom und Skalartransport. Druckkorrektur löst
`sum_faces rho A_f (U_f - d_f grad_face p_correction)·n=0`.
Nach Korrektur wird der konservierte Facefluss gespeichert; Zellgeschwindigkeit
allein ersetzt ihn nicht. Der Operator muss konstante/lineare Druckfelder exakt
rekonstruieren und einen alternierenden Druckmodus erkennen. Nichtorthogonale
Korrektur wird bis zu ihrem Residualgate iteriert; höchstens 20 Schritte.
Monotone Transportreferenz nach [NIST FiPy](https://pages.nist.gov/fipy/en/stable/numerical/scheme.html).

Äußere RANS-Iteration: mu_t aktualisieren, Impuls/SIMPLE konvergieren,
SA lösen, zulässigen relaxierten Zustand bilden, alle Residuen neu aufbauen.
Maximal 500 äußere Schritte, alpha_U=alpha_p=alpha_a=0.5; Werte konfigurierbar
und in Provenienz. Sparse-LU ist Referenz; vorhandene iterative Backends erst
nach separater Prüfung der unsymmetrischen Fluidblöcke. Kein CG auf dem
gekoppelten Druck-/Impulssystem. Fehlende Konvergenz ist ein Fehler.

Für U_ref>0, L_ref>0 gelten pro finalem Lauf (explizite Pflichtskalen):

| Größe | Norm und Standardgate |
| --- | --- |
| Kontinuität | max abs(sum Face-Massenströme)/V <= 1e-8 rho U_ref/L_ref |
| Globale Masse | abs(sum äußere Massenströme) <= 1e-10 rho U_ref A_ref |
| x-/y-Impuls | je max abs(R_i/V) <= 1e-8 rho U_ref²/L_ref |
| SA | max abs(R_a/V) <= 1e-8 U_ref nu/L_ref |
| Globale Impulsbilanz | je abs(sum rho U(U·n)A - sum sigma_f n A - sum f V) <= 1e-8 rho U_ref² A_ref |

A_ref=b L_ref; V=b Zellfläche. SA kann unmultipliziert mit rho integriert
werden; dann Einheit des integrierten R_a ist m⁵/s². Rohresiduen und Skalen
werden beide ausgegeben. Größere SA-Werte dürfen eine zusätzliche relative
Skala verwenden, aber das obige absolute Gate nicht ersetzen.
Leistung: Dissipation 2 mu_eff S:S, Körper-/Randarbeit und konvektive
Energieflüsse getrennt, numerische Transportdissipation explizit. Rest <=
1e-6 rho U_ref³ A_ref. Wirbelviskose Dissipation ist ein Modellterm, keine
vollständige turbulente Energiebilanz. Iterandenänderung genügt keinem Gate.

## 5. Stationäre FSI, Lastübertragung und Netzbewegung

Auf gemeinsamer Schnittstelle n_s=-n_f. Last auf Struktur:
`t_s=-sigma_f n_f=sigma_f n_s`. Positiver Fluiddruck drückt die Struktur
vom Fluid weg; viskose Scherlast wird mitübertragen. Absolute Differenz
p-p_external bestimmt den Drucklastanteil. Ein Gaugewechsel darf Lasten nur
unverändert lassen, wenn p_external um denselben Wert verschoben wird.
Kein Multiplikator 2/3 rho k; dieser ist in der gewählten SA-Schließung unbekannt.

Ebene Fluidrandsegmente werden über b zu 3D-Flächen extrudiert. Struktur-
Triangle3-Facetten müssen diese Fläche geometrisch abdecken und benannte
Sätze besitzen; keine Tet10-Schnittstelle im ersten Nachweis. Quadratur auf
der gemeinsamen Segment-/Dreiecksüberdeckung definiert H, die
Strukturverschiebung an Fluidquadraturpunkten interpoliert. Lasten sind
bereits integrierte Kräfte F_q in N, nicht Druckwerte:

```text
d_fluid = H d_structure
F_structure = H^T F_fluid
sum F_structure = sum F_fluid                 if H reproduces constants
sum x_structure cross F_structure = sum x_q cross F_q if H reproduces coordinates
d_structure^T F_structure = d_fluid^T F_fluid
```

H basiert auf baryzentrischen P1-Gewichten in den zugeordneten Dreiecken.
Keine nearest-neighbor-Lastverteilung, keine Extrapolation; fehlende Überdeckung,
Lücken oder Mehrdeutigkeit führen zum Fehler. Geometrische Zuordnungstoleranz
ist 1e-10 L_ref; Integrationstoleranzen sind unabhängig davon. Überdeckung
und Koordinatenreproduktion auf aktueller Geometrie erneut prüfen. Proben über
Breite werden flächengewichtet auf die ebene Verschiebung reduziert; Abweichung
vom gewichteten Mittel >1e-8 L_ref oder |d_z|>1e-8 L_ref ist nicht unterstützte
3D-Verformung. Die Reduktion und ihre transponierte Lastabbildung müssen
Bestandteil desselben H sein. [preCICE](https://precice.org/configuration-mapping.html)
unterscheidet entsprechend konservative Kräfte und konsistente Interpolation;
wir implementieren die obige überprüfbare Paarung selbst, ohne preCICE-Abhängigkeit.

Netzbewegung #90: Für beide Koordinaten Laplace(d_mesh)=0 im Referenzgebiet,
Schnittstelle d_mesh=H d_s, übrige feste Wände d_mesh=0. Bei verschobenen
Kanalecken werden Ein-/Auslassknoten linear zwischen ihren Wandendpunkten
verschoben; widersprüchliche Ecke-Vorgaben sind ein Fehler. Auf übrigen offenen
Rändern normaler Gradient null. Referenzpunkte plus absolute Verschiebung,
keine kumulierten Inkremente. Neue Flächen, Normalen, Volumina, H und Wandabstände
vor dem nächsten Fluidsolve. Stationäre Kopplungsiteration ist keine physikalische
Zeit: keine Wandgeschwindigkeit d/Delta_t und kein ALE-Zeitfluss.

Zulässigkeit: max |d_s|/L_ref<=0.02, max Frobeniusnorm grad(d_s)<=0.02,
alle Netzgates aus Abschnitt 4, minimale Zellfläche >=0.5 Referenzfläche.
Kleine rotationsfreie Verzerrung genügt allein nicht; Geometriegates bleiben
Pflicht. Überschreitung beendet den Versuch mit Diagnose; keine stille
Last-/Verschiebungsbegrenzung. Das sind Abnahmegrenzen für lineare Struktur auf
Referenzgeometrie und aktualisierte Fluidgeometrie, keine Großverformungs-FSI.

Partitionierte Iteration k:

1. Aus d^k aktuelle Fluidgeometrie bauen und prüfen, RANS vollständig lösen.
2. Aktuelle Fluidtraktion konservativ übertragen, Struktur auf unverändertem
   Referenznetz lösen; starre Moden und Strukturresiduen prüfen.
3. r_d=d_solved-d^k bilden, d^(k+1)=d^k+omega r_d, Standard omega=0.2.
4. Wiederholen, höchstens 200 Kopplungsschritte. Kein Erfolg bei bloß kleiner
   relaxierter Änderung. Prüfen des unrelaxierten r_d und der Laständerung.

Abbruch verlangt max |r_d|<=1e-8 L_ref UND max |F^k-F^(k-1)|<=
1e-6 rho U_ref² A_ref sowie sämtliche Fluid-, Netz- und Strukturgates.
Anschließend Fluid und Struktur auf dem akzeptierten Geometriestand nochmals
prüfen; bei Verletzung weiter iterieren statt alte Lasten als finales Ergebnis
zu deklarieren. Kraft-/Moment-/virtuelle Arbeitsfehler der Übertragung <=1e-10
relativ mit absoluten Floors 1e-12 N, 1e-12 N m und 1e-12 J für O(1)-Tests.
Virtuelle Arbeit ist ein Mappingtest; stationäre Wandgeschwindigkeit null
liefert keinen dynamischen Energieaustauschnachweis. Divergenz, Ressourcenlimit,
Netzfehler oder inneres Nichtkonvergieren beenden mit letzter vollständiger
Historie und klarer Fehlerart. Eine instationäre Lösung wird nicht ersatzweise
als stationäre FSI ausgegeben.

## 6. API, Provenienz, Adapter und Lizenzen

Geplante separate Python-Objekte `FluidMesh2D`, `RansProblem`, `RansResult`,
`StationaryFsiProblem`, `StationaryFsiResult` unter `cfd/` und `coupling/`.
Das sind Zielnamen, noch keine verfügbaren Imports. Model-Schema 1–3 und
FlowProblem-Schema 1 werden unverändert erhalten; neue deklarative Verträge
tragen eigene `rans_schema_version=1` und `fsi_schema_version=1`.
#93 liefert strikte TOML-Adapter mit Pflichtwerten und Fehlerpfaden.

Datenfluss: Meshadapter/BC-Preflight -> FV-Geometrie -> RANS -> Tractions ->
InterfaceMap -> vorhandener Struktursolver -> MeshMotion -> Kopplungskontrolle.
Physikmodule bleiben frei von Datei-I/O. Öffentlich dokumentierte Formen:
U[ncell,2] m/s, p[ncell] Pa, a/nu_t[ncell] m²/s,
mass_flux[nface] kg/s, wall_distance[ncell] m, traction[nwall,2] Pa;
Struktur d[nnode,3] m und F[nnode,3] N. Alle Rückgabedaten unveränderlich;
interne Eigentümer-/Normalenkonventionen und Mappingindizes sind nachvollziehbar.

Provenienz nennt vollständige Eingaben, Netz-/Interfacehashes, Fluidparameter,
SA-Variante und Konstanten, Druckreferenz/externe Last, Breite, Schemes,
Toleranzen/Skalen, Relaxation, Solver-/Toolbox-/Dependencyversionen.
Ergebnisartefakte enthalten rohe Feldresiduen, Bilanzen, y+, Netzgates,
Iterationstabellen und Erfolg/Fehler. Wiederverwendung nur mit passender
Geometrie und sämtlichen Eingabehashes; gespeicherte Konvergenz ist keine Validierung.

Keine neue Laufzeitabhängigkeit beschlossen. NumPy (BSD-3-Clause), SciPy
(BSD-3-Clause), meshio (MIT) bleiben bestehende Adapterbasis; die Fluid-JSON-
Geometrie wird eigenständig gelesen, VTU nutzt meshio. Gmsh ist optional zur
externen Geometrievorbereitung (GPL), kein eingebetteter neuer CFD-Generator.
preCICE, OpenFOAM und INL sind Methoden-/Referenzquellen, keine übernommenen
Codebestandteile. Fremde Daten/PDFs werden verlinkt, nicht ohne Rechtebeleg
redistribuiert. #88/#92 müssen Downloads mit Herkunft, Datum, SHA-256,
Parametern und Nutzungsrechten manifestieren; ungeklärte Weitergabe blockiert
Fixture-Veröffentlichung, nicht die eigenständig hergeleiteten Tests.

## 7. Unabhängige Referenzen und vorab gesetzte Gates

Alle Grenzen in diesem Abschnitt sind Projektentscheidungen vor Implementierung;
keines der neuen Gates ist bereits bestanden. Fehlende Datennutzungsrechte,
falsche Modellvariante oder ungedeckte Geometrie sind no-go. Ein gescheitertes
Budget wird untersucht; nachträgliche Änderung benötigt begründeten ADR und
neue Referenzplanung, keinen blind neu erzeugten Golden-Wert.

| Nachweis / Issue | Referenz und Messgrößen | Gate |
| --- | --- | --- |
| FV-Grenzfall #226/#88 | Poiseuille und nichtlineare MMS aus Doc/70, orthogonale Quad-Netze | bestehende vier Netzstufen und Fehlerbudgets unverändert; neuer Operator unabhängig gegen MAC |
| Variable Viskosität #87/#88 | U=(y(1-y),0), mu_eff=1+y Pa s; f_x=1+4y N/m³, p=0, Gebiet [0,1]²; analytischer voller Spannungstensor | Raster 16/32/64/128: feinster RMS-U-Fehler <1e-3 m/s, Ordnung >1.8; Bilanz-/Residualgates Abschnitt 4 |
| Netz-/Druckstabilität #85/#226 | affine Geschwindigkeiten/Druck, alternierender p, konstante Traction; orthogonale und 20-Grad-Scherung von Quad-Netzen | affine Gradienten/konstante Flüsse bis 1e-11 bei O(1)-SI-Daten; Checkerboard nicht im Nullraum, einziger Drucknullraum konstant |
| SA-Code #87/#88 | unabhängige symbolische Punktwerte aller Terme und Wand-/Null-/S-Verzweigungen; glattes positives a=nu[1+sin(pi x)sin(pi y)], U=(1,0), d=y+1 auf [0,1]², virtueller Wandlinie y=-1 | Quellen aus separater Differentiation, Dirichletdaten analytisch; 16/32/64/128: relativer RMS-a-Fehler <1e-2, Ordnung >0.8 mit begrenztem Upwind; keine negative a-Zelle |
| RANS-Code #88 | voll entwickelter periodischer SA-Kanal, delta=1 m, rho=1 kg/m³, nu=1/587.19 m²/s, f_x=1 N/m³; eigener unabhängiger 1D-BVP-Solver für volle Modellgleichungen | vier Wandnetze mit 32/64/128/256 Zellen pro Halbkanal; feinster RMS-U-Fehler <0.5% und Wandkraftfehler <0.1%; y+<=1 |
| Modellvergleich #88 | [Moser/Kim/Mansour Re_tau=590](https://turbulence.oden.utexas.edu/MKM_1999.html), mittleres U+ und Bulkgeschwindigkeit; keine Kalibrierung der SA-Konstanten | RMS-U+-Abweichung <=5% des DNS-RMS-U+, relative Bulkabweichung <=5%; Fehler der letzten Netzverdopplung <1% und höchstens 1/5 des Modellbudgets |
| Mapping #89/#92 | konstante/affine Verschiebung und konstante/lineare Traction auf nicht passenden Segment-/Dreiecksnetzen, Zwei-Punkt-Rechnung unten | Kräfte, Momente, virtuelle Arbeit gemäß Abschnitt 5; explizite Fehlfälle Überdeckung und falsche Breite |
| Struktur #92 | vorhandene analytische Zug-/Biegetests, extrudierte Struktur mit kleiner Belastung | vorhandene Struktur-V&V bleibt Pflicht; keine neue ungeprüfte Strukturformulierung |
| Bidirektionale FSI #92 | elastisch geführter paralleler Kanaldeckel, feste Q=1e-3 m³/s, H0=0.1 m, L=1 m, b=1 m, mu=1 Pa s, rho=1 kg/m³, k_s=1e4 N/m, p_external=p_out=0 | positive Wurzel d(H0+d)³=6 mu L² Q/k_s; feinster Deckelweg und Druckabfall je <1% Abweichung; vier Fluid-/Strukturnetzstufen, letzte Änderung <0.2%, alle Kopplungsgates |

Die SA-MMS verwendet einen expliziten Test-Wandabstand, nicht den Abstand zu
allen Dirichletgrenzen; sie ersetzt den physischen Wandabstandstest nicht.
Die DNS-Daten sind eine unabhängige Referenz zur Modellbeurteilung, keine
DNS-Implementierung und keine experimentelle Validierung. Re_tau=587.19 (Datensatzname 590) wird
über delta, nu und Druckantrieb festgelegt; u_tau=sqrt(f_x delta/rho)=1 m/s.
Die [konkrete Profilquelle](https://turbulence.oden.utexas.edu/data/MKM/chan590/profiles/chan590.means)
nennt Re_tau=587.19 und die Normierung U_tau/h; ihr Rundname 590 wird nicht
als exakter Parameter übernommen. Fluss ist Ergebnis, kein gleichzeitig
vorgeschriebener Zusatz. Der 1D-Solver
ist ein anderer numerischer Code, aber dieselbe Modellgleichung; deshalb
separater Codevergleich und separater DNS-Modellvergleich. Alle Daten werden
auf denselben Wandabstand interpoliert, ohne Fits an SA-Konstanten.

Der FSI-Deckel ist ein geführter starrer Lastverteiler auf einem linearen
elastischen Träger; die Struktursteifigkeit k_s wird zuvor unabhängig geprüft.
Der vollständig entwickelte laminare Kanal nutzt einen parabolischen Einlass
mit festem Gesamt-Q und Druckauslass. Gleichmäßige Deckelverschiebung verändert
H=H0+d. Druckkräfte liefern F=b L Delta_p/2, Delta_p=12 mu L Q/(b H³).
Damit folgt die angegebene nichtlineare Wurzelgleichung. Der Deckel verteilt
vertikale Lasten; horizontale Scherkräfte werden über Führungsreaktionen
bilanziert. Dies verifiziert bidirektionale Geometrie-/Lastkopplung im laminaren
Grenzfall. Zusätzlich verpflichtend für #92: ein SA-Deckel mit H0=0.1 m, L=b=1 m,
rho=1 kg/m³, nu=(0.05 m)(0.01 m/s)/587.19 und k_s=1 N/m. Festes Q ist
der vorab berechnete Fluss der unabhängigen 1D-SA-Lösung bei H0 und
f_x=0.002 N/m³ (Re_tau=587.19); Wert und Referenzhash vor dem gekoppelten
Lauf festhalten. Bei jeder neuen Höhe bestimmt der unabhängige BVP die
Druckdifferenz bei diesem festen Q. Weg-/Lastabweichung <2%, letzte
Netzänderung <0.4%; dieselben Deformations- und y+-Gates bleiben Pflicht. Ein laminarer Deckeltest allein gibt RANS-FSI nicht frei.
#93 liefert beide reproduzierbaren Rezepte samt tatsächlichen Gates.

## 8. Algebraische Plausibilität und Umsetzungsgates

Vorzeichenbeispiel: Fluid unter horizontaler Wand, n_f=(0,1), p=2 Pa,
S=0, Segmentlänge 3 m und b=4 m. Fluidkraft auf Wand ist
-sigma_f n_f A=(0,24) N. Auf Struktur nach oben; Gegenkraft auf Fluid
(0,-24) N. Ein ausgelassener Breitenfaktor würde fälschlich 6 N liefern.

Mappingbeispiel: H=((1,0),(1/2,1/2)), d_s=(0,2) m und
F_q=(3,5) N. H d_s=(0,1) m, H^T F_q=(5.5,2.5) N;
Gesamtkraft 8 N und virtuelle Arbeit 5 J auf beiden Seiten. Bei
x_s=(0,2) m und x_q=(0,1) m ist auch das Lastmoment 5 N m gleich.
Diese Rechnungen sind unabhängig von einem CFD-/Struktursolver.

| Gate | Freigabevoraussetzung |
| --- | --- |
| #84 Vertragsabnahme | Gleichungen/Einheiten, Quellen und algebraische Beispiele geprüft; unabhängiges fachliches Review und anwendbare PR-CI vor Integration |
| #85 Netz | Import-/Geometrie-/Wandabstand-TDD und vollständige Randgruppen |
| #226 Quad-FV-Operator | SIMPLE/Rhie-Chow-Stabilitäts- und Patchnachweise, variable Spannungsterme und Vier-Netz-Vergleiche; uniforme MAC-Regression erhalten |
| #87 Schließung | SA-TDD und variable Spannungsterme, Positivitäts-/Nichtkonvergenzdiagnosen, getrennte Residuen; kein Validierungsanspruch |
| #88 CFD-/Modellnachweis | alle obigen CFD-/SA-/DNS-Gates, gehashte Quellen und Nutzungsrechte; unabhängiger Referenzvergleich |
| #89/#90 Übertragung/Bewegung | Kräfte/Momente/Arbeit, vollständige Schnittstelle, stationäre Netzgates und Wiederaufbau der Geometrie |
| #91/#92 bidirektionale FSI | innerer und äußerer Abbruch, laminare und SA-Deckelreferenz, Ressourcen-/Fehlerdiagnosen; Einwegkopplung reicht nicht |
| #93/#94 Workflow/Release | dokumentierte API/Rezepte/Artefakte, aktuelle CI und Review; Daten-/Numerik-/Geometriegates sämtlich erfüllt |

#85 liefert Netz/Import ohne CFD-Solver. #226 liefert den gesonderten
Quad-FV-Operator mit Druckkopplung; reine Importdaten geben #87 nicht frei.

Umsetzungspräzisierung vom 05.10.2026 (ADR-051): #226 ist in Grundoperatoren
#229 und SIMPLE-/Strömungsabnahme #230 geteilt. Alle obigen Gates und die
Abhängigkeit von #87 bleiben bestehen; ein integrierter #229-Patch allein
schließt #226 nicht. [Operatorvertrag und separate Rekonstruktionsprüfung](73_QUAD_FV_OPERATORS.md).
#90 liefert Netzbewegung für dieselbe FV-Geometrie, ohne MAC-Netze
stillschweigend zu verformen. Reihenfolge: #84 -> #85 -> #226 -> #87 -> #88;
#89 nach #84/#85, #90 nach #85/#226, #91 nach #88/#89/#90,
#92 -> #93 -> #94. Kein M11-Optimierer als Voraussetzung.

Der langfristige M12-Milestone nennt 3D-Strömung. Die hier begrenzte ebene
Abnahme hebt dieses Ziel nicht auf: #94 darf damit keinen vollständigen
3D-Milestone behaupten. Verbindliche Fortsetzung vor einer solchen Freigabe:
strukturierte konvexe Hex8-Fluidzellen und Quad4-Außenflächen, dieselbe
kollokierte FV-/SA-Form in drei Raumrichtungen, Wandabstand zu triangulierten
Flächen und unbeschränkte 3D-P1-Schnittstelleninterpolation. Alle Gleichungen
und geometrischen Kriterien gelten dann dimensionsentsprechend; Face-Fläche,
Zellvolumen, drei Impulsresiduen und Momente um alle drei Achsen ersetzen die
Extrusionsreduktion. Kein Tet-CFD und keine neue Modellfamilie.
Die Netz-/Operator-/FSI-Issues müssen diese zweite Stufe vor #94 separat
abnehmen. Zusätzliche verpflichtende 3D-Referenz: divergencefreies
U=curl(0,0,psi) mit psi=sin²(pi x)sin²(pi y)sin²(pi z),
p=sin(2pi x)sin(2pi y)sin(2pi z) auf [0,1]³, rho=0.1 und mu=1 in SI;
analytisch hergeleitete Quellen, vier Netze 8/16/32/64, feinste RMS-U-/p-Fehler
<0.02 m/s bzw. Pa, Ordnung >0.8 mit begrenztem Upwind. Zusätzlich SA-MMS mit
a=nu[1+sin(pi x)sin(pi y)sin(pi z)], virtuelle Wandebene z=-1,
und denselben SA-Gates wie oben, sowie extrudierter SA-Deckel in echten
3D-Fluidelementen: Ergebnisse innerhalb 2% der unabhängigen ebenen Referenz,
alle drei Kräfte/Momente und Arbeit bilanziert. Diese sind noch offene
Implementierungs-/Nachweisgates, keine vom 2D-Code geerbte Freigabe.

Quellenprüfung am 04.10.2026: TMR ist von turbmodels.larc.nasa.gov nach
[tmbwg.github.io/turbmodels](https://tmbwg.github.io/turbmodels/) umgezogen.
[Die NASA bestätigt den Umzug](https://www.nasa.gov/nasa-turbulence-modeling-resource/).
Modellgleichungen und Quellen dienen der Methodenwahl; keine externe Quelle
behauptet die hier projektspezifisch gesetzten Fehlerbudgets oder deren Bestehen.


## 9. Ausgeführte Vertragsprüfung (04.10.2026)

- Primärquellen für SA-Variante/S-Modifikation, FV-Druckkopplung,
  konservative Übertragung und DNS-Metadaten gelesen; aktueller TMR-Umzug
  berücksichtigt. Weitergaberechte der Referenzdaten werden nicht behauptet.
- Dimensionen der Impuls-/SA-Terme und Rand-/Lastvorzeichen manuell geprüft.
  Unabhängige NumPy-Algebra bestätigt 8 N, 5 N m und 5 J im Mappingbeispiel.
  Druck-/Breitenbeispiel liefert 24 N mit korrekter Gegenkraft.
- Positive laminare Deckelwurzel unabhängig mit SciPy Brent berechnet:
  d=0.000589512662694 m und Delta_p=11.7902532539 Pa;
  Feder-/Druckkräfte stimmen mit absolutem Rest <1e-10 N.
- Lokale Dokumentlinks und Querverweise geprüft. Es wurde kein ausführbares
  Verhalten geändert; ein künstlicher Red-Test ist nicht anwendbar.
- 62 bestehende CFD-Unit-/Beispieltests bestanden mit `PYTHONPATH=src
  .venv/bin/python -m pytest tests/unit/test_cfd_channel.py
  tests/unit/test_cfd_stokes.py tests/unit/test_cfd_convection.py
  tests/unit/test_cfd_flow.py tests/unit/test_cfd_flow_config.py
  tests/regression/test_cfd_flow_example.py -q`. Ruff, Formatprüfung, mypy
  und `git diff --check` bestanden. Keine numerischen Operatoren geändert;
  keine erneute vollständige Konvergenz-/Gesamtsuite für diesen Dokumentationsdiff.
- Unabhängiges fachliches Review steht aus. Die hier berechneten algebraischen
  Beispiele sind keine bestandenen RANS-/FSI-Implementierungsgates.
