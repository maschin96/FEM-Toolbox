# M12 — SA-noft2-Implementierung und getrennte Abnahme (#87/#88)

Entwicklungsstand 20.10.2026, Paketkandidat `1.7.0.dev1`; keine M12-Freigabe.
Maßgeblich bleibt der [vorab festgelegte Vertrag](71_M12_RANS_FSI_CONTRACT.md).
Die aktuelle M12-Planung begrenzt die Freigabe auf ebene 2D-RANS und
spannweiteninvariante extrudierte Kleinverformungsstruktur. Die frühere
3D-Fortsetzungsplanung bleibt als spätere Erweiterung gekennzeichnet.

## Modell und Schnittstellen

`fem_toolbox.cfd.sa_model.sa_terms` liefert unveränderliche Zellfelder für
Wirbelviskosität, Diffusivität, Produktion, Destruktion, Gradientquadratterm,
implizite Destruktionsrate, modifizierte Wirbelstärke, r und f_w.
Eingaben sind a[ncell] in m²/s, grad(U)[ncell,2,2] in 1/s,
grad(a)[ncell,2] in m/s, positiver Zellwandabstand und positives molekulares nu.
Die Gradientenkonvention ist [Geschwindigkeitskomponente, Ableitungsrichtung].

SA-noft2 verwendet weder Trip noch f_t2. Die Allmaras-Modifikation ist
[TMR Note 1(c)](https://tmbwg.github.io/turbmodels/spalart.html) mit
c2=0.7, c3=0.9; S_tilde=0 setzt r=10. Das modellgemäße r-Limit ersetzt keine
Begrenzung von a oder nu_t. Negative, nichtendliche oder numerisch nicht
darstellbare Zustände werden konkret abgelehnt. Ein gemeinsamer relaxierter
Schritt wird höchstens zwanzigmal versucht und dabei halbiert; kein Zellwert
wird auf null geklemmt. Die molekulare Viskosität bleibt strikt positiv.

`solve_sa_transport` löst nur die skalare Gleichung für einen vorgegebenen
konservativen Face-Massenstrom und Geschwindigkeitsgradienten. Einlaufwerte
kommen aus FluidDomain; Wände setzen a=0, Symmetrie/Auslass haben Nullnormalgradient,
periodische Paare teilen denselben Transportfluss mit entgegengesetztem Vorzeichen.
Begrenztes lineares Upwind oder Donor-Cell wird explizit gewählt. Volle
Nichtorthogonalitätskorrektur, variable Diffusivität und alle nichtlinearen
Quellen werden nach jedem akzeptierten Schritt neu ausgewertet.

Optional explizite Rand-/Wandabstandsfelder dienen allein der Manufactured-Solution-
Verifikation und werden in der Provenienz so ausgewiesen. Sie ersetzen weder
den echten Abstand zu Wandsegmenten noch die y+-Prüfung eines RANS-Laufs.
Transport- und Diffusionsflüsse haben Einheit m⁵/s²; Zellresiduen m²/s².
Nichtkonvergenz bewahrt letzten Zustand, Minimum, Zell-ID und Residuenhistorie.

`solve_rans` koppelt U, p und a gemeinsam. Es verlangt ein explizites
nichtnegatives Startfeld mit mindestens einer positiven Zelle und positive
Uref/Lref-Skalen. Dieselbe konservative Rhie-Chow-Massenflussfolge transportiert
Impuls und SA. Die volle variable Newton-Spannung verwendet mu+rho nu_t.
Positivitäts-Backtracking verkleinert sämtliche Feldupdates zusammen.
Frische Residuen von Kontinuität, beiden Impulsen und SA sowie globale
Masse/Impuls/Leistung müssen gleichzeitig die unveränderten Vertragsgates
erfüllen. Erst danach wird y+ aus tatsächlicher tangentialer Wandtraktion
berechnet; y+>1 oder wandnormales Schichtwachstum>1.2 ergibt einen Fehler.
`RansConvergenceError` bewahrt die gemeinsame Historie und alle letzten Felder.

## Bisherige getrennte Prüfungen

- TDD: fehlende SA-Punkt-, Transport- und RANS-APIs zunächst Importfehler;
  danach Punkt-/Transport-/Kopplungstests bestanden. Analytische affine
  SA-Diffusion auf geschertem Netz, Einheiten, starre Rotation, Nullzustand,
  Allmaras-Zweig und Nullwirbelstärke sind getrennt geprüft.
- SA-MMS aus Doc/71: U=(1,0), a=nu[1+sin(pi x)sin(pi y)], d=y+1.
  Vor Auswertung festgelegte SI-Skalen Uref=Lref=1, nu=0.1 m²/s und b=0.7 m;
  analytisch unabhängig differenzierte Quellen, keine Aufrufe der produktiven
  Modellterme zur Erzeugung der Quelle. Vier Netze 16/32/64/128, lineares Upwind.
  Relative RMS-a-Fehler: 9.50617e-4 / 2.48245e-4 / 6.36585e-5 / 1.61304e-5.
  Ordnungen 1.937 / 1.963 / 1.981; unveränderte Gates <1e-2 und >0.8 erfüllt.
  Alle Zellen nichtnegativ, alle finalen SA-Residuen <=1e-9 m²/s².
- Der erste gemeinsame Kanaltest ist ein kontrollierter Low-Re-Grenzfall.
  Seine Konvergenz ist keine RANS-Modellvalidierung. Unabhängiger 1D-SA-Kanal,
  DNS-Modellvergleich, vier Wandnetze und Sensitivitäten gehören zu #88.

## Vorab fixierter Kanal-/Sensitivitätsplan für #88

Die folgenden Einstellungen werden vor der formalen Vier-Netz-Auswertung
fixiert. Die ersten 32-/64-Halbkanalläufe dienten der Solverdiagnose;
SA-Relaxation 0.5 erschöpfte auf dem 32er-Netz das 500er-Budget, während 0.1
konvergierte. Der Fehlversuch bleibt als No-go erhalten. Geändert wird die
explizite Iterationsrelaxation, keine Modellkonstante oder Abnahmegrenze.

- Voller Kanal 0<=y<=2 m, L=b=1 m, rho=1 kg/m³, nu=1/587.19 m²/s,
  Druckantrieb 1 N/m³. Sechs periodische x-Zellen; 32/64/128/256 Zellen
  je Halbkanal. Gemeinsame glatte Wandverteilung
  y(eta)=expm1(4.5 eta)/expm1(4.5), am Zentrum gespiegelt.
- Anfangsfeld a=3nu; U=p=0. Begrenztes lineares Upwind, alpha_U=alpha_p=0.5,
  alpha_a=0.1, Uref=20 m/s, Lref=1 m, maximal 500 äußere Schritte und
  500 SIMPLE-Schritte, maximal 50000 gemeinsame Unbekannte. Alle absoluten
  Residuen-, Bilanz-, y+- und Geometriegates aus Doc/71 bleiben bestehen.
- Unabhängiger kontinuierlicher 1D-SA-BVP, ohne Produktions-CFD-/SA-Importe,
  Kollokationstoleranz 1e-8; zweite Initialisierung und strengere Toleranz
  prüfen Referenzunabhängigkeit. Numerisches RMS-U relativ zu diesem BVP
  wird mit physischen Zellvolumina gewichtet. Feinster Fehler <0.5 %,
  Wandkraftfehler <0.1 %, wie im ursprünglichen Vertrag.
- DNS-Profilvergleich auf den ursprünglichen MKM-Wandabständen; ein kubischer
  Interpolant der numerischen Zellwerte verwendet U=0 an der Wand, die
  tatsächlich berechnete Wandschubspannung für die Wandableitung und die
  Symmetrieableitung null am Zentrum. Keine Fits von SA-Konstanten.
  Relative RMS-U+-Abweichung und relative Bulkabweichung jeweils <=5 %.
  Bulk wird aus physischer Zellquadratur bzw. Simpson-Integration des DNS-Profils
  berechnet. Letzte Netzverdopplung <1 % und höchstens 1/5 des Modellbudgets.
- Getrennte Wandsensitivität: bei 128 Zellen/Halbkanal beta=5.0 statt 4.5,
  alle unveränderten Residuen-/y+-/BVP-Gates; Größenänderungen separat berichten.
  Physische Parametersensitivität: nu +/-5 % bei 64 Zellen/Halbkanal, jeweils
  eigener unabhängiger BVP, unveränderte Modellkonstanten; kein zusätzlicher
  DNS-Validierungsanspruch außerhalb Re_tau=587.19.
- Wiederholung bei 64 Zellen/Halbkanal mit alpha_a=0.08 statt 0.1; dieselben
  physikalischen Gates, Vergleich der finalen Felder und Kennwerte.

Profiling des 32er-Kanals zeigte bei zehn Schritten 7.89 von 8.24 Solversekunden
im erneuten Aufbau geometrischer Maps. Wiederverwendung erfolgt deshalb nur
innerhalb derselben unveränderten FluidDomain. Sämtliche Viskositäts-, Spannungs-,
Impuls- und Residuenoperatoren werden weiter neu aufgebaut. TDD prüft den
unabhängig analytischen vollen Variable-mu-Spannungstensor und unveränderte alte
Iteranden. Der 64er-Lauf liefert vor/nach dieser Änderung identische U/a/p,
Wandtraktionen, y+ und komplette Residuenhistorien; Laufzeit 157 gegenüber
10.2 Sekunden. Netzbewegung verwendet stets eine neue Geometrie.

Der Viskositätstest mit nu-Faktor 0.95 zeigte anschließend einen weiteren
Solverfehler: Eine allein an der Impulsgleichung orientierte innere Toleranz
ließ das frisch rekonstruierte SA-Residuum bei ungefähr 1e-8 m²/s² stagnieren,
oberhalb des unveränderten Gates von ungefähr 3.2e-10 m²/s². Ein eigener
Regressionstest reproduziert die Nichtkonvergenz. Das innere Impulsbudget
berücksichtigt deshalb die SA-Amplitude:
`0.1 * äußeres Impulsgate * nu / max(nu,max(a))`.
Die Regel und ihre tatsächlichen Werte je Iteration werden protokolliert.
Dies ist ein konservatives inneres Fehlerbudget, keine allgemeine mathematische
Fehlerschranke; sämtliche äußeren Feld- und Bilanzgates bleiben unverändert.
Der Regressionstest besteht nach der Korrektur. Die fünf Kanal-/Wand-/Parameter-
prüfungen bestehen danach gemeinsam (265 s); die übrige lokale Gesamtsuite
liefert 1108 bestanden und 83 optionale PETSc-Tests übersprungen.

Der unabhängige BVP und dessen Skalierungs-/Initialisierungsprüfungen sowie
die vollständige Kanalnetzfolge sind versionierte Testquellen. Die eigene
DNS-Modellabnahme wird separat mit Originaldatenhash und Quellenmanifest
abgeschlossen; die BVP-Parität allein ersetzt diese Modellbeurteilung nicht.

Offen bleiben die vollständige CFD-/Modellabnahme #88, Schnittstellenübertragung,
Netzbewegung, bidirektionale FSI und deren Referenz-/Workflow-/Releasegates.
Unabhängiges menschliches Fachreview liegt nicht vor und wird nicht behauptet.
