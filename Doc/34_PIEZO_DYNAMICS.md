# Lineare Piezoelastodynamik — #31

Status: 31a–d (Modalanalyse, Newmark, Harmonik, Model/Export) in v1.1.0.
[Release-Nachweis und Nutzungsgrenzen](35_M6_RELEASE_VERIFICATION.md).
Keine Rückwirkung der Verformung auf die Temperatur, keine experimentelle Validierung.

## DAE- und Elektrodenvertrag

Mit der stationären Vorzeichenkonvention gilt

`M ü + K u + B phi = f_u`, `Bᵀ u - D phi = f_phi`.

Nur mechanische Freiheitsgrade tragen Masse. Nach Dirichlet-Elimination und
Zusammenfassung äquipotentieller Elektrodenknoten ist `D` positiv definit.
Der elektrische Anteil ist algebraisch: `phi = D⁻¹ (Bᵀ u - f_phi)`.
Damit ist `K_eff = K + B D⁻¹ Bᵀ`; es wird niemals ein Eigenwertproblem mit
singulärer elektrischer Masse gelöst. Die Tet4-Masse ist konsistent integriert.

Jede zusammenhängende Komponente benötigt eine explizite elektrische Referenz
und mechanische Lagerung gegen alle Starrkörpermoden. Es gibt keine automatische
Erdung. `Electrode(potential=0)` bezeichnet eine auf die Referenz kurzgeschlossene
Elektrode; `Electrode(charge=0)` eine offene äquipotentielle Elektrode mit
verschwindender **Gesamtladung**. Ungesetzte elektrische Flächen tragen `D.n=0`,
sind aber keine Leiter. Leiterladung ist minus die Summe der ursprünglichen
elektrischen Residuen. Die Kondensation verwendet denselben Elektrodenvertrag
wie der stationäre Solver; überlappende Elektroden werden abgelehnt.

Modalbedingungen sind homogene Störungen um den linearen Referenzzustand.
Nichtnullwerte werden ausdrücklich abgelehnt. Ein statischer Bias ist separat
zu lösen und verändert ohne Vorspannungs-/Nichtlinearitätsmodell die Moden nicht.
Die Zeitintegration nimmt nur mechanische Anfangsauslenkung und
Geschwindigkeit entgegen. Potential und Beschleunigung sind bei `t0` aus der
algebraischen Gleichung bzw. dem mechanischen Gleichgewicht zu bestimmen;
elektrische Anfangswerte sind keine unabhängigen Zustände. Bewegte mechanische
Dirichlet-Ränder müssen Verschiebung, Geschwindigkeit und Beschleunigung
konsistent vorgeben; verschobene Trägheitsterme gehören zur reduzierten Last.

## Python-API

`solve_piezoelectric_modes(mesh, material, density, conditions, count=6,
electrodes=None)` aus `fem_toolbox.solvers.piezo_dynamics` nimmt ein homogenes
Material oder eines pro Tet4 sowie positive Dichte in kg/m³ (Skalar oder pro
Zelle) entgegen. DOFs kommen aus `piezoelectric_dofs(mesh)`.

Die niedrigsten `count` Moden werden auf `uᵀ M u=1` normiert. Der größte
mechanische Eintrag wird positiv gewählt; bei mehrfachen Eigenwerten ist die
Basis des Eigenraums nicht eindeutig. Ergebnisarrays sind unveränderliche
Kopien: Verschiebungsformen `(nmode,nactive,3)`, Potentialformen
`(nmode,nactive)`, Kreisfrequenz in rad/s und Frequenz in Hz. Formen sind keine
physisch angeregten Amplituden. Elektrodenspannungen/-ladungen, feldweise
relative Rückwärtsfehler und deterministische Provenienz werden mitgeführt.

Der dichte Schur-Referenzsolver begrenzt die Größe **vor** der Matrixallokation
auf 2048 aktive Gesamt-DOFs. Das ist ein transparenter Forschungsgrenzwert,
keine Skalierbarkeitszusage. Größere Probleme benötigen einen späteren dünnbesetzten
Operator-/Backendpfad. Starrkörpermoden und rein elektrisch/fest gehaltene
Modelle gehören nicht zu dieser ersten Modalanalyse.

## Unabhängige Verifikation

- Einheits-Tet4 mit nur einem freien mechanischen und einem elektrischen DOF:
  `M=rho/60`, `K=C/6`, `B=e/6`, `D=kappa/6`. Kurzschluss:
  `omega²=(C/6)/(rho/60)`; offen: `omega²=((C+e²/kappa)/6)/(rho/60)`.
  Toleranz `2e-14` relativ für das gut konditionierte skalare Problem.
- Zwei entkoppelte mechanische DOFs: eigenständige diagonale Frequenzreferenz
  und Massennormierung/-orthogonalität, Toleranz `3e-14`.
- Leitfähiger d33-Stab, `L=0.02 m`, Querschnitt `0.01² m²`,
  `C=1e9 Pa`, `e=5 C/m²`, `kappa=1e-8 F/m`, `rho=7500 kg/m³`.
  **Synthetische Verifikationsparameter**, keine Hersteller-Materialdaten.
  Querverschiebungen sind fest; `u(0)=0`, `T(L)=0`. Aus `D'=0` und
  `rho ü=T'` folgt mit `C_D=C+e²/kappa` für offene Elektroden
  `omega=(pi/2L)*sqrt(C_D/rho)`. Kurzschluss setzt nur die gesamte Spannung
  auf null: `lambda*cot(lambda)=e²/(kappa*C_D)`,
  `omega=lambda/L*sqrt(C_D/rho)`. `C` allein wäre für diesen nichtaffinen
  dynamischen Modus die falsche Kurzschlussreferenz.
- Uniforme Netze mit 1, 2 und 4 Unterteilungen pro Richtung müssen für beide
  Fälle fallende Raumfehler und unter 2 % Fehler auf dem feinsten Netz zeigen.
  Elektrische Gesamtladung der offenen Elektrode schließt bis auf Rundungsfehler
  (relativ `1e-10` zur modalen Spannung mal geometrischer Kapazität).
- Negative Dichte, fehlende Referenz, Starrkörpermoden, unzulässige Modenzahl,
  nichtnull Randwerte, vollständig feste Mechanik und Größenlimit werden geprüft.

Die Modalanalyse besitzt kein Zeit- oder Frequenzraster; ihr Raumfehler wird
nicht mit Integrations- oder Abtastfehlern vermischt. Allgemeine Feldresiduen
müssen kleiner als `1e-9` sein. Diese sind ergänzende Kontrollen und ersetzen
keine unabhängige Frequenzreferenz.

## Newmark-Zeitintegration (31b)

`solve_transient_piezoelectric` aus `fem_toolbox.solvers.piezo_transient`
verwendet `beta=1/4`, `gamma=1/2`. Physikalische und numerische Dämpfung sind
jeweils null; es gibt keine versteckten Dämpfungsparameter. Auf einem festen
Raumnetz ist die Integration für glatte Daten zeitlich zweiter Ordnung.

Argumente entsprechen der Modalanalyse, ergänzt um `times` in Sekunden,
`initial_displacement`, `initial_velocity` (jeweils `(nactive,3)`), optional
`loads(t)` und `prescribed(t)`. Der Lastcallback liefert einen Vektor
`[f_u,f_phi]` der Länge `4*nactive` mit **bereits integrierten** Kräften (N)
und vorzeichenbehafteten elektrischen Lasten (C). Für eine Gesamtleiterladung
`Q` wird `-Q` auf einem Potential-DOF der offenen Elektrode eingetragen;
Beiträge auf verbundenen Elektrodenknoten summieren sich. Diese Schnittstelle
ersetzt keine Volumen-/Flächenquadratur. Der Model-Adapter wird das übernehmen.

`prescribed(t)` liefert `KinematicValues(u,v,a,phi)` in aktiver Knotenreihenfolge,
mit Nulleinträgen außerhalb der vorab festgelegten Rand-DOFs. Ableitungen müssen
vom Aufrufer konsistent vorgegeben werden. Verschiebungs- und Geschwindigkeitssprünge
sind nicht unterstützt; Beschleunigungs- und Spannungssprünge sind möglich.
Die Lastreduktion enthält ausdrücklich `-M_fc*a_c`. Fehlende Anfangswerte sind
null auf freien DOFs und entsprechen den vorgeschriebenen Werten am Rand;
explizite Anfangswerte müssen die Randdaten bis `1e-12 m` bzw. `m/s` erfüllen.

Schritt-, Rampen- und CSV-Kurven aus `fem_toolbox.time` können Lastcallbacks
speisen. Ihre Knoten/Ereignisse sind als `events=signal.knots` zu übergeben,
bei mehreren Signalen als zusammengeführte Liste. Ereignisse werden ins Zeitraster
aufgenommen. Ein Intervall endet beim linken Grenzwert des Ereignisses; danach
werden Potential und Beschleunigung mit den rechtsstetigen Daten neu bestimmt.
Verschiebung und Geschwindigkeit bleiben kontinuierlich. Nicht angegebene
Sprungzeiten können nicht zuverlässig erkannt werden.

Das Ergebnis enthält vollständige `u/v/a/phi`-Zeitverläufe, ursprüngliche
Residuen einschließlich mechanischer Trägheit, feldweise Rückwärtsfehler,
positive Gesamtenergie und Provenienz aus den tatsächlich ausgewerteten Daten.
Bei homogenen, ungetriebenen Rändern gilt
`E=0.5*vᵀMv+0.5*uᵀKuu*u+0.5*phiᵀD*phi=const`. Bei Anregung ist diese Energie
nicht konstant; es wird keine unberechnete Arbeitsbilanz behauptet. Die
Speichergrenze ist `ntime * 4*nactive <= 2_000_000`; Streaming folgt separat.

Nachweise: freie skalare Schwingung erhält Energie relativ auf `3e-12` und
zeigt bei 40/80/160 Schritten Fehlerquotienten zwischen 3.8 und 4.2.
Zeitverfeinerung erfolgt auf **demselben** Ein-DOF-Raummodell. Kraft- und
Spannungssprünge erzeugen keine vorzeitige Auslenkung/Geschwindigkeit.
Eine quadratische Starrtranslation mit passender Volumenlast prüft bewegte
Ränder, Anfangsbeschleunigung und Massenterme bis `3e-14` absolut in SI.
Ein vorbelasteter statischer Zustand bleibt stationär, einschließlich
nichtverschwindender elektrischer Ladung. Diese Kontrollen sind numerische
Verifikation, keine Messdatenvalidierung.

## Harmonische Analyse (31c)

`solve_harmonic_piezoelectric(..., frequency_hz, loads=..., prescribed=...)`
aus `fem_toolbox.solvers.piezo_harmonic` löst eine Frequenz einschließlich
DC (`frequency_hz=0`). Konvention: `x(t)=Re(x_hat*exp(+i*omega*t))`,
`omega=2*pi*frequency_hz`. Amplituden sind Spitzenwerte, keine RMS-Werte;
`v_hat=i*omega*u_hat`, `a_hat=-omega²*u_hat`.

Last- und Randwertvektoren haben Länge `4*nactive` in derselben Reihenfolge
wie Newmark, dürfen aber komplex sein. Die additive komplexe Infrastruktur
ändert **nicht** den reellen Vertrag von `LinearSystem`. Imaginärteile werden
weder durch Casts noch beim Rückrechnen elektrischer DOFs verworfen.
Die dynamische Steifigkeit ist `K_eff-omega²*M`; bewegte Ränder tragen
`+omega²*M_fc*u_c` zur rechten Seite bei.

Es gibt keine Dämpfung. Zur Resonanzdiagnose wird die Matrix beidseitig mit
`diag(K_eff+omega²*M)^(-1/2)` skaliert. Beträge ihrer dimensionslosen Eigenwerte
unter oder gleich `1e-10` führen zu einem ausdrücklichen Resonanz-/Konditionsfehler,
auch bei verschwindender Last. Das ist eine numerische Nutzungsgrenze, keine
physikalische Resonanzbreite. Die Funktion liefert keine künstlich begrenzten
Resonanzamplituden. Ein Frequenzsweep enthält nur die ausdrücklich berechneten
Frequenzen und findet dazwischenliegende Resonanzen nicht automatisch.

Verifikation: komplexer skalarer Oszillator ober- und unterhalb der Eigenfrequenz,
DC-Grenzfall, komplexe Spannung, bewegter Rand mit konsistenter Masse, negative/
unrepräsentierbare Frequenzen und explizite Resonanz. Mit harmonisch passenden
Anfangswerten konvergiert Newmark auf einem festen Raummodell quadratisch gegen
das harmonische Ergebnis (100/200 Schritte). Ein Frequenzrasterfehler ist davon
getrennt: hier wird eine einzelne, identische Frequenz in beiden Pfaden verglichen.

## Model-Schema 3 und Export (31d)

`Model(mesh=..., document=...)` und `Model.from_file(...)` akzeptieren additiv
`schema_version=3`, `units="SI"` und `physics` gleich `piezo_modes`,
`transient_piezoelectric` oder `harmonic_piezoelectric`. Schema 1 und 2 behalten
ihre Bedeutung. Ein vollständiges TOML-/CSV-Beispiel mit ausführbarem Skript
steht unter [examples/piezo_dynamics](../examples/piezo_dynamics/README.md).

Das Dokument enthält:

- `base`: ein vollständiges Schema-1-Modell mit `physics="piezoelectric"`,
  Material-/Regionszuordnung, **homogenen** mechanischen/elektrischen Randwerten
  und optionalen Volumen-/Flächenlasten. Ein Netzpfad steht nur im Schema-3-Root.
  Der stationäre Referenzpfad bereitet das integrierte Lastmuster vor. Dessen
  Löseroptionen gelten nur für diese Vorbereitung; dynamische Verfahren sind
  die hier festgelegten Referenzen. Modalanalyse verlangt ein lastfreies Basismodell.
- `density={value=..., unit="kg/m^3"}`: positive Dichte, homogen oder pro Tet4.
- `analysis` mit den nachfolgend beschriebenen Optionen; unbekannte Schlüssel
  und abweichende Einheiten werden zurückgewiesen.

Modal: `count` legt die Anzahl Moden fest. Transient: `times={values=[...],unit="s"}`
und optional `load_factor` als Zahl, Step/Ramp/Series oder CSV-Profil skalieren das
Basismuster. CSV-Proben werden bei `from_file` ins kanonische Dokument übernommen;
alle Profilknoten werden automatisch eingefügt. `initial_displacement` und
`initial_velocity` sind optional, in m bzw. m/s und aktiver Knotenreihenfolge.

`analysis.prescribed` ist eine Liste mit globalen `dofs` aus
`piezoelectric_dofs(mesh)`. Transiente Spannungen verwenden `unit="V"` und
`value` als Profil. Mechanische Bewegungen verwenden `unit="m"` und
`coefficients=[c0,c1,...]` für `u(t)=sum(cj*t^j)`; Geschwindigkeit und
Beschleunigung werden exakt abgeleitet. Die Koeffizienten haben entsprechend
Einheiten m/s^j. Die DOFs müssen in der Basistopologie bereits vorgeschrieben
sein. Die erste Modellversion erlaubt hier glatte Polynome; allgemeinere
kinematische Vorgaben stehen über die Python-Solver-API bereit.

Harmonisch: `frequencies={values=[...],unit="Hz"}` muss nichtnegativ und streng
aufsteigend sein. `load_factor={real=...,imag=...}` skaliert das Basismuster.
Vorgeschriebene Phasoren verwenden `dofs`, `unit` und `real`/`imag`. Für mehrere
unabhängige Lastmuster oder zeitlich veränderliche Leiter-Gesamtladung steht die
vollständige vorzeichenbehaftete Lastvektor-Schnittstelle der Solver zur Verfügung.

`export_piezo_series(path, model_result, probes={...})` aus
`fem_toolbox.io.piezo_series` schreibt in einen **neuen** Ordner:

- `piezo.pvd` und VTU-Dateien mit u/phi, bei transienter/harmonischer Analyse
  außerdem v/a und ursprünglichen mechanischen/elektrischen Residuen.
- Tet4-Zellfelder: Dehnung, Spannung, elektrisches Feld und elektrische
  Verschiebung, ohne Glättung. Inaktive Knoten und Dreiecks-Zellfelder sind
  nullgefüllt und mit Gültigkeitsmasken versehen.
- Harmonische Felder: `_real`, `_imag`, `_amplitude`, `_phase_rad`,
  `_phase_valid`. Amplitude/Phase sind komponentenweise; bei Amplitude null
  wird Phase null mit ungültiger Maske gespeichert. Reelle Zeit-/Modalfelder
  erhalten nur `_real`, ohne künstliche Phase.
- `probes.csv` für benannte Knotenindizes oder `(x,y,z)` in m im Tet4-Gebiet.
  Komplexe Felder werden **vor** der Amplituden-/Phasenbildung interpoliert.
- `manifest.json` mit Achsenbedeutung, Hz/rad/s, SI-Einheiten, Energie pro
  Zeitrahmen, feldweisen Residuen, Datei-Hashes und vollständiger Model-Provenienz.

Der ParaView-Slider ist je nach Analyse Zeit in s, Frequenz in Hz oder
Modenindex, kein impliziter Wechsel der Einheit. Die Daten bleiben ans
berechnete Netz und die effektiven Zellmaterialien gebunden. Vorhandene
Verzeichnisse werden nicht überschrieben; bei Schreibabbruch bleiben PVD und
Manifest mit abgeschlossenen Frames und Status `incomplete` lesbar.

Nachweise: alle drei Model-Pfade, Einheit-/Schemafehler, CSV-Auflösung,
analytisch abgeleitete Lagerbewegung, VTU-Roundtrip einschließlich komplexer
Zellfelder, Messpunktinterpolation und kontrollierter Schreibabbruch.
CI führt das Beispiel zusätzlich aus frisch installierten Wheels und sdists
unter Python 3.11–3.13 aus.

## Tet10-Modalanalyse (#35d)

Konsistente Masse, offene/kurzgeschlossene Elektroden, Model und Modenexport
sind integriert. [Umfang und Verifikation](51_TET10_MODES.md).
Tet10-Harmonik und Gesamt-Abnahme-/Kostenvergleiche bleiben offen.

## Tet10-Newmark (#35d)

Zeitintegration mit konsistenter Masse, bewegten Rändern, Last-/Spannungsereignissen,
Model und Zeitserienexport ist integriert. [Nachweise und Grenzen](52_TET10_NEWMARK.md).
Tet10-Harmonik und Gesamt-Abnahme-/Kostenvergleiche bleiben offen.
