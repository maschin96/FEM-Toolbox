# Lineare Piezoelastodynamik — #31

Entwicklungsstand: 31a (Modalanalyse), additive Forschungs-API, kein Release.
31b (Newmark), 31c (Harmonik) und 31d (Model/Export) folgen in separaten PRs.
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
Die zukünftige Zeitintegration nimmt nur mechanische Anfangsauslenkung und
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
