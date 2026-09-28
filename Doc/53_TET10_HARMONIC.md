# Tet10-Harmonik (#35d)

`solve_harmonic_piezoelectric` und `Model`-Schema 3 mit
`physics="harmonic_piezoelectric"` unterstützen geradkantige Tet10, einschließlich
DC, komplexer Kraft-/Ladungsanregung, vorgeschriebener Spannung und bewegter Ränder.
Die Konvention bleibt `Re(x_hat exp(+i omega t))`; Amplituden sind Spitzenwerte,
v=iωu und a=-ω²u. Es gibt keine physikalische oder numerische Dämpfung.

## Operator- und Ausgabe-Vertrag

Die konsistente mechanische Masse und elektrische Schur-Elimination stammen aus
dem gemeinsamen Dynamikaufbau. Zur Rekonstruktion unter elektrischen Lasten
nutzt Harmonik dessen gespeicherten dielektrischen Block, wie Newmark.
Damit gilt dieselbe Tet10-Symmetrisierung gegen Quadratur-Rundungsreste.
Originalmatrizen bleiben für physikalische Residuen unverändert; Backend-
Toleranzen werden nicht gelockert. Elektrische Freiheitsgrade tragen keine Masse.

Offene Elektroden sind äquipotentiell mit vorgegebener Gesamtleiterladung Q;
die rechte Seite enthält -Q. Alle Triangle6-Knoten werden zusammengefasst.
Kurzschluss setzt das gesamte Elektrodenpotential auf null. Bewegte mechanische
Ränder tragen +ω² M_fc u_c zur reduzierten rechten Seite bei.

Tet10-Frequenzserien verwenden Manifest-Schema 2 und P2-Messpunktinterpolation.
u/v/a/φ sowie rohe Schwerpunktfelder erhalten Realteil, Imaginärteil, Betrag,
Phase in rad und Gültigkeitsmasken. Die Koordinaten der Schwerpunktwerte stehen
im Manifest. Phasen werden erst nach komplexer Interpolation berechnet.
Keine Glättung oder Projektion; die VTU-Flächen erhalten die bestehenden Nullwerte
und Masken. Der Tet4-Export bleibt unverändert.

## TDD und unabhängige Referenzen

Fünf neue Tests scheiterten zunächst an den direkten bzw. Model-Tet10-Sperren.
Nach Freischaltung bestand der skalare Pfad; der Mehrzellen-Stabtest scheiterte
an der SPD-Symmetrieprüfung des separat rekonstruierten dielektrischen Blocks.
Die Wiederverwendung des gemeinsamen Blocks behebt diesen Fehler. Der Beispieltest
scheiterte vor dem Anlegen des Skripts erwartungsgemäß. Anschließend bestanden
die gezielten Tests; keine Referenztoleranzen wurden zum Bestehen aufgeweitet.

### Skalarer Eckoszillator

Einheits-Rechtstetraeder, freie u₁ₓ/φ₁, ρ=6, Cxx=10, exxx=2, κxx=3:
∫N₁²=1/420, ∫(∂N₁/∂x)²=1/10. Daher M=6/420, K=1, B=0,2, D=0,3 und

```text
q = (f_u + (2/3) f_phi) / ((10+4/3)/10 - omega² (6/420))
phi = (0.2 q - f_phi) / 0.3.
```

f_u=1+2i N und f_phi=0,01-0,02i C werden bei 0, 0,2 und 2 Hz geprüft, also
auch oberhalb der Eigenfrequenz. Relative Toleranz 3e-14 für die kleine
algebraische Referenz. Ein weiterer Fall kombiniert komplexe Spannung mit
bewegtem y-Achsen-Eckknoten in x-Richtung: dessen Kreuzmasse zum freien
x-Eckknoten ist ρ∫N₁N₂=6/2520; die Normalsteifigkeitskopplung ist null.
Resonanz wird auch bei Last null explizit abgelehnt.

Der komplexe DC-Stab verwendet die unabhängigen affinen Felder aus
[Doc 52](52_TET10_NEWMARK.md), multipliziert mit 1+2i. Die Leiterladung und
Triangle6-Traktion sind separat hergeleitet. Toleranzen: u 2e-16 m,
φ 5e-11 V und summierte freie Elektrodenresiduen <1e-18 C.

### Analytisch erzwungener Stab auf vier Netzen

L=0,02 m, Querschnitt A=1e-4 m², ρ=7500 kg/m³, C=1e9 Pa, e=5 C/m²,
κ=1e-8 F/m, CD=C+e²/κ. Unteres Ende fest/geerdet, Querverschiebungen blockiert.
Eine obere harmonische Normalspannung τ=1000+2000i Pa treibt den Stab.
Mit k=ω√(ρ/CD), λ=kL=0,6 gilt u(z)=a sin(kz):

```text
offen (Q=0): a = tau / (CD k cos(lambda))
kurzgeschlossen: a = tau / (CD k cos(lambda) - e²/(kappa L) sin(lambda)).
```

Dies folgt aus konstantem D, T=CD u'-(e/κ)D und φ'=(e u'-D)/κ.
Die obere Traktion wird unabhängig mit Triangle6-Gewichten integriert:
0 an Ecken und Flächeninhalt/3 an Kantenknoten. Dieselben geometrischen Gewichte
messen die gemittelte Endverschiebung, verglichen mit a sin(λ).

| Axiale Teilungen | relativer Fehler offen | kurzgeschlossen |
|---:|---:|---:|
| 1 | 5,47640e-5 | 2,95071e-4 |
| 2 | 4,94529e-6 | 2,66513e-5 |
| 4 | 3,57523e-7 | 1,92681e-6 |
| 8 | 2,32301e-8 | 1,25195e-7 |

Letzte Ordnung jeweils 3,944. Der Test fordert streng sinkende Fehler,
letzte zwei Ordnungen 3,5–4,5 und feinsten Fehler <2e-6. Es ist eine axiale
Verfeinerung und eine gemittelte Randantwort, keine allgemeine H1/L2-Zusage.
Äquipotentialität und freie Elektrodenladung werden separat geprüft; Ladung
relativ zur Kapazität κA/L mal Spannung mit 1e-10, Feldresiduen <1e-9.

### Zeit-, Export- und Backend-Prüfung

Mit harmonisch konsistenten Anfangsdaten konvergiert Newmark bei 100/200
Schritten auf unverändertem Raumnetz gegen die komplexe Antwort (Fehlerquotient
3,8–4,2, feinster Fehler <1e-6 m). Raum-, Zeit- und Frequenzraster sind getrennt;
der Solver berechnet nur angegebene Frequenzen.

Model-Sweep und Rundreise prüfen die negative Tet10-Eckkraft -0,05(1+2i) N
einer konstanten Volumenlast. Ein komplex vorgeschriebener Kantenknoten liefert
nichtverschwindende, unabhängig berechnete Schwerpunktspannungen/-felder.
Messpunkte prüfen P2-Werte sowie komplexe Amplitude/Phase; keine Interpolation
bereits getrennter Beträge und Phasen. Absolute Feldtoleranz 2e-12 in den
jeweiligen synthetischen SI-Skalen deckt Rundung ab.

Native PETSc-direkt/FieldSplit-Parität wird in real/complex geprüft:
Stabantworten offen/kurzgeschlossen (rtol=2e-8, u atol=1e-15 m,
φ atol=1e-9 V) und geladener komplexer DC-Stab (u atol=2e-14 m,
φ atol=2e-9 V). Beide feldweisen Residuen bleiben <1e-9.
Alle Referenzen sind Code-Verifikation, keine experimentelle Validierung.

## Reproduktion und Grenzen

```sh
PYTHONPATH=src python examples/tet10_harmonic/run.py
PYTHONPATH=src python -m pytest tests/unit/test_tet10_harmonic.py tests/convergence/test_tet10_harmonic_convergence.py
```

Das eigenständige Gmsh/TOML-Beispiel prüft drei Frequenzen und exportiert VTU/PVD,
Manifest und Messpunkt-CSV. Die CI führt es aus Wheel und sdist aus.
2048 aktive Gesamt-DOFs, geradkantige Tet10 und die unveränderte dimensionslose
Resonanzschwelle 1e-10 bleiben Grenzen. Kein gedämpfter Resonanzsolver und keine
automatische Resonanzsuche zwischen Frequenzpunkten. #35 bleibt bis zum
Gesamtabgleich und zu den Genauigkeits-/Kostenvergleichen mit #37 offen.
