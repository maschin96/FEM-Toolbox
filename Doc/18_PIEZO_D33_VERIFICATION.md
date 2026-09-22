# Direkter und inverser d₃₃-Nachweis (#22)

Status: Code- und Lösungsverifikation für einen bewusst eingeschränkten
quasistatischen d₃₃-Block; keine experimentelle Validierung und keine
allgemeine Floating Electrode.

## Modell, Einheiten und Randbedingungen

Ein Block mit Länge `L=0,02 m` und Querschnitt `A=10⁻⁴ m²` liegt in
Polarisationsrichtung `z`. Für den homogenen Benchmark gelten
`C₃₃=10⁹ Pa`, `e₃₃=5 C/m²`, `κ₃₃=10⁻⁸ F/m`. Alle übrigen
Kopplungseinträge sind null; die Steifigkeit und Permittivität sind diagonal
und positiv. Die analytische 1D-Reduktion lautet

```text
ε = du_z/dz,  E = -dφ/dz
σ = C₃₃ ε - e₃₃ E = C₃₃ u'_z + e₃₃ φ'
D = e₃₃ ε + κ₃₃ E = e₃₃ u'_z - κ₃₃ φ'
dσ/dz = 0,  dD/dz = 0
```

Im 3D-Modell werden `u_x=u_y=0` an allen Knoten gesetzt und `u_z=0`
unten. Für diese diagonal gewählten Tensoren entstehen dabei keine
transversalen Spannungen; die Seitenflächen sind tractionsfrei und haben
`D·n=0`. Die untere Elektrode hat `φ=0`. Kraft `F` bedeutet gleichmäßig
verteilte Traktion `F/A` oben. Ein Randpotential wird an allen Knoten der
oberen Fläche vorgegeben. `Q_top` ist die Ladung der angrenzenden Elektrode,
also das Negative der elektrischen Dirichlet-Reaktion aus dem vorzeichenbehafteten
Solver. Diese Vorzeichen folgen [Doc/12](12_ELECTROSTATICS.md) und
[Doc/17](17_PIEZOELECTRIC_SOLVER.md).

## Direkte Wirkung und Kurzschluss

Bei freier Oberseite und `φ_top=V=100 V` folgt aus `σ=0`

```text
d₃₃,eff = e₃₃/C₃₃ = 5·10⁻⁹ m/V
u_top = -d₃₃,eff V = -5·10⁻⁷ m
Q_top = A/L (κ₃₃ + e₃₃²/C₃₃) V = 1,75·10⁻⁸ C
H = -1/2 Q_top V = -8,75·10⁻⁷ J
```

Der negative Wert ist **elektrische Enthalpie** unter vorgegebenem Potential,
keine negative gespeicherte Energie. Beim inversen Effekt werden beide
Elektroden auf `0 V` gehalten und oben `F=10 N` aufgebracht. Dann gilt

```text
u_top = F L/(A C₃₃) = 2·10⁻⁶ m
Q_top = -e₃₃ F/C₃₃ = -5·10⁻⁸ C
H = F u_top/2 = 10⁻⁵ J
```

Die Reziprozität wird unabhängig an den zwei Lastfällen geprüft:
`u_top/V = Q_top/F = -e₃₃/C₃₃`. Zusätzlich prüfen die Tests die
getrennten Enthalpien, Kraft- und Ladungsbilanzen. Die relativen
Benchmarktoleranzen betragen `2·10⁻¹¹`, weit oberhalb des erwarteten
Rundungsfehlers der kleinen Sparse-LU-Systeme und weit unterhalb einer
relevanten physikalischen Abweichung.

## Einfacher Leerlauf

Für eine einzige Potentialreferenz unten und **punktweise** `D·n=0` auf
allen übrigen Flächen gilt in diesem 1D-invarianten Block `D=0`. Mit
`F=10 N` erhält man

```text
C_open = C₃₃ + e₃₃²/κ₃₃ = 3,5·10⁹ Pa
u_top = F L/(A C_open) ≈ 5,7142857·10⁻⁷ m
φ_top = e₃₃ F L/(A κ₃₃ C_open) ≈ 285,7142857 V
H = F u_top/2 ≈ 2,8571429·10⁻⁶ J
```

Die obere Fläche wird hier allein wegen der 1D-Symmetrie equipotentiell.
In allgemeinen Geometrien ersetzt die punktweise Nullflussbedingung **keine**
leitende Floating Electrode mit nur insgesamt null Ladung. Deren allgemeines
Modell ist Aufgabe von M5. Der Leerlauf-Benchmark deckte einen numerischen
Randfall auf: Wenn die physikalische Nettoladung null ist, muss die Toleranz
der globalen Ladungsbilanz auch die ungekürzten elektrischen Matrixprodukte
berücksichtigen. Die Korrektur bleibt auf Rundungsfehlerordnung begrenzt.

## Zweite unabhängige Referenz

`tests/reference/piezo_bar_1d.py` ist eine separate, dichte 1D-FE-Rechnung
mit **knotenweise verschachtelten** Unbekannten `[u₀,φ₀,u₁,φ₁,…]`. Sie
assembliert ihre Zweiknoten-Elemente direkt aus der unabhängig angegebenen
Enthalpiedichte `C u'²/2 + e u'φ' - κ φ'²/2` und löst mit NumPy. Sie importiert
keine Produktionsmodule für Geometrie, Material, Assembly, Randbedingungen
oder Solver. Der Vergleich verwendet vier Längenschichten mit zwei
Materialien `(C;e;κ)=(12;1,5;2,5)` unten und `(18;2;3)` oben, `F=0,6`
und `V=0,8` bei Einheitsgeometrie. Für die 3D-Tet4-Lösung werden nodale
Felder auf allen Ebenen, obere Elektrodenladung und gesamte Enthalpie mit
dieser zweiten FE-Implementierung verglichen. Der eigenständige 1D-Code
prüft Diskretisierung und Assembly unabhängig vom 3D-Rechenpfad; er ist
kein Vergleich mit einem etablierten Fremdsolver oder Messdaten.

## Vierstufige Konvergenz

`tests/convergence/test_piezo_d33_convergence.py` verwendet auf dem
Einheitswürfel die gekoppelten quadratischen Felder
`u_z=0,2 z(1-z)` und `φ=0,3 z(1-z)` mit `C=10`, `e=2`, `κ=3`.
Aus den starken Gleichungen folgen unabhängig `b_z=5,2 N/m³` und
`ρ=1 C/m³`. Randwerte unten/oben sind null; Seitenflächen sind natürlich.
Die Netze haben `n=2,4,6,8` Abschnitte je Raumrichtung. Für den regulären
Aufbau sind die erwarteten Fehler analytisch

```text
L2(u,φ) = sqrt((0,2² + 0,3²)/30) · n⁻²
H1_positiv = sqrt((10·0,2² + 3·0,3²)/3) · n⁻¹
H_exakt = (10·0,2²/2 + 2·0,2·0,3 - 3·0,3²/2)/3 = 0,061666… J
```

`H1_positiv` nutzt die positive entkoppelte Gradientennorm; die
indefinite elektrische Enthalpie wäre als Fehlernorm ungeeignet. Die Tests
integrieren beide Fehler mit einer unabhängigen Gauss-Regel, fordern für
die letzten zwei Verfeinerungsschritte Ordnungen 1,8–2,2 beziehungsweise
0,85–1,15 und abnehmenden Enthalpiefehler. Außerdem werden feldweise
Residuen, Kraft- und Ladungsbilanz kontrolliert. Vier Netzstufen und
asymptotische Raten belegen Lösungsverifikation, nicht nur ein kleines
algebraisches Residuum.
