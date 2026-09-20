# Tet4-Formulierung

Status: Entwurf

## 1. Geometrie und Formfunktionen

Für ein Tetraeder mit Knoten `a = 1..4` und Koordinaten `x_a` gilt

```text
N_a(x,y,z) = a_a + b_a x + c_a y + d_a z
```

Mit

```text
X = [[1, x_1, y_1, z_1],
     [1, x_2, y_2, z_2],
     [1, x_3, y_3, z_3],
     [1, x_4, y_4, z_4]]
```

stehen die Koeffizienten in den Spalten von `X^{-1}`. Die Gradienten der
Formfunktionen sind im Element konstant. Das Volumen ist

```text
V = abs(det([x_2-x_1, x_3-x_1, x_4-x_1])) / 6.
```

Die Implementierung soll keine explizite Matrixinverse bilden, sondern die
Koeffizienten mit einem linearen Löser bestimmen. Zusätzlich werden
Determinante, Kantenmaß und Konditionsindikator geprüft.

## 2. Operatoren

Für skalare Felder ist

```text
G = [grad N_1, grad N_2, grad N_3, grad N_4] in R^(3x4).
```

Die mechanische `B`-Matrix in der festgelegten Voigt-Konvention hat Größe
`6 x 12`. Für jeden Knoten mit `grad N_a = [N_x,N_y,N_z]^T` ist der Block

```text
B_a = [[N_x, 0,   0  ],
       [0,   N_y, 0  ],
       [0,   0,   N_z],
       [0,   N_z, N_y],
       [N_z, 0,   N_x],
       [N_y, N_x, 0  ]].
```

Alle Volumenintegranden aus konstanten Materialmatrizen, `B` und `G` sind
konstant und werden mit einem Integrationspunkt exakt integriert.

## 3. Elementmatrizen

Mit nodalen Vektoren `u_e`, `T_e`, `phi_e`:

```text
K_uu       = V B^T C^E B
K_uphi     = V B^T e^T G
K_phiu     = K_uphi^T
K_phiphi   = V G^T kappa^S G
K_TT       = V G^T k G
```

Das gekoppelte piezoelektrische Element wird als

```text
[ K_uu    K_uphi  ]
[ K_phiu -K_phiphi]
```

assembliert. Die thermische Expansionskopplung lautet

```text
K_uT = - integral_Omega_e B^T C^E alpha N dOmega.
```

Weil `B`, `C^E` und `alpha` konstant sind und
`integral N_a dOmega = V/4`, ist dieser Term analytisch exakt. Falls die
mechanische Dehnung in der elektrischen Verschiebung um die thermische Dehnung
reduziert wird, folgt zusätzlich konsistent

```text
K_phiT = - integral_Omega_e G^T e alpha N dOmega.
```

Dieser Block darf nicht ohne den zugehörigen Materialmodellschalter erscheinen.

## 4. Lastvektoren und Oberflächen

Für konstante Volumenlasten gilt

```text
f_b,a = (V/4) b.
```

Eine Tet4-Fläche ist ein linearer Dreiknoten-Dreieck. Bei konstanter
Oberflächenlast `t_bar`, Wärmestrom `q_bar` oder Oberflächenladung `omega_bar`
erhält jeder Flächenknoten ein Drittel des Gesamtbeitrags. Variable Lasten
werden mit einer Dreiecksquadratur passender Ordnung integriert.

Die Flächennormale wird aus der Elementtopologie bestimmt und gegen den
gegenüberliegenden Tetraederknoten nach außen orientiert. Dadurch ist die
Vorzeichenkonvention unabhängig von der Reihenfolge importierter Dreiecke.

## 5. Massen- und Kapazitätsmatrizen für spätere Transienten

Die konsistente skalare Tet4-Massenmatrix ist

```text
M_s = (V/20) [[2,1,1,1],
              [1,2,1,1],
              [1,1,2,1],
              [1,1,1,2]].
```

Daraus folgen die mechanische Massenmatrix `rho (M_s kron I_3)` und die
thermische Kapazitätsmatrix `rho c_p M_s`. Mass Lumping muss als eigene,
getestete Option implementiert werden.

## 6. Grenzen des Tet4

Der lineare Tetraeder besitzt konstante Verzerrung, Spannung, Feldstärke und
Temperaturgradient innerhalb eines Elements. Daraus folgen:

- niedrige Genauigkeit in Biegung und bei starken Gradienten,
- volumetrisches Locking nahe Inkompressibilität,
- diskontinuierliche rohe Fluss- und Spannungswerte über Elementgrenzen,
- hoher Netzbedarf für lokale Spannungsauswertung.

Nodale Mittelung verbessert nur die Darstellung und ist kein Genauigkeitsbeleg.
Ergebnisse werden primär elementweise und mit Netzkonvergenz berichtet.

## 7. Element-Invarianten

Jede Elementroutine muss automatisiert prüfen:

- Partition der Eins: `sum N_a = 1`
- `sum grad N_a = 0`
- Symmetrie von `K_uu`, `K_TT`, `K_phiphi`
- sechs Starrkörpermoden für ein unbelagertes mechanisches Einzelelement
- konstantes Potential erzeugt kein elektrisches Feld
- konstante Temperatur erzeugt keinen Wärmestrom
- `K_phiu = K_uphi^T` innerhalb numerischer Toleranz
- Invarianz unter starrer Translation und zulässiger Knotenumordnung
