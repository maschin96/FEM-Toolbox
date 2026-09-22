# Gekoppeltes piezoelektrisches Tet4-Element (#20)

Status: Lokale Operatoren und Patch-Verifikation implementiert. Globale
Assembly, Randbedingungen und Lösung des Sattelpunktproblems folgen in #21.

## Diskrete Form und Vorzeichen

`fem_toolbox.elements.piezoelectric` verwendet die Materialkonvention aus
[Doc/15](15_PIEZOELECTRIC_MATERIAL.md) und die ursprüngliche lokale
Knotenreihenfolge von `Tet4Geometry`. Der lokale Vektor hat die Reihenfolge
`[u_0x,u_0y,u_0z, ..., u_3z, phi_0, ..., phi_3]`.

Mit mechanischem Operator `B` (6×12), Potentialgradient `G` (3×4),
positivem Elementvolumen `V`, Steifigkeit `C^E`, Kopplung `e` und
Permittivität `kappa^S` ist

```text
K_e = V [[B^T C^E B,  B^T e^T G],
         [G^T e B,   -G^T kappa^S G]].
```

Der negative elektrische Block folgt der mit `-1` multiplizierten
Gauß-Zeile des [mathematischen Modells](02_MATHEMATICAL_MODEL.md). Die
Kopplungsblöcke sind Transponierte. Sie dürfen nicht getrennt mit
unterschiedlichen Vorzeichen assembliert werden. Ein konstantes Potential
hat `G phi = 0` und erzeugt deshalb keine Kopplungskraft.

```python
from fem_toolbox.elements.piezoelectric import piezoelectric_matrix

local_matrix = piezoelectric_matrix(geometry, material)  # (16, 16)
```

Die Einträge haben je nach Block die Einheiten N/m, N/V, C/m und C/V.
`Tet4Geometry` lehnt entartete Elemente bereits vor dem Aufruf ab.

## Mechanische Eigendehnung und elektrischer Modellschalter

Eine konstante thermische Eigendehnung `epsilon_th` (6,) erzeugt immer die
mechanische rechte Seite `V B^T C^E epsilon_th`. Standardmäßig verwendet
`D` die **gesamte** Dehnung `epsilon`, sodass die elektrische thermische
rechte Seite null ist. Der explizite Schalter `electric_thermal_strain=True`
verwendet stattdessen `epsilon-epsilon_th` in `D`; dann lautet der zusätzliche
elektrische RHS-Beitrag `V G^T e epsilon_th`. Der Tangentenblock bleibt in
beiden Fällen gleich. Diese Entscheidung folgt der projektweiten Regel,
den Term `-e alpha theta` niemals implizit zu aktivieren.

```python
from fem_toolbox.elements.piezoelectric import piezoelectric_thermal_load

local_rhs = piezoelectric_thermal_load(
    geometry,
    material,
    thermal_strain,
    electric_thermal_strain=True,
)  # (16,)
```

## Felder und elektrische Enthalpie

`piezoelectric_response` nimmt nodale Verschiebungen `(4,3)` in m und
Potentiale `(4,)` in V. Es gibt konstante Dehnung, elektrisches Feld,
Spannung und elektrische Verschiebung sowie integrierte Energieterme zurück.
Für `epsilon_m=epsilon-epsilon_th` und
`epsilon_c=epsilon_m` bei aktivem Schalter beziehungsweise `epsilon_c=epsilon`
sonst gilt

```text
H_mech =  V/2 epsilon_m^T C^E epsilon_m
H_cpl  = -V epsilon_c^T e^T E
H_diel = -V/2 E^T kappa^S E
H      = H_mech + H_cpl + H_diel.
```

`H_diel` ist der **negative Beitrag zur elektrischen Enthalpie**, keine
negative gespeicherte Feldenergie. Die erste Ableitung von `H` nach den
nodalen Variablen ist `K_e [u,phi] - f_th` (ohne äußere Lasten). Ohne
vorgegebenen Potentialbezug enthält der elektrische Block eine konstante
Nullmode; mechanische Starrkörpermoden bleiben ebenfalls bestehen.

## Unabhängige Nachweise und Grenzen

`tests/unit/test_piezo_element.py` prüft die Blockvorzeichen, Reziprozität,
den Grenzfall eines konstanten Potentials, analytische affine Felder und
Energieableitungen durch zentrale Differenzen. Die Handrechnung für den
Einheits-Tet4 mit `C_xxxx=10 Pa`, `e_xxx=2 C/m²`, `kappa_xx=3 F/m`,
`epsilon_xx=0.2` und `grad(phi)_x=0.3 V/m` ergibt
`sigma_xx=2.6 Pa` und `D_x=-0.5 C/m²`. In
`tests/patch/test_piezo_element_patch.py` werden die Elementresiduen eines
verzerrten Mehrzellen-Netzes mit unabhängig über Randdreiecke integrierten
konstanten Traktionen und D-Flüssen verglichen, auch nach Knotenumordnung.
Das sind Code- und Patch-Verifikation, keine experimentelle Validierung.
