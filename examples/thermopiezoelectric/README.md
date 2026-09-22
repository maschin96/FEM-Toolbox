# Deklarativer T → (u, φ)-Referenzfall (#27)

Dieses synthetische Ein-Tet4-Beispiel läuft von der
[Gmsh-Netzdatei](mesh.msh) über das [TOML-Modell](case.toml) bis zu
`fields.vtu` und `report.json`. Es demonstriert die öffentliche
`Model.from_file(...).solve()`-API. Die Zahlen sind bewusst klein und in SI;
sie beschreiben keinen konkreten Werkstoff und keine experimentelle Validierung.

## Ausführung

Nach Installation des Pakets aus dem Repository-Root:

```sh
python examples/thermopiezoelectric/run.py --output /tmp/fem-thermopiezo-example
```

`fields.vtu` enthält rohe, ungemittelte Knotenfelder `temperature_K`,
`displacement_m`, `potential_V` und zellweise `heat_flux_W_m2`, `stress_Pa`
in Engineering-Voigt-Reihenfolge sowie `electric_displacement_C_m2`. Die
Datei lässt sich in ParaView öffnen. `report.json` enthält Modell-/Toolbox-
Version, SHA-256 der effektiven Modelleingabe, beide Solverprovenienzen,
feldweise Residuen und Wärme-, Kraft-, Moment- und Ladungsbilanz. Der Bericht
ist ein **beispielspezifischer** Vertrag; der allgemeine Export wird in #140
festgelegt. Die Ausgabedateien werden bei erneutem Lauf im gewählten Ordner
ersetzt. `run.py` benötigt keine Änderungen an Paket- oder Netzdateien.

## Unabhängige Referenz

Das Netz hat die Knoten `(0,0,0)`, `(1,0,0)`, `(0,1,0)` und `(0,0,1)` m.
`T=300+10x` K ist vollständig vorgeschrieben. Mit `k=2 W/(m K)` folgt
`q=-k grad T=(-20,0,0) W/m²`. Die lineare Tet4-Interpolation ist hierfür
exakt; die mittlere Zelltemperatur ist 302,5 K. `alpha_x=0,01/K` und
`T_ref=300 K` ergeben `epsilon_th,x=0,025`.

Nur `u_x` und `phi` am Knoten `(1,0,0)` sind frei. Mit `C_xx=10 Pa`,
`e_xx=2 C/m²`, `kappa_xx=4 F/m` und aktivem elektrisch-thermischem Term
lauten die konstanten konstitutiven Größen

```text
sigma_x = 10*(u_x - 0,025) + 2*phi
D_x     = 2*(u_x - 0,025) - 4*phi
```

Die freie mechanische Gleichung ergibt `sigma_x=0`. Gemäß der negativen
elektrischen Zeile der Toolbox ergibt die positive Volumenladung
`rho=4 C/m³` am freien Knoten `D_x=-rho/4=-1 C/m²`. Das lineare
2×2-System liefert unabhängig vom FEM-Solver:

| Größe | Referenz |
|---|---:|
| `u_x(1,0,0)` | `0,025 - 1/22 = -0,0204545454545 m` |
| `phi(1,0,0)` | `1/4,4 = 0,227272727273 V` |
| `q_x` | `-20 W/m²` |
| `sigma_x` | `0 Pa` |
| `D_x` | `-1 C/m²` |

Der Regressionstest verwendet für Felder absolute Grenzen von `1e-12` in
jeweiligen SI-Einheiten und `rtol=0`; für globale Bilanzen `1e-10 W`,
`1e-10 N` beziehungsweise `1e-10 C`. Die Grenzen liegen weit unter den
Anregungsskalen 20 W/m², 4 C/m³ und 0,2 V und decken Rundungsfehler ab,
ohne Vorzeichen- oder Kopplungsfehler zu verschleiern.

Die Geometrie ist nur ein Element groß und alle Temperaturen sind
vorgeschrieben. Daraus folgt kein Nachweis für Netzkonvergenz,
Mehrmaterial-Grenzflächen oder reale Bauteilgrößen; diese Prüfungen werden
in #141 behandelt.
