# Stationäre sequenzielle Thermomechanik (#16)

`solve_thermomechanics` löst zuerst die stationäre Wärmeleitung und übergibt die
Knotentemperaturen auf demselben Tet4-Netz an `solve_elasticity`. Die Kopplung ist
einseitig. Wärmequellen aus mechanischer Arbeit und transiente Effekte gehören
nicht zu diesem Modell.

`ThermalExpansion(alpha, reference_temperature)` erwartet den globalen
Engineering-Voigt-Vektor `[xx, yy, zz, yz, xz, xy]` in `1/K` und eine absolute
Referenztemperatur in K. `ThermalExpansion.from_principal` rotiert drei
Hauptachsenkoeffizienten mit einer orthonormalen Matrix, deren Spalten die
Materialachsen im globalen System darstellen. `expansion` kann ein Wert für alle Zellen oder eine Liste
in Tet4-Eingabereihenfolge sein. Damit sind regionsweise Koeffizienten und
Referenztemperaturen möglich. `temperature` in `solve_elasticity` ist ein
Knotenfeld mit `mesh.n_points` absoluten Temperaturen in K.

Pro Zelle gilt `epsilon_th = alpha * (mean(T_nodes) - T_ref)`. Da `B` und `C`
auf Tet4 konstant sind, ist der thermische Lastvektor
`V * B.T @ C @ epsilon_th` für ein linear interpoliertes Temperaturfeld exakt.
Ausgegebenes `strain` ist die Gesamtdehnung; `stress` und `element_energy`
verwenden die mechanische Dehnung `strain - epsilon_th`. Reaktionen und
Lastvektor enthalten den äquivalenten thermischen Lastterm. Die Provenienz
enthält Temperatur, Expansionsparameter und den tatsächlich integrierten
Eigenverzerrungsvektor je Zelle.

```python
from fem_toolbox.materials import ThermalExpansion
from fem_toolbox.solvers.elasticity import solve_elasticity

result = solve_elasticity(
    mesh,
    elastic_material,
    mechanical_conditions,
    expansion=ThermalExpansion([1e-5, 1e-5, 1e-5, 0, 0, 0], 293.15),
    temperature=temperature_K,
)
```

Für den vollständigen Ablauf akzeptiert `solve_thermomechanics` zusätzlich
`thermal_material` und `thermal_conditions` und liefert getrennte `heat`- und
`elasticity`-Ergebnisse. Beide Randbedingungssätze verwenden die jeweiligen
feldspezifischen DOF-Indizes aus `thermal_dofs(mesh)` und
`mechanical_dofs(mesh)`.

## Verifikation zu #17

Der gekoppelte Manufactured-Solution-Test verwendet den Einheitswürfel mit
`T=300+x(1-x)` K und `u=(0.01x²,0,0)` m. Aus `k=2 W/(m K)` folgt unabhängig
`Q=4 W/m³`. Für `E=100 Pa`, `nu=0.25`, `alpha_x=0.002/K` und `T_ref=300 K`
ist `C11=120 Pa`; damit folgt aus `-div sigma=f` der Volumenlastterm
`f_x=-120[0.02-0.002(1-2x)] N/m³`. Die exakte Spannung ist
`sigma=(120,40,40,0,0,0)[0.02x-0.002x(1-x)] Pa`. Die Randwerte stammen
direkt aus diesen Funktionen. `tests/convergence/test_thermomechanics_convergence.py`
integriert die Fehlerquadrate unabhängig mit Gauss-Quadratur auf jedem Tet4.

| Zellen je Achse | `L2(T)` K | `L2(u)` m | `L2(sigma)` Pa |
|---:|---:|---:|---:|
| 2 | 0.0456435 | 0.000456435 | 0.383401 |
| 4 | 0.0114109 | 0.000114109 | 0.191659 |
| 6 | 0.00507151 | 0.0000507151 | 0.127767 |
| 8 | 0.00285272 | 0.0000285272 | 0.0958242 |

Die letzten beiden Verfeinerungsschritte ergeben jeweils Ordnung 2.00 für
Temperatur und Verschiebung sowie 1.00 für Spannung. Die Fenster `[1.8,2.2]`
und `[0.85,1.15]` trennen die erwarteten P1-Ordnungen von falschen Ordnungen,
ohne eine bitgenaue Solverausgabe zu verlangen. Die Bilanzgrenzen `1e-9 W`,
`1e-9 N` und `1e-9 N m` liegen weit unter den Lastskalen von 4 W und etwa
2 N. Feldweise relative algebraische Residuen bleiben unter `1e-12`.

`tests/patch/test_thermomechanics_patch.py` prüft zusätzlich eine Drehung der
Materialhauptachsen um 90°: `20 µstrain/K` in der ersten Materialachse,
`ΔT=25 K`, `E=100 GPa`, `nu=0.25`. Frei entsteht `u_y=500 µstrain*y` ohne
Spannung. Vollständig behindert sind die Spannungen
`(-20,-60,-20,0,0,0) MPa`, die gespeicherte Energie `15 kJ` und die
Reaktionskraft auf der Fläche `y=1` gleich `-60 MN`. Hier sind `1e-5 Pa`
für Spannungen, `1e-5 N` für Kräfte sowie `1e-8 J` für Energie Rundungsgrenzen bei
Skalen von `10^7` bis `10^8` in SI-Einheiten; sie entsprechen relativen
Abweichungen weit unter `1e-12`. Die Materialeingaben werden ausdrücklich
von GPa und Mikro-Dehnung in Pa und `1/K` umgerechnet. Diese Nachweise sind
Code- und Lösungsverifikation, keine experimentelle Validierung.
