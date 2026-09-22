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
