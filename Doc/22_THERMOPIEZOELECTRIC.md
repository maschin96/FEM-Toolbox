# Stationärer Ablauf T → (u, φ) (#26)

Status: Implementiert und durch unabhängige Ein-Tetraeder-Referenzen sowie
End-to-End-Modelltests verifiziert; keine experimentelle Validierung.

`solve_thermopiezoelectric` löst die stationäre Wärmegleichung und überträgt
die absoluten Knotentemperaturen (K) auf demselben Tet4-Netz in den
monolithischen piezoelektrischen Solver. Die thermische Eigenverzerrung je
Zelle ist `epsilon_th = alpha * (mean(T_nodes) - T_ref)` in
Engineering-Voigt-Reihenfolge. `alpha` und `T_ref` können regionsweise
zugeordnet werden. Mechanische Arbeit speist die Wärmegleichung nicht zurück;
transiente Terme sind nicht enthalten.

`electric_thermal_strain` steuert den elektrischen thermischen Term ausdrücklich.
Bei `false` bleibt er aus; bei `true` wird `-e epsilon_th` im elektrischen
Verschiebungsfeld berücksichtigt. Der Referenzsolver verwendet weiterhin die
Stress-Charge-Konvention aus [Doc/02](02_MATHEMATICAL_MODEL.md).

```python
from fem_toolbox.solvers.thermopiezoelectric import solve_thermopiezoelectric

result = solve_thermopiezoelectric(
    mesh,
    thermal_material,
    piezoelectric_material,
    expansion,
    thermal_conditions,
    piezoelectric_conditions,
    electric_thermal_strain=True,
)
print(result.heat.power_balance)
print(result.piezoelectric.force_balance, result.piezoelectric.charge_balance)
```

Das Ergebnis enthält Temperatur, Wärmefluss und Leistungsbilanz unter `heat`;
Verschiebung, Potential, Spannung, elektrisches Feld, E-/D-Felder, Energie,
Kraft- und Ladungsbilanz unter `piezoelectric`. Beide Teilergebnisse behalten
ihre getrennten algebraischen Residuen und Solverprovenienzen. Die
deklarative `Model`-API gibt dieses Paar unter `field_result` aus und hasht
beide Eingaben sowie alle drei Materialzuordnungen in der Laufprovenienz.

## Unabhängiger Grenzfall

Ein Einheits-Tetraeder wird mechanisch vollständig geklemmt. Mit
`C_xx=10 Pa`, `e_xx=2 C/m²`, `kappa_xx=4 F/m`, `alpha_x=0.01/K`,
`T_ref=300 K` und `T=310 K` folgt `epsilon_th,x=0.1`. Die x=0-Knoten
liegen auf `phi=0`; der x=1-Knoten ist elektrisch frei. Bei aktivem
elektrisch-thermischem Term erzwingt Ladungsfreiheit
`D_x=e_xx(0-0.1)+kappa_xx E_x=0`, also `E_x=0.05 V/m` und
`phi(x=1)=-0.05 V`. Der Test prüft diese Werte mit absolut `1e-12` in den
jeweiligen SI-Einheiten. Damit ist die T-Übergabe im elektrischen Feld
sichtbar; bei deaktiviertem Term ist das Potential null.

Ein zweiter Test fixiert alle mechanischen und elektrischen DOFs und prüft die
lineare Überlagerung thermischer und elektrischer Anregung gegen direkt aus
dem Stoffgesetz berechnete `sigma_x=1 Pa` und `D_x=-4.2 C/m²`.
Die Einzelfeldgrenzen werden gegen die bestehenden Wärme- und
Piezo-Solver geprüft. Leistung, Kraft, Moment und Ladung behalten getrennte
Bilanzen. Diese Tests sind Code- und Lösungsverifikation eines linearen
Grenzfalls. Der ergänzende 3D-Vier-Netze-Nachweis mit Materialgrenze und
leitenden Elektroden steht in [Doc/25](25_THERMOPIEZO_CONVERGENCE.md).
