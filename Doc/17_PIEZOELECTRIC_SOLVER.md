# Monolithischer piezoelektrischer Referenzsolver (#21)

Status: Implementiert und mit algebraischen sowie affinen Patchtests verifiziert;
keine experimentelle Validierung.

## Gleichungen und Lösungsweg

`fem_toolbox.solvers.piezoelectric` assembliert für Tet4-Zellen die
symmetrisch-indefinite Form aus [Doc/16](16_PIEZOELECTRIC_ELEMENT.md):

```text
[ K_uu       K_uφ ] [ u ] = [ f_mech ]
[ K_uφ^T  -K_φφ  ] [ φ ]   [ q_signed ]
```

Der elektrische Lastvektor ist `q_signed = ∫ N ω dA - ∫ N ρ dV`.
Die DOFs liegen feldweise in der Reihenfolge `[u_x,u_y,u_z je aktivem
Knoten, φ je aktivem Knoten]`; ungenutzte Netzknoten erhalten keine DOFs.
`piezoelectric_dofs(mesh)` liefert die Indizes für `Dirichlet`.

Dirichlet-Werte werden symmetrieerhaltend eliminiert. Für die freien DOFs
skaliert eine kongruente Diagonalmatrix die beiden Blöcke mit dem Kehrwert
der Quadratwurzel ihres jeweiligen medianen positiven Diagonaleintrags.
Die Skalen haben vor dieser Transformation N/m beziehungsweise C/V und
stehen im Ergebnis. Danach löst SciPy Sparse LU das Sattelpunktsystem.
`block_scaling=False` deaktiviert nur diese Transformation; `method="direct"`
ist die einzige erlaubte Solverwahl. CG passt nicht zum indefiniten System;
MINRES benötigt zuerst einen nachgewiesen geeigneten Präconditioner.

Vor der Faktorisierung werden Starrkörpermoden und der Potentialbezug je
zusammenhängender Netzkomponente geprüft. Eine fehlende elektrische Gauge
wirft `ElectricGaugeError`; fehlende mechanische Lagerung, singuläre oder
schlecht konditionierte Matrizen werfen `ValidationError`. Die LU-Pivots
werden relativ mit `1e-12` kontrolliert.

## API-Beispiel

```python
import numpy as np

from fem_toolbox.boundary import Dirichlet
from fem_toolbox.materials import PiezoelectricMaterial
from fem_toolbox.mesh import Mesh
from fem_toolbox.solvers.piezoelectric import piezoelectric_dofs, solve_piezoelectric

mesh = Mesh(
    points=[[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetra4=[[0, 1, 2, 3]],
)
c = np.diag([10.0, 20.0, 30.0, 4.0, 5.0, 6.0])
e = np.zeros((3, 6))
e[0, 0] = 2.0
kappa = np.diag([3.0, 4.0, 5.0])
material = PiezoelectricMaterial(c, e, kappa)
dofs = piezoelectric_dofs(mesh)
free = np.r_[dofs.indices("u", [1], ["x"]), dofs.indices("phi", [1])]
fixed = np.setdiff1d(np.arange(dofs.size), free)
result = solve_piezoelectric(
    mesh,
    material,
    [Dirichlet(fixed, np.zeros(len(fixed)))],
    point_loads={1: [1.0, 0.0, 0.0]},
    volume_charge=-2.4,
)
print(result.displacement[1], result.potential[1])
print(result.relative_residual_u, result.relative_residual_phi)
```

Materialdaten und Lasten sind in SI: `C^E` Pa, `e` C/m²,
`κ^S` F/m, Verschiebung m, Potential V, Kraft N, Ladung C.
`body_force` ist N/m³; `tractions` ordnet äußeren Triangle3-Zeilen
Traktionsvektoren in N/m² zu. `volume_charge` ist ρ in C/m³;
`surface_charges` ordnet denselben Flächenzeilen den äußeren Fluss
`ω = D·n` in C/m² zu. Materialien und Expansion dürfen einzeln oder
zellweise angegeben werden. Nodaltemperaturen sind absolute Kelvinwerte;
`electric_thermal_strain=True` aktiviert den in [Doc/16](16_PIEZOELECTRIC_ELEMENT.md)
definierten elektrischen thermischen Term. Der Standardschalter lässt ihn aus.

## Ergebnis und Nachweise

`PiezoelectricResult` enthält rohe nodale Felder, Zell-Dehnungen,
Spannungen, E-/D-Felder und die drei Enthalpieanteile. Arrays sind
unveränderlich. `residual_u` (N) und `residual_phi` (C) stammen vom
ursprünglichen vorzeichenbehafteten System; Reaktionen sind nur auf
Dirichlet-DOFs belegt. `electrode_charges = -reactions_phi` beschreibt
angrenzende ideale Elektroden gemäß [Doc/12](12_ELECTROSTATICS.md).
`force_balance`, `moment_balance` und `charge_balance` prüfen globale
Bilanzen. Beide freien Feldresiduen werden separat durch
`|| |A| |x| + |b| ||` ihrer Zeilen normiert und gegen `1e-9` geprüft.
Die Provenienz nennt die diskreten Eingaben, Solverparameter und Skalen.

Der 2×2-Referenzfall in `tests/unit/test_piezo_solver.py` prüft Lösung,
Vorzeichen, Blockskalierung und getrennte Residuen gegen eine unabhängig
aufgestellte dichte Matrix. Derselbe Testbestand prüft fehlende Lagerung,
Gauge, ungültige Solverwahl und den expliziten thermischen Modellschalter.
`tests/patch/test_piezo_solver_patch.py` rekonstruiert affine gekoppelte
Felder auf einem verzerrten Mehrzellen-Netz, einmal mit vorgegebenen
Randwerten und einmal mit unabhängig berechneten Traktionen und D-Flüssen.
Diese Nachweise decken den linearen stationären Tet4-Referenzumfang ab;
Floating Electrodes, transiente Effekte und große parallele Modelle folgen
in separaten Arbeitspaketen.

Der anschließende [d₃₃-Benchmark](18_PIEZO_D33_VERIFICATION.md) prüft
direkten/inversen Effekt, Kurzschluss, einfachen Leerlauf, Reziprozität und
vierstufige Konvergenz. Der Leerlauf-Benchmark präzisiert die
Rundungstoleranz der globalen Ladungsbilanz bei physikalisch verschwindender
Nettoladung.

## M7: zweite Elementfamilie

[Stationäre Tet10-Piezoelektrizität](49_TET10_PIEZOELECTRIC.md) verwendet
denselben Solver-/Elektrodenvertrag mit allen zehn Zell- und sechs Flächenknoten.
Zellfelder sind Schwerpunktwerte, Energiebeiträge sind Volumenintegrale.
Thermische und dynamische Tet10-Piezokopplung bleiben Folgearbeit.
