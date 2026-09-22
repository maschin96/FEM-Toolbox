# Piezoelektrische Materialdaten und Achsrotation (#19)

Status: Implementiert und durch Materialtests codeverifiziert; noch kein
gekoppelter Piezo-Solver und keine experimentelle Validierung.

## Konvention und öffentliche API

`fem_toolbox.materials.PiezoelectricMaterial` speichert räumlich konstante
Materialdaten in globalen kartesischen Achsen. Die Stress-Charge-Form lautet

```text
sigma = C^E epsilon - e^T E
D     = e epsilon + kappa^S E
E     = -grad(phi)
```

Dabei sind `C^E` (6×6) in Pa, `e` (3×6) in C/m² und `kappa^S` (3×3)
in F/m. Spannungen und technische Verzerrungen verwenden die Reihenfolge
`[xx, yy, zz, yz, xz, xy]`; die Verzerrungseinträge für Schub sind
`gamma_ij = 2 epsilon_ij`. Elektrische Komponenten folgen `[x, y, z]`.
Die API nimmt absolute Permittivitäten, keine relativen Materialzahlen, an.

```python
from fem_toolbox.materials import PiezoelectricMaterial

material = PiezoelectricMaterial(C_E, e, kappa_S, axes=local_axes)
```

`axes` ist eine orthonormale Matrix mit Determinante +1, deren Spalten die
lokalen Materialachsen im globalen Koordinatensystem darstellen. Ohne `axes`
werden die Eingabetensoren bereits in globalen Achsen verstanden. Alle
Eingabeformen, endlichen Werte, Symmetrien und die positive Definitheit von
`C^E` und `kappa^S` werden geprüft. Sehr schlecht aufgelöste Eigenwerte mit
`lambda_min/lambda_max <= 1e-12` werden abgelehnt. Akzeptierte
Symmetrieabweichungen bis `1e-12` relativ zum größten Koeffizienten werden
gemittelt. Eingaben werden kopiert und unveränderlich gespeichert; fehlerhafte
Daten lösen `ValidationError` aus. Die Steifigkeitsprüfung ersetzt keine
Anwendungseignung, etwa bei nahezu inkompressiblen Tet4-Materialien.

## Konvertierung aus Strain-Charge-Daten

Die alternative Darstellung lautet

```text
epsilon = s^E sigma + d^T E
D       = d sigma + kappa^T E
```

`s^E` hat die Einheit 1/Pa, `d` die Einheit C/N und `kappa^T` F/m. Die
Konvertierung erfolgt mit linearen Gleichungslösungen:

```text
C^E     = (s^E)^-1
e       = d C^E
kappa^S = kappa^T - d C^E d^T
```

Die letzte Differenz wird ausdrücklich auf positive Definitheit geprüft.
Eine positive `kappa^T` allein ist nicht ausreichend. Der umgekehrte Weg
`as_strain_charge()` liefert globale `s^E`, `d` und `kappa^T` als Kopien.

```python
material = PiezoelectricMaterial.from_strain_charge(
    compliance=s_E,
    strain_coupling=d,
    permittivity_at_stress=kappa_T,
    axes=local_axes,
)
s_global, d_global, kappa_T_global = material.as_strain_charge()
```

## Achsrotation

Aus einer Rotation `R` entstehen die Voigt-Transformationen `T_epsilon(R)`
für technische Verzerrungen und `T_sigma(R)` für Spannungen. Sie sind wegen
der Schubfaktoren verschieden. Der Code bildet sie aus rotierten symmetrischen
3×3-Basistensoren und berechnet

```text
C_global     = T_sigma(R) C_local T_epsilon(R)^-1
e_global     = R e_local T_epsilon(R)^-1
kappa_global = R kappa_local R^T.
```

Damit bleiben mechanische Arbeit `sigma:epsilon` und die Kopplungsrelation
erhalten. Eine Drehung um 90° um `z` vertauscht beispielsweise lokale
`C_xxxx=10` und `C_yyyy=20` zu global `C_xxxx=20`, `C_yyyy=10`; ein
lokaler Eintrag `e_zxx=2` wird zu global `e_zyy=2`.

## Provenienz und Nachweis

`to_provenance_dict()` liefert ein getrenntes, JSON-kompatibles Abbild der
ursprünglichen lokalen Eingabetensoren, der Eingabeform und Achsen sowie der
transformierten globalen Stress-Charge-Tensoren und ihrer Einheiten. Die
Struktur kann in `RunProvenance.solver_parameters` übernommen werden, sobald
der gekoppelte Solver (#21) diese Materialdaten verwendet.

Die Tests in `tests/unit/test_piezo_material.py` prüfen die skalare
Konvertierung unabhängig mit `C=1/s`, `e=d/s` und
`kappa^S=kappa^T-d²/s`, Hin-/Rücktransformation, 90°-Rotation, eine
45°-Dehnung mit technischer Scherung, Eingabekopien und ungültige Daten.
Die numerischen Testwerte sind skalierte algebraische SI-Referenzen, keine
Messwerte eines realen Piezomaterials. Toleranzen von etwa `1e-14` relativ
und `1e-13` absolut berücksichtigen wenige Matrixmultiplikationen und
lineare Lösevorgänge bei Koeffizienten der Größenordnung 1–30; reale
Materialdatensätze benötigen passende SI-Skalierung und eigene Nachweise.
