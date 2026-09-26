# Residualindikator für stationäre Tet4-Wärme (#36)

Status: implementierte Forschungs-API im M7-Issue-Branch. Code-Verifikation,
keine experimentelle Validierung und kein garantierter Fehlerschätzer.

## Vertrag und Norm

`fem_toolbox.post.heat_indicator(mesh, temperature, material, ...)` erhält
Temperaturen in K für **alle Netzknoten** (`mesh.n_points`), einschließlich
gegebenenfalls inaktiver Knoten. Ein kompaktes `HeatResult.temperature` muss der
Aufrufer über `HeatResult.node_ids` zuordnen. Leitfähigkeit ist homogen oder
zellweise als `ThermalMaterial` vorgegeben. Es gelten die stationären Vorzeichen
aus [dem mathematischen Modell](02_MATHEMATICAL_MODEL.md).

Zielnorm ist die thermische Energienorm

```text
||e||_E² = integral_Omega grad(e)^T k grad(e) dV  [W K].
q_h = -k grad(T_h), lambda_K = kleinster Eigenwert von k_K.
eta_K² = h_K²/lambda_K * integral_K Q² dV
       + 1/2 sum_F innen h_F/lambda_F * integral_F [q_h.n]² dA
       + sum_F Neumann/Robin h_F/lambda_K * integral_F r_F² dA.
lambda_F = min(lambda_K, lambda_Nachbar).
r_F = q_h.n - q_pres - h*(T_h-T_inf).
```

`h_K` und `h_F` sind maximale Kantenlängen. Innerhalb eines affinen Tet4 mit
konstanter Leitfähigkeit verschwindet `div(q_h)`, nicht aber die Quelle oder
der innere Flusssprung. Der Sprung ist die Summe der beiden nach außen gerichteten
Normalflüsse; jede innere Fläche wird einmal berechnet und hälftig verteilt.
An Materialgrenzen wird der **Fluss**, nicht der Gradient verglichen.
Nichtzusammenhängende externe Element-IDs verändern die zellweise Reihenfolge nicht.

`HeatIndicator.volume_squared`, `.interior_squared` und `.boundary_squared`
enthalten unveränderliche Arrays in W K, `.cell` enthält `eta_K` in sqrt(W K),
`.total` deren euklidische Norm. Die Gewichtung ist ein nachvollziehbarer
Residualindikator; keine Zuverlässigkeitskonstante, Effizienzkonstante oder
Robustheit gegenüber beliebigem Materialkontrast/Anisotropie wird zugesagt.
Eine Datenoszillation oder Dirichlet-Interpolationskorrektur ist nicht enthalten.

## Randdaten und Grenzen

Die Triangle3-Zeilen müssen den vollständigen äußeren Rand abdecken. Der
Topologie-Preflight lehnt fehlende, innere, fremde und nichtmannigfaltige
Randflächen ab. `dirichlet_facets` bezeichnet explizit die Flächen ohne
Randresidualterm. Alle übrigen Flächen sind ohne weitere Angaben isoliert.
`outward_fluxes` und `convection` verwenden dieselben Triangle3-Indizes und
SI-Konventionen wie `solve_heat`; Fluss und Konvektion addieren sich. Überlappung
mit Dirichlet-Flächen ist unzulässig. Der Aufrufer verantwortet die Übereinstimmung
mit dem tatsächlich gelösten stationären Problem; aus Temperaturen allein lässt
sich die Randart nicht rekonstruieren. Vorgegebene Dirichlet-Daten müssen auf den
markierten Flächen erfüllt sein; der Indikator prüft dies nicht selbst.

Volumenquadratur ist bis Grad 5, Flächenquadratur bis Grad 4 exakt, jeweils
bezogen auf das **quadrierte Residuum**. Nichtpolynomiale Quellen und Randdaten
brauchen eine eigene Quadraturbewertung. Tet10 wird gezielt abgelehnt;
Transienten, gekoppelte Felder, automatisches Remeshing und M10-Adaptivität sind
nicht Bestandteil dieser API. Eine Modellvalidierung oder Konvergenzstudie wird
nicht durch einen einzelnen Indikatorwert ersetzt.

## Anwendung und deterministische Markierung

Nach einem stationären Lauf, mit den tatsächlichen Problemranddaten:

```python
import numpy as np
from fem_toolbox.post import heat_indicator, mark_bulk

temperature = np.zeros(mesh.n_points)
temperature[result.node_ids] = result.temperature
indicator = heat_indicator(
    mesh,
    temperature,
    material,
    source=source,
    dirichlet_facets=dirichlet_triangle_rows,
    outward_fluxes=outward_fluxes,
    convection=convection,
)
marked = mark_bulk(indicator.cell, fraction=0.5)
```

`mark_bulk` wählt die kleinste Anzahl von Zellen, deren quadrierte Indikatoren
mindestens den Anteil `fraction` an der Gesamtsumme tragen. Gleichstände werden
nach aufsteigender Zellzeile aufgelöst; Ausgabe ebenfalls aufsteigend. `fraction`
muss in `(0, 1]` liegen. Exakt null ergibt eine leere Auswahl; es gibt keine
implizite dimensionsbehaftete Rauschschwelle. Ein anwendungsspezifisches globales
Abbruchkriterium muss der Aufrufer vor der Markierung anwenden. Normierung vor
dem Quadrieren schützt die Auswahl vor gemeinsamem Skalierungsüberlauf.

## Verifikation und vorab festgelegte Kriterien

Unit-Tests prüfen analytisch integrierte Einzeltetraeder-Quellen, innere Sprünge,
Neumann-/Robin-Terme, deterministische minimale Markierung, invalides Input,
Tet10-Ablehnung, verzerrte affine Patches und Permutationsinvarianz. Für Größen
von Ordnung eins dienen relative Toleranzen `1e-14` beziehungsweise absolute
Patchgrenzen `1e-13` bis `1e-12` ausschließlich der Rundungsfehlerkontrolle.
Der gelöste Materialtest verwendet `k=1` beziehungsweise `k=8 W/(m K)` in zwei Schichten, konstante Flussdichte und einen achtfachen
Gradientensprung. Der Indikator bleibt innerhalb der Rundungstoleranz null.

Der glatte Manufactured-Solution-Test verwendet dieselbe unabhängige quadratische
Temperatur und anisotrope Leitfähigkeit wie die bestehende Wärme-Konvergenzprüfung.
Vor dem ersten Lauf festgelegt: Indikatorordnung 0,8–1,2 auf den letzten beiden
Intervallen, streng sinkender Indikator, Effektivität 0,1–50 und Verhältnis
maximaler/minimaler Effektivität unter 1,3. Die breite absolute Spanne ist kein
Genauigkeitsversprechen; entscheidend sind stabiler Quotient und richtige Ordnung.

| Würfelunterteilungen n | Energienormfehler | Indikator eta | Effektivität eta/Fehler |
|---|---:|---:|---:|
| 2 | 0,500000 | 5,882661 | 11,765323 |
| 4 | 0,250000 | 2,998008 | 11,992033 |
| 6 | 0,166667 | 2,011110 | 12,066657 |
| 8 | 0,125000 | 1,512975 | 12,103796 |

Beobachtete Indikatorordnungen: 0,9725 / 0,9847 / 0,9893. Die Energienormfehler
werden unabhängig mit symmetrischer Grad-2-Quadratur exakt integriert.

Der lokalisierte Test setzt `T=exp(-80*(x-0.5)²)` und die analytische Quelle
`Q=(160-25600*(x-0.5)²)*T` auf einem n=8-Netz an. Vorab gilt: im mittleren
Halbvolumen liegen jeweils über 65 % der quadrierten Indikatoren und Fehler;
über 80 % der für 50 % Bulk-Anteil markierten Zellen liegen dort. Gemessen:
99,39 % Indikatoranteil, über 99 % Fehleranteil, 100 % der markierten Zellen.
Die nichtpolynomiale Fehlerquadratur wird unabhängig mit Gauss-Ordnungen 7 und 8
verglichen (relative Abweichung höchstens `1e-5`). Dies ist ein lokalisierter
Referenzfall, kein allgemeiner Beweis einer optimalen adaptiven Strategie.

Reproduktion im eingerichteten Repository:

```sh
PYTHONPATH=src .venv/bin/python -m pytest tests/unit/test_heat_indicator.py tests/convergence/test_heat_indicator_convergence.py -q -s
```

Der TDD-Red-Lauf scheiterte zunächst am fehlenden Modul `fem_toolbox.post`;
nach Implementierung bestehen die analytischen und numerischen Kriterien.
Issue-Abschluss und Integration setzen zusätzlich die geltenden CI-/Reviewregeln
voraus. PETSc ist für diesen reinen Postprocessing-Nachweis keine Voraussetzung.
