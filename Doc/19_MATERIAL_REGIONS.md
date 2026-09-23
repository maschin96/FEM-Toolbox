# Materialregionen und perfekte Grenzflächen (#137)

Status: Regionenzuordnung ist in die deklarative `Model`-API integriert und
an den stationären Tet4-Solvern sowie am vollständigen thermopiezoelektrischen
Modell verifiziert. `Model` prüft die Belegung vor der Assembly.

## Vertrag

`assign_region_materials(mesh, field, regions, materials)` ordnet einen
benannten Materialkatalog den `volume_sets` des Netzes zu. Unterstützte Felder
sind `T` (`ThermalMaterial`), `u` (`ElasticMaterial`), `phi`
(`DielectricMaterial`), `piezo` (`PiezoelectricMaterial`) und `expansion`
(`ThermalExpansion`). Die Materialobjekte enthalten validierte SI-Tensoren.
Jede Tet4-Zeile muss für das gewählte Feld genau einen Materialnamen erhalten.
Überlappende Sätze mit **demselben Namen** sind erlaubt und werden beide als
Herkunft gespeichert. Verschiedene Namen auf derselben Zelle gelten auch dann
als Konflikt, wenn ihre Zahlenwerte zufällig übereinstimmen. Unbekannte
Regionen oder Materialien, Lücken, ungenutzte Definitionen und leere
Volumennetze stoppen vor dem Solver. Diagnosen nennen die externe Tet4-ID,
sofern vorhanden.

`MaterialAssignment.cell_materials` folgt der Tet4-Zeilenreihenfolge und kann
direkt an den passenden Solver übergeben werden. `to_provenance_dict()` enthält
Zell-IDs, Materialnamen, alle zugehörigen Gruppen, SI-Tensoren und bei
Piezo-Materialien Eingabeform, lokale Eingaben, globale Tensoren und Achsen.
Bei `ThermalExpansion.from_principal` bleiben die verwendeten lokalen Achsen
zusätzlich erhalten; direkt eingegebene globale Voigt-Koeffizienten erhalten
keine behauptete lokale Orientierung. Die Option
`electric_thermal_strain` ist eine **globale Solveroption** und wird in den
Solverparametern des Piezo-Solvers im Ergebnisbericht festgehalten. Der
Bericht enthält auch die wirksamen Materialkataloge und Achsen pro Region.

## Konforme Grenzfläche

Die Solver verwenden dieselben globalen Knoten und Freiheitsgrade auf beiden
Seiten einer gemeinsamen Tet4-Facette. Damit sind `T`, `u` und `phi` dort
stetig. Die schwachen Formen assemblieren beidseitige Elementbeiträge auf
denselben Freiheitsgraden; ohne explizite Grenzflächenlast folgen daraus
entgegengesetzte Wärme- und elektrische Flüsse sowie Traktionen. Es gibt
keinen Kontaktwiderstand und keine Sprungfreiheitsgrade.

## Unabhängige Zweischicht-Referenz

Die Tests verwenden einen Quader `1 m × 1 m × 1 m` mit Grenzfläche bei
`x=0,5 m`, auf jeder Seite isotrope Koeffizienten `a₁` und `a₂=2a₁`. Für
einen Einheitsunterschied der Primärvariable ergibt die 1D-Reihenschaltung
den konstanten Fluss

`J = Δv / (0,5/a₁ + 0,5/a₂)`.

Links ist `v(x)=v(0)-J x/a₁`, rechts
`v(x)=v(0)-J(0,5/a₁+(x-0,5)/a₂)`. Diese Funktion ist auf jeder
Tet4-Zelle affin und auf der gemeinsamen Grenzfläche stetig; das FE-Modell
kann sie ohne Diskretisierungsfehler darstellen.

| Feld | Koeffizienten | Randwerte | unabhängige Grenzflächengröße |
|---|---:|---:|---:|
| Wärme | `k₁=2`, `k₂=4 W/(m K)` | `301/300 K` | `qₓ=8/3 W/m²` |
| Elektrostatik | `κ₁=2·10⁻⁹`, `κ₂=4·10⁻⁹ F/m` | `1/0 V` | `Dₓ=8/3·10⁻⁹ C/m²` |
| Mechanik | `Cₓₓₓₓ=2·10⁶/4·10⁶ Pa` | `uₓ=0/10⁻³ m` | `σₓₓ=8/3·10³ Pa` |

Für die Mechanik werden `u_y=u_z=0` vorgeschrieben und diagonale, positiv
definite Steifigkeitstensoren gewählt. Der Piezo-Regressionstest setzt
`e=0` und prüft, dass der monolithische Solver dieselben mechanischen und
elektrischen Grenzflächenwerte liefert; regionsweise Expansion mit
`T=T_ref` sowie `electric_thermal_strain=True` läuft durch denselben Pfad.
Dies ist ein entkoppelter Grenzfall und kein allgemeiner Nachweis für
inhomogene, voll gekoppelte Piezo-Grenzflächen.

Die Toleranzen wurden vor dem Lauf aus den SI-Größen und der exakten
stückweise affinen FE-Darstellbarkeit festgelegt: nodale Temperaturen
`3·10⁻¹¹ K` absolut und `2·10⁻¹³` relativ, Wärmestrom
`2·10⁻¹¹ W/m²` und `2·10⁻¹²`, Potential `3·10⁻¹² V` und
`2·10⁻¹²`, elektrische Flussdichte `2·10⁻²⁰ C/m²` und
`2·10⁻¹¹`, Verschiebung `2·10⁻¹⁴ m` und `2·10⁻¹¹` sowie
Spannung `2·10⁻⁷ Pa` und `2·10⁻¹¹`. Sie lassen lediglich Rundung
bei Sparse-LU und Elementrekonstruktion zu. Globale Kraft-, Wärme- und
Ladungsbilanzen werden zusätzlich in ihren jeweiligen SI-Einheiten geprüft.

Die Belege liegen in `tests/unit/test_material_regions.py`,
`tests/unit/test_thermopiezoelectric_model.py` und
`tests/patch/test_material_interfaces.py`. Die meshunabhängige 3D-
Konvergenzstudie des vollständigen `Model` liegt in
`tests/convergence/test_thermopiezo_model_convergence.py` (#141).
