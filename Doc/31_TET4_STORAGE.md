# Tet4-Masse und Wärmekapazität

Status: Implementierter erster Teil von [#29](https://github.com/maschin96/FEM-Toolbox/issues/29), M6.
Noch kein transienter Solver; Materialregionen, Zeitbasis und Serienvertrag aus
#29 sind weiterhin offen. Die stationären Materialkonstruktoren bleiben unverändert.

## Direkte Element-API

`fem_toolbox.elements.storage.mass_matrix(geometry, density, lumped=False)`
liefert eine mechanische 12×12-Matrix in kg, in ursprünglicher lokaler
Knotenreihenfolge mit xyz-Komponenten pro Knoten. Elektrische DOFs tragen keine Masse.
`capacity_matrix(geometry, density, specific_heat, lumped=False)` liefert die
thermische 4×4-Matrix in J/K. `lumped` ist jeweils ein Keyword-Argument.
Dichte in kg/m³ und spezifische Wärmekapazität in J/(kg K) sind positive,
endliche, elementweise konstante reelle Skalare. Ungültige Eingaben und nicht
repräsentierbare Speicherkoeffizienten lösen `ValidationError` aus.

Für baryzentrische Tet4-Formfunktionen gilt durch exakte Simplexintegration
`∫ Ni² dV = V/10` und `∫ Ni Nj dV = V/20` für `i != j`.
Daraus folgt `S = V/20 (I + 11ᵀ)`, `M = rho S ⊗ I3` und
`C = rho c S`. Explizites `lumped=True` verwendet die Zeilensummen:
`S_lumped = V/4 I`. Beide Varianten erhalten Masse pro Translationsrichtung
und Gesamtwärmekapazität. Die gelumpte Variante integriert allgemeine
quadratische Felder nicht exakt.

## Verifikation und Grenzen

`tests/unit/test_storage.py` prüft die analytischen Diagonal-/Nebendiagonalintegrale,
Gesamtwärmekapazität, kinetische Energie einer starren Translation,
Positivdefinitheit, Symmetrie, quadratische Feldintegration der konsistenten
Kapazität, Knotenumordnung, Volumenskalierung und ungültige Eingaben.
Relative Toleranz 1e-14 deckt Rundung der gut konditionierten Geometrie mit
Kantenlängen 1–3 m ab; algebraische Symmetrie wird ohne absolute Toleranz geprüft.
Raum-/Zeitkonvergenz und globale dynamische Bilanzen benötigen die folgenden
Solver-Arbeitspakete und werden hier nicht beansprucht.

TDD: Nach Setzen von `PYTHONPATH=src` scheiterte der erste Testlauf am fehlenden
Modul `fem_toolbox.elements.storage` (Red); nach Implementierung bestanden die
gezielten Tests (Green). Das lokale virtuelle Environment hatte das Paket
anfangs nicht installiert, daher verwenden die Prüfungen explizit den Quellpfad.
