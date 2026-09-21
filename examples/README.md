# Beispiele

Dieses Verzeichnis nimmt reproduzierbare Referenzmodelle auf. Ein Beispiel muss
seine Modellannahmen, Einheiten, Ausführung, erwarteten Kennwerte und deren
Toleranzen dokumentieren. Numerische Beispiele werden zusammen mit den jeweils
verifizierten Solverfähigkeiten ergänzt.

## Lineare Elastizität

`python examples/elastic_tension.py` berechnet nach der Paketinstallation einen
Zugwürfel und schreibt `elastic_tension.vtu`. Analytische Referenzen und
Toleranzen: [Lineare Elastizität](../Doc/10_LINEAR_ELASTICITY.md).

## Stationäre Wärmeleitung

`python examples/steady_heat.py` berechnet einen Quader mit Temperaturvorgabe
und Konvektion und schreibt `steady_heat.vtu`. Parameter, analytische
Wärmeleistungen und Toleranzen: [Stationäre Wärmeleitung](../Doc/11_STEADY_HEAT.md).

## Elektrostatik

`python examples/parallel_plate.py` berechnet einen idealen Plattenkondensator
und exportiert `parallel_plate.vtu`. [Referenz und Toleranzen](../Doc/12_ELECTROSTATICS.md).
