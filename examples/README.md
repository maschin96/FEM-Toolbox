# Beispiele

Dieses Verzeichnis nimmt reproduzierbare Referenzmodelle auf. Ein Beispiel muss
seine Modellannahmen, Einheiten, Ausführung, erwarteten Kennwerte und deren
Toleranzen dokumentieren. Numerische Beispiele werden zusammen mit den jeweils
verifizierten Solverfähigkeiten ergänzt.

## Lineare Elastizität

`python examples/elastic_tension.py` berechnet nach der Paketinstallation einen
Zugwürfel und schreibt `elastic_tension.vtu`. Analytische Referenzen und
Toleranzen: [Lineare Elastizität](../Doc/10_LINEAR_ELASTICITY.md).
