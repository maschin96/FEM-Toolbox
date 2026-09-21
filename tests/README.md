# Teststruktur

Die Tests spiegeln die Verifikationsstufen aus
[`Doc/05_VERIFICATION_AND_VALIDATION.md`](../Doc/05_VERIFICATION_AND_VALIDATION.md):

- `unit/`: Algebra, Datenmodelle und kleine deterministische Funktionen
- `patch/`: Element- und Patch-Tests physikalischer Felder
- `convergence/`: räumliche und später zeitliche Konvergenzstudien
- `regression/`: fachlich begründete Referenzkennwerte und Fehlerregressionen

Tests werden aus dem Repository-Root mit `python -m pytest` ausgeführt.

## Manuelle VTU-Sichtprüfung für M1

Nach Änderungen am Export wird zusätzlich ein repräsentatives Netz erzeugt und
in ParaView geöffnet:

```sh
python - <<'PY'
from fem_toolbox.io import read_gmsh, write_vtu

mesh = read_gmsh("tests/fixtures/overlapping_groups.msh")
write_vtu("/tmp/fem-toolbox-overlapping-groups.vtu", mesh, binary=False)
PY
```

Zu prüfen sind fünf Punkte, zwei positiv orientierte Tetraeder, eine
Dreiecksfläche und die vier getrennten Zellmasken `ceramic`, `active`, `ground`
und `electrode`. Diese Sichtprüfung ergänzt den automatischen
Gmsh-intern-VTU-Rundreisetest; sie gilt erst als ausgeführt, wenn der PR dies
ausdrücklich dokumentiert.
