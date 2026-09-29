# Tet10-Newmark

Vom Repository-Root mit installierter Toolbox:

```sh
python examples/tet10_newmark/run.py --output tet10_newmark_results
```

Ein synthetischer Tet10-Eckoszillator schwingt frei aus 0,01 m Anfangsauslenkung.
Das Skript prüft die diskrete Newmark-Phase, konstante Energie und algebraisches
Potential und exportiert VTU/PVD, Manifest und Messpunkt-CSV. Ein Zielverzeichnis
darf noch nicht existieren. Keine Gmsh-Installation nötig.

[Herleitung, Toleranzen und Grenzen](../../Doc/52_TET10_NEWMARK.md).
