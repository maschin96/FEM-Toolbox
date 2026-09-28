# Tet10-Harmonik

Mit installierter Toolbox vom Repository-Root:

```sh
python examples/tet10_harmonic/run.py --output tet10_harmonic_results
```

Ein synthetischer Tet10-Eckoszillator erhält eine komplexe Volumenlast. Die
Frequenzen 0, 0,2 und 2 Hz prüfen DC sowie Antworten unter-/oberhalb der Resonanz.
Das Skript vergleicht u/φ mit unabhängigen skalaren Formeln und exportiert
Real-/Imaginärteile, Amplituden, Phasen und einen quadratisch interpolierten
Messpunkt. Das Zielverzeichnis darf noch nicht existieren. Gmsh ist nicht nötig.

[Herleitung, Konvergenz und Grenzen](../../Doc/53_TET10_HARMONIC.md).
