# Modellweite Randbedingungs- und Nullraumprüfung (#138)

Status: Implementiert und mit getrennten Tet4-Komponenten, redundanten
Randwerten sowie festen und schwebenden Elektroden geprüft.

`Model(mesh=..., document=...)` und `Model.from_file(...)` validieren
Randbedingungen nach Material-, Topologie- und Schemaprüfung, aber vor
Element-Assembly oder LU-Faktorisierung. Ein Fehler ist ein
`ModelSchemaError` mit Dokumentpfad, Feld, zusammenhängender Komponente und
betroffenen Knoten-/DOF- oder Facetten-IDs. Die Konstruktion eines `Model`
ist damit ein Preflight für strukturelle Randbedingungsfehler.

Die Komponenten werden aus gemeinsam verwendeten Tet4-Knoten bestimmt;
ungenutzte Netzknoten erzeugen keinen Feld-DOF. Ihre Nummern folgen dem
kleinsten aktiven Knotenindex und stimmen für die geprüften Fälle mit den
direkten Solverdiagnosen überein.

| Feld | Erforderliche Verankerung je Komponente |
|---|---|
| `T` | Dirichlet-Temperatur oder Konvektion auf einer Außenfläche mit `h > 0` |
| `u` | Rang 6 der auf vorgeschriebene Komponenten eingeschränkten sechs Starrkörpermoden; dieselbe normierte Geometrie und Rangtoleranz `1e-10` wie im direkten Solver |
| `phi` | mindestens ein festes Potential durch Dirichlet oder feste Elektrode; eine Floating Electrode mit Ladung allein fixiert keine Gauge |

Gleiche Dirichlet-Vorgaben für denselben Feld-DOF sind zulässig. Abweichende
Werte sind ein Fehler mit Verweis auf beide Einträge. Bei Elektroden werden
doppelte Namen, gemeinsame Knoten verschiedener Elektroden, Floating-/Dirichlet-
Kollisionen und abweichende feste Potentiale erkannt. Eine feste Elektrode und
eine gleichwertige explizite Potentialvorgabe bleiben zulässig. Ein auf
derselben Elektrodenfläche vorgeschriebener `surface_charge`-Fluss ist eine
widersprüchliche zweite elektrische Randbedingung und wird sowohl im Model als
auch im direkten Piezo-Solver abgelehnt. Mehrfachbelegte Facetten innerhalb
einer Flächenlastart werden mit Auswahlname und Facetten-ID berichtet.

Die Prüfung bewertet **Struktur**, keine Lastkompatibilität oder vollständige
Matrixkondition. Ein unterbestimmtes reines Neumann-Problem kann im direkten
Solver zusätzlich als kompatibel oder inkompatibel anhand seiner Wärme- bzw.
Ladungsbilanz diagnostiziert werden. Die Modellprüfung ersetzt dessen
physikalische Residuen- und Bilanztests nicht. Sie setzt keine Referenzwerte
automatisch und verändert keine Randbedingungen.

`tests/unit/test_model_preflight.py` vergleicht die Feld-/Komponentendiagnosen
für Wärme, Elektrostatik und Mechanik mit den direkten Solvern. Gültige
Konvektion, identische Mehrfachvorgaben und eine ladungsfreie Floating
Electrode bleiben lösbar. Falsche Vorzeichen oder Bilanzen würden weiterhin
von den vorhandenen Patch-, Benchmark- und Residualtests erkannt.
