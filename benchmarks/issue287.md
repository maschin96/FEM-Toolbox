# Issue #287 — vektorisierte Randtopologie

Nachweis vom 08.11.2026 für [#287](https://github.com/maschin96/FEM-Toolbox/issues/287).
Basis: `8820e0075be9c0a3d6db0660cfa6cfa8949ab604` mit den hier geprüften Änderungen.

## Verfahren und Kompatibilität

`inspect_boundary_topology` prüft Tet4-Determinanten und skaliertes Volumen in
Blöcken von höchstens 65.536 Zellen. Die Entartungstoleranz bleibt `1e-12`.
Nichtendliche und grenznahe Ergebnisse werden in Eingabereihenfolge durch
`Tet4Geometry` geprüft; dessen Fehlerklasse, externe Zell-ID und Diagnosewerte
bleiben maßgeblich. Bei Tet10 bleibt die bisherige zellweise Reihenfolge von
Eckgeometrie- und Geradkantigkeitsprüfung erhalten.

NumPy sortiert die drei Eckknoten jeder Fläche und gruppiert die Inzidenzen
lexikographisch und stabil. Besitzer bleiben in Zellreihenfolge. Nur äußere
Flächen benötigen öffentliche `BoundaryFace`-Objekte; Flächeninhalt und äußere
Normale werden dort gesammelt berechnet. Normale bleiben unveränderlich.
Öffentliche Rückgabetypen und Oberflächengruppen bleiben erhalten.

Die eingefrorene Referenz in `tests/reference/legacy_boundary_topology.py`
entspricht der Implementierung von v1.4.0, lediglich ihr Modulkommentar wurde
angepasst. Der SHA-256 des ursprünglichen Moduls ist
`5c55b80e65405e2e65079423fb6191beda6b83461cde114b149d37db0b5c74f5`;
auch `elements/tet4.py` ist gegenüber v1.4.0 unverändert.

## Vorab gewählte Vergleichsgrenzen

Inzidenzen, Reihenfolge, Orientierung, Knoten, Besitzer, externe IDs, Gruppen,
Flächen-/Dreieckszuordnung und Fehlerdiagnosen müssen exakt übereinstimmen.
Für Flächeninhalte gilt relative Toleranz `2e-14`, ohne absolute Aufweichung;
für Normalen relative und absolute Toleranz `2e-14`. Diese Grenzen decken die
unterschiedliche Double-Summationsreihenfolge in den gut konditionierten
Referenznetzen ab. Sie sind kein allgemeines physikalisches Fehlerbudget.

Die Tests vergleichen vier zufällige Nummerierungen und Zellreihenfolgen für
verzerrte strukturierte Tet4-/Tet10-Netze, überlappende Gruppen und externe IDs.
Hinzu kommen leere/partielle Ränder, nichtmannigfaltige Inzidenz ohne
Oberflächenzeile und Entartung einschließlich benachbarter Double-Werte an der
Toleranzgrenze. Bestehende Tests für fremde/innere/fehlende Flächen,
nichtkonforme Triangle6 und gekrümmte Tet10 bleiben wirksam.
Der Test gegen zellweise `Tet4Geometry`-Erzeugung schlägt mit dem installierten
Altstand fehl und besteht mit `PYTHONPATH=src` am neuen Stand.

## Messung

[Messskript](boundary_topology.py) · [Rohdaten](evidence/issue287.json).
Aufruf aus dem Repository-Root:

```sh
PYTHONPATH=src:. .venv/bin/python benchmarks/boundary_topology.py \
  --sizes 12 26 55 --repeats 1 --output results/issue287/boundary_topology.json
```

Je Variante und Netz ein frischer Prozess, angeforderte BLAS-/OpenMP-Threadlimits
auf 1. Gemessen wird ausschließlich der vollständige Funktionsaufruf mit
`require_complete=True`, einschließlich Validierung und Rückgabeaufbau.
Netzaufbau, Imports und Prüfsummenbildung sind ausgeschlossen. Spitzenspeicher
ist hingegen das Prozessmaximum einschließlich Netzaufbau und Imports.
Einheitswürfel mit sechs Tet4 je kartesischer Zelle und vollständigem Außenrand.
Python 3.13.7, NumPy 2.4.6, macOS 27.2, arm64.

| Zellen | v1.4.0 skalar | Vektorisiert | Faktor | Prozess-RSS alt/neu |
| ---: | ---: | ---: | ---: | ---: |
| 10.368 | 1,048 s | 0,0445 s | 23,6 | 71,5 / 57,6 MiB |
| 105.456 | 11,640 s | 0,4857 s | 24,0 | 423,1 / 252,2 MiB |
| 998.250 | 115,195 s | 5,3222 s | 21,6 | 3.534,3 / 1.783,6 MiB |

Die Prüfsummen aller diskreten Inzidenz-/Außenflächenfelder stimmen pro Netz
überein; beide Varianten liefern Gesamtfläche 6. Nach der abschließenden
Absicherung der Toleranzgrenze wurden die vektorisierten Messungen wiederholt;
die skalaren Messungen verwenden unverändert dieselbe Referenz. Die Rohdaten
halten diese Zuordnung und beide Quellhashes fest. Einzelmessungen belegen
diesen lokalen Effekt, keine plattformübergreifende Zeitgarantie. Tet10 wird
auf Ergebnisparität geprüft; ein Tet10-Beschleunigungsfaktor wird nicht behauptet.

## Qualitätsprüfung

- 42 gezielte Topologie-/Paritätsprüfungen bestehen.
- Ruff-Prüfung, Formatprüfung, strikte Typprüfung aller 70 Produktionsmodule
  und `git diff --check` bestehen.
- Der vollständige erste Testlauf ergab 854 bestanden, 83 übersprungen und zwei
  Fehler in `test_handbook_build.py`. Beide Fehler sind auch ohne den neuen
  Produktionscode reproduzierbar: Die bereits vorhandenen lokalen Dateien
  `Doc/79_PRODUCT_VISION_AND_NEXT_TEN_MILESTONES.md` und `Doc/80_PRODUCT_BACKLOG.md`
  fehlen in `Doc/latex/coverage.json`. Die vorhandenen Planungsänderungen wurden
  durch dieses Issue nicht bearbeitet.
- Abschließender Regressionslauf am finalen Stand: **861 bestanden,
  83 übersprungen, 2 ausgeklammert** in 113,61 s. Ausgeklammert wurden exakt
  die beiden oben genannten, separat bestätigten bestehenden Handbuchfehler:

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q \
  --deselect tests/unit/test_handbook_build.py::test_repository_matrix_is_complete \
  --deselect tests/unit/test_handbook_build.py::test_planned_topics_have_existing_targets_and_owners
```

Die 83 übersprungenen Tests benötigen das optionale PETSc-Backend.
Keine Veröffentlichung und kein Release-Gate werden mit diesem lokalen Nachweis
behauptet.

## Prüfung des GitHub-PR-Stands auf aktuellem main

Für den PR wurde ausschließlich der Umfang von #287 auf
`20243494ffa9fd0317835f5d14d0f071dc9acd9a` (main, v1.6.0/M11) übertragen.
Der Produktionscode der Randtopologie war dort gegenüber der obigen Messbasis
unverändert; die archivierten Messungen beziehen sich weiterhin auf ihre
angegebenen Quellhashes und wurden nicht umgeschrieben.

In einer separaten Umgebung wurde das Wheel des PR-Stands mit
`pip install --no-deps --no-build-isolation .` installiert. Die vorhandenen
NumPy-/SciPy-/Prüfwerkzeuge wurden wiederverwendet, das installierte
`fem-toolbox`-Paket stammt aus dem neuen Wheel (Version 1.6.0).
Der vollständige Lauf `python -m pytest -q` gegen diese Installation ergab
**890 bestanden, 83 wegen PETSc übersprungen, keine ausgeklammerten Tests**
in 111,98 s. Auch die beiden zuvor lokal fehlschlagenden Handbuchprüfungen
bestehen auf diesem sauberen main-Stand.

Zusätzlich bestehen `python -m ruff check .`,
`python -m ruff format --check .`, `python -m mypy src/fem_toolbox`
(74 Produktionsmodule) und `git diff --check`. Die ursprüngliche lokale
Installation 1.4.0 wurde dabei nicht verändert.
