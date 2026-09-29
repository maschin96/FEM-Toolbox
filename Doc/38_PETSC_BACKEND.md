# Optionaler PETSc-Direktbackend (#34b)

Teilbeitrag zu [Issue #34](https://github.com/maschin96/FEM-Toolbox/issues/34)
auf M7, kein Release. [Dokument 39 (#34c)](39_PETSC_FIELDSPLIT.md) ergänzt
FieldSplit/Schur und iterative Verfahren. Der [interne Vertrag](37_LINEAR_BACKEND.md) und die bestehenden
physikalischen Einheiten, Randbedingungen und Residuen gelten unverändert.

## Auswahl und Installation

```python
from fem_toolbox.solvers.petsc_backend import PetscBackend
from fem_toolbox.solvers.heat import solve_heat

backend = PetscBackend()
# mesh, material und conditions wie im bestehenden Wärmebeispiel:
result = solve_heat(mesh, material, conditions, backend=backend)
```

Alle acht direkten Solver-Hooks aus Dokument 37 akzeptieren das Backend.
`None` bleibt SciPy. `Model`/TOML konfigurieren es noch nicht. Core und selbst
Adaptermodul importieren ohne PETSc; erst `PetscBackend()` lädt `petsc4py.PETSc`
und diagnostiziert fehlende oder inkompatible Installationen mit `ValidationError`.
Backendkennung und Provenienz enthalten PETSc-Version, Scalar-Typ und Verfahren.

PETSc ist **keine Core-Abhängigkeit**. Seine native Installation und Scalar-Auswahl
werden bewusst getrennt von `pip install fem-toolbox` verwaltet. Ein generisches
pip-Extra würde weder die Scalar-Konfiguration noch die ABI-Kompatibilität mit
NumPy oder einer vorhandenen PETSc-Installation festlegen. Getesteter Einstieg
ist deshalb die separate Linux-Umgebung aus `ci/petsc.Dockerfile`:

```sh
docker build -f ci/petsc.Dockerfile -t fem-toolbox-petsc:test .
docker run --rm -v "$PWD:/work:ro" fem-toolbox-petsc:test \
  python ci/run_petsc_tests.py --scalar real -q tests/petsc -p no:cacheprovider
docker run --rm -v "$PWD:/work:ro" fem-toolbox-petsc:test \
  python ci/run_petsc_tests.py --scalar complex -q tests/petsc -p no:cacheprovider
```

Das Image installiert die Toolbox regulär; der Runner verlangt ihren Import aus
der virtuellen Umgebung, ohne lokale Quellbaum-/egg-info-Metadaten zu verwenden.
Nach Quelländerungen das Image neu bauen. Der Runner wählt **vor dem ersten Import** genau einen der beiden Ubuntu-Builds
und prüft dessen tatsächlichen Scalar-Typ. Die native Umgebung darf nicht im
laufenden Python-Prozess gewechselt werden. Eigene Installationen nach der
[offiziellen petsc4py-Anleitung](https://petsc.org/release/petsc4py/install.html)
sind möglich, aber keine zusätzliche getestete Supportzusage.

| Komponente | Getestete Konfiguration |
|---|---|
| Betriebssystem | Ubuntu 24.04, lokaler Docker-Test Linux arm64 |
| Python | 3.12.3 |
| PETSc / petsc4py | 3.19.6 / 3.19.6, Ubuntu-Pakete real und complex |
| Scalar / Index | float64 bzw. complex128 / int32 |
| NumPy / SciPy | 1.26.4 / 1.14.1 |
| Verfahren | ein Prozess, COMM_SELF, KSP PREONLY + PC LU, PETSc-eigene SeqAIJ-LU |

Die CI führt dieselben beiden Scalar-Konfigurationen auf Linux x86_64 aus und
protokolliert die tatsächlich installierten Versionen. Das Basisimage ist per
Digest, die Python-Testpakete sind per Version fixiert; Ubuntu-Sicherheitspakete
werden weiter aktualisiert. Die bestehende Core-Matrix Python 3.11–3.13 bleibt
getrennt davon. Native macOS-/Windows-Installationen und andere Versionen werden
hier nicht als getestet ausgewiesen. PETSc und petsc4py verwenden BSD-2-Clause;
siehe [Lizenz des Projekts](https://petsc.org/release/install/license/).

## Numerik, Scalar-Vertrag und Lebensdauer

- Gemeinsame Validierung, Skalierung und unskalierte feldweise Residuen stammen
  aus `_LinearData`; die Faktorisierung und linearen Solves laufen in PETSc.
  Es gibt keinen SciPy-Fallback und keinen Wechsel auf externe LU-Pakete.
- Ein complex128-Build verarbeitet komplexe Daten direkt. Bei float64-PETSc
  löst ein reeller Operator getrennte Real-/Imaginärteile einer komplexen RHS.
  Ein komplexer Operator wird äquivalent in
  `[[Re A, -Im A], [Im A, Re A]]` eingebettet. Lösungen werden rekonstruiert und
  gegen das ursprüngliche komplexe System geprüft. Keine Scalar-Anteile gehen
  stillschweigend verloren; reelle Eingaben ergeben weiter reelle Ergebnisse.
- `kind="positive_definite"` prüft Hermiteschheit und positive Eigenwerte vor
  der PETSc-LU. Dieser dichte Nachweis ist für die bereits begrenzten M6-Blöcke
  vorgesehen und kostet kubischen Aufwand. Das bestehende 2048-DOF-Limit bleibt.
- RCM und eine Nichtnulldiagonal-Permutation unterstützen auch rein imaginäre
  Operatoren mit strukturell fehlender Diagonale. Explizite Nulleinträge erfüllen
  den SeqAIJ-Symbolikvertrag; sie ändern den Operator nicht. Es gibt keinen
  Diagonalshift. Native LU hat andere Pivotierung als SuperLU und kann schwierige
  Systeme ablehnen; identische Pivotfolgen oder Konditionsschätzungen werden
  nicht versprochen. Fehler, nichtendliche Pivots, überschrittene Pivotgrenzen
  und schlechte Originalresiduen werden diagnostiziert.
- `SolveReport.method` ist `preonly_lu` bzw. `preonly_lu_real_block`, `reason`
  enthält den positiven PETSc-Konvergenzcode; `iterations=None` kennzeichnet den
  Direktlöser. Ein leerer Operator braucht keine native Faktorisierung.
- Faktoren besitzen ihre Matrixdaten. Mehrere RHS und Wiederverwendung sind
  möglich; `close()` gibt native Ressourcen frei und sperrt weitere Solves.
  Ohne explizites Schließen erledigt der Finalizer dies bei Freigabe des Faktors.
  Der Solver liest keine globale PETSc-Optionsdatenbank ein.

Keine MPI-Verteilung, FieldSplit, iterative Präconditionierung, Performance-
oder experimentelle Validierungsbehauptung. Der modale Eigenwertsolver bleibt
SciPy `eigh`; nur seine lineare elektrische Elimination verwendet PETSc.

## Verifikationsnachweise

Vor Umbau: 428 Core-Tests grün. TDD-Red: Importtest ohne optionale Abhängigkeit
und elf Tests im echten reellen PETSc-Runtime scheiterten am fehlenden Adapter.
Anschließend bestanden sie. Ein zusätzlicher Red-Test mit rein imaginärem
Operator deckte fehlende strukturelle Diagonaleinträge auf; nach Ergänzung und
Nichtnulldiagonal-Permutation ist auch dieser Test grün.

Je **18 echte PETSc-Tests** bestehen lokal für float64 und complex128:
Skalierung, Faktorbesitz, Wiederverwendung, komplexe Matrix/RHS, leere Systeme,
Singularität, kleine Pivots, ungültige Eingaben und explizit gesperrter SciPy-LU-
Fallback. Physikalische Tests prüfen affine Wärme-/Verschiebungs-/Potentialfelder,
homogene Speichererwärmung, SciPy-Parität für feste/geladene schwebende Elektroden
samt Energie, Ladung und Kraft/Moment sowie einen unabhängig hergeleiteten
Tet4-Einmodenoszillator mit Moden, Newmark-Energieerhaltung und komplexer Harmonik.

Algebraische Referenzen O(1): `rtol=0, atol=1e-13`; physikalische Parität
`rtol=1e-11, atol=1e-12` für die verwendeten SI-Testgrößen. Wärme bei 300 K:
`atol=1e-10 K` berücksichtigt Subtraktion des Temperaturniveaus. Die unabhängige
Oszillatorformel und Energie verwenden `rtol=1e-12`. Feldresiduen müssen jeweils
unter `1e-10` liegen. Bestehende wissenschaftliche Tests/Toleranzen bleiben
unverändert; keine neue Diskretisierung und kein neues Zeitintegrationsverfahren.

Core ohne PETSc: **429 bestanden, 18 explizit übersprungen**. Die beiden
PETSc-CI-Jobs setzen `FEM_REQUIRE_PETSC=1`: fehlende native Laufzeit ist dort ein
Fehler, kein Skip. Ruff, Formatierung, mypy und Diff-Prüfung gehören zusätzlich
zur Abnahme. CI- und Reviewstatus werden am PR nachgewiesen; diese Dokumentation
behauptet keine unabhängige Reviewfreigabe oder Fertigstellung des Gesamtissues.
