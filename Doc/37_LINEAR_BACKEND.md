# Interner linearer Backend-Vertrag (#34a)

Status: implementiert und lokal geprüft im M7-Issue-Branch. Teilbeitrag zu
[Issue #34](https://github.com/maschin96/FEM-Toolbox/issues/34), Entscheidung
[ADR-024](07_RISKS_AND_DECISIONS.md#adr-024--interner-linearer-backend-scalarvertrag-34a).
PETSc-Direktlösen ergänzt [Dokument 38 (#34b)](38_PETSC_BACKEND.md).
FieldSplit bleibt #34c. Die SciPy-Referenz benötigt keine neue Abhängigkeit.

## Verantwortlichkeiten und Schnittstelle

`fem_toolbox.solvers.linear_backend` ist eine interne Erweiterungsgrenze.
`LinearBackend` beschreibt `name` und
`factorize(operator, *, fields, scale=None, kind="general",
pivot_tolerance=0, relative_pivot=False)`. Das resultierende
`FactorizedOperator.solve(rhs)` liefert `LinearSolution(values, report)`.
`solve_linear(..., backend=None)` fasst beide Aufrufe für einmaliges Lösen
zusammen. `None` wählt eine neue SciPy-Implementierung; keine globale Auswahl,
kein impliziter Cache und kein stiller Fallback bei Fehlern.

- Operator: endliche quadratische dichte oder SciPy-Sparse-Matrix `(n,n)`.
  Intern eigenes CSC-Snapshot in float64 oder complex128; der bisherige
  reelle `assembly.LinearSystem`-Vertrag bleibt unverändert.
- Rechte Seite: endlich, `(n,)` oder `(n,nrhs)`, float64/complex128. Eine
  komplexe RHS bei reellem LU-Faktor wird durch getrennte reelle und imaginäre
  Solves behandelt; komplexe Anteile werden nicht verworfen.
- `fields`: nichtleere Feldnamen auf eindimensionale ganzzahlige lokale
  Zeilenindizes. Alle Zeilen genau einmal, keine Überlappung oder Lücke.
  Leere Felder sind zulässig. Es handelt sich ausdrücklich nicht um originale
  globale DOF-Nummern vor Randbedingungen oder Elektrodenkondensation.
- `scale`: endlicher positiver reeller Vektor `(n,)`, standardmäßig Eins.
  Für `S=diag(scale)` wird `S A S y = S b` gelöst und `x=S y` zurückgegeben.
  Vorzeichenänderungen (Elektrostatik) werden auf beide Gleichungsseiten im
  Physikadapter angewendet, nicht im Backend erraten.
- `kind="general"`: Sparse-LU; absolute beziehungsweise relative Schranke
  für die U-Diagonalpivots, falls `pivot_tolerance>0`. Das ist ein numerischer
  Schutz, keine Konditionszahlschätzung. Bestehende stationäre Pivotgrenzen
  bleiben unverändert. `kind="positive_definite"`: dichte Cholesky nach
  Hermitesch-Prüfung; Pivotgrenze ist hier nicht anwendbar. Diese Variante ist
  für die bereits begrenzten dichten M6-Teilprobleme vorgesehen.

Der Faktor kopiert Matrix, Skala und Indexmengen. Änderungen an Eingaben
verändern ihn nicht. Wiederverwendung ist explizit möglich, wenn Operator,
Feldpartition und Skalierung unverändert sind; sonst neu faktorisieren.
Eine leere `(0,0)`-Matrix liefert eine leere Lösung ohne Bibliotheksfaktorisierung.
Ungültige Daten, Singularität, fehlende positive Definitheit, nichtendliche
Zwischenwerte oder fehlgeschlagene Residuenprüfung führen zu `ValidationError`.
Ein alternatives Backend muss diesen Vertrag einschließlich Validierung und
Residuenprüfung erfüllen; der Hook vertraut seiner Implementierung.

## Bericht und physikalische Residuen

`SolveReport` enthält Backendkennung, Verfahren (`sparse_lu`/`cholesky`),
`converged=True`, `reason="direct_solve"`, `iterations=None`, den unveränderlichen
Residualvektor beziehungsweise die Residualmatrix sowie Gesamt- und
feldweise relative Residuen. Fehler liefern keine scheinbar erfolgreiche
Lösung. Die Feldkennwerte sind über eine unveränderliche MappingProxy verfügbar.

Für jedes Feld I gilt, auf dem **unskalierten übergebenen System**:

```text
r = A x - b
eta_I = ||r[I]|| / || (|A| |x| + |b|)[I] ||
```

Matrix-RHS verwendet die Frobeniusnorm über Zeilen und Spalten. Normen werden
vor dem Quadrieren gemeinsam skaliert, damit große/kleine SI-Zahlen keinen
unnötigen Overflow/Underflow verursachen. Nenner und Zähler null ergeben null;
ein von null verschiedener Zähler bei Nenner null ergibt Unendlich. Gesamt-
und jeder Feldwert müssen höchstens `1e-9` sein. Gerade bei gemischten N-/C-
Größen darf ein kleiner Gesamtwert kein fehlerhaftes elektrisches Feld verdecken.
Der rohe Residualvektor behält die Gleichungseinheiten der Eingabe.

Das Backend sieht bereits eliminierte Randwerte. Reaktionen, physikalische
Bilanzen und Residuen des ursprünglichen vollständigen Systems werden deshalb
weiterhin nach Rekonstruktion in jedem Physiksolver berechnet und geprüft.
Der Backendbericht ersetzt diese vorhandenen Nachweise nicht. Er steht dem
Aufrufer des internen Backends zur Verfügung; bestehende öffentliche
Physikergebnisobjekte und Exportschemas erhalten keine neue Pflichtkomponente.

## Angebundene Pfade

| Pfad | Backend-System und lokale Felder | Zusätzlicher unveränderter Nachweis |
|---|---|---|
| Stationäre Wärme | freie Temperaturzeilen, `T` | Temperaturgradient, Wärmefluss und Leistungsbilanz |
| Elastizität | freie Verschiebungen, `u` | Starrkörperprüfung, Kräfte/Momente und Energie |
| Elektrostatik | freie vorzeicheninvertierte Zeilen, `phi` | Gauge, ursprüngliches elektrisches Residuum und Ladung |
| Piezo stationär | freie kondensierte Zeilen, `u`/`phi` | Feldresiduen, Elektrodenladung und Kraft-/Ladungsbilanz |
| Transiente Wärme | ein freies System pro Backward-Euler-Schritt, `T` | Speicher- und Reservoirbilanz |
| Piezo modal | elektrische Schur-Elimination, `phi` | Eigenproblem weiterhin `scipy.linalg.eigh`, Modalresiduen |
| Piezo Newmark | elektrisches Eliminationssystem `phi`, Masse und Newmark-System `u` | konsistente Anfangswerte, Trägheit, Energie und Feldresiduen |
| Piezo harmonisch | elektrisches Eliminationssystem `phi`, dynamisches Schur-System `u` | Resonanzprüfung und komplexe rekonstruierte Feldresiduen |

In der Harmonik ersetzt allgemeine Sparse-LU den vorherigen dichten symmetrischen
Direktaufruf. Cholesky bleibt für positive M6-Blöcke erhalten. Skalenmultiplikation
statt Division kann Rundung ändern; bitgleiche alte Lösungen werden nicht
versprochen. Physikalische Referenzen und festgelegte Toleranzen bleiben gleich.
Das 2048-DOF-Limit und der elektrische DAE-/Phasorvertrag bleiben erhalten.

Die acht direkten Solve-Einstiegspunkte besitzen den additiven experimentellen
Keyword-Hook `backend=None`; bestehende Aufrufe bleiben gültig. Ohne explizite
Auswahl werden weiterhin SciPy-Löser verwendet. Stationäre Provenienz behält
bei Default-Auswahl ihre bisherige Backendkennung; explizite Backends werden
mit `backend.name` aufgezeichnet. M6-Provenienz ergänzt `linear_backend`.
`Model` und die sequenziellen Kopplungs-Wrapper verwenden weiterhin den Default;
TOML-Backend-Auswahl und öffentliche Stabilitätszusage für diesen internen
Vertrag folgen dem späteren PETSc-Supportvertrag. Der Modensolver-Hook betrifft
nur die lineare Elimination, nicht den Eigenwertsolver.

## Beispiel der algebraischen Forschungsgrenze

```python
from fem_toolbox.solvers.linear_backend import solve_linear

result = solve_linear(
    [[4.0, 1.0], [1.0, -2.0]],
    [5.0, 8.0],
    fields={"u": [0], "phi": [1]},
    scale=[0.5, 2.0],
)
# x = [2, -3]; 4*2-3=5, 2-2*(-3)=8
print(result.values, result.report.field_relative_residuals)
```

## Prüfungen und Grenzen

Vor dem Umbau: 400 bestehende Tests bestanden. TDD-Red: neue Backend-Vertragstests
scheiterten am fehlenden Modul, anschließend grün; vier Integrationstests
scheiterten am fehlenden Backend-Hook, anschließend grün. Unabhängig aufgestellte
Mini-Systeme prüfen reelle/komplexe Lösung, Skalierung, mehrere RHS, Faktorbesitz,
Feldpartitionen, Singularität, ungültige Eingaben und leere Systeme.
Algebraische Toleranz `atol=1e-14, rtol=0` für gut konditionierte Größen O(1).
Ein gezielt falscher elektrischer Teilsolve an der SciPy-Infrastrukturgrenze
wird trotz kleinem Gesamtfehler abgelehnt; ein Backendfehler wird weitergegeben.

Integrationstests prüfen affine Einzelfelder, lokale Feldindizes nach festen
und schwebenden Elektroden sowie austauschbare lineare M6-Teilprobleme. Bestehende
Patch-, Bilanz-, Energie-, Raum-/Zeitkonvergenz-, Modal-, Harmonik- und Modeltests
bleiben verpflichtend. Das 37a-Profiling wurde an die neue LU-Grenze angepasst;
archivierte Messdaten werden nicht umgeschrieben. Keine Leistungsverbesserung,
keine PETSc-Parität, keine MPI-Skalierung oder experimentelle Validierung behauptet.

Prüfbefehle aus dem Repository-Root:

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

Lokales Ergebnis: **428 Tests bestanden**, Ruff-Lint und Formatprüfung
bestanden, mypy für 51 Quelldateien bestanden, Diff-Prüfung bestanden.

Lokaler Import weiterhin über `PYTHONPATH=src`; dies ist ein Quelltest.
CI/Artefakte und Reviewstatus werden am PR nachgewiesen. Kein unabhängiger
Reviewer beteiligt; keine Reviewfreigabe behauptet. #34 bleibt bis zur
vollständigen Abnahme einschließlich FieldSplit (#34c) offen.
