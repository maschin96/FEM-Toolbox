# Iterativer PETSc-Backend mit FieldSplit/Schur (#34c)

M7-Teilbeitrag zu [Issue #34](https://github.com/maschin96/FEM-Toolbox/issues/34).
Aufbau auf dem [Backend-Vertrag](37_LINEAR_BACKEND.md) und dem
[optionalen PETSc-Direktadapter](38_PETSC_BACKEND.md), Entscheidung ADR-026.
Kein Release und keine Aussage über verteilte Skalierung oder Laufzeitgewinn.

## Auswahl und Vertrag

```python
from fem_toolbox.solvers.petsc_backend import PetscFieldSplitBackend
from fem_toolbox.solvers.piezoelectric import solve_piezoelectric

backend = PetscFieldSplitBackend(
    rtol=1e-13,
    atol=0.0,
    max_iterations=500,
    restart=50,
)
# mesh, material, conditions und electrodes wie im bestehenden Piezo-Beispiel:
result = solve_piezoelectric(
    mesh,
    material,
    conditions,
    electrodes=electrodes,
    backend=backend,
)
```

`PetscBackend()` bleibt der direkte PREONLY/LU-Adapter; `backend=None` bleibt
SciPy. Der neue Adapter ist wie sein Vorgänger experimentell und explizit an
allen acht direkten Solve-Einstiegspunkten auswählbar. Core-Import bleibt ohne
PETSc möglich. `Model`/TOML wählen weiterhin SciPy; eine neue deklarative
Konfiguration ist nicht Teil dieses internen Backend-Schritts.

Der Adapter verwendet **FGMRES mit Rechtspräkonditionierung**. Das ursprüngliche
Piezo-System ist symmetrisch-indefinit, aber die approximative Blockfaktorisierung
muss weder symmetrisch noch positiv definit sein. Deshalb werden weder CG noch
MINRES auf diesen präkonditionierten Operator angewendet. Die frühere allgemeine
GMRES-Einschränkung im Architekturentwurf wird für diesen Fall präzisiert.

| Freies kondensiertes System | Verfahren |
|---|---|
| Zwei nichtleere Felder | FGMRES + PETSc FieldSplit, volle Schur-Blockfaktorisierung, A11 als Schur-Präkonditionierer |
| Ein nichtleeres Feld, auch nach vollständiger Vorgabe des zweiten | FGMRES + native LU-Präkonditionierung |
| Keine freien Zeilen | Leere Lösung, null Iterationen, keine native Faktorisierung |
| Mehr als zwei nichtleere Felder | Früher `ValidationError`; kein implizites Zusammenlegen |

## Blockstruktur und Scalar-Abbildung

Für `A=[[A00,A01],[A10,A11]]` ist der Schur-Operator
`S=A11-A10 A00⁻¹ A01`. Die volle Blockfaktorisierung verwendet die beiden
Kopplungsblöcke. Als Approximation von `S⁻¹` dient `A11⁻¹`. Beide Teilinversen
werden mit PETSc SeqAIJ-LU berechnet. Der volle Schur-Operator wird dabei nicht
dicht assembliert. Das ist ein nachvollziehbarer erster Blockpräkonditionierer,
keine AMG- oder Großmodell-Performancezusage. Für ein einzelnes Feld ist LU
bereits eine sehr genaue Präkonditionierung; die äußere Iteration braucht
normalerweise einen Schritt. Das ist ein expliziter Pfad, kein Fehler-Fallback.

Die Feldindizes stammen ausschließlich aus der übergebenen lokalen Partition
**nach** Randbedingungs- und Elektrodenkondensation. Nicht zusammenhängende
Indexmengen sind erlaubt. `u` wird zuerst angeordnet, andere Feldnamen danach
lexikographisch; die Einfügereihenfolge des Python-Mappings ändert dies nicht.
Leere Felder bleiben im Originalresiduenbericht enthalten.

Auf einem reellen PETSc-Build bleibt die komplexe Blockeinbettung aus #34b
erhalten. Für ein Feld mit Indizes `I` wird dessen native Menge `I ∪ (I+n)`;
Real-/Imaginärteile werden nicht versehentlich als verschiedene Physikfelder
behandelt. complex128-PETSc arbeitet direkt mit den ursprünglichen Indizes.
Beide Varianten erhalten komplexe Matrixeinträge und RHS vollständig.

Die LU-Blöcke müssen numerisch invertierbar sein. A11=0 oder singuläre
Diagonalblöcke werden nicht durch eine alternative Schur-Näherung oder
Diagonalverschiebung repariert. Pivotgrenzen werden pro LU-Block geprüft;
sie sind keine Konditionszahl des Gesamtoperators. `kind="positive_definite"`
behält die separate Prüfung und die bisherige Pivotgrenzen-Semantik.

## Nullräume und Fehler

Die bestehenden Physikadapter prüfen die Eindeutigkeit vor dem Backend-Aufruf:

- Wärme: fehlende Temperaturreferenz; kompatible Neumann-Bilanz wird von einer
  inkompatiblen Nettoleistung unterschieden (`ThermalNullspaceError`).
- Elektrostatik: fehlende Gauge; kompatible Ladungen werden von einer
  inkompatiblen Gesamtladung unterschieden (`ElectricGaugeError`).
- Mechanik/Piezo: verbleibende Starrkörpermoden pro verbundenem Gebiet und
  erforderliche elektrische Referenzen bleiben Fehler.

Es gibt keine automatische Fixierung, RHS-Projektion oder implizite Änderung
der Lastbilanz. Ein expliziter PETSc-Nullraumbasis-/Projektionsmodus wird nicht
hinzugefügt: Die unterstützten öffentlichen Analyseverträge fordern bereits
entfernte Nullräume. Das präzisiert die frühere Planung „Nullraumübergabe“ aus
ADR-025 zugunsten unveränderter, überprüfbarer physikalischer Randbedingungen.

An der internen algebraischen Grenze müssen Aufrufer ein eindeutiges Problem
liefern. Singuläre LU-Blöcke werden schon beim Aufbau erkannt, auch bei RHS=0.
Ein Krylov-Verfahren garantiert jedoch keine allgemeine Rangprüfung beliebiger
gekoppelter Matrizen: Bei einem singulären Gesamtoperator und kompatibler RHS
kann ein kleines Residuum allein die Eindeutigkeit nicht belegen. Daher ersetzt
der Backend-Status die physikalische Nullraumprüfung ausdrücklich nicht.

`max_iterations` und `restart` sind positive ganze Zahlen; `0<rtol<1` und
`atol>=0` müssen endlich sein. KSP prüft die **unpräkonditionierte Norm des
skalierten** Systems mit diesen Toleranzen. Der Divergenzfaktor ist `1e8`.
Ein negativer Konvergenzgrund liefert `ValidationError` mit PETSc-Grundname und
-code, ausgeführter Iterationszahl, Limit und skaliertem Restfehler. Native
Aufbau-/Solve-Fehler werden weitergereicht, kein SciPy- oder Direktfallback.

Erfolgreiches KSP allein genügt nicht: Danach gelten weiterhin Gesamt- und
feldweise Residuen `A x-b` am ursprünglichen unskalierten, gegebenenfalls
komplexen System mit der unveränderten Grenze `1e-9`. Die Physiksolver prüfen
zusätzlich rekonstruierte Residuen, Reaktionen, Ladungs-/Kraft-/Wärmebilanzen
und Energien. Eine absichtlich lockere KSP-Toleranz wird im Test von diesem
zweiten Kriterium zurückgewiesen.

## Berichte, Wiederverwendung und Isolation

`SolveReport.method` ist `fgmres_fieldsplit_schur` oder `fgmres_lu`, ergänzt um
`_real_block` bei komplexer Einbettung. `reason` nennt den PETSc-Konvergenzgrund
des letzten nativen RHS-Solves; sämtliche vorherigen müssen erfolgreich sein.
`iterations` summiert die **äußeren** Iterationen dieses `solve`-Aufrufs über alle
RHS-Spalten und gegebenenfalls getrennte Real-/Imaginärsolves. Das Limit gilt
für jeden nativen RHS-Solve separat. Innere PREONLY-Anwendungen sind keine
zusätzlichen Krylov-Iterationen. Backendkennung/Provenienz enthalten die vier
konfigurierbaren Parameter. `close()`, Datenbesitz und Wiederverwendung gelten
wie bei #34b; Iterationszähler und Nullstart werden für jeden Aufruf erneuert.

PETSc lädt beim Aufbau der FieldSplit-Unterlöser intern Optionen. Ein privates
Präfix je Solver isoliert deshalb diese Optionen von anderen KSPs im Prozess.
Der Adapter verändert die globale Optionsdatenbank nicht und bietet darüber
keine benutzerseitige Konfiguration an. Alle unterstützten Einstellungen sind
explizite Konstruktorparameter.

## Verifikation und Reproduktion

Die Installation und beiden Docker-Kommandos aus Dokument 38 bleiben gültig;
nach Codeänderungen das Image neu bauen. Die bestehenden echten PETSc-CI-Jobs
führen das gesamte Verzeichnis `tests/petsc` aus, inklusive FieldSplit.
Fehlende native Laufzeit darf dort nicht zu übersprungenen Tests führen.

TDD: 14 neue Tests scheiterten zunächst am fehlenden Backend. Der Test zum
Ein-Schritt-Limit verwendet zwei verschiedene Schur-Eigenwerte und erzeugt
den echten PETSc-Iterationsabbruch (3.19.6: `DIVERGED_MAX_IT`, Code -3).
Ein weiterer Red-Test zeigte, dass fremde `fieldsplit_u_pc_type`-Optionen die
Unterlöser beeinflussen; das private Präfix behebt dies ohne globale Mutation.

Zusätzliche Nachweise: komplexe/operatorseitige Scalar-Treue, nicht zusammen-
hängende Feldindizes, Skalierung und mehrere RHS, Wiederverwendung nach
Nichtkonvergenz, leere Felder/Systeme, ungültige Optionen, Nullraumdiagnostik
und die unverändert strengen Originalresiduen. Die Physikparität aus #34b läuft
nun sowohl direkt als auch mit FieldSplit: Einzelfelder, feste/geladene
schwebende Elektroden und modale/Newmark-/harmonische M6-Teilprobleme mit
Energien und Bilanzen. Der vorhandene unabhängige Wärme-Konvergenzfall läuft
zusätzlich mit dem iterativen Adapter auf vier Netzen (2, 4, 6, 8): L2-Ordnung
zwischen 1,8 und 2,2, Energieordnung zwischen 0,85 und 1,15.

Die Referenzwerte und Toleranzen aus Dokument 38 bleiben unverändert. Neue
algebraische Parität für O(1)-Daten: `rtol=1e-11, atol=1e-12`, Feldresiduen
unter `1e-11`; integrierte Neumann-Lasten: `atol=1e-14 W` bzw. `C`, `rtol=0`.
Keine neue Diskretisierung oder Referenzwertanpassung.

Lokaler Abschlussstand: **429 Core-Tests bestanden, 43 optionale Tests ohne
PETSc sichtbar übersprungen**. Im installierten Docker-Paket jeweils **43
PETSc-Tests bestanden** für float64 und complex128 (Python 3.12.3,
PETSc/petsc4py 3.19.6, Linux arm64). Ruff-Lint und Formatprüfung, mypy für 52
Quelldateien sowie `git diff --check` bestanden. CI auf Linux x86_64 und
Reviewstatus werden am PR nachgewiesen; keine unabhängige Freigabe behauptet.

## Primärquellen

- [PETSc PCFIELDSPLIT](https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/)
- [Schur-Präkonditionierer und A11](https://petsc.org/release/manualpages/PC/PCFieldSplitSetSchurPre/)
- [FGMRES und Rechtspräkonditionierung](https://petsc.org/release/manualpages/KSP/KSPFGMRES/)
