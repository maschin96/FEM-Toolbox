# M7-Release-Kandidat und Supportmatrix (#38)

Die Supportmatrix gilt auch für den unveränderten Funktionsumfang von 1.2.0.
Aktueller Status: [finaler Release-Nachweis](58_M7_RELEASE_VERIFICATION.md).
Die folgenden Kandidatenversionen und Testergebnisse sind historisch.

Stand: 28.09.2026. Kandidat **1.2.0rc1**, Zielversion **v1.2.0**.
Stabile veröffentlichte Version bleibt v1.1.0. Der Kandidat enthält die in M7
integrierten PETSc-, Tet10-, Indikator- und Performanceänderungen. Er ist kein
finales Release und setzt keinen Release-Tag.

## Verbindliche Supportmatrix des Kandidaten

„Ja“ bedeutet Unterstützung im unten beschriebenen Umfang und Nachweis an den
Referenzfällen. PETSc wird ausschließlich über direkte Python-Solveraufrufe
gewählt; `Model`/TOML verwendet SciPy und bietet keinen Backend-Schlüssel.
Daher bedeutet eine PETSc-Spalte keine deklarative PETSc-Unterstützung.

| Analyse | Tet4 SciPy | Tet10 SciPy | Tet4 PETSc direkt / FieldSplit | Tet10 PETSc direkt / FieldSplit | Model |
|---|---|---|---|---|---|
| Stationäre Wärme | Ja | Ja | Ja / Ja | Ja / Ja | Schema 1 |
| Transiente Wärme (Backward Euler) | Ja | Ja | Ja / Ja | Ja / Ja | Schema 2 |
| Statische Elastizität | Ja | Ja | Ja / Ja | Ja / Ja | Schema 1 |
| Elektrostatik | Ja | Ja | Ja / Ja | Ja / Ja | Schema 1 |
| Stationäre Thermomechanik T → u | Ja | Ja | Ja / Ja | Ja / Ja | Nur Python-API |
| Stationäre Piezoelektrizität u–φ | Ja | Ja | Ja / Ja | Ja / Ja | Schema 1 |
| Stationäre Thermopiezoelektrizität T → (u, φ) | Ja | Ja | Ja / Ja | Ja / Ja | Schema 1 |
| Piezo-Moden | Ja | Ja | Lineare Elimination / ebenso | Lineare Elimination / ebenso | Schema 3 |
| Piezo-Newmark | Ja | Ja | Ja / Ja | Ja / Ja | Schema 3 |
| Piezo-Harmonik einschließlich komplexer Lasten | Ja | Ja | Ja / Ja | Ja / Ja | Schema 3 |

Tet10 bedeutet ausschließlich konforme geradkantige Tet10/Triangle6 mit mittigen
Kantenknoten. Die Modalanalyse nutzt unabhängig vom linearen Backend SciPy `eigh`;
PETSc ersetzt nur lineare Teilschritte. Die dichte Dynamikreferenz bleibt auf
2048 aktive Gesamt-DOFs und bei Zeitverläufen 2.000.000 skalare History-Werte
begrenzt. FieldSplit-Konvergenz ist problemabhängig und wird diagnostiziert;
kein stiller Wechsel zum Direktlöser.

Die nativen PETSc-Spalten werden sowohl im float64- als auch complex128-Build
geprüft. Komplexe Systeme werden bei reellem PETSc gemäß
[Scalar-Vertrag](38_PETSC_BACKEND.md) real eingebettet beziehungsweise in zwei
RHS zerlegt; sie werden nicht still reell abgeschnitten. PETSc ist sequentiell
(`COMM_SELF`), kein MPI-Support. Thermokopplungen bleiben stationär und einseitig.

### Nachweise je Matrixzeile

Die [Tet10-Abnahme](55_TET10_ACCEPTANCE.md) ordnet jede Analyse den Patch-,
Konvergenz-, Bilanz-, Model- und Exporttests zu. Native Parität beider Backend-
Verfahren steht in [test_physics.py](../tests/petsc/test_physics.py), die echte
Laufzeit wird durch [conftest.py](../tests/petsc/conftest.py) erzwungen.
Tet4- und Tet10-Einzelfelder, Thermokopplungen und Dynamik sowie feste und
schwebende Elektroden werden dort ausgeführt. Die expliziten nativen
Thermomechanik-/Thermopiezo-Paritätstests sind nun nach Elementordnung
parametrisiert; der Tet4-Thermomechanikfall nutzt ein exakt darstellbares affines
Temperaturfeld, der Tet10-Fall weiterhin sein quadratisches Feld. Die
Tet4-Thermopiezo-Prüfung vergleicht denselben diskreten Fall mit SciPy und
kontrolliert zusätzlich Wärme-, Kraft-, Momenten- und Ladungsbilanz.

Die [Performance-Abnahme](56_PERFORMANCE_ACCEPTANCE.md) enthält reproduzierbare
Laufzeit-/Speichernachweise; eine unterstützte Kombination ist keine
Geschwindigkeitszusage. [Wärmeindikatoren](43_HEAT_INDICATOR.md) sind ausschließlich
für stationäre Tet4-Wärme verfügbar, unabhängig vom zur Lösung gewählten Backend.
Sie sind Forschungsindikatoren, keine garantierte Fehlerschranke.

## Nicht unterstützte Kombinationen und Diagnose

| Kombination | Verhalten / Nachweis |
|---|---|
| Gekrümmte Tet10, nichtmittige Kanten, hybride Tet4/Tet10-Netze | Geometrie-/Topologiefehler vor der Lösung; [Topologietests](../tests/unit/test_tet10_topology.py), [Geometrietests](../tests/unit/test_tet10_geometry.py) |
| Tet10 mit Tet4-Zeilensummenlumping | Explizite Validierung; konsistente Kapazität oder ausdrücklich HRZ; [Wärmetests](../tests/unit/test_tet10_heat_integration.py) |
| Tet10-Wärmeindikator | `ValidationError`; [Indikatortests](../tests/unit/test_heat_indicator.py) |
| PETSc gewählt, native Laufzeit fehlt/inkompatibel | `ValidationError` beim Backendaufbau, kein SciPy-Fallback; [Optionalitätsprüfung](../tests/unit/test_petsc_optional.py) |
| PETSc-Backend in TOML | Unbekannter Schlüssel wird als `ModelSchemaError` abgelehnt; keine stille Backendwahl |
| Nicht unterstützter Model-Analysemodus / Schema | `ModelSchemaError`; bestehende Model-Schematests |
| Dynamik oberhalb der dokumentierten Grenzen, singuläre/resonante Systeme | Explizite Größen-/Lösbarkeitsdiagnostik; keine Freigabe durch kleinen linearen Residualwert allein |

Keine neue transiente Thermokopplung, reziproke Thermoelastizität, automatische
Adaptivität, nahezu inkompressible Formulierung oder experimentelle Validierung.

## Installation und Artefakte

CPython 3.11–3.13, MIT, NumPy/SciPy/meshio bleiben unverändert. PETSc/petsc4py
sind keine Core-Abhängigkeit. Die Core-CI installiert Wheel und sdist in frischen
Umgebungen ohne Checkout und ohne PETSc und führt alle alten und neuen Beispiele
aus. Zusätzlich werden Kandidatenversion, MIT-Metadaten und mitgelieferte
Lizenzdatei sowie die explizite Diagnose bei angefordertem fehlendem PETSc geprüft.

Die [optionale PETSc-Installation](38_PETSC_BACKEND.md) nutzt das gepinnte
Ubuntu-24.04-Containerrezept, Python 3.12 und beide nativen Scalar-Builds.
Core Linux/Python 3.11–3.13 und lokale macOS-arm64-Nachweise sind getrennt von
PETSc Linux. Keine zusätzliche Windows-/macOS-PETSc- oder MPI-Zusage.

```sh
python -m pip install ./fem_toolbox-1.2.0rc1-py3-none-any.whl
# Alternativ das sdist mit isolierten Build-Abhängigkeiten:
python -m pip install ./fem_toolbox-1.2.0rc1.tar.gz
```

Diese Namen bezeichnen Kandidatenartefakte aus der CI, noch keine veröffentlichten
Release-Downloads. Kein PyPI-Upload wird vorgenommen.

## Verifikation dieses Kandidatenschnitts

Am 28.09.2026 lokal ausgeführt:

- Core-Suite: **605 bestanden, 83 optionale PETSc-Tests übersprungen**.
- Neu gebauter Kandidatencontainer, Python 3.12.3/PETSc 3.19.6:
  **83 bestanden je real- und complex-Build**, keine übersprungenen Tests.
- Frische getrennte Wheel-/sdist-Installationen ohne PETSc: Kandidatenversion,
  Import aus `site-packages`, MIT-Metadaten/Lizenzdatei und explizite Fehlermeldung
  bei Backendwahl geprüft. Je **17 Beispiele** aus v1.0/v1.1 und M7 mit `-I`
  aus separaten temporären Arbeitsverzeichnissen erfolgreich ausgeführt.
  Lokale Artefaktumgebung: macOS arm64, Python 3.13.5, NumPy 2.5.3,
  SciPy 1.18.1, meshio 5.3.5. Beispielskripte stammen aus dem Checkout,
  Toolbox-Import ausschließlich aus dem installierten Archiv.
- Ruff, Formatprüfung, Mypy, `git diff --check` und lokale Markdown-Linkziele
  bestanden. Die CI wiederholt die Linux-/Python- und Artefaktmatrix am PR-Stand.

TDD zur Versionsgrenze: Der aktualisierte Pakettest scheiterte zunächst mit
`1.1.0 != 1.2.0rc1`; nach Metadatenänderung und isolierter Neuinstallation bestand
er. Die ergänzten Tests zur Backendwahl sichern bestehende Ablehnungen ab;
Produktionssolver und numerische Toleranzen wurden nicht verändert. Die zuerst
fehlenden benannten Flächensätze der erweiterten Tet4-Testfixture wurden vor den
erfolgreichen nativen Gesamtläufen ergänzt.

```sh
PYTHONPATH=src python -m pytest -q
uv build --out-dir /tmp/fem38-dist
# Native Containerbefehle: siehe Doc/38_PETSC_BACKEND.md
```

Ein unabhängiger Reviewer dieses Kandidatenschnitts ist nicht nachgewiesen;
daraus wird keine unabhängige Freigabe abgeleitet.

## Verbleibende Freigabeschritte

1. Kandidaten-PR gegen M7 prüfen und integrieren; CI und unabhängigen Reviewstatus
   am aktuellen Stand offenlegen. Fehlender unabhängiger Reviewer ist keine Freigabe.
2. Integrierte Funktionsissues explizit abschließen. Finale Version 1.2.0,
   Changelog und Release-Datum konsistent setzen und erneut als Artefakte prüfen.
3. Aktuelles `main` in M7 übernehmen und den resultierenden Stand prüfen.
4. Freigabefähigen M7-Release-PR nach `main` integrieren; CI am tatsächlichen
   `main`-Zielcommit abwarten.
5. Erst dann annotiertes `v1.2.0`-Tag, GitHub Release und geprüfte Artefakte
   veröffentlichen; #38 und den Milestone abschließen.

Dieser Kandidat erfüllt nicht bereits die letzten Veröffentlichungsschritte.
Details bleiben im [Release-Prozess](09_RELEASE_PROCESS.md) verbindlich.
