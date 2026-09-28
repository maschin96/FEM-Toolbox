# Tet10-Gesamtabnahme (#35)

Stand: 28.09.2026, integrierter M7-Funktionsstand
`e56b429611bd653832fff70a9b5524e7142c133b` (einschließlich PR #198).
Die Kriterien aus #35 sind durch die unten zugeordneten Implementierungen und
Verifikationsnachweise abgedeckt. Dieser Bericht ergänzt deren Gesamtabgleich;
er ist keine Freigabe von v1.2.0. Review und Integration dieses Abnahmeberichts
sowie das übergreifende Release-Gate #38 bleiben gesonderte Schritte.
Ein unabhängiges Review dieses Gesamtabgleichs liegt noch nicht vor.

## Kriterien und prüfbare Nachweise

| Kriterium aus #35 | Integrierter Umfang | Ausführbare Nachweise |
|---|---|---|
| 35a: Topologie, I/O, Mengen, Regionen, Kompatibilität | [Topologievertrag](40_TET10_TOPOLOGY.md), PR #183; Gmsh 2.2/4.1, VTU, Kanten-/Flächenzuordnung, unveränderliche Arrays | [Topologietests](../tests/unit/test_tet10_topology.py), [I/O-Tests](../tests/unit/test_tet10_io.py); Tet4-Regression in der Gesamtsuite |
| 35b: P2-Basis, Quadratur, Ableitungen, Orientierung | [Geometrie](41_TET10_GEOMETRY.md), PR #184; mittige Kantenknoten, affine Geometrie | [Geometrietests](../tests/unit/test_tet10_geometry.py): Partition der Eins, Ableitungssummen, Polynomreproduktion, unabhängige Monomintegrale, alle 24 Eckpermutationen, Skalierung und Negativfälle |
| 35c: Feldoperatoren, Speicher, Kopplung, Randintegrale | [Operatorvertrag](42_TET10_OPERATORS.md), PR #185; konsistente mechanische Masse und Wärmekapazität; Triangle6-Elektroden mit allen sechs Knoten | [Operatortests](../tests/unit/test_tet10_operators.py): unabhängige baryzentrische Speicherintegrale, Energie, Starrkörpermoden, Piezo-Patch, Lastvorzeichen, Robin-Bilanz und Elektrodenkondensation |
| 35d: Solver, Preflight, Model, Backends, Export | PRs #188, #190–#197; Analyseübersicht unten | Integrations-/Model-/Exporttests in den unten verlinkten Teilberichten; [native Backend-Parität](../tests/petsc/test_physics.py) |
| Abnahme: glatte Lösungen auf mindestens vier Netzen | Stationäre Wärme, Elektrostatik, Mechanik und Kopplungen; P2-Interpolation unabhängig zusätzlich geprüft | [Konvergenznachweise](#konvergenz-und-bilanzen) |
| Abnahme: Biegung/Gradient bei Fehler- und Kostenbudgets | PR #198, [Messvertrag und Ergebnisse](54_ACCURACY_COST_COMPARISON.md): Tet4/Tet10 × SciPy/PETSc, drei Größen, drei Wiederholungen, 72 Läufe | [Unveränderte Rohdaten](../benchmarks/baselines/m7-37c-linux-arm64.json), [Norm-/Budgettests](../tests/regression/test_accuracy_cost.py), [nativer Worker-Vergleich](../tests/petsc/test_native_accuracy_cost.py) |

Die PR-Nummern beziehen sich auf das Repository
[maschin96/FEM-Toolbox](https://github.com/maschin96/FEM-Toolbox/pulls?q=is%3Apr+is%3Amerged+base%3Acodex%2Fmilestone-m7).
Die Integration erfolgte ausschließlich in `codex/milestone-m7`.

## Analyse- und Ergebnisumfang

Alle folgenden Pfade unterstützen geradkantige Tet10. SciPy ist der Standard;
die nativen Paritätstests verwenden PETSc direkt sowie FieldSplit, jeweils mit
reellem und komplexem PETSc. Das ist ein Nachweis für die dokumentierten kleinen
Referenzfälle, keine uneingeschränkte Konvergenzzusage des iterativen Solvers.

| Analyse | Einstieg und Model-Vertrag | Ergebnis-/Verifikationsvertrag |
|---|---|---|
| Stationäre Wärme | `solve_heat`, Model-Schema 1 | [Wärmeintegration](45_TET10_HEAT_INTEGRATION.md) |
| Transiente Wärme | `solve_transient_heat`, Schema 2 | [Kapazität, Zeitprofile und Export](45_TET10_HEAT_INTEGRATION.md) |
| Elektrostatik | `solve_electrostatics`, Schema 1 | [Elektrische Felder und Energie](46_TET10_ELECTROSTATICS.md) |
| Statische Elastizität | `solve_elasticity`, Schema 1 | [Mechanik, Biegung und Energie](47_TET10_ELASTICITY.md) |
| Stationäre Thermomechanik | `solve_thermomechanics`, Python-API; kein eigener deklarativer Analysemodus | [Einseitiges T → u](48_TET10_THERMOMECHANICS.md) |
| Stationäre Piezoelektrizität | `solve_piezoelectric`, Schema 1 | [Elektroden und Enthalpie](49_TET10_PIEZOELECTRIC.md) |
| Stationäre Thermopiezoelektrizität | `solve_thermopiezoelectric`, Schema 1 | [T → (u, φ)](50_TET10_THERMOPIEZOELECTRIC.md) |
| Piezo-Modalanalyse | `solve_piezoelectric_modes`, Schema 3 | [Masse und Moden](51_TET10_MODES.md) |
| Piezo-Newmark | `solve_transient_piezoelectric`, Schema 3 | [Zeitintegration und Ereignisse](52_TET10_NEWMARK.md) |
| Piezo-Harmonik | `solve_harmonic_piezoelectric`, Schema 3 | [Komplexe Antwort und DC](53_TET10_HARMONIC.md) |

Nodale Primärfelder verwenden die quadratische Basis auch für Messpunkte.
Abgeleitete Felder wie Spannung, Fluss, E und D sind rohe **Schwerpunktwerte**,
keine elementweise konstanten Felder, Knotenprojektionen oder vollständigen
Quadraturdatensätze. Auswerteorte und Exportnamen/Metadaten kennzeichnen dies.
Energie und Enthalpie werden integriert und dürfen nicht aus einem einzelnen
Schwerpunktwert mal Volumen rekonstruiert werden. Die jeweiligen Exporttests
prüfen diese Unterscheidung einschließlich komplexer Felder und Zeitserien.
Tet4 behält seine bisherigen Feldnamen und Verträge.

## Konvergenz und Bilanzen

Die stationären Manufactured Solutions verwenden jeweils vier Netze mit
N=1/2/3/4. Analytische Lasten und glatte Referenzfelder sind im Test angegeben.
Für die letzten beiden Verfeinerungsintervalle sind vorgegebene Fenster
L2: 2,7–3,3 und H1: 1,8–2,2 aktiv. Diese bestätigen die erwarteten Ordnungen
3 beziehungsweise 2; sie sind keine Zusage für Singularitäten oder beliebige
Netze/Materialkontraste.

- [Wärme](../tests/convergence/test_tet10_heat_operator.py): quartische Lösung,
  Produktions-/Testassemblerparität, unabhängige Überintegration, Leistungsbilanz.
- [Elektrostatik](../tests/convergence/test_tet10_electrostatics.py): Potential,
  H1-Fehler, unabhängig integrierte Energie und Ladungsbilanz.
- [Mechanik](../tests/convergence/test_tet10_elasticity.py): Verschiebung,
  Gradientenfehler, unabhängige Energie sowie Kraft- und Momentenbilanz.
- [Thermomechanik](../tests/convergence/test_tet10_thermomechanics.py): T und u,
  Spannungsfehler, thermische Eigenverzerrung, Energie und getrennte Bilanzen.
- [Piezo](../tests/convergence/test_tet10_piezoelectric.py) und
  [Thermopiezo](../tests/convergence/test_tet10_thermopiezoelectric.py): feldweise
  Normen, integrierte Enthalpie sowie Wärme-/Kraft-/Momenten-/Ladungsbilanzen.
- [Interpolation](../tests/convergence/test_tet10_interpolation.py): kubisches
  Polynom, unabhängige Quadratur; L2/H1-Ordnungen 3/2 bis 1e-9.

Erneute lokale Ausführung am 28.09.2026: **11 Konvergenztests bestanden**.
Beobachtete Ordnungen im letzten Intervall (N=3 → 4):

| Feld | L2 | H1 |
|---|---:|---:|
| Wärme / Elektrostatik | 2,985 | 1,978 |
| Mechanik | 2,989 | 1,980 |
| Thermomechanik u | 2,989 | 1,980 |
| Piezo u / φ | 3,007 / 2,985 | 1,994 / 1,978 |
| Thermopiezo u / φ / T | 3,008 / 2,985 / 2,985 | 1,994 / 1,978 / 1,978 |

Die Zahlen sind gerundete Ausgaben der oben verlinkten unveränderten Tests.

M6-Nachweise ergänzen diese räumlichen Feldnormen: Backward Euler auf vier
Zeitgittern, Newmark-Zeitordnung und Energie sowie vierstufige analytische
Stabnachweise für Eigenfrequenzen und harmonische Endverschiebungen. Die dort
beobachtete Ordnung ungefähr 4 betrifft diese **skalaren Zielgrößen**, nicht
eine H1-Feldordnung 4. Quellen:
[transiente Wärme](../tests/convergence/test_tet10_transient_heat.py),
[Newmark](../tests/unit/test_tet10_newmark.py),
[Modalanalyse](../tests/convergence/test_tet10_modal.py),
[Harmonik](../tests/convergence/test_tet10_harmonic_convergence.py).

## Grenzen und bewusste Nichtziele

- Nur konforme geradkantige Tet10/Triangle6 mit mittigen Kantenknoten.
  Gekrümmte Geometrie, falsche Kanten-/Randzuordnung und hybride Tet4/Tet10-
  Übergänge werden abgelehnt; dies ist kein isoparametrischer Geometriesolver.
- Kleine Verzerrungen und lineare Materialmodelle bleiben Voraussetzung.
  Tet10 behebt Inkompressibilitätslocking nicht automatisch; kein Nachweis für
  nahezu inkompressible Materialien oder Singularitäten.
- Wärmeindikatoren bleiben auf stationäre Tet4-Wärme beschränkt. Keine Tet10-
  Adaptivität, kein automatisches Remeshing, kein MPI-Skalierungsnachweis.
- Einseitige stationäre Thermokopplung; keine reziproke Thermoelastizität und
  keine neue transiente thermoelektromechanische Kopplung.
- Tet10-Wärme nutzt konsistente Kapazität oder ausdrücklich gewähltes HRZ;
  Tet4-Zeilensummenlumping wird nicht still auf Tet10 übertragen.
- Die bestehenden Dynamikgrenzen (2048 aktive Gesamt-DOFs,
  2.000.000 skalare History-Werte) und Harmonik-Resonanzdiagnosen bleiben aktiv.
  Die dichten dynamischen Teilschritte werden durch PETSc nicht unbegrenzt skalierbar.
- Verifikation des Codes, keine experimentelle Validierung oder pauschale
  Leistungsüberlegenheit. Im gemessenen Raster ist Tet4 bei groben Fehlerzielen
  günstiger; Tet10 erreicht kleinere Fehler. Nicht erreichte Ziele bleiben sichtbar.

## Reproduktion und Abnahmestatus

```sh
PYTHONPATH=src python -m pytest -q -s tests/convergence/test_tet10*
PYTHONPATH=src python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy src/fem_toolbox
git diff --check
```

Native PETSc- und isolierte Wheel-/sdist-Beispieltests sind in der
[CI](../.github/workflows/ci.yml) konfiguriert. Die native Fixture behandelt
fehlendes PETSc bei `FEM_REQUIRE_PETSC=1` als Fehler; diese CI-Nachweise sind keine
lokalen Skip-Ergebnisse. PETSc-Installation und Scalar-Auswahl stehen im
[Backendbericht](38_PETSC_BACKEND.md).

Am integrierten Vorgänger [PR #198](https://github.com/maschin96/FEM-Toolbox/pull/198)
bestanden alle 19 CI-Prüfungen am Head `188a43dc6ad7d4b17b6d7f7f95a0641d35606f4d`,
einschließlich nativer PETSc-real/complex-Ausführung und Wheel-/sdist-Beispielen
auf Python 3.11–3.13. Der Mergecommit `e56b429` ist die Basis dieses Berichts;
diese PR-Evidenz wird nicht als separater CI-Lauf auf dem Mergecommit ausgegeben.

Lokale Gesamtprüfung dieses Abnahmeschnitts: **603 bestanden, 77 optionale
Tests übersprungen**. Ruff, Formatprüfung, Mypy und `git diff --check` bestanden;
219 lokale Markdown-Linkziele der geänderten Dokumente wurden auf Existenz
geprüft. Native Laufzeiten werden durch die CI nachgewiesen.

Dieser Abnahmeschnitt ändert ausschließlich Dokumentation. Daher ist kein
künstlicher TDD-Red-Test erforderlich. Bestehende fachliche Tests werden erneut
ausgeführt; Toleranzen und Referenzdaten bleiben unverändert. Die frühere TDD-
Historie ist in den jeweiligen Implementierungsberichten dokumentiert.

#35 kann nach geprüftem Merge dieses Berichts in den M7-Branch geschlossen
werden. Die übergreifende Release-Supportmatrix, finale Paket-/Releaseversion,
main-Integration, Zielcommit-CI, Tag und Veröffentlichung gehören weiterhin zu
#38; der stabile Hauptbranch bleibt bis dahin v1.1.0.
