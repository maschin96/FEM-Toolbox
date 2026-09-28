# Statische Tet10-Elastizität durch Model, Solver und Export (#35d)

**Aktueller Gesamtstand:** [Tet10-Abnahme](55_TET10_ACCEPTANCE.md).
Die folgenden Status- und Restumfangsangaben dokumentieren den jeweiligen
Implementierungsschnitt; später integrierte Analysen und Vergleiche sind im
Gesamtbericht zusammengeführt.

Historischer Umfang dieses Integrationsschnitts: statische lineare Elastizität nach Wärme und
Elektrostatik. #35 bleibt für thermische Kopplungen, Piezo, Dynamik und
verbleibende Kosten-/Vergleichsgates offen. Keine M7-Releasefreigabe oder
experimentelle Validierungsbehauptung. Die nachfolgende
[Tet10-Thermomechanik](48_TET10_THERMOMECHANICS.md) hebt die unten
beschriebene damalige Sperre thermischer Eigendehnungen inzwischen auf.

## Unterstützter Vertrag

`solve_elasticity`, `mechanical_dofs` und Model-Schema 1 unterstützen geradkantige
Tet10 mit mittigen Kantenknoten. Jeder aktive Knoten trägt u=(x,y,z), ein Element
30 DOFs. Inaktive Netzknoten erhalten keine DOFs. Benannte Triangle6-Flächen
beziehen alle sechs Knoten in Lager und Lasten ein. Materialregionen und
anisotrope Steifigkeit verwenden unverändert die bestehenden SI-/Voigt-Verträge.

Volumenkräfte, globale Flächentraktionen und Punktlasten werden konsistent
assembliert; Lasten auf festgehaltenen DOFs bleiben in der Reaktionsbilanz.
Tet10-Volumenquadratur hat Grad 5, Triangle6-Flächenquadratur Grad 4, jeweils
bezogen auf den gesamten Integranden. Keine Druck- oder Folgelasten. Der
Model-Preflight und die direkte Lagerprüfung erfassen sechs Starrkörpermoden
pro zusammenhängender Komponente einschließlich Kantenknoten. Die Backend-Grenze
bleibt unverändert: SciPy als Referenz, optional PETSc über die direkte Python-API.

Gemischte Elementfamilien, gekrümmte Kanten und thermische Tet10-Eigendehnungen
werden ausdrücklich abgelehnt. Auch `expansion`/`temperature` am direkten Solver
führen bei Tet10 zu einem konkreten Fehler, statt die Tet4-Zellmittelung ungeprüft
zu übernehmen. Die vorhandene Tet4-Thermomechanik bleibt erhalten. Tet10 behebt
Inkompressibilitätslocking nicht automatisch; die bestehenden Materialgrenzen
gelten weiterhin. Piezo-/Thermokopplungen erhalten durch diesen Schnitt keine Freigabe.

## Auswerteorte, Energie und Export

`ElasticResult.strain/stress` bleiben `(ncells,6)`-Arrays in Eingabezellreihenfolge.
`field_points` ist bei Tet4 `None`, bei Tet10 ein unveränderliches `(ncells,3)`-
Array physischer Schwerpunkte in Metern. Tet10-Spannungen und -Dehnungen sind
rohe Auswertungen dort, keine konstanten Zellfelder und keine Knotenprojektion.
Voigt: `[xx, yy, zz, yz, xz, xy]`, Dehnungen mit technischen Schubkomponenten.

`element_energy` bleibt integrierte positive Verzerrungsenergie pro Zelle in J.
Für Tet10 wird `1/2 integral epsilon(u)^T C epsilon(u)` mit punktabhängigem B
integriert. Ein Schwerpunktwert allein liefert keine exakte Energie.
Die bestehenden Tet4-Formeln und thermischen Energiebeiträge bleiben erhalten.

Schema 2 schreibt `strain_1_at_centroid` und `stress_Pa_at_centroid` als
`raw`/`cell_centroid`, mit expliziten Auswerteorten und korrekten Voigt-Komponenten.
`mechanical_energy_J` bleibt ohne Schwerpunktsuffix und ist als
`cell_integral`/`cell_volume` gekennzeichnet. Tet4 behält Schema 1 und seine
bisherigen Feldnamen. Fehlende oder netzfremde Tet10-Auswerteorte werden abgelehnt.

## Verifikation

Lokal: `PYTHONPATH=src .venv/bin/python -m pytest -q` mit 535 bestandenen und
49 optional übersprungenen Tests. Ruff, Formatprüfung, Mypy und
`git diff --check` bestanden; native PETSc-/Artefaktprüfungen erfolgen in CI.

TDD: drei initiale Integrationstests scheiterten vor Umsetzung an der
Model-Sperre, fehlenden Tet10-DOFs und der Tet4-Solverbeschränkung. Anschließend
bestanden quadratischer Patch, Biegefall und Lagerprüfung. Ein weiterer Test
scheiterte vor Ergänzung des ausführbaren Gmsh-/TOML-Beispiels an dessen fehlender Datei.

Analytische Fälle und Negativnachweise:

- Quadratischer Patch u=(0,01 x²,0,0), E=100 Pa, ν=0, b=(-2,0,0) N/m³,
  rechte Traktion (2,0,0) N/m²: Energie 0,02/3 J, Schwerpunktspannungen sowie
  Kraft-/Momentenbilanz. Model, benannte Ränder, Dokumentreproduktion und
  unabhängiges VTU-Einlesen sind eingeschlossen.
- Reine Biegung auf dem Einheitswürfel: u=(-κx(y-1/2), κx²/2,0), κ=0,01 1/m,
  E=100 Pa und ν=0. Links u=0, rechts t=(-Eκ(y-1/2),0,0), sonst lastfrei.
  εxx=-κ(y-1/2), sämtliche Schubdehnungen null, W=Eκ²/24 J und äußeres
  Moment Mz=Eκ/12 N m. Bereits sechs Tet10 repräsentieren die Lösung exakt.
  Dies ist ein analytischer Biegepatch, noch kein Vergleich bei gleichem Kostenbudget.
- Anisotroper affiner Patch mit sechs unabhängig erwarteten Voigt-Komponenten,
  auch nach Orientierungsumkehr und mit inaktivem Netzknoten.
- Vorgegebene Starrkörpertranslation/-rotation: keine Verzerrungsenergie;
  fehlende Lager werden im Model und im direkten Solver abgelehnt.
- Punktlast: Kraft-/Momentenbilanz und Arbeitsidentität W=1/2 f·u.
- Zwei Model-Materialregionen mit E=100/800 Pa, ν=0 und σxx=1 Pa:
  kontinuierliche Verschiebung/Normalspannung, achtfacher Dehnungssprung,
  Energie 0,0028125 J.
- Gekrümmte Geometrie, unzulässige thermische Tet10-Eingaben, fehlende
  Exportpunkte und unveränderliche Schwerpunktkoordinaten sind geprüft.

Die kleinen synthetischen Fälle verwenden absolute Toleranzen: Verschiebung
2e-14 m, Dehnung 2e-13, Spannung 2e-11 Pa (anisotroper Patch/Starrkörpermotion
1e-10 Pa), Kraft/Moment 2e-12 N bzw. N m, Energie 2e-15 J (anisotroper Patch
2e-14 J). Der Biegefall hat 2e-16 J Energietoleranz; Starrkörperenergie <1e-24 J.
Sie decken Double-Rundung und Assemblierung ab, keine allgemeine Fehlergrenze.

## Vier-Netze-Nachweis

Manufactured Solution: u=(0,01(x⁴+y⁴/2+z⁴/4),0,0), E=100 Pa, ν=0,25
(λ=μ=40 Pa). Daraus folgt b=(-0,01(1440x²+240y²+120z²),0,0) N/m³.
Verschiebungen werden am gesamten Rand exakt vorgegeben. Unabhängig formulierte
baryzentrische P2-Basis/Ableitungen und eine Duffy-Produktregel mit sechs
Gauss-Punkten pro Richtung integrieren L2- und vollständigen H1-Seminormfehler.
Jedes Zellenergieintegral wird zusätzlich über die unabhängige Tensorform
`lambda/2 trace(grad u)^2 + mu sym(grad u):sym(grad u)` geprüft
(`rtol=2e-12`, `atol=1e-16 J`).

| Würfelteilung n | L2-Fehler | H1-Seminormfehler |
|---|---:|---:|
| 1 | 8,0234e-4 | 5,20219e-3 |
| 2 | 1,1196e-4 | 1,43559e-3 |
| 3 | 3,35718e-5 | 6,48902e-4 |
| 4 | 1,42095e-5 | 3,67120e-4 |

Die letzten beiden Intervalle liefern L2-Ordnungen 2,971/2,989 und
H1-Ordnungen 1,958/1,980. Akzeptanzfenster: 2,7–3,3 bzw. 1,8–2,2.

## Backend und ausführbares Beispiel

Neue native PETSc-Paritätstests werden durch die bestehenden real/complex-Jobs
jeweils direkt und mit FieldSplit ausgeführt. Toleranzen: 1e-11 m, 1e-8 Pa,
Gesamtenergie 1e-11 J, Kraft/Moment 1e-9 N bzw. N m; ohne PETSc sind diese
Tests lokal ausdrücklich übersprungen. Das
[Gmsh-/TOML-Beispiel](../examples/tet10_elasticity/README.md) läuft zusätzlich
in der Artefakt-CI aus frischen Wheel-/sdist-Installationen auf Python 3.11–3.13.
