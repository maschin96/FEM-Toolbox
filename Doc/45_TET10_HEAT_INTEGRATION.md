# Tet10-Wärmepfad: Model, Solver, Backend und Export (#35d)

Status: implementierter und lokal geprüfter thermischer Teil von #35d im
M7-Issue-Branch. #35 bleibt für die übrigen Physiken, Kopplungen und Vergleiche
offen. Der spätere [elektrische Integrationsschnitt](46_TET10_ELECTROSTATICS.md)
ergänzt inzwischen Elektrostatik; die folgende Übersicht beschreibt den
thermischen Teil zum Zeitpunkt seiner Umsetzung. Keine experimentelle Validierung oder M7-Releasefreigabe.

## Unterstützter Umfang

| Pfad | Tet4 | Geradkantige Tet10 |
|---|---|---|
| Stationäre Wärme, Python und Model-Schema 1 | weiterhin unterstützt | unterstützt |
| Backward Euler, Python und Model-Schema 2 | weiterhin unterstützt | unterstützt |
| Dirichlet, Quelle, Fluss, Konvektion, Materialregionen | Triangle3 | Triangle6 mit allen sechs Knoten |
| Lineares Backend an direkten Wärme-APIs | SciPy / optional PETSc | gleiche Backend-Grenze |
| Masse/Kapazität | konsistent oder boolesches `lumped` | konsistent oder explizites `lumping="hrz"` |
| Ergebnis-/Serienexport | bestehendes Schema 1 | Schema 2 mit Schwerpunktkennzeichnung |
| Mechanik, Elektrostatik, Piezo und Kopplungen im Model | bisheriger Umfang | noch ausdrücklich abgelehnt |

Die Elementfamilie folgt aus `Mesh.cells` beziehungsweise dem Gmsh-Zelltyp;
keine implizite Hochstufung und kein neuer TOML-Elementschalter. Gemischte
Elementordnungen, gekrümmte Kanten und nichtmittige Kantenknoten bleiben
abgelehnt. Der gemeinsame Preflight verwendet die aktive Zell-/Flächentopologie
und erfasst auch die Kantenknoten bei Zusammenhang, Referenzen und Randvorgaben.
Fehlende Temperaturreferenzen werden stationär erkannt; isolierte transiente
Modelle benötigen bei positiver Kapazität keine künstliche Referenz.

Die internen thermischen Dispatchfunktionen verwenden die vorhandenen
[quadratischen Operatoren](42_TET10_OPERATORS.md); es gibt keinen zweiten
Zeitintegrator, Backend-Solver oder neuen Materialtyp. Tet4-Operatorformeln
bleiben die bestehenden Referenzfunktionen. Quellen und Randlasten behalten
alle Vorzeichen und SI-Einheiten. Tet10-Volumenquadratur hat Grad 5,
Triangle6-Flächenquadratur Grad 4 (jeweils Gesamtintegrand, nicht Last allein).

## Punktabhängige Ergebnisse und versionierter Export

`HeatResult.gradient/heat_flux` und `HeatFrame.gradient/heat_flux` behalten die
Form `(ncells, 3)` und die Einheiten K/m bzw. W/m². Der zusätzliche optionale
Eintrag `field_points` macht den Unterschied explizit:

- Tet4: `None`, da der Gradient in jeder Zelle konstant ist.
- Tet10: unveränderliches `(ncells, 3)`-Array der physischen Zellschwerpunkte in m.
  Die Feldwerte sind rohe Auswertungen an genau diesen Punkten. Sie sind **nicht**
  elementweise konstant, an Knoten gemittelt oder projiziert.

Die Schwerpunktwerte beschreiben nicht den gesamten linearen Gradientenverlauf.
Für andere Punkte sind die nodalen Temperaturen mit
`Tet10Geometry.shape_gradients(point).T @ T_cell` auszuwerten; für
Energienormfehler ist die punktabhängige Integration erforderlich. Aus einem
Schwerpunktwert allein darf kein konstanter Wärmefluss angenommen werden.

Stationäres Ergebnisbündel und HeatSeries schreiben für Tet10 Schema 2:

- Zelltyp `tetra10`, Randtyp `triangle6`, tatsächliche Zell-IDs.
- Feldnamen mit `_at_centroid`, etwa `heat_flux_W_m2_at_centroid`.
- Schwerpunktkoordinaten in `field_evaluation_points_m` und Kennzeichnung
  `location` beziehungsweise `field_location = "cell_centroid"`.
- Unveränderte nodale Temperaturen, Reaktionen, Residuen und SI-Einheiten.
- Randzell-Padding bleibt null; Zeitserien besitzen weiterhin Gültigkeitsmasken.

Export lehnt fehlende oder zum Netz unpassende Tet10-Auswertepunkte ab.
Tet4 bleibt bei Schema 1 und bisherigen Feldnamen. Der eigenständige
[Manifestleser](../scripts/read_result_bundle.py) unterstützt beide Versionen
und prüft weiterhin die VTU-Prüfsumme. Ältere Leser dürfen Schema 2 ablehnen;
sie sollen Tet10-Schwerpunktwerte nicht als Tet4-Konstantfelder interpretieren.
Model-Eingabeschemas bleiben 1 beziehungsweise 2: die Erweiterung betrifft
unterstützte Netze und ausdrücklich ausgewählte Kapazität, keine neue Physik.

## Transiente Besonderheiten

Quadratische Formfunktionen können innerhalb eines Elements negativ sein.
Sonden werden deshalb mit den affinen baryzentrischen Eckkoordinaten lokalisiert
und anschließend mit allen zehn Formfunktionen interpoliert. Sonden außerhalb
des Netzes werden weiterhin abgelehnt. Zeitgitter, Ereignisse, linke Lastgrenzen,
rechte Dirichletgrenzen und Speicherschranken folgen dem [M6-Vertrag](33_TRANSIENT_HEAT.md).
Die History-Abschätzung umfasst die zusätzlichen Schwerpunktkoordinaten.

Standard ist konsistente Kapazität. Für Tet10 lehnen Python und Model
`lumped=True` ab: quadratische Zeilensummen an Ecken wären negativ. Stattdessen
ist `lumping="hrz"` ausdrücklich wählbar, auch als TOML-`solver.lumping`.
Ungültige Werte werden nicht still auf den Standard zurückgesetzt. Tet4
verwendet weiter ausschließlich den bisherigen booleschen Vertrag.

HRZ erhält die Gesamtkapazität, nicht die nodalen Momente des konsistenten
Quellvektors. Bei räumlich gleichförmiger Quelle kann eine HRZ-Lösung daher auf
einem groben Tet10-Netz räumlich variieren. Der Nachweis prüft hier Gesamtenergie,
Bilanz und den ungestörten konstanten Zustand. Der exakt gleichförmige
Aufheizfall und die Zeitordnung verwenden konsistente Kapazität. Es wird keine
Quellenumlagerung erfunden, um beide Diskretisierungen gleich erscheinen zu lassen.

## Verifikation

TDD-Red: fünf neue Integrationstests scheiterten zuerst an fehlenden Tet10-DOFs,
Model-Sperre und fehlendem transientem Dispatch. Separat scheiterten der
Schema-2-Manifestleser und das ausführbare Beispiel vor ihrer Ergänzung.
Negative Lumping-Fälle deckten stilles Defaulting von `""`/`False` auf; jetzt
werden beide gezielt abgelehnt.

Lokaler Gesamtstand: 517 Tests bestanden, 45 optionale Tests übersprungen;
Ruff, Formatprüfung, Mypy und `git diff --check` bestanden. Testaufruf im
Checkout: `PYTHONPATH=src .venv/bin/python -m pytest -q`. Native PETSc-Prüfungen
benötigen die unten genannte CI-Umgebung.

Die Nachweise umfassen:

- Quadratischer Patch mit 300-K-Offset, Q=-6 W/m³, exakten nodalen Temperaturen,
  analytischen Schwerpunktgradienten/-flüssen und globaler Leistungsbilanz.
- Triangle6-Fluss und Robin-Rand über benannte Model-Auswahlsätze; portable
  Dokumentreproduktion, Anfangszustände und zeitabhängige Kantenknoten-Randwerte.
- Materialgrenze mit k=1 beziehungsweise k=8 W/(m K), konstantem
  Normalfluss und achtfachem Gradientensprung. Dieser Test verwendet
  Temperaturabweichungen der Größenordnung 1, um den Kontrast unabhängig von
  Auslöschung durch große Temperaturreferenzen zu prüfen; der 300-K-Patch ist separat.
- Isolierte Aufheizung, Gesamtenergie, nichtuniformer Anfangszustand mit exakter
  quadratischer Sondeninterpolation, HRZ-Erhaltung und explizite Invalidfälle.
- Vollständiger Gmsh-/TOML-Pfad mit unabhängiger VTU-Leserprüfung, unveränderlichen
  Auswerteorten, Schema-2-Manifeste, falsche Exportpunkte und History-Budget.
- Bestehender unabhängiger Testassembler **und** Produktionssolver auf vier
  Netzen n=1/2/3/4 für T=x⁴+y⁴/2+z⁴/4; Produktions-/Referenzparität mit `2e-12`
  absolut und unabhängiger überintegrierter Fehlernorm. Die letzten beiden
  Intervalle müssen L2-Ordnung 2,7–3,3 und H1-Ordnung 1,8–2,2 erfüllen.
- Transiente Manufactured Solution T=300+t², Q=24t bei rho*c=12 auf vier
  Zeitgittern mit 2/4/8/16 Schritten: Endfehler genau dt bis `5e-11 K`,
  erwartete Ordnung eins im Fenster 0,99–1,01.
- Native SciPy/PETSc-Parität für stationäre und transiente Tet10-Wärme wird in
  den vorhandenen reellen/komplexen PETSc-CI-Jobs geprüft, jeweils direkt und
  mit FieldSplit. Lokal ohne PETSc sind diese Tests ausdrücklich übersprungen.

Toleranzen: analytische Knotentemperaturen bei 300 K bis `5e-12 K`,
Schwerpunktflüsse bis `2e-11 W/m²`, stationäre Bilanzen bis `2e-11 W`,
transiente Bilanzen bis `2e-10 W` (Materialpreset mit großer Speicherkapazität:
`1e-7 W`). Sie decken Double-Rundung und Assemblierung der kleinen
Referenzfälle ab, keine allgemeine Diskretisierungsgenauigkeit.

Ausführbare [stationäre und transiente Beispiele](../examples/tet10_heat/README.md).
Die Artefakt-CI führt sie zusätzlich aus frischen Wheel-/sdist-Installationen
auf Python 3.11–3.13 aus. Das ist ein zusätzliches Release-Artefaktgate,
kein lokal bereits behaupteter Installationsnachweis.

## Restumfang

Die thermische Integration ist ein eigenständig prüfbarer PR-Schnitt innerhalb
von #35d. Noch offen sind die produktive Tet10-Anbindung von Mechanik,
Piezo-/Thermokopplungen und Piezo-Dynamik (Elektrostatik ist inzwischen in
[Dokument 46](46_TET10_ELECTROSTATICS.md) integriert) sowie deren
ortsabhängige Ergebnisverträge. Biegungs-/Fehlerbudgetvergleiche und #37c
bleiben weitere Gates; #35 wird durch diesen Teilbeitrag nicht geschlossen.
