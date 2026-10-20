# Verifizierter generischer Unsicherheitsworkflow (#64)

Stand: M9-Kandidat, 10.10.2026. Keine experimentelle Validierung. Die
[genehmigte Neuplanung](66_M9_DATA_SEARCH_AND_RELEASE_OPTIONS.md) trennt diese
Softwareverifikation von den offenen Messdaten-/Anwendungsissues #63/#65.

## Annahmen und mathematischer Vertrag

Der Anwender liefert Wahrscheinlichkeitsgesetze und begründet sie ausdrücklich.
`Normal(mean, standard_deviation)` fordert endlichen Mittelwert und positive
Standardabweichung; `Uniform(lower, upper)` endliche wachsende Grenzen. Alle
Parameter liegen in expliziten SI-Einheiten. Ein bekanntes Intervall allein
begründet keine Gleichverteilung. Es gibt keine geschätzten Verteilungen,
Abschneidung negativer Werte, Optimierung oder inverse Identifikation.

Die Eingabematrix X hat n gemeinsame Stichprobenzeilen und p Parameter-Spalten.
Default: unabhängige Eingaben. Eine explizite Korrelation R[p,p] beschreibt nur
bei ausschließlich normalen Eingaben eine multivariate Normalverteilung.
R muss symmetrisch, positiv semidefinit und auf der Diagonale eins sein
(absolute Rundungstoleranz 1e-12). Singuläre R sind erlaubt; Eigenwerte zwischen
-1e-12 und null werden auf null gerundet. Abhängige uniforme/mischverteilte
Eingaben werden abgelehnt; es wird keine unbenannte Copula unterstellt.

Ein lokaler PCG64-Generator erzeugt standardnormale gemeinsame Zeilen Z.
Für normale Abhängigkeit werden sie mit dem Eigenwertfaktor von R transformiert;
Normalwerte sind mu+sigma*Z. Unabhängige Uniformwerte verwenden die Normal-CDF
Phi(Z) und a+(b-a)*Phi(Z). Seed, n, Parameterreihenfolge, R, Gesetze,
Begründungen sowie NumPy-/SciPy-Versionen stehen in `sampling.json`.
Identische Spezifikation in derselben Laufzeit reproduziert die Werte; keine
Byteidentität über verschiedene BLAS-/Plattform-/Bibliotheksstände versprechen.
Der SHA-256 bezieht sich auf die tatsächlich verwendeten float64-Zeilen in
Little-Endian-Reihenfolge.

Für Y=f(X) werden Stichprobenmittel und unverzerrte Stichprobenvarianz mit
Nenner n-1 sowie linear interpolierte empirische Quantile berichtet. Varianz
hat die quadrierte Ausgabe-SI-Einheit. Standardquantile sind 0.025, 0.5, 0.975.
Das Quantilband schätzt Ausgangswahrscheinlichkeiten unter den Eingabegesetzen;
es ist kein Konfidenzintervall, keine garantierte Fehlerschranke und keine
experimentelle Validierung. Präfixgrößen 2,4,8,...,n zeigen Mittel-/Varianzänderung;
ein flaches Diagramm ist kein allgemeiner Konvergenzbeweis. Für kritische
Quantile sind zusätzliche unabhängige Seeds und größere Stichproben nötig.

Numerischer FEM-Fehler, Eingabeunsicherheit und Modellabweichung bleiben getrennt.
Jede Zielgröße verlangt ein explizites deterministisches numerisches Budget
in ihrer SI-Einheit oder `None`, wenn unbekannt. Das Budget wird separat
berichtet und niemals zur statistischen Varianz addiert. Die Begründung eines
Budgets muss aus einem unabhängigen Netz-/Zeit-/Solvernachweis kommen; diese
API ermittelt ihn nicht. Ein kleiner algebraischer Rest belegt keinen Netzfehler.

Die zentrale lokale Sensitivität nutzt f(x+h*ej)-f(x-h*ej), geteilt durch 2*h.
Zentrum und positive SI-Schritte sind explizit; Schrittverfeinerung und zulässige
Perturbationen bleiben Aufgabe des Anwenders. Gradienteneinheit:
Ausgabeeinheit/Eingabeeinheit. Diese Ableitung ist keine globale
Varianzzerlegung oder Sobol-Analyse.

Methodischer Bezug: [JCGM 101:2008, insbesondere Abschnitte 5–7](https://www.bipm.org/documents/20126/2071204/JCGM_101_2008_E.pdf)
beschreibt die Fortpflanzung expliziter Eingabeverteilungen. Hier wird ein
begrenzter Softwareworkflow umgesetzt, keine Konformitätsbewertung oder
vollständige adaptive JCGM-Prozedur behauptet.

## Öffentliche Python-API

Modul `fem_toolbox.uncertainty`; unveränderte Root-Exporte.

| API | Eingaben und Ergebnisse | Fehler/Mutabilität |
|---|---|---|
| `Normal(mean, standard_deviation)` | reelle SI-Skalare | frozen; ValueError für bool/nichtendlich/sigma<=0 |
| `Uniform(lower, upper)` | reelle SI-Skalare | frozen; ValueError für ungültige Grenzen/überlaufende Breite |
| `UncertainParameter(name, target, unit, distribution, justification)` | eindeutiger ASCII-Name, vorhandener skalarer M8-Zielpfad, SI-Einheit, Gesetz, Begründung | frozen; ValueError für leere/ungültige Metadaten |
| `SamplingPlan(parameters, *, size, seed, correlation=None)` | n>=2, nichtnegativer Integer-Seed, R[p,p] | detached; ValueError für Duplikate, Formen, unzulässige Abhängigkeit |
| `SamplingPlan.sample()` | float64[n,p] in deklarierter Reihenfolge | bytebasierter unveränderlicher Array; ValueError bei Überlauf |
| `SamplingPlan.to_document()` | neue JSON-kompatible Schema-1-Spezifikation | Rückgabe ist unabhängig vom Plan veränderbar |
| `summarize_samples(values, *, unit, numerical_error_bound, quantile_probabilities=(.025,.5,.975))` | endlicher reeller Vektor[n>=2]; geordnete Wahrscheinlichkeiten in [0,1] | neue JSON-Statistik; ValueError für bool/complex/NaN/Form/Budget/Quantile |
| `local_sensitivity(evaluate, *, center, steps)` | Vektoren[p], genau 2*p reelle Callback-Ausgaben | detached Callback-Arrays; Callbackfehler propagieren; ValueError für ungültige/zu kleine Schritte |
| `UncertaintyStudy(plan, *, model, root)` | SamplingPlan, relativer stationärer M8-TOML-Pfad | kein Solve im Konstruktor |
| `UncertaintyStudy.run(output, *, quantities, numerical_error_bounds, sensitivity_center=None, sensitivity_steps=None)` | eindeutige M8-Quantity-Definitionen, je Ziel Budget/None, optional Zentrum/Schritte[p] | M8-Schemafehler propagieren, physische Fehler persistieren; FileExistsError bei Lock; IO-/Hashfehler propagieren |

## M8-Anbindung und Ergebnisvertrag

Jede gemeinsame Zeile wird als genau ein M8-Study-Fall mit skalaren Parametern
aufgelöst; kein kartesisches Produkt. Stationäre Wärme, Elastizität,
Elektrostatik, Piezo und Thermopiezo folgen den bestehenden Model-/SI-Verträgen.
Keine neuen Netzgeneratoren, transienten Unsicherheitsläufe oder Tensor-
komponentensweeps. Ungültige Model-Schemas werden vor dem ersten Solve abgelehnt.
Physisch unzulässige Normal-/Uniformwerte werden als M8-Fehlerfälle gespeichert.

Das Verzeichnis darf wiederverwendet werden, wenn `sampling.json` genau zur
Spezifikation passt. Geänderte Modell-/Netzinputs invalidieren M8-Caches.
`sample-NNNNNN/catalog.json` und result bundles behalten Status, Versuche,
Fingerprints, Laufprovenienz, Hashes und Einheiten. Quantity-Tabellen liegen
unter `table-<SHA256 aus Katalog und Definitionen>`; geänderte Fälle oder
Definitionen erhalten eigene Tabellen. Fehlende/manipulierte Dateien sind
Fehler, keine erfolgreichen Messwerte. Ein exklusiver `.uncertainty.lock`
verhindert konkurrierende Writer. Bei hartem Prozessabbruch vor manuellem
Lock-Entfernen prüfen, dass kein Writer mehr läuft. Erfolgreiche M8-Fälle
werden bei erneuter Ausführung übernommen; alte Versuche bleiben erhalten.

`uncertainty.json` wird atomar geschrieben und enthält:

- `uncertainty_schema_version=1`, `sampling`, `sample_sha256`, Quantity-Verträge.
- `samples`: geordnete n Zeilen mit Eingabe-/Ergebnisprovenienz, Katalog- und
  Tabellenpfaden/Hashes, Status, Werten und fehlenden Größen samt Grund.
- `complete`: alle Stichproben erfolgreich und alle angefragten Größen vorhanden.
  Bei false bleiben `summaries` leer; keine Statistik einer unbenannten Teilmenge.
- `summaries`: je Ziel Mittel, Varianz, Quantile, Präfixkonvergenz und separates
  numerisches Budget; keine Diskretisierungs- oder Modellfehler-Schätzung.
- `sensitivity`, `sensitivity_samples`, `sensitivity_status`: zusätzliche
  2*p Perturbationen, Zentrum, Schritte, Gradienten und Einheiten. Misslungene
  Perturbationen unterdrücken Sensitivität, nicht bereits vollständige MC-Statistik.
- `scope`: Eingabeunsicherheit unter angenommenen Gesetzen; keine Messvalidierung.

## Reproduktion und Verifikation

[Ausführbares Wärme-Lehrbeispiel](../examples/uncertainty_heat/README.md):
128 gemeinsame Zeilen, Seed 64, Quelle Uniform(4,12) W/m³ und Temperatur
Normal(300,1) K. Für den diskreten Tet4-Fall gilt Tmax=T0+Q/8,
E=301 K, Var=13/12 K². Jeder FEM-Wert wird mit 1e-11 K, Gradienten [1/8,1]
mit 2e-10 in jeweiligen SI-Einheiten geprüft. Das ist ein diskreter Referenzfall,
kein Kontinuums-Netzfehlernachweis; das numerische Netzbudget bleibt unbekannt.

Vorab definierte statistische Tests: korrelierter Normalfall
X~N(1,2²), Y~N(-2,3²), rho=.5, F=4+2X-Y: E=8, Var=13.
100000 Samples/Seed641; Mittelgrenze 6*sqrt(13/n), Varianzgrenze
6*13*sqrt(2/(n-1)). Uniformer linearer Fall X~U(-1,1), F=3+2X:
E=3, Var=4/3. 50000 Samples/Seed64; Grenzen aus zweitem und viertem
unabhängigen Uniformmoment. Diese festen Tests dienen Softwareverifikation,
keiner nachträglichen Verteilungsanpassung.

```sh
PYTHONPATH=src python -m pytest tests/unit/test_uncertainty.py tests/regression/test_uncertainty_study.py
python examples/uncertainty_heat/run.py --output /tmp/fem-uncertainty-heat
```

Zusätzlich: ungültige Gesetze/Korrelationen/Metadaten, singuläre Normalgesetze,
Konstantgrößen, Seeds/Prefix/Mutationsschutz, sichtbare fehlende Größen,
Materialfehler, Lock, identischer Resume, Modeländerung und zwei
Sensitivitätsschritte. TDD-Red: fehlendes Modul, danach fehlende M8-Anbindung;
Green: 38 gezielte Kern-/FEM-Tests, zusätzlich zwei API-Driftprüfungen. Kein unabhängiges Fachreview wird behauptet.

## Kandidatenprüfung am 10.10.2026

Vollständige lokale Suite nach Wiederaufnahme-Korrekturen: 813 passed, 83 skipped (PETSc lokal nicht installiert),
98.82 s. Ruff, Formatprüfung (341 Python-Dateien), mypy (68 Quelldateien) und
Diffprüfung bestanden. Zusätzlich zwei AST-Vertragstests mit vorherigem Red:
fehlendes Prüfwerkzeug; danach Signatur-/Feldänderung und aktueller Modulvertrag.
Quellenmatrix: 138 Einträge. Der isolierte TeX-Build erzeugt fünf PDFs, verwirft
drei Negativfälle und baut das Quellenarchiv erfolgreich erneut.

Das neue Rezept wurde zusätzlich mit einem frisch installierten Wheel aus dem
entpackten Quellenarchiv unter `python -I` ausgeführt. Import aus dessen
site-packages, 128 vollständige Fälle; größter diskreter Fehler
1.1368683772161603e-13 K. Der Kontinuumsfehler bleibt unbekannt.
Lokale PDFs zeigen den ehrlichen dirty-Kandidatenstand mit Inputhashes;
CI und finale Artefakte müssen den tatsächlichen PR-/Releasecommit identifizieren.

Bei erneuter Prüfung wurden zwei Wiederaufnahme-Risiken mit zuerst roten
Regressionstests geschlossen: Tabellen müssen den aktuellen Kataloghash
belegen; ein abgebrochener Neubau darf keinen älteren Erfolg als aktuellen
`uncertainty.json` stehen lassen. Die CI prüft das neue Rezept einschließlich
Resume aus Wheel und sdist auf Python 3.11–3.13 ohne Checkout im Artefaktjob;
auch älteste/neuste Abhängigkeitskombinationen erhalten UQ-Verifikation.
Die neuen Benutzerseiten und die API-Ergänzung wurden gerendert und visuell
geprüft; lange API-Zeilen sind korrigiert, API-PDF ohne Overfull-Warnungen.
