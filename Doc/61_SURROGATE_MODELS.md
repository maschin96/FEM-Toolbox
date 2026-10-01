# Polynomiale Surrogatmodelle aus FEM-Kennwerten (#205)

`fem_toolbox.surrogates` approximiert reelle skalare Zielgrößen stationärer
FEM-Modelle in einem ausdrücklich begrenzten Parameterraum. Mehrere Ausgänge
sind möglich; vollständige räumliche Felder, neuronale Netze, POD/ROM, inverse
Probleme, Optimierung und Unsicherheitsquantifizierung gehören nicht zum Umfang.
NumPy genügt; keine neue ML-Abhängigkeit. Kein Training/Vorhersage-Aufruf löst
selbst ein FEM-Modell. Das eigene Beispiel erzeugt die FEM-Daten zuvor explizit.

## Öffentliche API und Datenvertrag

```python
from fem_toolbox.surrogates import (
    DataSplit,
    PolynomialSurrogate,
    SurrogateData,
    Variable,
    fit_surrogate,
)

data = SurrogateData(
    [[0.0], [1.0], [2.0], [0.5]],  # X[n,p]
    [[300.0], [302.0], [304.0], [301.0]],  # Y[n,q]
    ["a", "b", "c", "test"],
    [Variable("load", "W", "applied heat power")],
    [Variable("temperature", "K", "temperature at the selected node")],
    provenance={"source": "independently specified example"},
)
model = fit_surrogate(data, DataSplit(("a", "b", "c"), ("test",)), degree=1, ridge=0)
prediction = model.predict([0.75])
model.save("surrogate.json")
loaded = PolynomialSurrogate.load("surrogate.json")
```

- `SurrogateData(X, Y, case_ids, parameters, targets, *, provenance)` verlangt
  nichtleere endliche reelle Matrizen `X[n,p]`, `Y[n,q]` und eindeutige nichtleere
  Fall-IDs. Komplexe, boolesche, Objekt-, NaN-/Inf-Daten sind unzulässig. Metadaten
  müssen genau zu den Spalten passen. Namen sind je Liste eindeutig und zwischen
  Parametern/Zielen verschieden. Quelle/Provenienz ist eine nichtleere endliche
  JSON-Tabelle. Arrays sind eigene unveränderliche Speicherabbilder; zurückgegebene
  Provenienz ist eine Kopie. Neue Daten statt Änderung privater Attribute verwenden.
- `Variable(name, unit, definition)` verlangt drei nichtleere Strings. `unit`
  deklariert die SI-Einheit; die Array-API rechnet keine Einheiten um und kann
  die physikalische Wahrheit einer vom Aufrufer angegebenen Einheit nicht prüfen.
  Der Tabellenadapter übernimmt die vom FEM-Export geprüften Einheiten.
- `DataSplit(train_ids, test_ids, seed=None)` bezeichnet eine **vollständige**,
  disjunkte Partition der Fall-IDs. Beide Mengen müssen nichtleer sein.
  Wiederholte IDs und identische Parametervektoren über die Splitgrenze sind Fehler.
  Wiederholungen innerhalb einer Gruppe bleiben erlaubt, zählen aber nicht als
  zusätzliche unabhängige Stützstellen. Identität ist exakte Gleichheit der
  numerischen Parametervektoren; keine versteckte Nähe-/Rundungstoleranz.
- Optional `split_cases(data, test_fraction=..., seed=...)`: deterministische
  Permutation eindeutiger Parametergruppen, Anteil aufgerundet und auf mindestens
  eine Gruppe je Seite begrenzt. Duplikate bleiben gemeinsam; Reihenfolge je Teil
  folgt dem Datensatz. Seed und endgültige Fall-IDs werden gespeichert. Eine kleine
  Trainingsmenge kann anschließend am Rang-/Stützstellengate scheitern.

Testdaten werden ausschließlich zur Auswertung benutzt. Skalierung, Koeffizienten,
Rang und Modellwahl verwenden sie nicht. Grad/Regularisierung bestimmt der Aufrufer
vorher; es gibt keine automatische Hyperparameteroptimierung. Wer wiederholt am
Testfehler auswählt, benötigt für eine belastbare Güteaussage neue unabhängige
Testdaten.

## Adapter für echte Kennwerttabellen aus #60

```python
from fem_toolbox.surrogates import SurrogateData

data = SurrogateData.from_table(
    "study_table",
    parameter_columns=["parameter:source"],
    target_columns=["temperature", "source_power"],
    exclude_invalid=False,
)
```

Der [Tabellenleser](60_STUDY_REPORTS.md) prüft Schema, CSV-Hash und vollständige
CSV-/JSON-Konsistenz. Spalten müssen existieren, eindeutige Rollen haben und
explizit in der gewünschten Modellreihenfolge angegeben werden. Nur skalare
endliche Zahlen sind als Parameter/Ziel zulässig; Vektoren, Strings und Boolesche
Werte werden nicht still kodiert oder imputiert.

Fehlende, nichtendliche oder ungültige ausgewählte Werte, nicht erfolgreiche Fälle,
nicht konvergierte Solver sowie ausdrücklich fehlgeschlagene Bilanzprüfungen werden
standardmäßig abgelehnt. `exclude_invalid=True` erlaubt ihren Ausschluss; die
Provenienz enthält **jede ausgeschlossene Fall-ID mit allen festgestellten Gründen**.
Eine leere Restmenge ist ein Fehler. Nicht geprüfte Bilanzen/Netzkonvergenz werden
nicht zu bestanden umgedeutet; Status und Herkunft bleiben im Datensatz erhalten.

Fall-IDs, Einheiten, Zieldefinitionen, Feldherkunft, Model-/Netz-/Solver-Fingerprints,
Versionen, Ergebnispfade, Katalog-/Tabellen-/CSV-Hashes und Auswahl sind im
Quellnachweis enthalten. Das Training braucht keine ursprünglichen Model-/Netzdateien.
Ungültige globale Tabellen-/Spaltenmetadaten sind keine ausschließbaren Einzelzeilen.

## Referenzverfahren

`fit_surrogate(data, split, degree=1, ridge=0.0)` baut eine vollständige
Polynombasis vom **Gesamtgrad** 1 oder 2, einschließlich Interaktionen.
Die Reihenfolge ist konstant, lineare Terme in Parameterreihenfolge und bei Grad 2
alle `z_i z_j` mit `i <= j` in lexikografischer Indexreihenfolge. Exponentenvektoren
werden im Artefakt gespeichert und beim Laden gegen diese Reihenfolge geprüft.

Aus **Trainingsdaten** werden Minimum `l` und Maximum `h` bestimmt:

```text
center = l/2 + h/2
scale  = h/2 - l/2
z      = (X - center)/scale
```

Konstante oder numerisch nicht skalierbare Trainingsparameter werden ausdrücklich
abgelehnt; sie müssen vom Aufrufer entfernt und im Quellenvertrag dokumentiert
werden. Konstante Zielgrößen sind erlaubt. Mindestens so viele **verschiedene**
Trainingspunkte wie Basisterme sind erforderlich, auch mit Regularisierung:
`1+p` für Grad 1, `1+p+p(p+1)/2` für Grad 2. Unregularisierter Rangmangel ist ein
Fehler; NumPy bestimmt den numerischen Rang nach seiner Standard-SVD-Toleranz.

Gelöst wird stabil mit `numpy.linalg.lstsq`, ohne Normalengleichungen:

```text
min_B ||A B - Y_train||_F² + ridge ||B[1:,:]||_F²
```

Bei `ridge > 0` werden `sqrt(ridge)*diag(0,1,...,1)` und Null-Zielzeilen an das
System angehängt. **Der konstante Term wird nicht bestraft.** `ridge` ist eine
explizite nichtnegative endliche Zahl für dieses skalierte, nicht durch Fallzahl
normalisierte Problem. Ridge darf kollineare Designs stabilisieren; numerischer
Rangverlust auch im erweiterten System wird abgelehnt, etwa bei extremem Gewicht.
Kein automatisches Nachjustieren, keine Optimierungs-/Unsicherheitsbehauptung.

## Vorhersage und Bericht

`predict(X, allow_extrapolation=False)` akzeptiert `[p]` oder `[n,p]`; Ergebnis
`Prediction.values` ist **immer `[n,q]`**, in Zielspaltenreihenfolge und deren
Einheiten. `parameters` und `targets` enthalten die vollständigen Metadaten.
Der Aufrufer muss Eingaben in deklarierter Reihenfolge/Einheit liefern.
`Prediction.extrapolated[n]` markiert jede Zeile außerhalb der gespeicherten
Minimum-/Maximum-Box. Standardmäßig wird eine solche Anfrage abgelehnt;
`allow_extrapolation=True` erlaubt sie ausdrücklich und erhält die Markierung.
Die Grenzen sind inklusive, ohne stilles Clipping oder Toleranzerweiterung.

Die Box ist **keine** konvexe Hülle der Trainingspunkte und kein Genauigkeitsbeweis
innerhalb des Bereichs. Insbesondere Rand-, Zwischenpunkt- und mehrdimensionale
Datenlücken werden durch die Box allein nicht erkannt. Überlauf/nichtendliche
Polynome oder Vorhersagen sind Fehler. Keine Konfidenzintervalle oder garantierten
Fehlerschranken werden ausgegeben.

`model.report` ist eine unabhängige JSON-Kopie mit Datensatz-/Fall-Fingerprints,
Train/Test-IDs/Seed, Grad, Ridge, Designrang, Parameter-/Zielmetadaten,
Trainingsbox, Quellnachweis und Toolbox-/NumPy-/Schema-Version. Getrennt für
Training/Test enthält er je Zielgröße in ihrer Einheit:

```text
MAE     = mean(abs(prediction - reference))
RMSE    = sqrt(mean((prediction - reference)²))
max_abs = max(abs(prediction - reference))
```

Testpunkte außerhalb der Trainingsbox werden im Prüfbericht ausdrücklich als
`test_extrapolated` markiert; sie verändern die Box nicht. Es gibt keine relativen
Fehlerquotienten, die bei Nullzielen undefiniert wären. Empirischer Surrogatfehler,
FEM-Diskretisierungsfehler und experimentelle Modellabweichung bleiben verschieden.

## Inspizierbares JSON-Artefakt, Schema 1

`save(path, overwrite=False)` verlangt `.json` und ein vorhandenes Elternverzeichnis.
Es schreibt zuerst eine temporäre Datei, synchronisiert deren Inhalt und veröffentlicht
sie exklusiv; bestehende Dateien führen zu `FileExistsError`. `overwrite=True`
ersetzt ausdrücklich und atomar. Ein Fehler hinterlässt kein halbes neues Ziel.
Keine Pickle-/Code-Serialisierung, keine neue Bibliotheksabhängigkeit.

Die äußere Tabelle enthält `model` und `artifact_sha256` (SHA-256 über kanonisches
JSON des Modells). Das Modell enthält genau:

- `surrogate_schema_version=1`, Methode `scaled_polynomial_lstsq`, Grad und Ridge;
- Exponentenbasis, Center/Scale, Minimum/Maximum und Koeffizientenmatrix;
- Parameter-/Zieldefinitionen einschließlich SI-Einheiten;
- vollständigen Prüf-/Quellbericht.

`PolynomialSurrogate.load(path)` prüft Hash, Version, Typen, endliche Zahlen,
Formen, Basisreihenfolge, positive Skalierung, Box-/Skalierungskonsistenz und
Bericht-/Metadatenübereinstimmung. Fremde Schemata werden nicht still migriert.
Der öffentliche Konstruktor aus einem Mapping verwendet dieselbe Strukturprüfung;
`to_document()` liefert eine unabhängige Kopie ohne äußere Hashhülle.
Ein Hash schützt gegen unbeabsichtigte Beschädigung, nicht gegen absichtliche
Manipulation. Laden prüft strukturelle Integrität, rekonstruiert aber keine
historischen Trainingsnachweise. Ursprüngliche CSV-/FEM-Dateien sind für spätere
Vorhersagen nicht erforderlich.

## Verifikation und Beispiel

[Wärmebeispiel](../examples/surrogate_heat/README.md) erzeugt acht echte FEM-Fälle,
exportiert eine #60-Tabelle, trainiert auf sechs Fällen und prüft zwei festgelegte
Holdouts. Es speichert/lädt das Artefakt und vergleicht geladene Vorhersagen gegen
FEM und unabhängig hergeleitete diskrete Erwartungen.

`tests/unit/test_surrogates.py` prüft lineare Mehrfachausgänge/konstante Ziele,
quadratische Interaktionen an ungesehenen Punkten, Ridge ohne Intercept-Strafe,
handberechnete MAE/RMSE/Maximalfehler, Train-only-Skalierung, Duplikat-Leckage,
Seeds, Rang/Unterbestimmtheit, Form-/Metadaten-/NaN-Fehler, reale Tabellen mit
Ausschlussprotokoll, Extrapolation sowie Persistenz/Artefaktkorruption. Das bewusst
unzureichende lineare Modell für `y=x²` meldet die Abweichung statt sie zu kaschieren.
Rundungstoleranzen der analytischen Tests: 1e-12 für Größenordnung 1–10,
2e-13 für das kleine gut konditionierte quadratische Gitter; jeweils absolut,
`rtol=0`. Bestehende FEM-Operatoren werden nicht geändert.
