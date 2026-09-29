# Wärme-Surrogat: FEM → Tabelle → Training → geprüftes Modellartefakt

Mit installiertem Paket:

```sh
python examples/surrogate_heat/run.py --output surrogate_heat_results
```

Das Ausgabeziel muss neu sein. Der vollständige Lauf enthält acht echte Tet4-
Wärmefälle, CSV-/JSON-Kennwerte, ein gespeichertes `surrogate.json` und
`report.json`. Keine optionale ML-Bibliothek oder Gmsh-Installation nötig.
Das mitgelieferte Netz ist der Einheitstetraeder aus dem stationären Referenzbeispiel.

Konstantes k=2 W/(m K), drei Knoten auf 300 K, freier Knoten bei x=1 m;
Q variiert in W/m³. **Vorab festgelegt**: Training bei Q=0,2,4,6,8,10;
Test bei Q=3,7. Ein lineares Modell bildet zwei Ausgänge ab: maximale aktive
Knotentemperatur in K und integrierte Quellenleistung in W. Die positive
Quellenrichtung macht die freie Knotentemperatur zum Maximum.

Unabhängige diskrete Referenz: V=1/6 m³, K11=k/6, f1=Q/24,
also T1=300+Q/(4k)=300+Q/8 und P=Q/6. Vor der Auswertung festgelegte
Grenzen: 1e-10 K und 1e-12 W, jeweils absolut; sie berücksichtigen Rundung
bei 300-K-Grundtemperatur bzw. Leistungsgrößen der Ordnung 1 W. Die FEM-
Leistungsbilanz wird mit 1e-11 W geprüft. Dies ist ein synthetischer
Verifikationsfall, keine experimentelle Validierung oder exakte
Kontinuumslösung des volumetrisch beheizten Tetraeders.

Die End-to-End-Prüfung vergleicht geladene Modellvorhersagen mit den beiden
zurückgehaltenen FEM-Läufen. Die Originaldateien werden später nicht benötigt:

```python
from fem_toolbox.surrogates import PolynomialSurrogate

model = PolynomialSurrogate.load("surrogate_heat_results/surrogate.json")
print(model.predict([5.0]).values)  # [[300.625, 0.833333...]]
```

Außerhalb Q∈[0,10] ist eine ausdrückliche Extrapolationsfreigabe nötig; auch
innerhalb gilt keine allgemeine Genauigkeitsgarantie. Trainingsfehler und
Testfehler stehen getrennt im Bericht. [Vollständiger Vertrag](../../Doc/61_SURROGATE_MODELS.md).
