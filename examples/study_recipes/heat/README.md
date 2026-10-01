# Wärme: zwei Leitfähigkeiten

[Installation, Ausführung und Lebenszyklus](../README.md).
Aufruf: `python -I study_recipes/run.py heat --output heat-results` aus der
sauberen Umgebung außerhalb des Checkouts.

`mesh.msh` ist ein Einheitstetraeder mit Knoten 0=(0,0,0), 1=(1,0,0),
2=(0,1,0), 3=(0,0,1) m, Volumen 1/6 m³, einer Tet4-Zelle und Region
`specimen`. Drei Knoten 0,2,3 haben T=300 K, Knoten 1 ist frei.
Homogene isotrope Leitfähigkeit k=2 bzw. 4 W/(m K), konstante Quelle
Q=8 W/m³, keine Konvektion und keine vorgeschriebenen Flächenflüsse.
Die Materialien sind synthetische Referenzparameter, keine Messdaten.

Die unabhängig integrierten Tet4-Beiträge ergeben K11=k/6 und f1=Q/24.
Daraus folgt **diskret** T1=300+Q/(4k). Das ist keine exakte kontinuierliche
Lösung für eine verteilte Quelle auf einem randarmen Einzeltetraeder.

| k [W/(m K)] | T1 [K] | Quellenleistung [W] | Dirichlet-Inflow [W] |
|---|---|---|---|
| 2 | 301 | 4/3 | −4/3 |
| 4 | 300,5 | 4/3 | −4/3 |

Größeres k senkt den Temperaturanstieg. Negative Dirichlet-Inflow bedeutet
Abfuhr der Quellenleistung an den temperaturgehaltenen Knoten; die globale
Bilanz `source + inflow - prescribed_outflow - convection_outflow` ist null.
Alle vier Knotenfelder werden zusätzlich mit separaten direkten `Model`-Läufen
verglichen, Quellenleistung und Bilanz ebenfalls. Tabelle: Messort Knoten 1,
Leistungen über die gesamte Region. Fall-ID/Fingerprints stehen in jeder Zeile.

Vorab festgelegte Toleranzen: `rtol=2e-11`, `atol=1e-11 K` für Temperaturen,
`atol=1e-11 W` für Leistungskennwerte; absolute Bilanzgrenze 1e-11 W.
Das kleine gut konditionierte lineare System erwartet Rundungsfehler weit
unter diesen Grenzen; Temperaturbasis 300 K und Leistung Größenordnung 1 W
motivieren die absoluten Grenzen. Keine nachträgliche Referenzanpassung.

Budget: zwei Fälle, je vier Knoten/eine Tet4; sequenziell SciPy-LU.
Gemessene Zeiten und Umgebung stehen im [Protokoll](../../../Doc/63_STUDY_RECIPE_ACCEPTANCE.md).
Der automatische Abnahmelauf rechnet außerdem zwei direkte Modelle, zwei
Fehlfallversuche und zweimal zwei Fälle zur Änderung/Wiederherstellung.

Annahmen: stationär, linear, konstante Stoffwerte, SI und ideale Punktvorgaben.
Tet4 stellt T stückweise linear und den Wärmestrom zellweise konstant dar.
Das grobe Netz löst keine allgemeinen Quellen-/Randgradienten auf. Keine
Netzkonvergenz, Unsicherheitsabschätzung oder experimentelle Validierung
wird hier beansprucht. [Wärmevertrag](../../../Doc/11_STEADY_HEAT.md),
[Studienvertrag](../../../Doc/59_PARAMETER_STUDIES.md).

Kontrollierter Stopp/Resume, negative Leitfähigkeit und Verdopplung der
Quelle werden im [gemeinsamen Ablauf](../README.md) reproduzierbar erklärt.
Nach Quellenverdopplung sind T1=302 bzw. 301 K und Quellenleistung 8/3 W;
gültige unveränderte Ergebnisse werden beim Resume nicht neu gerechnet.
