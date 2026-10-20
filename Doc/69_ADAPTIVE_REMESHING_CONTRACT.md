# ADR-055 — Konforme Tet4-Verfeinerung und Adaptivitätsvertrag (#67)

2026-10-11 · M10 / v1.5.0 · **Go für den unten begrenzten Umfang**

## Entscheidung und Machbarkeit

Wir verwenden deterministische Kantenstern-Bisektion im eigenen MIT-lizenzierten
Python-Code; kein zusätzlicher Remeshing-Adapter und keine neue Lizenzabhängigkeit.
Für jede markierte Eingabezelle wird ihre längste Kante gewählt (Gleichstand:
lexikographische Knotenindizes). Die sortierte Kantenliste wird abgearbeitet;
jeweils alle inzidenten Tetraeder und Randdreiecke teilen denselben Mittelpunkt.
Auch unmarkierte Nachbarzellen werden dadurch verfeinert. Die Beschränkung auf
geradkantige Tet4 bewahrt die Geometrie exakt.

Der Mehrgebietsfall in `tests/unit/test_refinement.py` prüft vollständige
Randtopologie, unverändertes Volumen beider Regionen, Tags, benannte Rand- und
Volumengruppen, deterministische Konnektivität und lokale statt globale Verfeinerung.
Der Red-Test scheiterte am fehlenden Modul `mesh.refinement`; anschließend bestanden
beide Machbarkeitstests. Manuelle Vertragsprüfung: Gemeinsame Kanten werden im
vollständigen Stern geteilt; jede gemeinsame Fläche erfährt dieselben Teilungen.
Damit gibt es keine hängenden Knoten. Kinder erben ihren ursprünglichen Elternindex;
Materialgrenzen werden nicht überschritten. Go bedeutet Machbarkeit dieser
Strategie, keine universelle Effizienz- oder Netzqualitätsgarantie.

`validate_refinement` prüft Elternindizes, Tags, auch überlappende Gruppen,
Kindlage im Elternsimplex, summierte Elternvolumina/-flächen und vollständige
konforme Randtopologie. Es ist ein Vertrag für diesen Bisektionsadapter, kein
allgemeiner Schnittdetektor für beliebige geometrisch überlappende Netze.
Externe Element-/Knoten-IDs werden bei Verfeinerung verworfen; neue Identität
ist Stufe + Netz-Hash, Abstammung sind explizite Elternzeilen.

## Fehlermaß, Budgets und Abbruch

Stationäre Wärme, affine Tet4, zellweise konstante SPD-Wärmeleitfähigkeit.
Das vorhandene Residualmaß `heat_indicator` liefert eta in sqrt(W K).
`eta <= tolerance` ist eine **heuristische Indikatortoleranz**, keine garantierte
Schranke des tatsächlichen Energiefehlers. Dörfler-Markierung verwendet standardmäßig
50 % der quadrierten Indikatorsumme. Zeitabhängigkeit und Kopplung sind ausgeschlossen.

`AdaptiveOptions`: nichtnegative endliche Toleranz; positive ganze maximale
Anzahl gelöster Netze (`max_steps`) und aktiver Temperatur-DOFs (`max_dofs`);
Markierungsanteil `(0,1]`; minimale relative Indikatorabnahme `[0,1)` und positive
Stagnationsfensterlänge. Standards: 1e-3, 10 Netze, 100000 DOFs, 0.5, 0.01, 3.
Ein überschrittenes Verfeinerungsbudget akzeptiert kein Teilnetz. Endstatus:
`tolerance`, `dof_budget`, `step_budget`, `stagnation`, `error`.
Bereits gelöste Stufen bleiben bei Fehlern erhalten. Priorität nach dem Solve:
Toleranz, Schrittbudget, Stagnation, Verfeinerung/DOF-Budget. Ungültige Anfangsnetze
werden vor dem Solve diagnostiziert. Softwarefehler außerhalb erwarteter numerischer
oder Validierungsfehler dürfen nicht als erfolgreiche Rechnung erscheinen.

Randwerte werden auf Eingabe-Facetten definiert: Dirichlet-K, Auswärtsfluss W/m²,
Robin W/(m² K) und Umgebung K. Nicht angegebene Facetten bleiben isoliert.
Kinder erben Facettenzuordnung; räumliche Funktionen werden neu ausgewertet.
Punktbedingungen, Koarsening, Glättung, gekrümmte Flächen und Zustandsübertragung
sind nicht enthalten. Eine formale gleichmäßige Formregularitätsgarantie wird
nicht behauptet; degenerierte Zellen werden mit vorhandenen Geometriegates abgelehnt.

## Nachweisvertrag

Jede Stufe speichert unveränderliches Netz, Elternzeilen, Rohtemperatur/-gradient,
Fluss, Markierung, eta, DOFs, Wärmebilanz und Solverprovenienz. Striktes JSON enthält
Netz/Tags/Gruppen und SHA-256 sowie die Abbruchkonfiguration. Identische Eingaben
liefern identische Records; Laufzeiten gehören ausschließlich in den Benchmark.
Ein Hash ersetzt keinen Funktionsquelltext: für erneute Ausführung von Callables
ist das ursprüngliche Problemskript erforderlich. Das ausführbare Beispiel ist
Teil des Reproduktionsnachweises.

#68 implementiert die Schleife; #69 prüft glatte und lokalisierte analytische
Referenzen unabhängig vom Indikator über mindestens vier Stufen und vergleicht
mit uniformer Kantenverfeinerung. #70 bleibt bis zur tatsächlichen Release-Abnahme
mit Installation, CI und Veröffentlichungsnachweis offen.
