# Synthetische Unsicherheitsfortpflanzung

Aus dem Repository-Root nach Installation:

```sh
python examples/uncertainty_heat/run.py --output /tmp/fem-uncertainty-heat
```

128 gemeinsame Stichproben, Seed 64: Quelle Q uniform auf [4,12] W/m³ und
Randtemperatur normal mit Mittel 300 K und Standardabweichung 1 K. Beide Gesetze
sind ausdrücklich synthetische Lernannahmen, keine gemessenen Unsicherheiten.
Konstante Leitfähigkeit 2 W/(m K), ein Tet4 mit drei temperaturfixierten Knoten.
Für diesen diskreten Fall gilt exakt Tmax=T0+Q/8, E[Tmax]=301 K und
Var[Tmax]=13/12 K². Das Skript prüft jeden diskreten FEM-Wert mit 1e-11 K
und die zentralen Ableitungen [1/8,1] mit 2e-10 in den jeweiligen SI-Einheiten.
Der Kontinuums-Netzfehler ist damit nicht bestimmt; sein Budget bleibt `None`.

`sampling.json` hält Gesetze, Begründungen, Korrelation, RNG, Versionen und Seed
fest; `uncertainty.json` enthält Quantile, Präfixkonvergenz, Sensitivitäten und
vollständige M8-Belege. `verification.json` berichtet die analytische Referenz.
Unter `sample-*` und `sensitivity-*` liegen Kataloge, Ergebnisse und unveränderliche
Mengentabellen. Erneuter identischer Aufruf übernimmt erfolgreiche FEM-Läufe;
eine andere Sampling-Spezifikation verlangt ein neues Ausgabeverzeichnis.

Quantilband: geschätzte Ausgangswahrscheinlichkeit unter den angenommenen
Gesetzen, kein Konfidenzintervall und keine garantierte Fehlergrenze. Eine
128-Punkte-Stichprobe belegt keine allgemeine Quantilkonvergenz. Normalgesetze
sind unbeschränkt; unphysikalische Fälle werden sichtbar abgelehnt und verhindern
Statistiken über eine still ausgewählte Restmenge. Experimentelle Validierung
bleibt außerhalb dieses Beispiels offen.
