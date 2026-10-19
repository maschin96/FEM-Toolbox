# SIMP-Optimierung: Material dort einsetzen, wo es trägt

Eine Einführung zur Struktur-Topologieoptimierung in M11.

## 1. Welche Frage beantwortet SIMP?

Bei einer gewöhnlichen FEM-Rechnung stehen Geometrie, Material, Lagerung und
Lasten bereits fest. Gesucht sind die daraus entstehenden Verschiebungen und
Spannungen. Die Topologieoptimierung stellt eine zusätzliche Frage: **Wo sollte
innerhalb eines erlaubten Bauraums Material liegen, damit eine vorgegebene
Aufgabe möglichst gut erfüllt wird?**

Ein Beispiel ist ein Kragarm: Eine Seite ist eingespannt, am freien Ende wirkt
eine Kraft. Ein vollständig ausgefüllter Quader kann diese Last tragen, benötigt
aber viel Material. Die Optimierung sucht innerhalb desselben Bauraums einen
steifen Entwurf mit begrenztem Materialeinsatz. Bereiche an Klemme und
Lasteinleitung können dabei ausdrücklich als Vollmaterial vorgeschrieben bleiben.

SIMP steht für **Solid Isotropic Material with Penalization**. Die Methode
beschreibt eine mögliche Materialverteilung durch kontinuierliche Dichtewerte.
Dadurch kann sie mit Gradienten arbeiten, statt jede mögliche Kombination aus
„Material vorhanden“ und „Material entfernt“ einzeln durchzurechnen.

## 2. Von einer Geometrie zu einem Dichtefeld

Der Bauraum wird zunächst mit finiten Elementen vernetzt. In M11 sind das
Tetraeder mit vier Knoten, also Tet4. Jedes Element bekommt eine dimensionslose
Dichte zwischen 0 und 1:

- **1** bedeutet Vollmaterial.
- **0** bedeutet einen rechnerisch weitgehend leeren Bereich.
- **Zwischenwerte** sind mathematische Entwurfszustände, häufig „Graubereiche“ genannt.

Diese Dichte ist keine aus einer Messung bestimmte Massendichte in kg/m³. Sie
beschreibt den Materialanteil im Optimierungsmodell. Ein Wert von 0,5 ist auch
keine unmittelbare Bauanweisung für einen halbporösen Werkstoff. Dazu wäre ein
zusätzliches, physikalisch begründetes Materialmodell nötig.

Die gesamte Materialmenge ergibt sich aus den Elementvolumina und Dichten:

$$
v = \frac{\sum_e V_e\rho_e}{\sum_e V_e}.
$$

Bei einer Volumengrenze von 0,6 darf das Dichtefeld höchstens 60 Prozent des
ursprünglichen Bauraumvolumens belegen. Feste Vollmaterialbereiche zählen dabei
mit. Lassen sie bereits zu wenig Spielraum, kann die Aufgabenstellung unzulässig
sein; die Software muss das melden.

## 3. Warum werden mittlere Dichten „bestraft“?

Ein gewöhnlicher elastischer Körper besitzt eine Steifigkeitsmatrix. SIMP
skaliert den Beitrag jedes Elements abhängig von seiner physischen Dichte:

$$
K_e(\rho_e)
= \left[\varepsilon+(1-\varepsilon)\rho_e^p\right]K_e^0.
$$

Dabei ist $K_e^0$ die Elementsteifigkeit bei Vollmaterial. Der Exponent $p$
steuert die Bestrafung von Zwischenmaterial. Die kleine positive Zahl
$\varepsilon$ verhindert, dass leere Elemente exakt ihre gesamte Steifigkeit
verlieren und dadurch das Gleichungssystem unlösbar wird.

Für die M11-Standardwerte $p=3$ und $\varepsilon=0{,}001$ ergibt sich:

| Dichte | Angerechnetes Materialvolumen der Zelle | Anteil der Vollmaterialsteifigkeit |
|---:|---:|---:|
| 1,00 | 100 % | 100 % |
| 0,50 | 50 % | etwa 12,59 % |
| 0,25 | 25 % | etwa 1,66 % |

Mittlere Dichten verbrauchen damit vergleichsweise viel Material für wenig
Steifigkeit. Das begünstigt eine Verteilung mit tragenden Vollmaterialbereichen
und weitgehend leeren Bereichen. Es erzwingt jedoch keinen vollständig binären
Entwurf. Insbesondere nach der Filterung können erhebliche Graubereiche bleiben.

Auch wichtig: Eine Zelle mit Dichte null ist wegen der Reststeifigkeit noch kein
geometrisch entferntes Loch. Dichteentwurf und tatsächliche Bauteilgeometrie
müssen deshalb getrennt beurteilt werden.

## 4. Was wird minimiert – und wie läuft die Suche ab?

M11 minimiert die **Nachgiebigkeit**. Bei einer festen äußeren Kraft ist ein
Körper umso nachgiebiger, je weiter er sich unter dieser Kraft verformt. Die
Zielfunktion lautet für einen Lastfall:

$$
C=f^T u.
$$

Der Kraftvektor $f$ und der Verschiebungsvektor $u$ ergeben eine Größe mit der
Einheit Joule. Bei linearer Elastizität und homogenen Verschiebungsvorgaben ist
sie doppelt so groß wie die gespeicherte Formänderungsenergie. Unter denselben
festen Kräften bedeutet kleineres $C$ einen steiferen Entwurf.

Die Optimierung wiederholt im Wesentlichen fünf Schritte:

1. Ein zulässiges Dichtefeld festlegen und daraus die Steifigkeiten aufbauen.
2. Mit der FEM die Verschiebungen für jeden Lastfall berechnen.
3. Nachgiebigkeit, Materialvolumen und Gleichgewicht prüfen.
4. Berechnen, wie empfindlich Ziel und Volumen auf jede Entwurfsvariable reagieren.
5. Die Variablen innerhalb ihrer Schranken und des Volumenbudgets verändern.

Die Empfindlichkeiten heißen **Sensitivitäten** oder Gradienten. Sie zeigen,
wo eine Änderung des Materialeinsatzes besonders wirksam wäre. M11 berücksichtigt
in dieser Ableitung auch den Dichtefilter. Unabhängige zentrale Differenzen
prüfen die Gradienten, indem leicht veränderte Entwürfe tatsächlich neu gerechnet
werden. Ein optisch plausibles Ergebnis allein wäre dafür kein ausreichender Test.

M11 verwendet SLSQP als Optimierungsverfahren. Es berücksichtigt Zielgradient,
Volumenrestriktion und Variablenschranken. Ein erfolgreicher Stopp bezeichnet eine
lokale numerische Konvergenz. Er beweist nicht, dass der weltweit beste Entwurf
gefunden wurde. Ein ausgeschöpftes Iterationsbudget bleibt ausdrücklich ein
Budgetstopp.

## 5. Wozu dient der Dichtefilter?

Unabhängig veränderbare Elementdichten können sehr feine, netzabhängige Muster
erzeugen. Der Filter mittelt deshalb Entwurfswerte aus einer räumlichen
Nachbarschaft. Nahe Elemente erhalten mehr Gewicht als weiter entfernte;
ungleich große Elemente werden über ihr Volumen berücksichtigt.

Der Filterradius besitzt die Einheit Meter. Er ist damit eine räumliche
Modellentscheidung, keine bloße Anzahl von Nachbarzellen. Ein größerer Radius
führt typischerweise zu glatteren Dichteverteilungen, kann aber auch mehr
Zwischenmaterial hinterlassen. Er garantiert in diesem Verfahren keine bestimmte
Mindestwandstärke des späteren Bauteils.

Die Software unterscheidet deshalb zwischen den **Entwurfsvariablen** und den
**gefilterten physischen Dichten**. Steifigkeit und Volumen beziehen sich auf
Letztere. Fest vorgeschriebene Voll- und Leergebiete müssen dabei auch nach der
Filterung ihre festgelegten Werte behalten.

## 6. Mehrere Lastfälle und die Grenze des Optimierungsziels

Ein Bauteil kann in mehreren Richtungen belastet werden. M11 berechnet diese
Lastfälle getrennt am selben Entwurf und minimiert ihre positiv gewichtete Summe:

$$
C_{\mathrm{gesamt}} = \sum_l w_l C_l.
$$

Die Gewichte legen fest, wie stark jeder Fall die Entscheidung beeinflusst.
Eine gewichtete Summe ist keine Garantie, dass der ungünstigste Lastfall allein
optimal behandelt wird. Ebenso wenig werden dadurch alle Lasten automatisch
zu einer gleichzeitigen physikalischen Belastung überlagert.

Nachgiebigkeitsminimierung ist außerdem keine Spannungsoptimierung. Ein steifer
Entwurf kann örtlich hohe Spannungen besitzen. Zulässige Spannungen, Ermüdung,
Beulen und Fertigbarkeit sind zusätzliche Fragestellungen und werden durch das
hier verwendete Ziel nicht automatisch abgesichert.

## 7. Vom Dichtebild zum erneut geprüften Bauteil

Aus einem Grauwertfeld entsteht noch keine freigegebene Geometrie. M11 wählt
zunächst über eine dokumentierte Schwelle aus, welche Tetraeder erhalten bleiben.
Dabei können getrennte Teile entstehen oder tragende Verbindungen verloren gehen.
Solche Ergebnisse werden als ungültig ausgewiesen.

Ein zweiter, ausdrücklich gekennzeichneter Weg ergänzt deterministische
Verbindungszellen. Diese Ergänzung fügt Material hinzu und wird vollständig
protokolliert. Sie darf nicht als unveränderte Schwellengeometrie ausgegeben werden.

Anschließend rechnet der vorhandene Mechaniksolver die diskrete Geometrie mit
Vollmaterial, denselben Lasten und den erhaltenen Lagerungen erneut. Erst dieser
Schritt zeigt, welche Nachgiebigkeit und welches Volumen das abgeleitete Bauteil
besitzt. Seine gesonderte Geometrievolumengrenze ersetzt nicht das ursprüngliche
Dichtevolumenbudget. Rohe Spannungsspitzen werden berichtet, begründen aber noch
keine Festigkeitsfreigabe.

Der M11-Workflow ist damit ein nachvollziehbarer Forschungsworkflow für
stationäre, lineare Strukturmechanik mit kleinen Verformungen und festen Lasten.
Er unterstützt eine begründete Entwurfsentscheidung; eine globale Optimalität,
experimentelle Bauteilvalidierung oder Fertigungsfreigabe folgt daraus nicht.

**Weiterlesen:** [Ausführbares Beispiel und API](84_M11_TOPOLOGY_WORKFLOW.md),
[mathematischer Vertrag](83_M11_TOPOLOGY_CONTRACT.md) und
[Verifikationsnachweise](85_M11_VERIFICATION.md).
