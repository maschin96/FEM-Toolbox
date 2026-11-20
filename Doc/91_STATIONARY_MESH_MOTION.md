# M12 — stationäre harmonische Fluidnetzbewegung (#90)

Stand 21.10.2026, Entwicklungskandidat 1.7.0.dev1.
Maßgeblich ist der [M12-Vertrag](71_M12_RANS_FSI_CONTRACT.md).
Die [konservative Schnittstelle](90_CONSERVATIVE_FSI_INTERFACE.md) liefert
die ebene strukturelle Spur; die bidirektionale Iteration bleibt #91.

## API und Bezugskonfiguration

`fem_toolbox.cfd.move_fluid_domain(reference, interface_displacements_m,
length_scale_m=..., preserve_inlet_flow=True, max_nodes=50000)` erhält eine
validierte `FluidDomain`. Für genau alle ausdrücklich mit `fsi=true`
markierten Wände wird ein Array `[nsegment,2 endpoints,2]` in m übergeben.
Die Reihenfolge folgt `mesh.boundary_sets[name]` und den Endpunkten in `line2`.
Die Strukturspur kommt aus `InterfaceMap.segment_endpoint_displacements`;
dessen Prüfung auf ebene spannweiteninvariante Darstellung bleibt Pflicht.

Für beide Komponenten wird auf demselben unveränderten Referenznetz gelöst:

```text
Laplace(d_mesh) = 0
d_mesh = d_interface             auf FSI-Wänden
d_mesh = 0                       auf übrigen festen Wänden
x_current = x_reference + d_mesh
```

Standard-Bilinearformen für Quad4 und 2x2-Gaußquadratur ergeben die skalare
Laplace-Matrix. Die Elementkoordinaten werden vor Jacobianbildung lokal
verschoben, um große absolute Ursprünge nicht in Ableitungen einzumischen.
Die zwei rechten Seiten werden gemeinsam mit Sparse-LU gelöst. Kein neuer
Laufzeitbaustein und keine neue Abhängigkeit. Der rohe freie Rest hat Einheit m;
sein Rückwärtsfehler, bezogen auf die absoluten Matrixbeiträge mit Floor
1e-14 Lref, muss <=1e-12 sein. Dirichletreaktionen sind keine freien Residuen.

Die Funktion addiert keine vergangenen Inkremente. Wiederholte gleiche
absolute Wege am gleichen Referenznetz liefern dieselben Punkte und Hashes.
Bei jeder stationären Kopplungsiteration muss derselbe Referenzbereich
übergeben werden. Iterationsnummern sind keine physikalische Zeit;
es entstehen weder Wandgeschwindigkeiten noch ALE-Zeitflüsse.

## Ränder, Periodizität und konstante Einlassmenge

Mehrfach zugewiesene Eckknoten müssen innerhalb 1e-12 Lref übereinstimmen.
Widerspruch zwischen FSI und fester Wand wird diagnostiziert, nicht gemittelt.
Jeder Ein-/Auslass muss eine zusammenhängende offene Segmentkette mit zwei
bekannten Wandecken bilden. Dazwischen werden die absoluten Eckwege linear
über die Referenzbogenlänge interpoliert. Andere offene Ränder erhalten die
natürliche Nullnormalgradient-Bedingung.

Explizite periodische Facepartner liefern passende Endpunktpaare. Ihre
Verschiebungen teilen nach Reduktion denselben Freiheitsgrad. Widersprüchliche
periodische Eckvorgaben werden abgelehnt; die neue FluidDomain prüft auch die
aktuelle geometrische Kongruenz erneut. Die vorhandenen Drucksprünge bleiben
physikalische Eingaben; ein periodischer Kanal hat hier keinen zusätzlich
festgeschriebenen Einlassvolumenstrom.

Im Standardmodus erhält jeder echte Einlass seinen ursprünglichen integrierten
Volumenstrom: Das ursprüngliche Geschwindigkeitsprofil wird mit dem Verhältnis
alter zu neuer signierter Flächenintegration skaliert. Breite, Dichte und
SA-Einlasswerte bleiben explizite unveränderte Eingaben. Bei einem gleichmäßig
geöffneten Parallelkanal skaliert U mit H0/H. Das entspricht der unabhängig
geprüften SA-Feststromähnlichkeit aus Doc/89; Druckkräfte werden daraus nicht
analytisch in den Fluidsolver eingesetzt. Ein ausdrücklich gewähltes
`preserve_inlet_flow=False` erhält stattdessen die ursprünglichen
Geschwindigkeitswerte. Die Wahl gehört zur Eingabeidentität und Provenienz.

## Gates und Fehler

Nach der Bewegung werden Punkte, Flächen, Face-Normalen, Qualitätsmaße,
Wandabstände und Randbedingungen neu aufgebaut. `quad4`, `line2`, benannte
Rand- und Volumengruppen behalten ihre Identität. Es gibt keine Wiederverwendung
alter geometrieabhängiger Fluidoperatoren. Die aktuelle Überdeckung zur
verformten Struktur muss vor einer Lastübertragung erneut aufgebaut werden.

Unveränderte Vertragsgrenzen:

- max|d_mesh|/Lref <=0.02; zusätzlich bleiben die getrennten Strukturgrenzen
  für max|u_s|/Lref und max||grad u_s||_F aus Doc/71 verpflichtend.
- Positive konvexe Zellen und positive Quad-Jacobians an Ecken/Gaußpunkten,
  Nichtorthogonalität <=60 Grad, Skewness <=0.25.
- Jede aktuelle Zellfläche mindestens 0.5 ihrer Referenzfläche.
- Wandnormales Schichtwachstum <=1.2. Das finale y+<=1 erfordert anschließend
  die tatsächliche tangentiale Wandschubspannung des neu gelösten Fluids.

`MeshMotionError` bewahrt Fehlerstufe, Referenznetz-SHA-256 und den letzten
absoluten Wegvektor. Stufen trennen Eingaben, Ressourcen, Rand-/periodische
Konflikte, lineare Lösung, Weggrenze, Geometrie, Flächenverhältnis und
Domain-/Wandqualität. Kein Clipping, kein automatisches Remeshing und kein
akzeptiertes Ergebnis bei ausgeschöpftem Knotenbudget.

`MeshMotionResult` enthält das neue unveränderliche Gebiet, absolute Wege,
freie Laplace-Rohresiduen, relativen Rest und minimales Flächenverhältnis.
SHA-256 unterscheiden Referenz-/aktuelles Netz, vollständige numerische
Gebiete, vorgeschriebene Schnittstellenspuren und Gesamteingaben inklusive
unveränderlicher Einstellungen. Dateinamen werden nicht als numerischer
Gebietsinhalt ausgegeben; der spätere Dateiworkflow protokolliert zusätzlich
die tatsächlichen Quelldateien. Auch gleiche Netze mit anderer Einlasspolitik
bleiben anhand der Einstellungen unterscheidbar.

## Nachweise

`tests/unit/test_mesh_motion.py` enthält 13 deterministische Tests:

- Affine Laplace-Patches auf orthogonalem und geschertem Quadnetz, absoluter
  statt kumulierter Bezug, unveränderte Konnektivität und feste Einlassmenge.
- Verzerrter isoparametrischer Patch: genaue affine Wege, erhaltene benannte
  Volumengruppe und neu berechneter Wandabstand. Ein Red-Test reproduzierte
  zunächst den Verlust dieser Gruppe; Übernahme der Referenzgruppen behebt ihn.
- Nicht passende Struktur-/Fluidnetze: konservative P1-Spur, harmonische
  Bewegung und neue vollständige Überdeckung auf aktueller Geometrie.
- Kongruente periodische Bewegung sowie abgelehnte widersprüchliche Partner.
- Getrennte No-go-Fälle für Weg >0.02 Lref, Zellinversion, Fläche <0.5,
  widersprüchliche feste Ecke, fehlende Schnittstelle und Knotenbudget.
- Eine glatte Wandabsenkung erzeugt positive ausreichend große Zellen,
  verletzt aber die Schichtwachstumsgrenze und wird separat abgelehnt.
- Eingabeidentität umfasst die Einlasspolitik; unskalierte Einlasswerte
  verändern den Volumenstrom erwartungsgemäß, ohne die Netzgeometrie zu ändern.

TDD begann mit dem fehlenden Netzbewegungsimport; zusätzlich beobachtete
Red-Fälle für Gruppenverlust und fehlende Einstellungsidentität. Alle 13
Tests sowie alle 21 Schnittstellentests bestehen gemeinsam. M12 bleibt
ein Entwicklungskandidat. Unabhängiges menschliches Fachreview liegt nicht vor.
