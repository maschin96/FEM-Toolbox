# Quadratische Feldoperatoren und Randintegrale (#35c)

M7 ergänzt auf der [affinen Tet10-Geometrie](41_TET10_GEOMETRY.md) lokale
Operatoren für Wärme, Elektrostatik, Mechanik und Piezoelektrizität einschließlich
Masse, Wärmekapazität, thermischer Ausdehnung und Triangle6-Randlasten.
Der additive API-Einstieg lautet `fem_toolbox.elements.quadratic` bzw.
`fem_toolbox.physics.quadratic`. Bestehende Tet4-Einstiege bleiben unverändert.
Die Auswahl im produktiven Model-/Solverablauf und quadratisches Postprocessing
folgen in #35d; Tet10-Modelle werden dort weiterhin ausdrücklich abgelehnt.

## Lokaler Vertrag

Alle Koordinaten und Materialien verwenden SI-Einheiten und die bisherigen
[Vorzeichen und Voigt-Konventionen](02_MATHEMATICAL_MODEL.md). Die Eingabe-
Knotenordnung bleibt erhalten. Eine skalare Größe hat zehn lokale DOFs,
Verschiebung 30 in knotenweiser xyz-Reihenfolge. Piezo verwendet `[u(30), phi(10)]`.
Materialtensoren und Expansionsdaten sind pro Element konstant; Geometrie ist
geradkantig mit mittigen Kantenknoten. Vorhandene validierte Materialklassen
werden wiederverwendet, ohne neue Bibliotheksabhängigkeit.

| Funktion in elements.quadratic | Ergebnis und Einheit |
|---|---|
| `strain_matrix(geometry, point)` | punktabhängiges B (6,30), 1/m |
| `conductance_matrix(geometry, material)` | Integral Gᵀ k G (10,10), W/K |
| `dielectric_matrix(geometry, material)` | positives Integral Gᵀ κ G (10,10), F |
| `stiffness_matrix(geometry, material)` | Integral Bᵀ C B (30,30), N/m |
| `piezoelectric_matrix(geometry, material)` | symmetrischer Blockoperator (40,40), unten beschrieben |
| `mass_matrix(geometry, density, lumping=...)` | mechanische Masse (30,30), kg |
| `capacity_matrix(geometry, density, specific_heat, lumping=...)` | Kapazität (10,10), J/K |
| `thermal_load(geometry, material, temperature, expansion, ...)` | thermischer RHS (30,) N oder (40,) [N,C] |

Der Piezo-Operator integriert

```text
[ Bᵀ C B      Bᵀ eᵀ G ]
[ Gᵀ e B     -Gᵀ κ G  ]
```

Die Blöcke haben Einheiten N/m, N/V bzw. C/m und C/V. Das negative elektrische
Vorzeichen gehört zur Enthalpieformulierung mit E=-grad(phi), nicht zu einem
negativen Materialtensor. Die transponierten Kopplungsblöcke werden aus derselben
Integration gebildet. Quadratische Gradienten sind ortsabhängig; die konstante
Tet4-Gradientenformel wird nicht übernommen.

`thermal_load` erwartet zehn nodale Temperaturen in K und `ThermalExpansion`.
Es integriert Bᵀ C α (N·T - Tref). Bei piezoelektrischem Material fügt
`electric_thermal_strain=True` zusätzlich Gᵀ e α (N·T - Tref) hinzu; standardmäßig
bleibt der elektrische thermische RHS null, gemäß der bisherigen Wahl zwischen
totaler und mechanischer Dehnung in D. Der Schalter ist für reine Elastizität
ungültig. Keine Pyroelektrizität und keine thermische Rückkopplung werden ergänzt.
Temperaturabweichungen werden vor der Interpolation gebildet, um Auslöschung
mit einem großen konstanten Referenzwert zu verringern.

## Masse und Kapazität

`lumping="consistent"` ist der Standard: exakt integriertes rho Ni Nj bzw.
rho c Ni Nj; mechanisch mit I3 tensoriert. Elektrische DOFs tragen keine Masse.
Dichte und spezifische Wärme müssen endlich und positiv sein. Nicht darstellbare
Gesamtspeicher und unterlaufende Diagonalkoeffizienten werden abgelehnt.

Für Tet10 gilt Integral Ni = -V/20 an Ecken und V/5 an Kantenknoten.
Einfaches Zeilensummen-Lumping wäre daher nicht positiv und ist nicht angeboten.
Explizites `lumping="hrz"` skaliert die positive konsistente Diagonale:

```text
diag_lumped_i = total * diag_consistent_i / sum(diag_consistent)
```

Dies erhält Gesamtmasse je Raumrichtung bzw. Gesamtkapazität. Es reproduziert
nicht sämtliche Momente der konsistenten Matrix. Die Diagonalskalierung folgt
[Oasys GSA: Element Mass / HRZ](https://docs.oasys-software.com/structural/gsa/references-theory/element-mass/).
Ein boolesches Tet4-`lumped=True` wird hier nicht still auf ein anderes Verfahren
umgedeutet. Die Auswahl muss ausdrücklich über den neuen Parameter erfolgen.

## Rand- und Volumenlasten

`physics.quadratic` verwendet `Tet10Geometry` für Volumen und
`Triangle6Geometry` für Randflächen. Lasten sind endliche Skalare/Vektoren oder
Callables mit einem physischen (3,)-Punkt als Eingabe:

| Funktion | Vorzeichen und Größe |
|---|---|
| `source_load` | +Integral N Q, zehn Leistungen in W |
| `body_load` | +Integral N b, 30 Kräfte in N |
| `volume_charge_load` | -Integral N rho, zehn Ladungen in C |
| `flux_load` | -Integral N q_out, sechs Leistungen in W |
| `surface_charge_load` | +Integral N (D·n_out), sechs Ladungen in C |
| `traction_load` | +Integral N t, 18 Kräfte in N |
| `convection_terms` | Integral h Ni Nj und Integral h N T_inf, W/K bzw. W |

Vorgegebene skalare Randdaten wechseln ihr Vorzeichen nicht mit der lokalen
Knotenorientierung. Leiterladung hat weiterhin das Gegenzeichen zu D·n_out.
`Convection` enthält konstantes h>=0 und eine konstante oder räumliche
Umgebungstemperatur. Die Grad-5-Volumenregel ist für Lastpolynome bis Grad 3 exakt,
die Grad-4-Flächenregel für Lastpolynome bis Grad 2. Die Robinmatrix ist exakt.
Nichtpolynomiale Callables benötigen einen eigenen Quadraturfehlernachweis.
Konstante Triangle6-Lasten wirken konsistent mit null Ecklast und A/3 je
Kantenknoten; das ist kein Anlass, die Eckknoten aus einer Elektrode zu entfernen.

Die gemeinsame interne Elektrodenreduktion verwendet nun `mesh.facets` und
`mesh.n_facets`. Feste Potentiale und schwebende Elektroden erfassen alle sechs
Triangle6-Knoten. Gesamtladung, Reduktion und Überlappungsprüfung behalten ihren
Vertrag; die globale Tet10-Solverfreigabe bleibt Bestandteil von #35d.

## Verifikation

- Unabhängige Fakultätsintegrale aller baryzentrischen Basisprodukte prüfen die
  vollständige konsistente Masse/Kapazität. Positive Eigenwerte, Gesamtmasse,
  affine kinetische Energie und HRZ-Erhaltung werden separat geprüft.
- Quadratische skalare Energie, sechs mechanische Starrkörpermoden, Symmetrie
  und die korrekte mechanische Rangzahl 24.
- Affiner Piezo-Patch mit analytischer Spannung, D und Enthalpie; quadratischer
  gekoppelter Patch mit anisotropen Tensoren, analytischer Divergenz und sämtlichen
  Volumen-/Randlasten. Eine negative Orientierung permutiert nur die lokalen DOFs.
- Freie thermische Expansion mit linearer Temperatur und quadratischer
  Verschiebung, sowohl elastisch als auch piezoelektrisch mit explizitem Schalter.
- Lastvorzeichen, konsistente Eck-/Kantenlasten, Robin-Gleichgewicht und feste/
  schwebende Sechsknoten-Elektroden. Ungültige Felder, Parameter und Schalter.

Zusätzlich löst ein ausschließlich im Test enthaltener Sparse-Assembler die
Wärmegleichung mit k=2 auf vier Würfelnetzen n=1,2,3,4. Die glatte Referenz ist
`T=x⁴+y⁴/2+z⁴/4`, mit Q=-24(x²+y²/2+z²/4) und exakten Dirichletwerten. Die
Fehlernormen verwenden eine unabhängige überintegrierte 6-Punkt-Produktregel.
Dies ist ein Code-Verifikationsproblem für Temperaturabweichungen, kein
absolutes thermisches Anwendungsszenario.

| n | L2-Fehler | H1-Seminormfehler |
|---|---:|---:|
| 1 | 0.08030158 | 0.52014538 |
| 2 | 0.01114274 | 0.14324531 |
| 3 | 0.00334559 | 0.06476286 |
| 4 | 0.00141754 | 0.03665694 |

Die letzten Raten sind 2.985 und 1.978 (erwartet 3 und 2); Quellintegral -14 und
Reaktionsbilanz werden auf jedem Netz geprüft. Elementidentitäten verwenden
absolute Toleranzen 2e-16 bis 6e-14 bei den dokumentierten Größenordnungen,
relative Toleranzen bis 3e-14 für anisotrope Operatoren. Globale Bilanz und freie
Residuen verwenden 2e-12 bzw. 1e-12. Ratenfenster 2.7–3.3 und 1.8–2.2 gelten
für die letzten beiden Verfeinerungen. Keine Toleranzaufweitung kaschiert Fehler.

M6-Zeitintegratoren und dynamische Tet10-Ende-zu-Ende-Nachweise sind hier noch
nicht angebunden. Biegungs-/Kostenvergleiche, Model, Backend-Auswahl und
ortsabhängige Ergebnisdarstellung bleiben offen. Insbesondere behauptet dieses
Arbeitspaket keinen vollständigen Abschluss von #35 oder experimentelle Validierung.
