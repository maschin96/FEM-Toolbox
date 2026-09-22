# Vier-Netze-Nachweis des stationären Mehrmaterialmodells (#141)

Status: **Code- und Lösungsverifikation** des deklarierten
`T → (u, φ)`-Workflows. Es gibt keine experimentelle Validierung. Der
ausführbare Nachweis liegt in
[`tests/convergence/test_thermopiezo_model_convergence.py`](../tests/convergence/test_thermopiezo_model_convergence.py),
die unabhängige Kontinuumsreferenz in
[`tests/reference/thermopiezo_layer_1d.py`](../tests/reference/thermopiezo_layer_1d.py).

## Geometrie, Netz und Modellvertrag

Ein Quader `[0,1]³ m` mit Einheitsquerschnitt besteht aus zwei perfekt
gebundenen Schichten `0≤z<0,5` und `0,5<z≤1`. Die regulären, unverzerrten
Netze teilen jeden der `n³` Würfel in sechs Tet4. Für `n=2,4,6,8` entstehen
48, 384, 1296 und 3072 Elemente. Der Generator in `tests/conftest.py`
verwendet Koordinatengitter und alle sechs Achsenpermutationen; er ist von der
Produktions-Netzverwaltung unabhängig. Die Grenzfläche `z=0,5` ist auf jedem
Netz konform. `volume_sets` heißen `lower`/`upper`, `facet_sets`
`bottom`/`top`.

Jeder Lauf wird mit `Model(mesh=..., document=...)` aus demselben
versionierten Schema 1 konstruiert. Das Dokument enthält getrennte
Materialzuordnungen für `T`, `piezo` und `expansion`, konstante Volumenlasten,
mechanische Lagerungen, thermische Dirichlet-Werte und zwei benannte leitende
Elektroden. `electric_thermal_strain=true` ist explizit gesetzt. Unten liegt
`φ=0`; oben wird entweder `φ=0,7 V` oder die Gesamtladung `Q=0,3 C`
vorgegeben. Weitere Kurzschluss- (`φ=0`) und Leerlaufläufe (`Q=0`) prüfen die
physikalisch verschiedenen Grenzen.

| Größe | Untere Schicht | Obere Schicht | Einheit |
|---|---:|---:|---|
| Wärmeleitfähigkeit `k` | 2 | 4 | W/(m K) |
| `C^E_zz` | 12 | 18 | Pa |
| `e_zz` | 1,5 | 2 | C/m² |
| `κ^S_zz` | 2,5 | 3 | F/m |
| `α_zz` | 0,001 | 0,0015 | 1/K |

Alle anderen Steifigkeitsdiagonalen sind `20,20,4,4,4 Pa`, alle anderen
Permittivitätsdiagonalen `5 F/m`, die sonstige Kopplung ist null. Die
Engineering-Voigt-Reihenfolge ist `xx,yy,zz,yz,xz,xy`. `T_ref=300 K`,
`T(0)=310 K`, `T(1)=300 K`, Wärmeerzeugung `8 W/m³`, Volumenkraft
`b_z=0,5 N/m³` und freie Ladungsdichte `ρ=0,2 C/m³`. `u_x=u_y=0` an allen
Knoten und `u_z=0` an Ober- und Unterseite; die seitlichen Flächen sind
fluss- beziehungsweise traktionsfrei. Die Werte sind synthetische SI-Daten.

## Unabhängige Kontinuumsreduktion

Die Symmetrie und das Materialmodell liefern in jeder Schicht Felder, die
nur von `z` abhängen. Wärme-, Kraft- und Ladungsbilanz ergeben

```text
q(z)     = q₀ + 8z
σ_zz(z)  = σ₀ - 0,5z
D_z(z)   = D₀ + 0,2z
```

Aus `T(0)-T(1)=10 K` und der Reihenschaltung
`∫₀¹(q₀+8z)/k(z) dz=10 K` folgt `q₀=70/3 W/m²`. Damit sind
`q(0,5)=82/3 W/m²`, `T(0,5)=303⅔ K` und `q(1)=94/3 W/m²`.
Die Kontinuumslösung ist je Schicht quadratisch in `T` und stetig an der
Grenzfläche.

Für `ε_m=u'_z-α(T-T_ref)` und `E_z=-φ'` gilt

```text
[σ_zz, D_z]ᵀ = [[C, e], [e, -κ]] [ε_m, φ']ᵀ.
```

Die Referenz invertiert diese **2×2-Stoffmatrix** analytisch und integriert
`u'` und `φ'` über beide Schichten. Die zwei Integrationskonstanten `σ₀,D₀`
folgen aus `u(1)=0` und entweder `φ(1)=0,7 V` oder
`Q_top=-D(1)·1m²=0,3 C`. Das Vorzeichen von `Q_top` bezeichnet die Ladung
des angrenzenden Leiters. Fünfpunkt-Gaußintegration liefert für die hier
auftretenden Polynome exakte Integrale. Es werden keine FEM-Matrizen oder
Produktions-Solver im Referenzmodul verwendet.

| Referenzgröße | Spannungselektrode | Ladungselektrode |
|---|---:|---:|
| `Q_top` [C] | `1,813745944993` | `0,3` |
| `φ_top` [V] | `0,7` | `0,145064180455` |
| `σ₀` [Pa] | `1,353734661495` | `0,405099663909` |
| Gesamtenthalpie [J] | `-0,674201464141` | `-0,029980511931` |

Der getrennt integrierte Wärme-Energieterm `∫ k|∇T|² dΩ` beträgt
`268,5` in den hier aus SI folgenden Einheiten `W·K`; er ist keine
gespeicherte thermische Energie.

## Fehlermaße und Gates

Eine Fünfpunkt-Duffy-Quadratur je Referenzkoordinate integriert die
quadratischen Fehler der stückweise kubischen Primärfelder auf jedem Tet4.
Für `T`, `u_z` und `φ` werden jeweils ein absolutes 3D-`L²`-Maß und ein
`H¹`-Seminormfehler berechnet. Letzterer umfasst auch die numerisch
auftretenden Quergradienten. Die Energiefehler vergleichen
`∑_K k_K|∇T_h|²|K|` und die vorzeichenbehaftete elektrische Enthalpie mit
den unabhängig integrierten Referenzen. Bei `h=1/n` sind für glatte
Teilbereiche die Ordnungen `L²≈2`, `H¹≈1` und Energiefehler `≈2` erwartet.
Die drei aufeinanderfolgenden Schätzer werden auf allen vier Netzen geprüft,
mit festen Intervallen `1,75–2,25`, `0,75–1,25` beziehungsweise `1,5–2,5`.

| Netz `n` | `L²(T)` [K] | `L²(u_z)` [m] | `L²(φ)` [V] | Wärme-Energiedifferenz [W·K] | Enthalpiedifferenz [J], Spannung |
|---:|---:|---:|---:|---:|---:|
| 2 | `7,217e-2` | `8,350e-4` | `2,030e-3` | `5,000e-1` | `1,510e-4` |
| 4 | `1,804e-2` | `2,086e-4` | `5,075e-4` | `1,250e-1` | `3,777e-5` |
| 6 | `8,019e-3` | `9,271e-5` | `2,255e-4` | `5,556e-2` | `1,679e-5` |
| 8 | `4,511e-3` | `5,215e-5` | `1,269e-4` | `3,125e-2` | `9,444e-6` |

Die Spannungselektrode und die Ladungselektrode erfüllen dieselben
Ordnungsfenster. Die integrierte Elektrodenantwort wird gegen die
Kontinuumsreferenz mit absolut `10⁻³h²` geprüft. Für die Grenzfläche
werden die Zellflüsse auf `z=0,5` mit den analytischen linearen
Quellsteigungen extrapoliert und von beiden Seiten getrennt verglichen:
`10⁻⁹ W/m²` für `q`, `0,02h² Pa` für `σ_zz` und
`0,004h² C/m²` für `D_z`. Die Fenster folgen aus der glatten
stückweisen Lösung und sind klein gegen die Flussgrößen `27,3 W/m²`,
`0,1–1,1 Pa` und `0,2–1,9 C/m²`.

Globale Wärme-, Kraft-, Moment- und Ladungsbilanzen bleiben jeweils unter
`10⁻⁸` in ihrer SI-Einheit; die feldweisen relativen algebraischen
Residuen unter `10⁻¹⁰`. Unabhängig von diesen Residuen werden
`source_power=8 W`, `dirichlet_inflow=-8 W`, die gesamte Lagerreaktion
`-0,5 N` und die Summe der Elektrodenladungen `-0,2 C` geprüft. Eine
Variation von oberer Verschiebung und Elektrodenspannung bestätigt
`∂F_top/∂V = -∂Q_top/∂u_top` bis `2·10⁻¹⁰` in den jeweiligen
Quotienteneinheiten. Ein absichtlich umgekehrtes Vorzeichen der
vorgegebenen Elektrodengesamtladung erzeugt ein anderes Potential und
wird als Negativfall erkannt.

Der Test läuft als Teil von `python -m pytest`; er kann gezielt mit
`python -m pytest tests/convergence/test_thermopiezo_model_convergence.py`
ausgeführt werden. Auf der lokalen Python-3.13-Umgebung benötigt er etwa
7 Sekunden. Er gehört zu den für numerische Änderungen und vor dem
v1.0-Release erforderlichen Konvergenzprüfungen.

Die analytische Reduktion prüft den vollständigen gekoppelten Material- und
Elektrodenpfad auf einem echten 3D-Netz. Sie deckt keine allgemeinen
Quergradienten des **exakten** Feldes, geometrisch gekrümmte Grenzflächen,
nichtkonforme Netze oder reale Materialdaten ab. Dafür sind gesonderte
Patch- und Validierungsfälle notwendig.
