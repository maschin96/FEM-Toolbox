# Stationäre Tet10-Thermomechanik (#35d)

Dieser Teil von #35d verbindet die verifizierten Tet10-Wärme- und Mechanikpfade.
`solve_elasticity(..., temperature=..., expansion=...)` und der bestehende
Python-Einstieg `solve_thermomechanics` unterstützen geradkantige Tet10.
Die Kopplung bleibt stationär und einseitig: T → u. Keine mechanische
Wärmerückwirkung, keine reziproke Thermoelastizität und keine transiente
Strukturkopplung. Piezo-/Thermopiezointegration, Dynamik und Kostenvergleiche
bleiben offen; #35 ist damit nicht abgeschlossen.

Der nachfolgende [stationäre Piezoschnitt](49_TET10_PIEZOELECTRIC.md)
ergänzt inzwischen Tet10-Piezo; die obige Abgrenzung beschreibt den damaligen Stand.

## Thermischer Last- und Ergebnisvertrag

Temperatur ist weiterhin ein absolutes Knotenfeld der Länge `mesh.n_points` in K.
Expansionsdaten sind konstante globale Engineering-Voigt-Vektoren pro Element
in 1/K mit regionsweiser Referenztemperatur. Die vorhandene Rotation mittels
`ThermalExpansion.from_principal` bleibt verwendbar. Beide thermischen Angaben
müssen gemeinsam vorliegen; falsche Arrayformen, negative Knotentemperaturen
und unvollständige Materiallisten werden weiterhin abgelehnt.

Für Tet10 gilt an jedem Integrationspunkt:

```text
theta(x) = sum_i N_i(x) * (T_i - T_ref)
epsilon_th(x) = alpha * theta(x)
f_th = integral B(x)^T C epsilon_th(x) dV
sigma(x) = C (B(x) u - epsilon_th(x))
W = 1/2 integral (B(x)u - epsilon_th(x))^T C (B(x)u - epsilon_th(x)) dV
```

Die vorhandene quadratische Operatorroutine integriert die thermische Last.
Die Grad-5-Volumenregel integriert sowohl den Lastintegranden (maximal Grad 3)
als auch die mechanische Energiedichte (maximal Grad 4) exakt für die
unterstützten konstanten Materialdaten und affine Geometrie. Der Temperatur-
Referenzoffset wird vor der Interpolation subtrahiert.

`ElasticResult.strain` bleibt Gesamtdehnung. `stress` zieht die thermische
Eigendehnung **am Schwerpunkt** ab. Die Energie integriert dagegen die
punktabhängige mechanische Dehnung über die gesamte Zelle. Das arithmetische
Mittel der zehn Knotentemperaturen ist hierfür kein Ersatz. Die Hashgrundlage
der Provenienz benennt den berichteten Tet10-Eigenverzerrungswert als
`thermal_eigenstrain_at_centroid` und enthält außerdem das vollständige
Temperaturfeld, Expansionsparameter und den assemblierten Lastvektor.
Tet4 behält seinen bestehenden Zellmittelungs-/Energievertrag unverändert.

## Sequenzieller Ablauf, Backend und Export

Die Kopplung überträgt alle aktiven Knotentemperaturen ohne Projektion auf
demselben Netz. Die optionalen Keyword-Argumente `thermal_backend` und
`mechanical_backend` erlauben getrennte lineare Backends; `None` bleibt SciPy.
Beide Teilresultate behalten ihre eigene Provenienz, Residuen und Bilanzgrößen.

Das gemeinsame Ergebnisbündel verwendet Schema 2 und enthält Wärme- und
Mechanikfelder sowie `cell_integral`-Energien. Vor dem Schreiben werden die
Schwerpunktkoordinaten **beider** Teilresultate gegen das Netz geprüft.
Felder bleiben rohe Schwerpunktwerte; keine Knotenprojektion. Tet4 behält
Schema 1. Das [ausführbare Gmsh-Beispiel](../examples/tet10_thermomechanics/README.md)
verwendet die Python-API: Ein eigenständiger `thermomechanics`-TOML-Modus
existiert weiterhin nicht. Tet10-Thermopiezo-Modelle bleiben ausdrücklich gesperrt.

## Verifikation

Lokal: `PYTHONPATH=src .venv/bin/python -m pytest -q` mit 543 bestandenen und
51 optional übersprungenen Tests. Ruff, Formatprüfung, Mypy und
`git diff --check` bestanden. Native PETSc-/Artefaktprüfungen laufen in CI.

TDD: drei Integrationstests scheiterten zuerst an der bisher ausdrücklichen
Tet10-Thermalsperre. Nach punktabhängiger Last-/Energieintegration und
gekoppeltem Export sind sie grün. Das ausführbare Beispiel erhielt vor seiner
Implementierung einen wegen fehlender Datei fehlschlagenden Test.

- Vollständig behinderte quadratische Temperatur T=300+x² K, E=100 Pa, ν=0,25,
  αxx=0,01/K: σ=-(1,2;0,4;0,4;0;0;0)x² Pa und W=0,0012 J aus Integral x⁴=1/5.
  Dies erkennt eine falsche Mittelung oder reine Schwerpunktsenergie.
- Freie Expansion bei T=300+x K: u=(0,005x²,0,0) m, Spannung und Energie null.
- Kopplung T=1+x², u=(0,01x²,0,0), E=100 Pa, ν=0 und αxx=0,01/K:
  σxx=2x-x² Pa und W=0,008/3 J; gemischte thermische/mechanische Ränder,
  Wärme-, Kraft-/Momentenbilanz und unabhängiges VTU-Einlesen.
- Rotierte anisotrope Expansion und regionsweise Tref=300/310 K bei T=325 K;
  unabhängige behinderte Spannungen/Energie. Fehlende Expansionsdaten und
  negative Knotentemperaturen werden geprüft.
- Verschieben von T und Tref um 300 K erhält quadratische Spannungen und
  Energie; getrennte Backend-Instanzen erhalten jeweils nur ihr eigenes Feld.
- Native PETSc-Parität für beide Teilfelder, direkt/FieldSplit und real/complex,
  läuft in CI; lokal ohne PETSc sind diese Tests übersprungen.
- Artefakt-CI führt das Beispiel auf Python 3.11–3.13 aus Wheel/sdist aus.

Toleranzen der kleinen synthetischen Fälle: Temperatur 4e-12 K,
Verschiebung 2e-14 m direkt bzw. 1e-13 m gekoppelt, Spannung 2e-11 Pa
(behinderte analytische Fälle 2e-13 Pa), Energie 2e-15 J direkt bzw. 2e-14 J
gekoppelt. Freie Expansionsenergie <1e-24 J. Bilanzen 2e-11 W und 2e-12 N/N m.
Diese Grenzen sind Rundungs-/Assemblierungstoleranzen, keine allgemeine
Diskretisierungsfehlergarantie. Bei der ersten gemischten Kopplungsprüfung mit
300-K-Offset überschritt ein Temperaturwert die vorgegebene 4e-12-K-Grenze
knapp (4,09e-12 K). Die Grenze wurde nicht aufgeweitet: Der Kopplungsfall nutzt
nun Tref=1 K zur Trennung der Feldintegration vom großen Temperaturoffset;
die direkten 300-K-Fälle und der Offset-Invarianztest bleiben separate Nachweise.
Die synthetischen Materialdaten behaupten keine Tieftemperaturvalidierung.

## Vier-Netze-Manufactured-Solution

Setze P=x⁴+y⁴/2+z⁴/4, T=1+P K und u=(0,01P,0,0) m. Für k=2 W/(m K),
E=100 Pa, ν=0,25 (λ=μ=40 Pa), αxx=0,002/K und Tref=1 K folgen:

```text
Q = -24 (x²+y²/2+z²/4)
b = (-0,01(1440x²+240y²+120z²)+0,96x³, 0,16y³, 0,08z³)
```

Auf allen Randknoten werden die exakten T/u-Werte vorgegeben. Unabhängig
formulierte P2-Basis/Ableitungen und eine Duffy-Produktregel mit sechs
Gauss-Punkten pro Richtung integrieren die Fehlernormen. Die unabhängige
Tensorenergie `lambda/2 trace(epsilon_m)^2 + mu epsilon_m:epsilon_m`
prüft zusätzlich jedes Zellintegral (`rtol=2e-12`, `atol=1e-16 J`).

| Würfelteilung n | L2(T) | H1(T) | L2(u) | H1(u) | L2(Spannungstensor) |
|---|---:|---:|---:|---:|---:|
| 1 | 8,03016e-2 | 5,20145e-1 | 8,0179e-4 | 5,20429e-3 | 6,19225e-1 |
| 2 | 1,11427e-2 | 1,43245e-1 | 1,12009e-4 | 1,43559e-3 | 1,70130e-1 |
| 3 | 3,34559e-3 | 6,47629e-2 | 3,35811e-5 | 6,48886e-4 | 7,69473e-2 |
| 4 | 1,41754e-3 | 3,66569e-2 | 1,42120e-5 | 3,67113e-4 | 4,35752e-2 |

Letzte Ordnungen: 2,985 / 1,978 / 2,989 / 1,980 / 1,977. Die letzten beiden
Intervalle müssen für L2(T/u) im Fenster 2,7–3,3 und für H1(T/u)/Spannung
im Fenster 1,8–2,2 liegen. Es werden vollständige Gradienten und Spannungstensoren
an Integrationspunkten verglichen; Schwerpunktwerte allein reichen nicht aus.
