# Stationäre Tet10-Thermopiezoelektrizität (#35d)

Der bestehende einseitige Ablauf T → (u, φ) unterstützt geradkantige Tet10
über die Python-API und `Model`-Schema 1. Temperatur wird auf demselben Netz
an allen aktiven Knoten übertragen; Materialregionen, Referenztemperaturen,
anisotrope Expansion und feste/schwebende Triangle6-Elektroden bleiben nutzbar.
Keine thermische Rückwirkung, Pyroelektrizität oder transiente Kopplung.
Tet10-Harmonik und die Gesamt-Abnahme-/Kostenvergleiche von #35 bleiben offen.

## Thermischer Integrationsvertrag

Für jedes Element und jeden Integrationspunkt gilt mit θ=N·(T-Tref):

```text
epsilon_th = alpha theta
epsilon_m = B u - epsilon_th
sigma = C epsilon_m - e^T E
D = e epsilon_c + kappa E
```

Bei `electric_thermal_strain=False` (Standard) ist εc=Bu, bei `True` ist
εc=εm. Dieser vorhandene Schalter wählt ausdrücklich, ob D den thermischen
Eigenverzerrungsterm enthält; er aktiviert keine Pyroelektrizität.
Die quadratische Lastintegration verwendet ∫BᵀCεth und bei aktivem Schalter
zusätzlich ∫Gᵀeεth. Tet4 behält seinen bisherigen Zellmittelvertrag.

Spannung und D ziehen die am Schwerpunkt interpolierte Eigenverzerrung ab.
Die drei Enthalpiebeiträge werden separat integriert:

```text
Wm = integral 1/2 epsilon_m^T C epsilon_m dV
Wc = integral -epsilon_c^T e^T E dV
He = integral -1/2 E^T kappa E dV
```

Quadratische Knotentemperatur macht εth quadratisch; die Grad-5-Regel integriert
die maximal quartischen Energieterme exakt. Ein Temperaturmittel oder eine
Schwerpunktsenergie wäre hierfür unzureichend. `strain` bleibt Gesamtdehnung.
Die Hashgrundlage der Tet10-Provenienz enthält alle Knotentemperaturen,
Expansionsdaten und Referenztemperaturen sowie den gekennzeichneten
Eigenverzerrungswert am Schwerpunkt. Tet4-Provenienz bleibt unverändert.

`solve_thermopiezoelectric` erhält die optionalen Keywords `thermal_backend`
und `piezoelectric_backend`; `None` wählt jeweils SciPy. Die Teillösungen behalten
getrennte Residuen, Bilanzen und Provenienzen. Der gekoppelte Export verwendet
Schema 2 mit Schwerpunktfeldern und Zellintegralen; die Auswertekoordinaten
beider Teilergebnisse müssen zum Netz passen. Model benötigt keine neue Schemaversion.

## TDD und analytische Nachweise

Vor der Umsetzung scheiterten zwei Tests an der Tet10-Thermalsperre und zwei
an den fehlenden Backend-Keywords. Der Model-Test scheiterte an der bisherigen
Tet10-Thermopiezosperre, der Beispieltest am fehlenden Skript. Danach bestanden
die gezielten Tests. Tet10-Harmonik bleibt gesperrt; Modal-/Newmark-Folgearbeit siehe Doc 51/52.

Synthetisches Material für die kleinen Würfelfälle: E=100 Pa, ν=0, κ=3 I F/m,
eₓₓₓ=2 C/m² und αxx=0,01/K. Beide elektrischen Schalterstellungen werden geprüft.

- Behinderte quadratische Temperatur T=300+x² K bei φ=x² V und u=0:
  σxx=4x-x² Pa, Dₓ=-6x-s·0,02x² C/m²; Wm=0,001 J, Wc=-s·0,01 J, He=-2 J,
  wobei s=0/1 der elektrische Schalter ist.
- Sequenzieller gemischter Randfall T=1+x² K, uₓ=0,01x² m, φ=x² V:
  σxx=6x-x², Dₓ=-5,96x-s·0,02x². Quellen Q=-4 W/m³, bₓ=-6+2x N/m³ und
  ρ=-5,96-s·0,04x C/m³; passende rechte Flächenlasten und linke Dirichletwerte.
  Wm=0,008/3 J, Wc=0,08/3-s·0,01 J, He=-2 J.
- Mechanisch behinderter Model-Fall T=1+x mit schwebender rechter Elektrode Q=0:
  φ=-s·0,01x²/3 V, D=0; bei s=1 ist σxx=-76x/75 Pa. Model-Rundreise und Export.
- Regionsweise Tref, rotierte anisotrope Expansion, gemeinsame Verschiebung
  von T/Tref um 100 K, ungültige Temperatur-/Materialdaten und fehlende Paare.
- Freie Expansion bei T=300+x: uₓ=0,005x², σ=E=0 für s=1 und Wm<1e-24 J.
- Getrennte Backend-Instanzen erhalten nur T bzw. u/φ. Native PETSc-Parität
  für beide Schalterstellungen läuft mit direkt/FieldSplit und real/complex in CI.

Absolute Toleranzen der kleinen SI-Fälle: T 4e-12 K, u/φ 2e-13 bis 2e-12 m/V,
Spannung und Bilanzen 2e-11 in der jeweiligen Einheit, Energien bis 2e-12 J.
Diese Grenzen decken Rundung/Assembly ab; keine allgemeine Fehlergarantie.
PETSc-Parität: rtol=2e-8, atol=2e-10 für Piezo-Felder; Temperatur/Bilanzen bis
2e-9 wegen iterativer Lösung. Die 1-K-Referenz der synthetischen gekoppelten
Fälle trennt Temperaturgradienten von großen Offsets; keine Tieftemperaturvalidierung.

## Vier-Netze-Nachweis

P=x⁴+y⁴/2+z⁴/4, T=1+P, u=(0,01P,0,0), φ=P, αxx=0,002/K, Tref=1 K,
E=100 Pa, ν=0,25, eₓₓₓ=2 C/m², κ=3 I F/m und k=2 W/(m K), s=1.
Exakte Außenrandwerte und unabhängig hergeleitete Quellen:

```text
Q = -24 (x²+y²/2+z²/4)
b = (-38,4x²-2,4y²-1,2z²+0,96x³, 0,16y³, 0,08z³)
rho = -35,76x²-18y²-9z²-0,016x³
```

Unabhängige P2-Basis und Duffy-Quadratur (6 Punkte je Richtung) messen vollständige
L2-/H1-Fehler. Alle Zellenthalpien werden gegen unabhängige Tensorformeln geprüft
(rtol=2e-12, atol=1e-16 J), ebenso die vier getrennten Bilanzen.

| Teilung n | L2(u) | H1(u) | L2(φ) | H1(φ) | L2(T) | H1(T) |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 8,0122e-4 | 5,21255e-3 | 8,03015e-2 | 5,20145e-1 | 8,03016e-2 | 5,20145e-1 |
| 2 | 1,14912e-4 | 1,46402e-3 | 1,11431e-2 | 1,43245e-1 | 1,11427e-2 | 1,43245e-1 |
| 3 | 3,43718e-5 | 6,60672e-4 | 3,34566e-3 | 6,47629e-2 | 3,34559e-3 | 6,47629e-2 |
| 4 | 1,44673e-5 | 3,72249e-4 | 1,41756e-3 | 3,66569e-2 | 1,41754e-3 | 3,66569e-2 |

Letzte Ordnungen: 3,008/1,994 (u), 2,985/1,978 (φ und T). Die letzten beiden
Intervalle müssen für L2 im Fenster 2,7–3,3 und für H1 in 1,8–2,2 liegen.
Das [Gmsh/TOML-Beispiel](../examples/tet10_thermopiezoelectric/README.md) wird
zusätzlich aus Wheel/sdist unter Python 3.11–3.13 ausgeführt.

## Lokale Abschlussprüfungen

`PYTHONPATH=src .venv/bin/python -m pytest -q`: 565 bestanden, 57 optional
übersprungen. Ruff, Formatprüfung, Mypy und `git diff --check` bestanden.
Native PETSc- und Artefaktprüfungen erfolgen in CI; lokale Skips sind kein
PETSc-Nachweis.

## Tet10-Modalanalyse (#35d)

Konsistente Masse, offene/kurzgeschlossene Elektroden, Model und Modenexport
sind integriert. [Umfang und Verifikation](51_TET10_MODES.md).
Tet10-Harmonik und Gesamt-Abnahme-/Kostenvergleiche bleiben offen.

## Tet10-Newmark (#35d)

Zeitintegration mit konsistenter Masse, bewegten Rändern, Last-/Spannungsereignissen,
Model und Zeitserienexport ist integriert. [Nachweise und Grenzen](52_TET10_NEWMARK.md).
Tet10-Harmonik und Gesamt-Abnahme-/Kostenvergleiche bleiben offen.
