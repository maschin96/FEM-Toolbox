# Experimenteller CFD-Einstieg (#84/#85/#86)

Dieser zweite Experimental-Branch beginnt unabhängig vom SPICE-Branch auf
stabilem main. M12 bleibt RANS/FSI; der laminare Kanal ist ausschließlich ein
Verifikationsschritt des viskosen Operators, kein laminares Produktrelease.
Die Gesamtakzeptanz von #84–#86 bleibt offen.

## Verbindlicher Vertrag dieses ersten Schnitts

Voll entwickelte stationäre inkompressible isotherme Strömung zwischen zwei
unendlichen parallelen Wänden: u=(u(y),0,0), konstante dynamische Viskosität
mu>0 Pa s, vorgegebener konstanter Druckgradient dp/dx in Pa/m. Kein
Druckfeld wird gelöst; dessen additive Gauge ist für diesen Operator irrelevant.
G=-dp/dx und -(mu u')'=G. Wandgeschwindigkeiten sind Dirichletwerte in m/s;
positive Geschwindigkeit und Druckkraft zeigen in +x. Kontinuität ist durch
x-Unabhängigkeit und verschwindende Quergeschwindigkeit identisch erfüllt.
Damit wird keine allgemeine Druck-Geschwindigkeits-Kopplung verifiziert.

Zellzentrierte konservative Finite Volumen auf streng steigenden ebenen
Flächenkoordinaten y in m. Ein gemeinsamer Schubfluss tau=mu delta(u)/delta(y)
pro Innenfläche garantiert diskrete Impulserhaltung. Wandabstände sind echte
halbe Randzellbreiten, keine Ghost-Zell-Erfindung. Sparse symmetrisch positiv
bestimmter Diffusionsoperator; keine Übernahme der Festkörper-Tet4-Operatoren.

Pro Einheitsbreite und -länge: Impulsbilanz tau_oben-tau_unten+G H=0 (Pa).
Diskrete Leistungsbilanz G sum(u_i h_i)+tau_oben U_oben-tau_unten U_unten
minus sum(mu delta(u)^2/delta(y))=0 (W/m²). Dissipation enthält die
Randsegmente. Ergebnisse besitzen unabhängige rekonstruierte Flüsse,
Zellresiduen und unveränderliche Arrays. Keine Dichte nötig, weil kein
Konvektions-/Massenterm ausgewertet wird.

## Referenzen und Gate

Aus zweimaliger Integration: Poiseuille u=G y(H-y)/(2mu) bei ruhenden Wänden;
Couette u=U_unten+(U_oben-U_unten)y/H bei G=0. Uniforme Raster 8/16/32/64:
maximaler Geschwindigkeitsfehler muss Ordnung 2 (1.99–2.01) zeigen.
Mu=2, G=4, H=1: Randtau=(2,-2) Pa, Q=1/6 m²/s. Bei 64 Zellen:
Geschwindigkeitsfehler <0.00025 m/s, Q-Fehler <0.0002 m²/s.
Couette auf gestrecktem Raster ist affine-exakt. Bilanz-/Schubtoleranz
1e-12 in den angegebenen SI-Einheiten deckt float64-Rundung kleiner Systeme
ab; sie ist kein allgemeines Genauigkeitsbudget. Ungültige Netze,
Parameter, Wertebereiche und nicht endliche Lösungen werden abgelehnt.

## Nächste Entscheidungsgates

Vor allgemeinen CFD-Netzen müssen #84 und #85 vollständig bearbeitet werden:
Netzadapter, Konvektionsschema, Druckkopplung, Druck-Gauge, nichtlineare
Iteration und offene Randbedingungen. Als RANS-Kandidat wird Spalart-Allmaras
bewertet; Variantenwahl, Schutzmaßnahmen und Wandauflösung sind vor #87
festzulegen. Primärquelle: [Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/spalart.html).
Das ist noch keine abgeschlossene Modellentscheidung. FSI-Vertrag,
konservative Übertragung und Referenzdaten bleiben ebenfalls offen.

Keine Turbulenztransportgleichung, Wandfunktion, dreidimensionaler CFD-Solver,
FSI oder Validierungsbehauptung in diesem Schnitt. Vorhandenes NumPy/SciPy;
keine neue Abhängigkeit oder externer Solveradapter.

## Ausgeführter Nachweis (01.10.2026)

TDD: 14 Tests scheiterten zuerst an der fehlenden API. Nach Implementierung
bestanden 13; das 32-Zellen-Volumenstromergebnis überschritt das vorgegebene
Fehlerbudget. Verfeinerung auf 64 Zellen erfüllt das unveränderte Budget;
die zusätzliche Vier-Netz-Prüfung weist zweite Ordnung nach. Danach alle
14 gezielten Tests grün. Fehler beim Beispiel: 6.10352e-05 m/s;
Impulsrest 8.88178e-16 Pa, Leistungsrest 2.22045e-16 W/m².

```sh
PYTHONPATH=src .venv/bin/python -m pytest -q
PYTHONPATH=src .venv/bin/python examples/cfd_channel/run.py
.venv/bin/python -m ruff check .
.venv/bin/python -m ruff format --check .
.venv/bin/python -m mypy src/fem_toolbox
git diff --check
```

Gesamtsuite: 710 bestanden, 83 übersprungen (93.68 s). Ruff, Formatprüfung,
mypy und Beispiel bestanden. Übersprungene Prüfungen gelten nicht als
bestanden. Keine unabhängige Reviewfreigabe und kein CFD-Release.

## Erweiterung: Druckkopplung

Der historische reine Kanal-Schnitt wurde um einen geschlossenen 2D-Stokes-
Solver mit echtem Druckfeld ergänzt. [Vertrag und Nachweise](69_CFD_PRESSURE_COUPLING.md).
Die oben offenen allgemeinen CFD-/RANS-/FSI-Gates bleiben offen.

## Vollständiger stationärer Referenzkern

Die historischen Grenzen des ersten Kanal-Schnitts sind im neuen
[Vertrag #86](70_STATIONARY_FLOW_CORE.md) um Konvektion, nichtlineare Iteration
und offene Kanalgrenzen ergänzt. Die allgemeine RANS-/FSI-Planung bleibt M12.
