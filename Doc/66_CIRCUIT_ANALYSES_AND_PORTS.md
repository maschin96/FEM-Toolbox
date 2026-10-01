# Schaltungskern: Analysen, Ports und Kopplungsentscheidung

Experimenteller M13-Vertrag zu #114–#119. Ergänzt den
[DC-Vertrag](65_EXPERIMENTAL_CIRCUIT_CORE.md). Kein M13-Release und keine
experimentelle Validierungsbehauptung. Portimplementierung und gekoppelte
Solver bleiben #120–#128; hier werden ihre Schnittstellen vorab festgelegt.

## Modell und benannte Ports

Zusätzlich zur begrenzten Textnetzliste enthält TOML-Schema 1 explizite
SI-Komponenteneinheiten und benannte Portzuordnungen. `ElectricalPort`
verbindet einen vorhandenen Schaltungsknoten mit genau einem benannten
FEM-Elektrodenleiter. Ports liefern Potential V und Leiterladung Q in C.
Der Strom **in die FEM-Elektrode** ist `i_port=dQ/dt`; er ist ein abfließender
Strom im KCL des Schaltungsknotens. Für einen zweipoligen isolierten Piezo
sind zwei Elektrodenports nötig, deren Ströme sich zu null summieren.
Leiterladung bleibt `Q=-sum(r_phi)` wie in Doc/20; kein Wechsel der
FEM-Residualkonvention. Im Frequenzbereich gilt `i_port=i*omega*Q_hat`.

`ThermalPort` benennt eine FEM-Kontakt-/Bauteilregion und die zugeordneten
Widerstände. Positive `P=sum(v_R*i_R)` fließt **in** die FEM-Wärmedomäne.
Porttemperatur ist absolute T in K; bei späteren thermischen Netzwerken
ist sie ein explizit gewichteter Schnittstellenwert. Der eigenständige Kern
validiert Namen und eindeutige Zuordnung; Existenz, Maße und Gewichte der
FEM-Region/Elektrode prüft erst der gekoppelte Adapter. Derselbe Widerstand
darf nicht mehrfach vollständig als Wärmequelle eingetragen werden.
Quellenleistung und reaktive Energie sind keine impliziten Wärmequellen.

TOML: `circuit_schema_version=1`, `[[components]]` mit `name`, `positive`,
`negative`, `value`, `unit`; optionale `[[electrical_ports]]` mit `name`,
`node`, `electrode`, und `[[thermal_ports]]` mit `name`, `region`,
`components`. Einheiten sind exakt `Ohm`, `F`, `H`, `V`, `A`.
Unbekannte Schlüssel/Versionen, falsche Einheiten und ungültige
Zuordnungen sind Fehler. Textnetzlisten haben keine Portsyntax;
porthaltige Modelle werden nicht verlustbehaftet als Netzliste exportiert.
Modell und Provenienz besitzen alle numerischen und Port-Eingaben.

## Transiente Referenz

MNA-Zustand `x=[v_nodes, i_V, i_L]`, mit E-Kapazitätsblock und negativer
L-Diagonale im Induktorstromblock:

```text
E xdot + A x = b(t)
(A + E/dt) x_new = b(t_new,left) + E x_old/dt
```

Backward Euler ist deterministisch und erste Ordnung. Quellprofile sind
explizite V-/I-Zuordnungen auf vorhandene Quellen: SI-Skalar, `Step`, `Ramp`
oder `TimeSeries`; keine Callables/Ausdrücke. Ereignisknoten werden in die
Zeitliste eingefügt. Zustände am Ereigniszeitpunkt verwenden den linken
Quellwert, der nächste Schritt den neuen Wert. Anfangszeit verwendet ebenfalls
den linken Wert. Keine versteckten Schritte vor dem Anfangszeitpunkt.
Timeseries-Einheiten müssen V/A entsprechen.

Anfangswerte: benannte Kondensatorspannungen in V und Induktorströme in A,
fehlende Einträge ausdrücklich null. Die konsistente Anfangslösung ersetzt
C durch ideale Spannungsquellen und L durch ideale Stromquellen und berechnet
die algebraischen Größen. Unzulässige C-/V-Schleifen und L-/I-Schnittmengen
werden dabei als Initialisierungsfehler zurückgewiesen. Der erste Umfang
behandelt reguläre Index-1-Netze; keine Ableitung idealer Sprünge, Impulsströme
oder automatische Reparatur widersprüchlicher Speicherzustände. Eine direkt
parallel zur idealen Spannungsquelle liegende C-Zustandsvorgabe ist auch dann
außerhalb, wenn die Momentanwerte übereinstimmen (unbestimmter Anfangsstrom).

Ausgabe: vollständige Zeitliste, Knoten-/Zweigverläufe, KCL, Q_C=C*v_C,
Psi_L=L*i_L, Speicherenergie `E_s=sum(C*v_C²+L*i_L²)/2` und pro Schritt:

```text
delta(E_s) - source_work + resistor_loss + numerical_dissipation = 0
source_work = -dt * sum(p_V+p_I)
resistor_loss = dt * sum(p_R)
numerical_dissipation = sum(C*delta(v_C)^2+L*delta(i_L)^2)/2
```

Physikalische Dissipation und BE-Dämpfung bleiben getrennt. KCL, Spannung/
Induktorzustand, Ladungsfortschreibung und Energie sind eigene Prüfgates.
Statische Matrizen, explizite Zeitschritte, History-Grenzen und konkrete
Schritt-/Zeitdiagnosen; keine adaptive Schrittsteuerung oder nichtlineare
Iteration im Referenzkern. Zu grobe Schritte können trotz konsistenter Bilanz
ungenau sein; drei Schrittweiten und beobachtete Ordnung prüfen das separat.

## Harmonik

Peak-Phasoren: `x(t)=Re(x_hat exp(+i*omega*t))`, `omega=2*pi*f`, f in Hz.
`(A+i*omega*E)x_hat=b_hat`. R: 1/R, C: i*omega*C, L: v=i*omega*L*i_L.
Quellenamplituden können explizite komplexe Peak-Werte sein; ein nicht
überschriebener V-/I-Wert wird als reelle Peak-Amplitude interpretiert.
`S=0.5*v_hat*conj(i_hat)` in VA, P=Re(S) in W, Q=Im(S) in var.
Induktive Blindleistung ist positiv, kapazitive negativ. Gemittelte
Speicherenergie ist `sum(C*abs(v_hat)^2+L*abs(i_hat)^2)/4` in J.

Bei f=0 gilt die DC-Lösung; Leistung wird ohne Faktor 1/2 berechnet und
Speicherenergie mit Faktor 1/2. Imaginary Quellenwerte bei DC werden abgelehnt.
Frequenzen müssen endlich und nichtnegativ sein. Floating-Netze, redundante
V-Schleifen und singuläre oder schlecht konditionierte verlustlose Resonanzen
werden vor Ergebnisausgabe diagnostiziert, auch bei Nullanregung. Keine
künstliche Dämpfung. Grenz-/Resonanzdiagnosen und Skalierung sind explizit;
nahe, aber reguläre gedämpfte Resonanz darf berechnet werden.

## Künftige gekoppelte Strategie und Abnahmefälle (#114)

Elektrisch wird nach Elektrodenkondensation **monolithisch** gekoppelt:
gemeinsame Knotenpotentiale und KCL mit Q-FEM-Port ersetzen getrennte
Dirichlet-/Floating-Vorgaben an angeschlossenen Elektroden. Im DC-Fall ist
der Kondensatorstrom null; offene Anfangsladung beziehungsweise isolierte
Leiterpotentiale erfordern explizite Ladungs-/Gauge-Bedingungen. Für
Transienten werden die M6-Mechanikzustände und BE-Schaltungsschritte mit
gemeinsamen Ereignisknoten gekoppelt; unterschiedliche Integratoren verlangen
einen ausgewiesenen diskreten Energienachweis in #122. Kein stillschweigendes
Weiterreichen von Spannungen nach einem getrennten Solve.

Wärmeverlustübertragung startet **einseitig**: räumliche Gewichte w mit
`integral(w)dOmega=1` beziehungsweise Kontaktflächen-Äquivalent verteilen
P_R konservativ. Wärmeausgang/Netzwerkport und FEM tauschen dieselbe Leistung
mit entgegengesetztem Vorzeichen. Temperatur-Rückwirkung in #126 erfolgt
iterativ mit vorab festgelegten Relaxations-/Residuen-/Maximaliterationsgrenzen;
Nichtkonvergenz stoppt den Lauf. Grenzen und Jacobian-/Monotonieprüfung werden
vor #126 festgelegt, keine freie nichtlineare Bauteilbibliothek.

Unabhängige Referenzfälle vor der späteren Portimplementierung:

- Piezo: Doc/20-Quader, kurzgeschlossen 0 V sowie Leerlauf 0 C; externe
  Kapazität ergänzt C_ext*V in die Gesamtladungsbedingung. Bei fest vorgegebener
  mechanischer Dehnung aus Doc/20 gelten C_FEM=2 F und Q_mech=-0.05 C,
  damit `V=0.05/(2+C_ext)` V im ladungsneutralen Fall. Widerstandsentladung
  folgt `tau=R*(2+C_ext)` bei festgehaltener Dehnung. Verglichen werden V,
  Q, I, mechanische Reaktion und alle Energieanteile; #123 legt feinere
  gekoppelte Zeit-/Netzbudgets fest.
- Kühlkörper: 10 V an 100 Ohm ergeben 1 W Joule-Leistung. Homogener Stab,
  k=10 W/(m K), A=0.01 m², L=0.1 m, T_base=293.15 K und gleichmäßig
  zugeführtem Stirnflächenstrom: T_hot=294.15 K, Base-Abfluss=1 W.
  Verlusttransfer, Quellennormierung und FEM-Netzkonvergenz getrennt prüfen.

Beide sind analytische Code-Verifikation, keine Experimente. Norm-/Vorzeichen-
Toleranzen für lineare kleine Referenzen liegen oberhalb float64-Rundung und
unterhalb des Diskretisierungsbudgets; vollständige gekoppelte Abnahmen
stehen in #123/#127. Kein externer Schaltungssimulator als Laufzeitabhängigkeit.
