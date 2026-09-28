# Tet10-Newmark-Zeitintegration (#35d)

`solve_transient_piezoelectric` und `Model`-Schema 3 mit
`physics="transient_piezoelectric"` unterstützen geradkantige Tet10.
Der bestehende Newmark-Vertrag bleibt unverändert: β=1/4, γ=1/2, konsistente
mechanische Masse, algebraisches Potential, keine physikalische oder numerische
Dämpfung. Dichte und Material können zellweise vorgegeben werden.

Der gemeinsame Dynamikaufbau hält den für Tet10 symmetrisierten dielektrischen
Block vor. Modal- und Newmark-Solve verwenden damit denselben elektrischen
Operator; Quadratur-Rundungsasymmetrien an theoretischen Nulleinträgen gelangen
nicht in den SPD-Solve. Ursprüngliche Matrizen bleiben für Residuen und Energie
erhalten. Backend-Toleranzen werden nicht gelockert. Tet4 behält seinen Operator.

## Randwerte, Lasten und Ausgabe

Anfangsverschiebung und -geschwindigkeit haben Form `(nactive,3)`; Potential und
Beschleunigung werden konsistent berechnet. Bewegte Ränder liefern u, v und a
mit `KinematicValues`; die Reduktion berücksichtigt `-M_fc a_c`.
Callbacks liefern bereits integrierte Kräfte und elektrische rechte Seiten.
Model integriert seine Volumen-/Flächenlasten mit den vorhandenen Tet10-/
Triangle6-Operatoren. Negative P2-Eckgewichte sind korrekt: ∫N_corner=-V/20,
∫N_edge=V/5. Eine Ersetzung durch gleiche positive Knotengewichte ist falsch.

Eine schwebende Elektrode erfasst alle sechs Knoten jeder Triangle6-Fläche.
Leiterladung Q entspricht `-∫D·n dA` mit der Normalen aus dem Festkörper;
der elektrische Lastvektor enthält -Q. Beiträge an verbundenen Knoten summieren
sich. Die elektrische Bilanz wird für die gesamte Elektrode ausgewertet.

Ereignisse werden eingefügt, Intervalle linksseitig integriert und Potential/
Beschleunigung rechtsstetig gespeichert. u/v bleiben kontinuierlich. Für direkte
Callbacks müssen Sprungzeiten über `events` angegeben werden; Model sammelt die
Profilknoten. Details: [M6-Zeitvertrag](34_PIEZO_DYNAMICS.md).

Zeitserien verwenden den bestehenden Tet10-Export mit Manifest-Schema 2,
rohen Schwerpunktfeldern und P2-Interpolation von u/v/a/φ an Messpunkten.
VTU/PVD und CSV enthalten den tatsächlichen Zeitraster einschließlich Ereignissen.

## TDD und Verifikation

Die ersten vier Tests scheiterten an den direkten bzw. Model-Tet10-Sperren.
Nach Freischaltung scheiterten zwei Mehrzellen-Elektrodentests an der strengen
SPD-Symmetrieprüfung des separat aufgebauten elektrischen Newmark-Blocks.
Die gemeinsame Operatorhaltung behebt diesen Fehler; alle gezielten Tests
bestehen. Der Gmsh/TOML-Beispieltest wurde vor Anlegen des Skripts rot ausgeführt.

- Eckoszillator des Einheits-Rechtstetraeders mit ρ=6, Cxx=10, exxx=2, κxx=3:
  m=6/420, k=(10+4/3)/10 bei freiem Eckpotential, k=1 bei φ=0.
  q(0)=0,01 m, v(0)=0,02 m/s, ω=√(k/m),
  q(t)=0,01 cos(ωt)+(0,02/ω)sin(ωt). Energie bleibt
  0,5 k q(0)²+0,5 m v(0)² (rtol=3e-12). Potential und Beschleunigung
  erfüllen φ=2q/3 bzw. a=-ω²q.
- Bewegte starre Translation uₓ=0,003t² mit unabhängigen P2-Volumenlastgewichten
  und bₓ=ρ·0,006 prüft die nichtdiagonale Masse und vorgeschriebene Ableitungen.
  Fehlergrenzen: u/v/φ 3e-15 SI, a und Kraftresiduen 3e-14 SI.
- Kraft- und Spannungssprünge erzeugen keinen vorzeitigen Impuls. Model prüft
  insbesondere die negative Eckkraft -0,05 N einer konstanten Volumenlast.
- Der d33-Stab aus [Doc 51](51_TET10_MODES.md) erhält offen/kurzgeschlossen eine
  einzelne Anfangsmode. Die diskrete Newmark-Phase ist 2 atan(ω Δt/2) je Schritt.
  Modenform, Energie und Äquipotentialbedingungen werden mit rtol=2e-9 und
  mechanischem atol=1e-14 m geprüft; offene Restladung <1e-18 C.
- Unabhängiger statischer Stabzustand u_z=1e-4 z, φ=z: σzz=100005 Pa,
  D_z=0,00049999 C/m². Obere Traktion und Leiterladung Q=-D_z·1e-4 m² halten
  diesen Zustand auf ungleichmäßigen Zeitschritten konstant. Triangle6-Traktion
  hat keine Eckkräfte, sondern A/3 je Kantenknoten. Toleranzen: u 2e-16 m,
  φ 2e-11 V, v 2e-10 m/s, summierte elektrische Residuen 1e-18 C.
- Model-Rundreise, Schwerpunktfelder, quadratische u/v/a-Messpunktwerte,
  unzulässige Anfangszustände, bewegte Randsprünge und Speichergrenze.
- Native PETSc-Parität direkt/FieldSplit in real/complex: freie Elektrodenmoden
  (rtol=2e-8, u atol=1e-14 m, φ atol=1e-9 V) und geladener Stab
  (u atol=2e-14 m, φ atol=2e-9 V). Feldresiduen jeweils <1e-9.

### Zeitkonvergenz auf unverändertem Raumnetz

Fehler von q(0,4 s) gegenüber der analytischen Schwingung:

| Schritte | freies Eckpotential | φ=0 |
|---:|---:|---:|
| 40 | 4,76947e-6 m | 6,12867e-7 m |
| 80 | 1,19865e-6 m | 1,49671e-7 m |
| 160 | 3,00055e-7 m | 3,71958e-8 m |

Letzte Ordnungen: 1,998 und 2,009. Raumverfeinerung der unveränderten
Masse-/Steifigkeitsoperatoren ist separat über vier Stabnetze in Doc 51 geprüft.
Die Einmoden-Zeitprüfung verwendet eine diskrete Eigenform und ist kein neuer
unabhängiger Raumkonvergenznachweis. Synthetische Referenzen sind Code-Verifikation,
keine experimentelle Validierung. Absolute Toleranzen beziehen sich auf die
angegebenen SI-Skalen und decken Rundung der jeweiligen kleinen Systeme ab.

## Beispiel und Grenzen

```sh
PYTHONPATH=src python examples/tet10_newmark/run.py
PYTHONPATH=src python -m pytest tests/unit/test_tet10_newmark.py
```

Das eigenständige Gmsh/TOML-Beispiel prüft freie Schwingung und Energie und schreibt
41 Zeitschritte mit einem inneren Messpunkt. CI führt es aus Wheel und sdist aus.
Die Grenzen bleiben 2048 aktive Gesamt-DOFs und 2.000.000 skalare Zustandswerte.
Keine gekrümmten Tet10, kein Streaming, keine transiente thermische Kopplung.
Gesamt-Abnahme-/Kostenvergleiche für #35 bleiben offen.

## Tet10-Harmonik (#35d)

Komplexe Frequenzantworten einschließlich DC, Elektroden, Model und P2-Export
sind integriert. [Nachweise und Grenzen](53_TET10_HARMONIC.md).
Gesamt-Abnahme und Genauigkeits-/Kostenvergleiche bleiben offen.
