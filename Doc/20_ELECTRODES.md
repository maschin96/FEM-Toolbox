# Equipotential- und Floating-Electrode-Randbedingungen (#24)

Status: Für den quasistatischen piezoelektrischen Tet4-Referenzsolver
implementiert und mit einem unabhängigen 1D-Grenzfall verifiziert.

## API und Vorzeichen

`Electrode(name, facets, potential=V)` setzt alle Knoten einer leitenden
Elektrode auf dasselbe Potential. `Electrode(name, facets, charge=Q)` lässt
dieses Potential unbekannt und gibt die gesamte **Leiterladung** in Coulomb
vor. `facets` ist ein Name aus `mesh.facet_sets` oder eine Folge von
Triangle3-Zeilenindizes. `charge=0` ist der Leerlauf; ein festes Potential an
beiden Elektroden stellt den Kurzschluss dar. Genau eine der Angaben
`potential` oder `charge` ist erforderlich.

Der Solver prüft die Außenflächen mit `inspect_boundary_topology`, lehnt
leere oder doppelte Facetten und sich an einem Knoten berührende Elektroden
ab. Ein Floating-Electrode-Knoten darf keine separate Dirichlet-Vorgabe
tragen. Eine feste Elektrode darf eine gleiche explizite Potentialvorgabe
überlagern; unterschiedliche Werte sind ein Fehler. Eine rein schwebende
elektrische Komponente benötigt weiterhin einen Potentialbezug.
Ein vorgeschriebener `surface_charge`-Fluss auf derselben Facette wie eine
Elektrode wird als widersprüchliche zweite elektrische Randbedingung abgelehnt.

Im bestehenden negativen elektrischen Block ist das ursprüngliche Residuum
`r=A x-b` auf einem Potential-Dirichlet-DOF die nach außen gerichtete
`D·n`-Ladung. Die angrenzende Leiterladung ist daher
`Q=-sum(r_i)` über die Elektrode. `PiezoelectricResult.electrode_states`
berichtet Name, gemeinsames Potential, Gesamtladung und gegebenenfalls
vorgegebene Gesamtladung. `electrode_charges` enthält die negativen
elektrischen Residuen auf festen und schwebenden Elektrodenknoten;
`reactions_phi` bleibt der Reaktionsvektor ausschließlich fester DOFs.
Elektrodennamen, Knoten und Vorgaben gehen in den kanonischen Eingabe-Hash
(`piezoelectric-discrete-v2`) ein. Läufe ohne Elektroden behalten das bisherige
`v1`-Schema und dessen Hash.

## Algebraische Formulierung

Für eine schwebende Elektrode bildet die Inzidenzmatrix `P` alle ihre
elektrischen Knotenpotentiale auf **einen** gemeinsamen Freiheitsgrad ab.
Der reduzierte Operator ist `Pᵀ A P` und die rechte Seite `Pᵀ b-Q` auf
diesem Freiheitsgrad. Damit folgt aus der reduzierten Gleichung
`sum(r_i)=-Q`. Mechanische Freiheitsgrade bleiben unverändert. Diese
Kondensation erhält die Symmetrie des monolithischen Operators. Der Solver
rekonstruiert anschließend die ursprünglichen nodalen Potentiale und
berechnet Reaktionen, Gesamtladung sowie feldweise Residuen am ursprünglichen
System. Das elektrische freie Residuum schließt individuelle
Elektrodenknoten aus, deren **Summe** durch die Ladungsbedingung geprüft
wird. Es erfolgt keine willkürliche Erdung.

## Unabhängiger Benchmark

Ein 1-m-Quader mit `Cₓₓₓₓ=10 Pa`, `eₓₓₓ=0,5 C/m²` und
`κₓₓ=2 F/m` erhält die vorgegebene Dehnung `εₓₓ=0,1` und
`u_y=u_z=0`. Die linke Elektrode liegt auf `0 V`. Aus
`Dₓ=eₓₓₓ εₓₓ-κₓₓ ∂φ/∂x` folgt:

| Fall | Rechter Rand | `φ_rechts` | `Dₓ` | `Q_rechts=-Dₓ A` |
|---|---|---:|---:|---:|
| Leerlauf | `Q=0 C` | `0,025 V` | `0 C/m²` | `0 C` |
| Kurzschluss | `φ=0 V` | `0 V` | `0,05 C/m²` | `−0,05 C` |
| Vorgegebene Ladung | `Q=1,95 C` | `1 V` | `−1,95 C/m²` | `1,95 C` |

Die affinen Primärfelder sind auf dem Tet4-Netz exakt darstellbar.
Die Tests setzen für Potential und Ladung absolute Toleranzen von
`10⁻¹² V` und `10⁻¹² C`, für `D` `10⁻¹² C/m²`, jeweils `rtol=0`.
Diese Grenzen berücksichtigen ausschließlich Sparse-LU- und
Elementrundung. Ein fester Elektrodenfall wird zusätzlich mit einzeln
vorgegebenen nodalen Dirichlet-Werten verglichen. Globale Ladungsbilanz und
feldweise Solverresiduen werden geprüft. Tests:
`tests/patch/test_electrodes.py`.

Die allgemeine Modellkonfiguration und Konfliktdiagnostik stehen in
[Doc/21](21_MODEL_SCHEMA.md) und [Doc/23](23_MODEL_PREFLIGHT.md).
Externe Schaltungen und zeitabhängige Leitfähigkeit sind nicht
Teil dieser Formulierung.
