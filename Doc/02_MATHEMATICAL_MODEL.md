# Mathematisches Modell

Status: Entwurf

## 1. Kinematik und Konventionen

Das Gebiet sei `Omega` mit äußerem Normalenvektor `n`. Die Unbekannten sind

- Verschiebung `u: Omega -> R^3`,
- Temperaturabweichung `theta = T - T_ref`,
- elektrisches Potential `phi`.

Es gelten

```text
epsilon(u) = sym(grad u)
E            = -grad phi
```

Die Voigt-Reihenfolge ist projektweit verbindlich:

```text
epsilon_v = [epsilon_xx, epsilon_yy, epsilon_zz,
             gamma_yz, gamma_xz, gamma_xy]^T
sigma_v   = [sigma_xx, sigma_yy, sigma_zz,
             tau_yz, tau_xz, tau_xy]^T
```

Dabei sind `gamma_ij = 2 epsilon_ij` technische Schubverzerrungen. Die
elektrischen Komponenten folgen `[x, y, z]`. Materialdaten aus Quellen mit
anderer Reihenfolge oder tensorialen Schubkomponenten müssen beim Import
explizit transformiert werden.

## 2. Bilanzgleichungen

Für den ersten quasistatischen Ausbau gelten im Volumen:

```text
-div(sigma)              = b
-div(k grad T)           = Q
 div(D)                  = rho_f
```

mit Volumenkraft `b`, Wärmequelle `Q` und freier Ladungsdichte `rho_f`.

Randbedingungen werden in Dirichlet- und Neumann-Anteile zerlegt:

```text
u = u_bar                 auf Gamma_u
sigma n = t_bar           auf Gamma_t

T = T_bar                 auf Gamma_T
-k grad(T) dot n = q_bar  auf Gamma_q

phi = phi_bar             auf Gamma_phi
D dot n = omega_bar       auf Gamma_omega
```

Die Vorzeichen von Wärmestrom und Oberflächenladung werden zusammen mit der
äußeren Normalenrichtung in der Eingabespezifikation festgelegt.

## 3. Konstitutive Gleichungen

Die lineare piezoelektrische Stress-Charge-Form mit thermischer Eigendehnung ist

```text
epsilon_m = epsilon - alpha theta
sigma     = C^E epsilon_m - e^T E
D         = e epsilon_m + kappa^S E + p theta
q         = -k grad T
```

`C^E` ist die Steifigkeit bei konstantem elektrischen Feld, `e` die
piezoelektrische Spannungs-Ladungs-Matrix, `kappa^S` die Permittivität bei
konstanter Verzerrung und `p` ein optionaler pyroelektrischer Vektor. Im
stationären Basismodell bis v1.0 gilt `p = 0`. Auch der Beitrag
`-e alpha theta` in `D` wird nur aktiviert,
wenn das Materialmodell ihn explizit vorsieht; diese Entscheidung wird durch
einen Materialmodelltyp sichtbar gemacht und niemals implizit getroffen.

Für isotrope Elastizität:

```text
lambda = E_young nu / ((1 + nu)(1 - 2 nu))
mu     = E_young / (2(1 + nu))
```

mit `E_young > 0` und `-1 < nu < 0.5`. Nahezu inkompressible Materialien sind
für den reinen Verschiebungs-Tet4 ungeeignet. Isotrope Eingaben mit `nu >= 0.49`
werden im implementierten Einzelfeldsolver ab v0.3.0 abgelehnt (ADR-007).

## 4. Schwache Formen

Für Testfunktionen `v`, `w`, `psi`, die auf den jeweiligen Dirichlet-Rändern
verschwinden, lautet die mechanische schwache Form:

```text
integral_Omega epsilon(v):sigma dOmega
= integral_Omega v dot b dOmega + integral_Gamma_t v dot t_bar dGamma
```

Für Wärme:

```text
integral_Omega grad(w) dot k grad(T) dOmega
= integral_Omega w Q dOmega - integral_Gamma_q w q_bar dGamma
```

Hier ist `q_bar = q dot n` positiv nach außen. Partielle Integration von
`-div(k grad T)=Q` liefert den Randterm `k grad T dot n = -q_bar`.
Für Konvektion `q dot n = h (T - T_inf)` steht daher
`+ integral_Gamma_h w h T` links und `+ integral_Gamma_h w h T_inf` rechts.
Diese Korrektur des früheren Pluszeichens ist in ADR-008 dokumentiert.

Für Elektrostatik verwenden wir die mit `-1` multiplizierte Gauß-Gleichung. Das
liefert zusammen mit der Mechanik eine symmetrische Sattelpunktstruktur:

```text
integral_Omega grad(psi) dot D dOmega
= integral_Gamma_omega psi omega_bar dGamma
 - integral_Omega psi rho_f dOmega
```

Der genaue Ladungsvektor wird durch einen Einflächen-Kondensatortest geprüft.
Diese Festlegung ist wesentlich: Eine isolierte Änderung des elektrischen
Residualvorzeichens verändert Kopplungsblöcke und Lastvorzeichen.

## 5. Diskrete Blockstruktur

Mit der DOF-Reihenfolge `[u, T, phi]` entsteht stationär schematisch

```text
[ K_uu   K_uT    K_uphi ] [u    ]   [f_u  ]
[  0      K_TT      0   ] [theta] = [f_T  ]
[ K_phiu K_phiT -K_phiphi] [phi  ]   [f_phi]
```

Für reine Piezoelektrizität gilt bei konsistenter Konvention
`K_phiu = K_uphi^T`. Der negative dielektrische Block macht das System
symmetrisch-indefinit. Deshalb ist Conjugate Gradient für das monolithische
Gesamtsystem nicht zulässig. Ein direkter Sparse-Solver ist die Referenz;
MINRES mit geeignetem Blockpräconditioner ist eine spätere Option.

Die thermische Zeile ist im stationären Ausbau M3–M5 unabhängig. Das entspricht
stationärer Wärmeleitung mit thermischer Ausdehnung, nicht vollständig reziproker
Thermoelastizität. Für letztere muss im transienten Modell der thermomechanische
Entropieterm konsistent ergänzt werden.

## 6. Eindeutigkeit und Nullräume

- Mechanik benötigt ausreichende Lagerung zur Entfernung aller Starrkörpermoden.
- Elektrostatik benötigt mindestens ein Referenzpotential oder eine explizite
  Nullraumbehandlung bei reinen Ladungsrandbedingungen.
- Wärmeleitung benötigt eine Temperaturreferenz, sofern ausschließlich
  Neumann-Randbedingungen wirken und die Nettowärmequelle kompatibel ist.

Der Solver soll diese Situationen vor der Faktorisierung diagnostizieren und
nicht nur eine generische Singularitätsmeldung ausgeben.

## 7. Einheiten

Intern gilt SI:

| Größe | Einheit |
|---|---|
| Länge, Verschiebung | m |
| Spannung, Elastizitätsmodul | Pa |
| Temperaturdifferenz | K |
| Wärmeleitfähigkeit | W/(m K) |
| Wärmekapazität | J/(kg K) |
| Potential | V |
| elektrisches Feld | V/m |
| elektrische Verschiebung | C/m^2 |
| Permittivität | F/m |
| piezoelektrische Matrix `e` | C/m^2, äquivalent N/(V m) |

Andere Einheitensysteme sind nur über einen vorgeschalteten, getesteten
Konverter zulässig.

## Experimenteller CFD-Operator

Der reduzierte Kanaloperator, SI-Vorzeichen und Erhaltungsbilanzen sind in
[Doc/68](68_EXPERIMENTAL_CFD_CORE.md) festgelegt; kein allgemeiner RANS-Solver.

## Experimentelle Stokes-Druckkopplung

Für den geschlossenen 2D-MAC-Schnitt gelten -mu Laplace(u)+grad(p)=f und
div(u)=0 mit G=-D^T und explizitem Druck-Gauge.
[Vorzeichen, Randwerte und Bilanzkonventionen](69_CFD_PRESSURE_COUPLING.md).

## Stationärer inkompressibler MAC-Kern

Konservative Konvektion, Picard/Oseen-Schritte, gewichtete Druckadjungiertheit
und der Geschwindigkeitseinlass/Druckauslass ergänzen den Stokes-Vertrag.
[Mathematik und Akzeptanz #86](70_STATIONARY_FLOW_CORE.md).

## M12: RANS-/FSI-Zielvertrag (#84, ADR-049)

Die nächste Schließung ist stationäres inkompressibles SA-noft2 mit expliziter
S-Modifikation und aufgelöster Wand. Für variable Wirbelviskosität gilt
`sigma_f=-p I+2(mu+rho nu_t)S(U)` und
`rho div(U tensor U)-div(sigma_f)=f`, `div(U)=0`; ein skalares
Viskositäts-Laplaceprodukt ersetzt diesen Spannungstensor nicht.
Auf der FSI-Schnittstelle ist die Strukturtraktion `-sigma_f n_f` mit
entgegengesetzten Fluid-/Strukturnormalen. Eine transponierte Lastabbildung
zur Verschiebungsinterpolation erhält Kräfte, Momente und virtuelle Arbeit.
[Der vollständige Vertrag](71_M12_RANS_FSI_CONTRACT.md) enthält schwache Formen,
SA-Gleichung/Konstanten, Wandabstand, Druckreferenz, externe Drucklast,
SI-/Breitenkonventionen, Netzgrenzen und unabhängige Abnahmegates.
Dies ist eine Entwicklungsentscheidung; RANS und FSI sind noch nicht implementiert.
