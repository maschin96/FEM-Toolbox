# M12 — Verifikation der stationären bidirektionalen FSI (#92)

Stand 19.11.2026, Entwicklungskandidat 1.7.0.dev1. Maßgeblich sind die
#92-Zeilen und Abschnitt 5 des [M12-Vertrags](71_M12_RANS_FSI_CONTRACT.md);
geprüft wird die [Kopplungsiteration](92_STATIONARY_FSI_COUPLING.md) (#91).
Benchmark, Referenzen, Auswertegrößen und Toleranzen standen vor der
Ergebnisrechnung fest. Dies ist eine Codeverifikation des gekoppelten
Grenzfalls, keine experimentelle Validierung und keine M12-Freigabe.

## Benchmark: elastisch geführter Kanaldeckel

Paralleler Kanal, Höhe H0 = 0.1 m, Länge L = 1 m, Breite b = 1 m, fester
Volumenstrom Q über ein entwickeltes Einlassprofil, Druckauslass p_out = 0,
p_external = 0. Die obere Wand ist ein geführter Deckel auf einer Feder k_s.

**Strukturmodell** (`tests/reference/fsi_lid.py`): Tet4-Deckelplatte
(Dicke 0.1 m) auf einer Federschicht (Höhe 0.1 m), beide ν = 0. Alle Knoten
sind in x/z geführt, die Oberseite in y fest. Die Federschicht ist daher
einachsig und k_s = E_s b L / h_s gilt exakt. Der Plattenmodul beträgt
1e9 E_s; als Balken auf elastischer Bettung ist β L < 1e-2, der Deckel verteilt
den linearen Druck also praktisch gleichmäßig (gemessene Spreizung ≤ 3.4e-11 m).
Ein Verhältnis 1e6 verbiegt den Deckel bereits so weit, dass die
Überdeckungsprüfung auf aktueller Geometrie (1e-10 Lref) korrekt scheitert.
Die Kondition des Verhältnisses 1e9 begrenzt die Rekonstruktion der
Lagerreaktionen K u − f auf etwa 1e-7 der Last; geprüft wird daher
Gleichgewicht ≤ 1e-6 der Deckellast, zusätzlich zum unveränderten
Rückwärtsfehlergate ≤ 1e-10 des Struktursolvers.

**Analytische Referenz.** Starrer Deckel, entwickelte Strömung:
F = b L Δp/2. Laminar gilt Δp = 12 μ L Q/(b H³), also
d (H0+d)³ = 6 μ L² Q / k_s (positive Wurzel, `brentq`).

**SA-Referenz.** Der unabhängige 1D-SA-Randwertlöser aus #88 liefert bei H0 und
f_x = 0.002 N/m³ (Re_τ = 587.19) den festen Volumenstrom
Q = 0.018576697 m³/s und Δp0 = f_x L. Bei festem Q und ν gilt
Q/(b ν) = 2 Re_τ U_b⁺(Re_τ); Re_τ hängt also nicht von H ab und
Δp(H) = Δp0 (H0/H)³ — dieselbe Wurzelform wie laminar. Der Test bestimmt
die Druckdifferenz zusätzlich mit dem BVP bei der finalen Höhe; Abweichung
von dieser Ähnlichkeit ≤ 1e-6. Einlass: flächengemittelte BVP-Profile von U
und ν̃ (64-Punkt-Gauß). Die Netzbewegung skaliert U mit H0/H und behält ν̃;
genau das ist das entwickelte Profil bei derselben Re_τ.

**Skalen** (vor dem Lauf festgelegt): Lref = H0, damit Deformations- und
Kopplungsgates auf die Kanalhöhe bezogen sind. Laminar Uref = 1 m/s:
ρ Uref² = 1 Pa ist mit der Stokes-Drucklast von 12 Pa vergleichbar; mit
ρ U_bulk² = 1e-4 Pa läge das Lastgate bei 1e-12 der Last. SA Uref = 0.2 m/s
= 20 u_τ wie in der bestehenden kleinen SI-RANS-Prüfung. Relaxation ω = 1
(zulässig nach Doc/71: ω ∈ (0,1]); die Kopplung ist mit |g′| ≈ 3 d/H stark
kontrahierend.

## Toleranzen (Doc/71)

| Fall | Weg/Druckabfall fein | Letzte Netzänderung | weitere Gates |
| --- | --- | --- | --- |
| Laminar | < 1 % | < 0.2 % | alle Kopplungsgates |
| SA-noft2 | < 2 % | < 0.4 % | y⁺ ≤ 1, ν̃ ≥ 0, alle Kopplungsgates |

Pro Netzpaar zusätzlich: Konvergenz beider Kopplungsgates, Transferaudit
≤ 1e-10 in jedem Schritt, virtuelle Schnittstellenarbeit mit der finalen
Strukturverschiebung ≤ 1e-10, Kleinverformung ≤ 0.02, Flächenverhältnis ≥ 0.5,
fester Einlassstrom und Auslassstrom innerhalb des Massengates, globale
Massenbilanz, Lagerreaktion = −Fluidkraft, monoton fallender Diskretisierungsfehler,
Einwegweg innerhalb derselben Toleranz und gemessener bidirektionaler Effekt
innerhalb 5 % des analytischen.

Druckabfall: lineare Kleinste-Quadrate-Anpassung p(x) über alle Zellen,
Δp = p(0) − p(L), verglichen mit der Referenz bei der tatsächlichen
Geometriehöhe H0 + d.

## Ergebnisse

`tests/convergence/test_stationary_fsi_lid.py`, lokal Python 3.13 / NumPy 2.4.6 /
SciPy 1.18.1, macOS arm64. CI führt beide Fälle als getrennte Pflichtjobs
`fsi-laminar-lid` und `fsi-sa-lid` auf Python 3.11–3.13 aus.

### Laminar (Q = 1e-3 m³/s, μ = ρ = 1, k_s = 1e4 N/m)

Analytisch: Einweg d1 = 6.000e-4 m, gekoppelt d* = 5.895127e-4 m.

| Fluid nx×ny | Struktur nx | Deckelweg [m] | Fehler Weg | Fehler Δp | Schritte |
| --- | --- | --- | --- | --- | --- |
| 10×4 | 5 | 5.71964e-4 | −2.977 % | −3.028 % | 5 |
| 20×8 | 10 | 5.85022e-4 | −0.762 % | −0.777 % | 5 |
| 40×16 | 20 | 5.88382e-4 | −0.192 % | −0.195 % | 5 |
| 80×32 | 40 | 5.89229e-4 | −0.048 % | −0.049 % | 5 |

Ordnung 2.0; letzte Änderung 0.144 % (Weg) und 0.144 % (Δp).
Bidirektionaler Effekt 1 − d*/d1: gemessen 1.7471 %, analytisch 1.7479 %.
Transferaudit ≤ 2.4e-15, virtuelle Arbeit ≤ 2.1e-15, finale
Laständerung ≤ 8.3e-8 N (Gate 1e-7 N). Laufzeit 58 s.

### SA-noft2 (ν = 0.05·0.01/587.19 m²/s, k_s = 1 N/m)

Analytisch: Einweg d1 = 1.0000e-3 m, gekoppelt d* = 9.714147e-4 m.

| Zellen/Halbkanal × nx | Deckelweg [m] | Fehler Weg | Fehler Δp | max y⁺ |
| --- | --- | --- | --- | --- |
| 32 × 4 | 9.64275e-4 | −0.735 % | −0.993 % | 0.498 |
| 48 × 6 | 9.67547e-4 | −0.398 % | −0.442 % | 0.324 |
| 64 × 8 | 9.69322e-4 | −0.215 % | −0.240 % | 0.240 |
| 96 × 12 | 9.70258e-4 | −0.119 % | −0.126 % | 0.158 |

Letzte Änderung 0.097 % (Weg) und 0.112 % (Δp). Bidirektionaler Effekt:
gemessen 2.8553 %, analytisch 2.8585 %. Jeweils 5 Kopplungsschritte,
Transferaudit ≤ 1e-16, finale Laständerung ≤ 1.4e-10 N (Gate 4e-9 N).
Laufzeit 9.4 min.

Die Netzstufen bleiben an die Schichtwachstumsgrenze 1.2 der SA-Wandauflösung
gebunden; 24 Zellen je Halbkanal würden sie mit Exponentialgradierung β = 4.5
verletzen, deshalb beginnt die Folge bei 32.

## Abgrenzung

- Der Deckel ist ein geführter, praktisch starrer Lastverteiler; ein
  ungleichmäßig nachgiebiger Tet4-Körper erzeugt spannweitenabhängige Spuren,
  die das Schnittstellengate korrekt ablehnt (vgl. Doc/92).
- Voll entwickelte Kanalströmung; Einlauf-, Ablösungs- oder
  Krümmungseffekte werden nicht nachgewiesen.
- Der SA-Fall verifiziert Geometrie-/Lastkopplung gegen die Modellgleichung
  (unabhängiger BVP), nicht das Turbulenzmodell gegen Messungen.
- Keine dynamische FSI, keine große Verformung, keine allgemeine 3D-Freigabe.
- Wiederholung in sauberer Umgebung erfolgt über die CI-Jobs; unabhängiges
  menschliches Fachreview liegt nicht vor.
