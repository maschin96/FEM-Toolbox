# M12 — stationäre bidirektionale FSI-Kopplung (#91)

Stand 18.11.2026, Entwicklungskandidat 1.7.0.dev1, experimentelle API in
`fem_toolbox.coupling`. Maßgeblich ist Abschnitt 5 des
[M12-Vertrags](71_M12_RANS_FSI_CONTRACT.md). Die Bausteine sind die
[konservative Schnittstelle](90_CONSERVATIVE_FSI_INTERFACE.md) (#89) und die
[harmonische Netzbewegung](91_STATIONARY_MESH_MOTION.md) (#90). Dies ist keine
M12-Freigabe: Die verpflichtenden laminaren und SA-Deckelreferenzen auf vier
Netzstufen bleiben #92, die deklarativen Rezepte #93.

## API

```python
from fem_toolbox.coupling import StationaryFsiProblem, solve_stationary_fsi

result = solve_stationary_fsi(
    StationaryFsiProblem(
        fluid_domain=reference_domain,  # unveränderte FluidDomain
        wall_name="top",  # genau eine Wand mit fsi=true
        structure=model,  # validiertes statisches Elastizitätsmodell
        facet_set="wet",  # benannter Triangle3-Satz
        length_scale_m=1.0,
        velocity_scale_m_s=1.0,
        fluid_model="laminar",  # oder "sa-noft2" mit initial_nu_tilde_m2_s
        fluid_options={...},  # unverändert an solve_quad_flow/solve_rans
        external_pressure_Pa=0.0,
        relaxation=0.2,
        max_coupling_iterations=200,
    )
)
```

`fluid_options` werden unverändert an `solve_quad_flow` oder `solve_rans`
übergeben. Skalen, Modellwahl und SA-Startfeld sind Kopplungseingaben und als
Optionen verboten. Die Strukturkoordinaten müssen den Referenzkoordinaten der
Schnittstelle entsprechen. Vorhandene Dirichletwerte, Körperlasten und
Traktionen des Modells bleiben erhalten; übertragene Knotenkräfte kommen hinzu.
Genau eine FSI-Wand wird unterstützt; mehrere markierte Wände sind ein
Eingabefehler statt einer stillen Teilkopplung.

## Iteration

Für k = 0, 1, … mit d⁰ = 0 auf dem unveränderten Struktur-Referenznetz:

1. Spur `InterfaceMap.segment_endpoint_displacements(dᵏ)` der Referenzabbildung
   bilden und `move_fluid_domain` auf das Referenzgebiet anwenden (absolute
   Wege, keine aufsummierten Inkremente).
2. Fluid auf der neuen Geometrie vollständig lösen.
3. Die Überdeckung zwischen aktueller Fluidwand und verformter Struktur
   (`x + dᵏ`) neu aufbauen und die aktuelle Fluidtraktion konservativ als
   Knotenkräfte Fᵏ übertragen; das Erhaltungsaudit 1e-10 bleibt Pflicht.
4. Lineare Struktur auf dem Referenznetz mit Fᵏ lösen; Gleichgewicht
   ≤1e-10, max|u|/Lref ≤ 0.02 und max‖∇u‖_F ≤ 0.02 prüfen.
5. r_d = d_solved − dᵏ, dᵏ⁺¹ = dᵏ + ω r_d mit konstantem ω ∈ (0,1].

Erfolg verlangt im selben Schritt **beide** Kriterien:
max_Knoten |r_d| ≤ 1e-8 Lref und max_Knoten |Fᵏ − Fᵏ⁻¹| ≤ 1e-6 ρ U_ref² b Lref.
Im ersten Schritt ist die Laständerung ∞; deshalb gibt es frühestens im
zweiten Schritt einen Erfolg. Eine kleine relaxierte Änderung allein genügt nie.
Fluid, Lasten und Struktur des akzeptierten Schritts stammen alle frisch aus
der Endgeometrie dᵏ; es werden keine alten Lasten als Endergebnis ausgegeben.
Das Ergebnis enthält sowohl die Geometrieverschiebung dᵏ als auch die
Strukturlösung d_solved; ihr Abstand ist durch das Weggate begrenzt.

## Ergebnis und Historie

`StationaryFsiResult` enthält Fluidresultat (`QuadFlowResult` oder
`RansResult`), `MeshMotionResult`, aktuelle `InterfaceMap`, `InterfaceLoads`,
`ElasticResult`, beide Verschiebungen [nnode,3] in m, die vollständige Historie,
unveränderliche Einstellungen mit `fsi_schema_version=1`, Gates in SI sowie
`input_sha256`. Dieser Hash umfasst Einstellungen, numerischen
Referenzgebiets-Hash der Netzbewegung sowie Strukturdokument und -netz.

Jede `FsiIteration` protokolliert unrelaxiertes max|r_d| absolut und relativ,
Laständerung, Fluid- und Außendruckkraft, Summe der Lagerreaktionen, rohe
Kraft-/Momentfehler und relatives Transferaudit, beide Kleinverformungsmaße,
minimales Zellflächenverhältnis und Fluiditerationen.

## Kontrollierte Abbrüche

`FsiCouplingError.kind` ∈ `input`, `interface`, `mesh`, `fluid`, `structure`,
`divergence`, `budget`. Der Fehler enthält die vollständige Historie bis zum
letzten abgeschlossenen Schritt und das letzte Iterat dᵏ (unveränderlich).

- `mesh`: `MeshMotionError`, etwa Weggrenze, Geometrie oder Eckkonflikt
  zwischen FSI- und fester Wand.
- `fluid`: Nichtkonvergenz, Ressourcen-, Rückström- oder y⁺-Fehler des Solvers.
- `interface`: fehlende Überdeckung, spannweitenabhängige oder nicht ebene
  Spur, Erhaltungsaudit.
- `structure`: Kleinverformungs- oder Gleichgewichtsgates.
- `divergence`: nicht endliche Werte oder `divergence_window` (Standard 10)
  aufeinanderfolgende Zunahmen des unrelaxierten max|r_d|.
- `budget`: `max_coupling_iterations` (Standard 200) erschöpft.

Kein Clipping, keine adaptive Relaxation, kein Remeshing und kein instationärer
Ersatz. Iterationsnummern sind keine physikalische Zeit: keine
Wandgeschwindigkeit und kein ALE-Zeitfluss.

## Nachweise

`tests/unit/test_stationary_fsi.py` (18 Tests):

- Ruhendes periodisches Fluid bei 7 Pa, außen 5 Pa, Block E = 1000 Pa,
  ν = 0, Höhe 1 m: unabhängig d_y = 0.002 (2 − y) m. Konstante Last liefert
  exakt r_k = 0.8ᵏ d* und die kleinste Schrittzahl mit r ≤ 1e-8 Lref;
  Fluid-, Außen- und Reaktionskraft (0, 84, 0), (0, −60, 0), (0, −24, 0) N.
- Hydrostatischer Deckel mit ρg = 100 N/m³: Anheben senkt den Deckeldruck.
  Harmonische Bewegung ist affin, daher gilt unabhängig vom Solver
  d* = (p₀ − p_ext) h / (E + ρ g h (1 − y_c/H₀)) und die konstante
  Kontraktion 1 − ω (1 + 0.0875). Lösung innerhalb der aus dem Weggate
  folgenden Schranke 2.4e-9 m, rund 8 % unter der Einweglösung.
- Reproduzierbare Historie und Identität; geänderte Relaxation ändert den Hash.
- Budgetabbruch mit Historie und letztem Iterat (1 − 0.8³) d*; Struktur-,
  Fluid- und Netzfehler getrennt klassifiziert; Divergenzkriterium; elf
  Eingabefehler vor jedem Solve.
- SA-noft2-Pfad: wandaufgelöster periodischer Kanal mit y⁺ ≤ 1 und
  gleichmäßigem Druck reproduziert die analytische Strukturverschiebung.

TDD: Zunächst fehlender Kopplungsimport. Ein erster Kanalfall mit
Druckverlauf entlang x scheiterte korrekt am Spannweitengate: Ein grober
Tet4-Block antwortet auf ungleichmäßige Last nicht exakt spannweiteninvariant.
Solche Strukturen benötigen den geführten Deckel aus #92; der Nachweis nutzt
deshalb den exakt linearen hydrostatischen Rückkopplungsfall.

Nicht abgedeckt: laminare und SA-Deckelreferenz mit vier Netzstufen (#92),
TOML-Workflow und Provenienzartefakte (#93), unabhängiges menschliches
Fachreview.
