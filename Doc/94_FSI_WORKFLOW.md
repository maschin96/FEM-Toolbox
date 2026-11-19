# M12 — deklarativer stationärer FSI-Workflow (#93)

Stand 19.11.2026, Entwicklungskandidat 1.7.0.dev1, experimentelle API
`fem_toolbox.coupling.fsi_config`. Grundlage sind der
[M12-Vertrag](71_M12_RANS_FSI_CONTRACT.md) (Abschnitt 6: eigene
`fsi_schema_version=1`, strikte TOML mit Pflichtwerten),
die [Kopplung](92_STATIONARY_FSI_COUPLING.md) und der
[Nachweis](93_STATIONARY_FSI_VERIFICATION.md). Model-Schema 1–3,
FluidDomain-/FluidMesh-Schema 1 und FlowProblem-Schema 1 bleiben unverändert;
die FSI-Datei verweist nur auf diese bestehenden Formate.

## Schema

```toml
fsi_schema_version = 1

[fluid]
domain = "domain.json"           # FluidDomain-JSON, relativ, mit eigenem Netzpfad
model = "laminar"                # oder "sa-noft2"
initial_nu_tilde = "nu_tilde.json"  # nur sa-noft2: Zahl oder JSON {"nu_tilde_m2_s": [...]}

[fluid.options]                  # optional, nur freigegebene Solveroptionen
sa_relaxation = 0.1

[structure]
model = "structure.toml"         # bestehendes Elastizitätsmodell mit Gmsh-Netz
facet_set = "wet"

[coupling]                       # alle Werte Pflicht, keine Defaults
wall = "top"
length_scale_m = 0.1
velocity_scale_m_s = 1.0
external_pressure_Pa = 0.0
relaxation = 1.0
max_coupling_iterations = 200
divergence_window = 10
span_origin_m = 0.0
preserve_inlet_flow = true
```

Regeln: exakte Schlüsselmengen; ganze Zahlen und Wahrheitswerte typstreng
(`2.0` ist kein Budget, `1` kein Boolean), Gleitkommawerte endlich. Pfade sind
relativ, bleiben im Konfigurationsverzeichnis und müssen die erwartete Endung
haben. `fluid.options` erlaubt für beide Modelle `force_N_m3`,
`pressure_reference_Pa`, `convection_scheme`, Relaxationen und
Ressourcenbudgets, für SA zusätzlich `sa_relaxation`; Skalen, Modell und
SA-Startfeld sind keine Optionen. Fehler nennen den Dokumentpfad
(`ModelSchemaError`). JSON und TOML laden nur Daten, keinen Code.

## API

```python
from fem_toolbox.coupling.fsi_config import load_fsi_config, run_fsi_config

config = load_fsi_config("fsi.toml")  # FsiConfig: problem, document, input_files
result = run_fsi_config("fsi.toml", "out")  # löst und schreibt Artefakte
```

`load_fsi_config` liest und validiert, baut `StationaryFsiProblem` und
berechnet SHA-256 von Konfiguration, FluidDomain, optionalem Startfeld,
Strukturmodell und Gmsh-Netz. `run_fsi_config` lehnt ein vorhandenes
Ausgabeverzeichnis vor dem Laden ab, löst, schreibt immer Artefakte und
wirft einen `FsiCouplingError` nach dem Schreiben der Abbruchartefakte erneut.
`write_fsi_artifacts` ist separat verwendbar. Der Koppler selbst liest keine
Dateien.

## Artefakte

| Datei | Inhalt (SI-Einheit im Schlüssel) |
| --- | --- |
| `report.json` | Status, Konfiguration, Eingabe-SHA-256, Versionen (fem_toolbox, NumPy, SciPy, Python), Einstellungen und Gates, `input_sha256`, Netz-/Gebiets-/Spurhashes, Laplace-Rest, Kopplungs-, Fluid- und Strukturbilanzen, Fluidprovenienz; SA: Konstanten, max y⁺, SA-Rest |
| `history.csv` | je Kopplungsschritt: Residuen, Laständerung, Kräfte in y, Transferaudit, Kleinverformung, Flächenverhältnis, Fluiditerationen |
| `fluid_fields.json` | aktuelles Netz und Netzverschiebung, U, p, Face-Massenfluss, Wandtraktion samt Wandzeilen; SA: ν̃, Wirbelviskosität, y⁺ |
| `structure_fields.json` | Referenznetz, Geometrie- und Strukturverschiebung, Schnittstellenkräfte, Reaktionen, Spannungen |

Abbruch: `report.json` mit `status = "failed"`, Fehlerart, Grund und letztem
Iterat sowie `history.csv`; keine Feldartefakte, da kein akzeptierter Zustand
existiert. Die Versionen stammen aus dem tatsächlich importierten Code.

## Beispiele

`examples/fsi_lid/` enthält beide Deckelfälle als vorbereitete Dateien
(gröbstes #92-Netzpaar). Sie werden mit `tests/reference/fsi_lid.py` erzeugt;
ein Test prüft Byte-Gleichheit (JSON-Zahlen relativ 1e-12, weil das
SA-Profil aus einer libm-abhängigen Randwertlösung stammt), damit Netz,
Randbedingungen und Parameter rekonstruierbar bleiben.

```sh
python examples/fsi_lid/run.py laminar --output fsi_laminar_results
python examples/fsi_lid/run.py sa --output fsi_sa_results
```

Exitcode 0 bei Konvergenz, 2 bei kontrolliertem Abbruch. CI kopiert die
Beispiele in die Release-Artefakte und führt beide aus frisch installiertem
Wheel und sdist mit `python -I` außerhalb des Checkouts aus.

## Wiederholung und Toleranzen

- Identische Dateien auf derselben Plattform: `history.csv` byte-gleich,
  gleiche `input_sha256` und Netzhashes (deterministische direkte Löser).
- Gegen die dokumentierten #92-Werte des gröbsten Netzpaars (Deckelweg
  laminar 5.71964e-4 m, SA 9.64275e-4 m): Abweichung ≤ 2e-9 m. Begründung:
  Das Kopplungsgate begrenzt max|r_d| auf 1e-8 Lref = 1e-9 m; die
  Gmsh-/JSON-Rundreise ordnet Zellen anders, verändert aber nicht die Lösung.

Nachweise: `tests/unit/test_fsi_config.py` (Problemaufbau, Hashes,
14 Fehlerdokumente mit konkretem Pfad, vorhandene Ausgabe, Abbruchartefakte)
und `tests/regression/test_fsi_lid_example.py` (beide Beispiele als
Unterprozess, Artefakte und Einheiten, Überschreibschutz, exakte
Wiederholung, reproduzierbare Eingaben).

Nicht abgedeckt: unabhängiges menschliches Fachreview, M12-Freigabe (#94).
