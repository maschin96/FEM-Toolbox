# Stationäre FSI am geführten Kanaldeckel (#93)

Zwei vorbereitete Fälle aus dem Nachweis [Doc/93](../../Doc/93_STATIONARY_FSI_VERIFICATION.md),
jeweils das gröbste Netzpaar: laminar (Fluid 10×4, 0,9 s) und SA-noft2
(32 wandaufgelöste Zellen je Halbkanal × 4, ca. 25 s). Schnittstelle:
[Doc/94](../../Doc/94_FSI_WORKFLOW.md).

Nach Installation des Pakets, auch aus Wheel oder sdist:

```sh
python examples/fsi_lid/run.py laminar --output fsi_laminar_results
python examples/fsi_lid/run.py sa --output fsi_sa_results
```

Jeder Fall besteht aus `fsi.toml` (fsi_schema_version 1, alle Kopplungswerte
Pflicht), `domain.json`/`mesh.json` (FluidDomain), `structure.toml`/`structure.msh`
(Elastizitätsmodell mit Gmsh-Gruppen lid/spring und wet/top/side) und im
SA-Fall `nu_tilde.json` (Startfeld, nur Anfangswert).

Ausgabe, ein vorhandenes Verzeichnis wird abgelehnt:

- `report.json`: Status, Konfiguration, SHA-256 aller Eingabedateien,
  Versionen, Einstellungen und Gates, Netz-/Gebietshashes, Kopplungs-,
  Fluid- und Strukturbilanzen, Fluidprovenienz, SA-Konstanten.
- `history.csv`: eine Zeile je Kopplungsschritt, Einheiten im Spaltennamen.
- `fluid_fields.json`: aktuelles Netz, U, p, Face-Massenfluss, Wandtraktion,
  im SA-Fall ν̃, Wirbelviskosität und y⁺.
- `structure_fields.json`: Referenznetz, Geometrie- und Strukturverschiebung,
  Schnittstellenkräfte, Reaktionen, Spannungen.

Bei einem kontrollierten Abbruch schreibt der Lauf `report.json` mit
Fehlerart und letztem Iterat sowie `history.csv` und endet mit Exitcode 2.
Die Eingaben werden mit `tests/reference/fsi_lid.py` erzeugt; ein
Regressionstest prüft, dass sie bytegleich reproduzierbar sind.
