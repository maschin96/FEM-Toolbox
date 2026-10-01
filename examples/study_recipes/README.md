# Zwei stationäre Studienrezepte (#61)

[Wärme: Leitfähigkeit](heat/README.md) · [Piezo: Geometrie und Material](piezo/README.md)

Diese Rezepte benötigen den M8-Kandidaten mit `Study` und `study_reports`.
Das veröffentlichte stabile v1.2.0 enthält diese APIs noch nicht. Die interne
Paketversion des Kandidaten lautet weiterhin 1.2.0; deshalb **Artefakthash und
Quellcommit**, nicht allein Versionsnummer, zur Identifikation verwenden.
Keine Visualisierungs-/Surrogat-Extras, kein Gmsh-Programm und kein Netzservice nötig.

## Saubere Installation und Ausführung

Aus dem Issue-/Kandidatencheckout (Python 3.11–3.13, Buildwerkzeug `uv`):

```sh
uv build --out-dir /tmp/fem61-artifacts
git rev-parse HEAD
shasum -a 256 /tmp/fem61-artifacts/*
tar -czf /tmp/fem61-recipes.tar.gz -C examples study_recipes
shasum -a 256 /tmp/fem61-recipes.tar.gz
```

Wheel enthält die Python-API; die Beispiele werden als separates Rezeptpaket
weitergegeben. Auch das bestehende sdist enthält keine Beispielnetze. Deshalb
**beide** Artefakte (Paketarchiv und Rezeptpaket) mit ihren Hashes aufbewahren.
Die Anleitung benötigt keinen Import aus einem Entwicklungscheckout.
In einem neuen Arbeitsverzeichnis außerhalb des Checkouts:

```sh
mkdir /tmp/fem61-user
cd /tmp/fem61-user
python3 -m venv env
env/bin/python -m pip install /tmp/fem61-artifacts/fem_toolbox-1.2.0-py3-none-any.whl
tar -xzf /tmp/fem61-recipes.tar.gz
env/bin/python -I study_recipes/run.py heat --output heat-results
env/bin/python -I study_recipes/run.py piezo --output piezo-results
```

Windows: neue venv, `env\Scripts\python.exe` statt `env/bin/python`, Pfade
entsprechend anpassen. Die Plattformprüfung hier erfolgte unter macOS.
Für die sdist-Reproduktion eine **weitere** neue venv erstellen und stattdessen
`fem_toolbox-1.2.0.tar.gz` installieren. `-I` ignoriert `PYTHONPATH` und den
Benutzer-Sitepfad. `verification.json.package_path` muss in der neuen venv liegen.
Für die exakt geprüfte Abhängigkeitskombination vor der Paketinstallation
`numpy==2.5.3 scipy==1.18.1 meshio==5.3.5` installieren (Python 3.13).

Jeder Aufruf verlangt ein neues Ausgabeziel und kopiert die Eingaben dorthin.
Er führt kontrolliert einen Fall aus, exportiert `partial/` mit sichtbaren
Fehlwerten, setzt fort und verifiziert unveränderte Wiederverwendung. `table/`
enthält `values.csv`, `metadata.json`, `report.md`; `runs/` enthält Katalog,
VTU und Bündelmanifest. Separate direkte Model-Läufe und die analytische
Referenz prüfen Kennwerte und Felder. `failure_table/` hält einen bewusst
fehlerhaften Fall sichtbar. Anschließend ändert der Runner eine Basiseingabe,
prüft Invalidierung und stellt die ursprüngliche Eingabe samt Ergebnissen wieder her.
`verification.json` enthält Soll/Ist, direkte Werte, Toleranzen, Paketpfad und Laufzeit.

## Fortsetzen und Fehler selbst nachvollziehen

Alternativ zur automatischen Abnahme in einer frischen Kopie des Rezeptes:

```python
from fem_toolbox.studies import Study
from fem_toolbox.study_reports import export_study_table

# Im Verzeichnis study_recipes; der Runner importiert nur öffentliche APIs.
from run import quantities

s = Study.from_file("heat/study.toml")  # alternativ piezo/study.toml
first = s.run("manual-runs", max_cases=1)
export_study_table("manual-runs/catalog.json", "manual-partial", quantities("heat"))
complete = s.run("manual-runs")
reused = s.run("manual-runs")
assert complete == reused
```

`max_cases=1` begrenzt neue Versuche und bildet einen kontrollierten Stopp
zwischen Fällen ab; es simuliert keinen harten Prozessabsturz während des Solves.
Bei tatsächlichem Ctrl-C erneut `run` am selben Ziel aufrufen; ein als `running`
verbliebener Fall wird nach dem bestehenden Studienvertrag erneut versucht.
In `heat/model.toml` Quelle 8→16 W/m³ ändern oder in `piezo/model.toml`
obere Spannung 100→50 V ändern, dann `s.run("manual-runs")` ausführen.
IDs bleiben gleich, Inputhashes ändern sich und `attempts` steigt; die Änderung
zurücknehmen und erneut ausführen. Eingaben nicht während eines Solves ändern.

Für einen sichtbaren Fehlfall in einer **separaten** Wärmestudienkopie
`values = [2.0, -1.0]` im Leitfähigkeitsparameter setzen. Der negative Fall
ist `failed`, der gültige `succeeded`; seine Tabellenwerte bleiben leer mit
Fehlergrund. Der automatische Piezo-Fehlfall verwendet eine negative definite
Permittivität neben einer positiven Variante. Die Fehlversuche dienen nur der
Lebenszyklusprüfung und gehören nicht zum fachlichen Parametervergleich.

## Interpretation und Troubleshooting

Fall-ID, Parameterwerte, Fingerprints und ursprüngliche Model-Provenienz in
`metadata.json` verbinden jede Zeile mit `runs/catalog.json` und dem Bündel.
`solver_status`, explizite `balance_checks` und `mesh_convergence=not_checked`
sind getrennte Aussagen. Kleine Bilanz/Residuen ersetzen keine Netzkonvergenz.
Die Rezepte sind **verifiziert**, nicht experimentell validiert.

- `ModuleNotFoundError`: Kandidatenarchiv in derselben venv installieren;
  Importpfad im Protokoll kontrollieren. Stabiles v1.2.0 genügt nicht.
- `FileExistsError`: neues Ziel für den Abnahmerunner wählen; normale Fortsetzung
  dagegen mit `Study.run` im bestehenden Studienziel ausführen.
- `failed`: Fehler im Katalog lesen. Leitfähigkeit/Permittivität müssen positiv
  definit sein; SI-Einheiten, Tensorform und vorhandene Zielpfade prüfen.
- Netz-/Setfehler: Netze zusammen mit TOML kopieren; Pfade innerhalb des
  Modellverzeichnisses halten, benannte Regionen/Elektroden nicht umbenennen.
- Singularität: mechanische Lagerung und Potentialreferenz kontrollieren;
  vorgegebene `x/y`-Sperren sind Teil des Piezo-Referenzmodells.
- Fehlwerte/`invalid`: CSV/JSON/Bündel nicht von Hand ändern. Hashfehler beheben
  durch erneuten Studienlauf; Bericht in ein neues Verzeichnis exportieren.
- Unerwartete Neuberechnung: geänderte Netzbytes, Basiseingaben oder
  Dependency-Versionen invalidieren bewusst den konservativen Fingerprint.
- Falsche Ladungs-/Energievorzeichen: Elektrode und Körpernormalen unterscheiden;
  die dielektrische Enthalpie unter festem Potential ist negativ.

[Abnahmeprotokoll und Laufbudget](../../Doc/63_STUDY_RECIPE_ACCEPTANCE.md).
