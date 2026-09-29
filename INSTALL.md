# Installation und erster Lauf

Stand: 29. September 2026. Diese Anleitung richtet eine **Entwicklungsumgebung
aus dem Repository** ein. Die stabile Version ist **v1.2.0**. Ein späterer
Checkout kann bereits weitere Entwicklungsänderungen enthalten.
[Supportmatrix und Grenzen](Doc/57_M7_RELEASE_CANDIDATE.md).
Ein PyPI-Release des aktuellen Checkouts wird hier nicht vorausgesetzt.

## Voraussetzungen und geprüfte Plattformen

- Reguläres **CPython 3.11, 3.12 oder 3.13** mit `venv` und `pip` sowie Git.
  Die Projektmetadaten erlauben `>=3.11,<3.14`.
- Die CI prüft Installation und Tests unter **Linux** auf allen drei
  Python-Versionen. **macOS/ARM64 mit Python 3.13** wurde zusätzlich lokal
  geprüft. Windows und andere Plattform-/Python-Kombinationen haben derzeit
  keine eigene Projekt-CI-Zusage. Details stehen in der
  [Kompatibilitätsmatrix](Doc/28_DEPENDENCY_COMPATIBILITY.md).
- Wheel und sdist werden zusätzlich ohne Repository-Checkout installiert und
  mit stationären und transienten/dynamischen Referenzbeispielen geprüft; Umfang und Freigabegrenzen
  stehen im [Artefaktnachweis](Doc/29_RELEASE_ARTIFACTS_AND_PROTECTION.md).
- Internetzugriff für den ersten Paketdownload; ungefähr 1 GB freier Platz
  für Python-Umgebung, numerische Wheels und optional Gmsh ist sinnvoll.
  C/C++-Compiler und separate Gmsh-Installation sind für die bereitgestellten
  `.msh`-Dateien nicht nötig.

Das Paket installiert **NumPy, SciPy und meshio** als Laufzeitabhängigkeiten.
Die Entwicklungsoption `dev` ergänzt pytest, Ruff und mypy. Die aufgelösten
Versionen bestimmt `pip` gemäß [pyproject.toml](pyproject.toml); die in
`ci/oldest-constraints.txt` enthaltenen Pins dienen ausschließlich der
CI-Prüfung der ältesten unterstützten Kombination. Gmsh ist ein **separates,
optionales Werkzeug zum Erzeugen neuer Netze** und keine Paketabhängigkeit.
Seine GPL-Lizenz betrifft die externe Installation; die Toolbox bleibt MIT.

## Projekt und Python installieren

In einem Terminal das Repository klonen und hineingehen:

```sh
git clone https://github.com/maschin96/FEM-Toolbox.git
cd FEM-Toolbox
```

Unter macOS und Linux eine virtuelle Umgebung mit einem unterstützten Python
anlegen und aktivieren (gegebenenfalls `python3.13` statt `python3` wählen):

```sh
python3 --version
python3 -m venv .venv
source .venv/bin/activate
```

Unter Windows in **PowerShell**:

```powershell
py -3.13 --version
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Wenn nur 3.11 oder 3.12 installiert ist, die Zahl im `py`-Aufruf entsprechend
ändern. In PowerShell mit gesperrter Skriptausführung kann man ohne Aktivierung
direkt `.\.venv\Scripts\python.exe` anstelle von `python` für die folgenden
Befehle verwenden; eine dauerhafte Änderung der Execution Policy ist unnötig.

Nach Aktivierung gelten auf allen Plattformen dieselben Befehle:

```sh
python --version
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pip check
python -c "import fem_toolbox, meshio, numpy, scipy; print(fem_toolbox.__version__)"
```

`-e` installiert den lokalen Checkout editierbar; Änderungen am Python-Code
sind danach ohne Neuinstallation sichtbar. Das Terminal muss im Repository-Root
stehen. Für eine reine Laufzeitumgebung ohne Testwerkzeuge kann stattdessen
`python -m pip install -e .` verwendet werden.

## Gmsh installieren und prüfen

Die Toolbox liest vorhandene Gmsh-2.2/4.1-`.msh`-Dateien mit `meshio` **ohne**
Gmsh-Programm. Wer neue Netze erzeugt, installiert Gmsh getrennt in derselben
virtuellen Umgebung:

```sh
python -m pip install gmsh
gmsh -version
python -c "import gmsh; print(gmsh.__version__)"
```

Das offizielle Gmsh-Paket enthält Anwendung und SDK. Alternativ können die
passenden Binärdateien für Windows, Linux, macOS x86 oder ARM von der
[offiziellen Downloadseite](https://gmsh.info/#Download) bezogen werden.
Bei manueller Installation muss das Verzeichnis mit `gmsh` (`gmsh.exe` unter
Windows) im `PATH` liegen; `gmsh -version` prüft das. Für die Toolbox selbst
sind weder `GMSH_HOME` noch `PYTHONPATH` nötig. Das Python-Paket aus `pip` und
eine manuell installierte Gmsh-Version sollten nicht versehentlich gemischt
werden: mit `python -m pip show gmsh` und `gmsh -version` die gewählten
Versionen vergleichen. Die [Gmsh-Dokumentation](https://gmsh.info/doc/texinfo/)
beschreibt CLI, SDK und unterstützte Binärplattformen.

## Erster Lauf und Installationsprüfung

Vom Repository-Root aus nach der Installation:

```sh
python -m pytest -q
python examples/thermopiezoelectric/run.py --output install-check
python -c "import json; print(json.load(open('install-check/report.json'))['balances'])"
```

Der erste Befehl prüft das Projekt einschließlich numerischer Referenzfälle.
Der zweite importiert das mitgelieferte Gmsh-Netz, löst `T → (u, φ)` und
erzeugt `install-check/fields.vtu` und `install-check/report.json`. Der dritte
Befehl zeigt die Bilanzen des beispielspezifischen Berichts; für das allgemeine
Ergebnisbündel gilt [Doc/24](Doc/24_RESULT_BUNDLE.md). Für einen erneuten
Beispiellauf `--output` auf ein anderes Verzeichnis setzen, wenn die alten
Dateien erhalten bleiben sollen. Die in
[der Beispielbeschreibung](examples/thermopiezoelectric/README.md)
angegebenen Referenzwerte erlauben eine unabhängige Plausibilitätsprüfung.

Gmsh kann zusätzlich ohne GUI an einem Minimalnetz geprüft werden. In einem
temporären Verzeichnis `box.geo` mit folgendem Inhalt anlegen:

```geo
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 1};
Physical Volume("body") = {1};
Mesh.CharacteristicLengthMax = 1;
```

Dann `gmsh box.geo -3 -format msh2 -o box.msh` ausführen. Die erzeugte Datei
ist ein neues Tet4-Netz; mit
`python -c "from fem_toolbox.io import read_gmsh; print(read_gmsh('box.msh').n_tetra4)"`
muss eine positive Zellenzahl erscheinen. Dieses zusätzliche Beispiel braucht
Gmsh; der obige Solver-Schnelltest nicht.

## Häufige Installationsprobleme

| Symptom | Prüfung und Abhilfe |
|---|---|
| `python` meldet eine Version außerhalb 3.11–3.13 | Mit `python3.13` bzw. `py -3.13` die Umgebung neu anlegen. `python --version` **nach** Aktivierung prüfen. |
| `No module named fem_toolbox` | `python -m pip show fem-toolbox` und `python -m pip --version` im selben Terminal prüfen; Aktivierung oder Installation aus dem Repository-Root wiederholen. `PYTHONPATH` ist bei editierbarer Installation nicht nötig. |
| `pip` kann NumPy/SciPy nicht installieren oder versucht einen Quellbuild | Eine unterstützte reguläre CPython-Version und aktuelle `pip`-Version verwenden; prüfen, ob für Plattform/Architektur Wheels vorhanden sind. Unter Linux ist nur der in der CI geprüfte Plattformumfang zugesagt. |
| `gmsh: command not found` beziehungsweise `gmsh` wird in PowerShell nicht erkannt | Virtuelle Umgebung aktivieren und `python -m pip show gmsh` prüfen; bei manueller Installation den Gmsh-Programmordner zum `PATH` hinzufügen und das Terminal neu öffnen. |
| Gmsh-Import scheitert an Grafikbibliotheken auf einem Linux-Server | Der Toolbox-Import vorhandener `.msh`-Dateien benötigt Gmsh nicht. Für neue Netze eine passende headless Gmsh-Installation oder die offiziellen Binär-/SDK-Hinweise verwenden; die numerischen Toolbox-Tests getrennt ausführen. |
| `ModelSchemaError` bei eigenen `.msh`-Dateien | Physische Volumengruppen und explizite Längenskalierung prüfen. Das [Modellschema](Doc/21_MODEL_SCHEMA.md) und der [Gmsh-Adapter](Doc/04_SOFTWARE_ARCHITECTURE.md) erläutern die Anforderungen. |

Nach der Arbeit `deactivate` ausführen. Eine neue Shell benötigt die Aktivierung
der vorhandenen `.venv` erneut; `python -m pip install -e ".[dev]"` ist erst
nach Änderung der Projektabhängigkeiten erneut nötig.

## Nachweis dieser Anleitung

Am 23. September 2026 wurde der Checkout in einer **neu angelegten, zuvor
leeren** virtuellen Umgebung mit regulärem CPython 3.13 auf macOS/ARM64
installiert. Die Umgebung bezog alle Python-Abhängigkeiten per `pip` neu und
enthielt keine System-Site-Packages. `pip check` meldete keine Konflikte;
`python -m pytest -q` bestand mit **286 Tests**. Der oben genannte
thermopiezoelektrische Beispiellauf erzeugte VTU und JSON; seine Wärme-,
Kraft-, Moment- und Ladungsbilanz lag innerhalb der dokumentierten
Toleranzen. Die separate Gmsh-Installation lieferte Version 4.15.2; der
`box.geo`-Aufruf erzeugte ein Netz mit 1101 Tet4-Zellen, das `read_gmsh`
einlas. Dieser lokale Nachweis ergänzt die Linux-CI-Matrix und ist keine
Windows-Verifikation.

## Optionales PETSc-Backend

Die normale Wheel-/sdist-Installation funktioniert ohne PETSc. Native PETSc-
Backends werden separat eingerichtet; die getestete Linux-Umgebung und die
Auswahl real/complex stehen im [PETSc-Installationsvertrag](Doc/38_PETSC_BACKEND.md).
`Model`/TOML verwendet SciPy; PETSc wird an der direkten Python-Solver-API
ausgewählt. Fehlendes PETSc wird beim Backendaufbau ausdrücklich diagnostiziert.
