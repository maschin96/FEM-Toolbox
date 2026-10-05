# Femantis — kanonischer Namespace und Altimporte

05.10.2026 · [Issue #237](https://github.com/maschin96/FEM-Toolbox/issues/237) · Paket B

Status: auf dem M14-Issue-Branch implementiert und lokal geprüft; keine
veröffentlichte Distribution und keine vollständige Namensmigration.
[Vertrag/ADR-053](77_FEMANTIS_MIGRATION_CONTRACT.md) ·
[Gesamtplan](76_FEMANTIS_RENAMING_PLAN.md)

## Öffentlicher Importvertrag

Die Implementierung liegt ausschließlich unter `src/femantis/`.
Der aktuelle Entwicklungsstand unterstützt:

```python
from femantis.mesh import Mesh
from femantis.model import Model
from femantis.errors import ValidationError
```

Vorhandene Imports bleiben während aller v2.x-Versionen erreichbar:

```python
from fem_toolbox.mesh import Mesh
from fem_toolbox.model import Model
from fem_toolbox.errors import ValidationError
```

Die beiden Pfade liefern dieselben Module, Klassen und Exceptions. Auch
`fem_toolbox` selbst ist nach dem Import dasselbe Modul wie `femantis`;
`__name__`, `__package__` und `__spec__` bleiben kanonisch. Daher zeigt
`fem_toolbox.__name__` künftig `femantis`. Die Klasse `FEMToolboxError`
behält ihren fachlichen API-Namen, Typen und Konstruktoren bleiben kompatibel.

Der Legacy-Rootimport erzeugt einmal pro Prozess eine `DeprecationWarning`
mit Verweis auf `femantis` und die Unterstützung während v2.x. Python blendet
diese Warnung nach seinen normalen Regeln gegebenenfalls aus; die Bibliothek
ändert keine Warnungsfilter. Der kanonische Import erzeugt diese Warnung nicht.

Rootimporte laden keine NumPy-/SciPy-/meshio-/PETSc-/Matplotlib-Module.
Fachmodule laden ihre erforderlichen Abhängigkeiten erst bei Bedarf. Der
optionale PETSc-Adapter bleibt ohne installierte PETSc-Laufzeit importierbar.

## Umsetzung und Typdaten

Ein eng auf `fem_toolbox.*` begrenzter Importfinder leitet Untermodule auf das
entsprechende `femantis.*`-Modul um. Er lädt nur das angeforderte Modul und
führt dessen Code nicht ein zweites Mal aus. Die Import-Spec wird nach der
Aliasinitialisierung auf die kanonische Spec zurückgesetzt, damit Reload und
Importmetadaten korrekt bleiben. Andere Namespaces werden nicht behandelt.

Die alten Pfade enthalten ausschließlich den Rootadapter, `py.typed` und
`.pyi`-Weiterleitungen. Es gibt keine zweite Solverimplementierung. Die
Weiterleitungen verweisen auf kanonische Typdefinitionen; beide Namespaces
werden mit mypy einschließlich eines gemischten Anwenderbeispiels geprüft.
Neue Module aus der parallelen SPICE-Entwicklung erhalten bei späterer
Integration ebenfalls passende Typweiterleitungen und Importtests.

Verbleibende Altzeichenketten im kanonischen Code sind bestehende
Datenschnittstellen, Solver-Provenienzkennungen oder backendinterne Kennungen.
Diese bleiben in Paket B unverändert. Distribution und Versionsmetadaten
heißen vorerst `fem-toolbox` / `1.3.0`; Paket C behandelt den Wechsel der
Distribution. Die lokal gebauten 1.3.0-Artefakte sind ausschließlich Prüfarbeits-
produkte und werden nicht als Ersatz für das veröffentlichte Release hochgeladen.

## Nachweise

- Red: Der kanonische Rootimport scheiterte nachweislich mit
  `ModuleNotFoundError: No module named 'femantis'`.
- Green: Nach Verschiebung und Importumstellung bestand der Rootimporttest.
- Red für Adapter: Beide Importreihenfolgen scheiterten am fehlenden Altpackage.
- Green: Alle 67 inventarisierten Modulpaare sind in beiden Importreihenfolgen
  identisch; kanonische Namen, Spec und Elternattribute stimmen.
- Einmalige Warnung, unbekannte Module, Lazy-Rootimporte und Reload sind geprüft.
- Analytischer Wärmelauf über beide Namespaces: Unit-Axis-Tet4 mit k=2 W/(m K),
  q=8 W/m³ und drei Knoten bei 300 K ergibt am freien Knoten 301 K,
  Quellleistung 4/3 W und entgegengesetzte Reservoirleistung. Für den kleinen,
  gut konditionierten SI-Fall gilt atol=1e-12 K bzw. W, rtol=0.
- VTU-Ausgabe behält `fem_toolbox:external_id` und `fem_toolbox:set:body`;
  IDs bleiben beim Roundtrip erhalten.
- Vollständige lokale Suite: **742 bestanden, 83 übersprungen** in 95,71 s
  auf CPython 3.13.7. Die übersprungenen Tests benötigen die lokal fehlende
  PETSc-Laufzeit; sie werden durch gesonderte reale/komplexe PETSc-CI geprüft.
- Ruff-Linting und Formatprüfung: bestanden; mypy für kanonischen Code,
  Legacy-Stubs und Anwendervertrag: 136 Quelldateien ohne Befund.
- Wheel und sdist isoliert gebaut, getrennt ohne Änderung der gemeinsamen
  Umgebung installiert: alle 67 Modulpaare, installierte Pfade, Version,
  `py.typed` und Typweiterleitungen vorhanden.

Befehle aus dem Repository-Root mit eingerichteter Umgebung:

```sh
python -m pytest tests/unit/test_namespace_compatibility.py tests/unit/test_package.py tests/unit/test_petsc_optional.py -q
python -m pytest -q
python -m ruff check .
python -m ruff format --check .
python -m mypy src/femantis src/fem_toolbox tests/typing
git diff --check
uv build --out-dir /tmp/femantis-237-artifacts
```

Die lokale vorhandene Entwicklungsumgebung wurde für pytest mit `PYTHONPATH=src`
auf diesen Worktree gerichtet. Artefaktprüfungen verwendeten den Interpreter
im isolierten Modus und prüften explizit, dass Module aus den temporären
Installationen stammen. Laufzeitabhängigkeiten wurden dabei aus der bestehenden
Umgebung genutzt; vollständige frische Plattform-/Backendabnahme erfolgt in CI.

## Parallele Integration

Am 05.10.2026 hat der Auftraggeber die Weiterarbeit parallel zum SPICE-
Meilenstein beauftragt. Damit wird der M14-Entwicklungsbeginn gegenüber der
ursprünglichen Reihenfolge vorgezogen: `codex/milestone-m14` basiert auf dem
aktuellen stabilen main; #237 integriert per PR ausschließlich dorthin.
Die SPICE-Branches bleiben eigenständig. Benötigte Quellverschiebungen und
Imports werden beim späteren Zusammenführen ausdrücklich abgeglichen;
SPICE-Code wird nicht vorab übernommen oder als freigegeben dargestellt.

Die vollständigen M14-/v2.0-Releasegates und die Nachweise für Registry-/Bild-/
Markenrechte bleiben erforderlich. Dieser Beginn ist keine Freigabe zum
Release-Merge, Taggen oder Veröffentlichen. Unabhängiges Review und aktuelle
PR-CI werden separat am tatsächlichen Commit beurteilt.
