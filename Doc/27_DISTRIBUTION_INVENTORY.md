# Weitergabe und Drittbestandteile (#143)

Stand: 23.09.2026. Diese Bestandsaufnahme betrifft den geplanten v1.0-Release.
Der Repository-Eigentümer hat für den aktuellen Quellstand externe Weitergabe
unter der [MIT-Lizenz](../LICENSE) gewählt (Entscheidung zu
[Issue #143](https://github.com/maschin96/FEM-Toolbox/issues/143)). Das
Repository bleibt durch diese Lizenzentscheidung privat, bis seine Sichtbarkeit
separat geändert wird; auch ein Paket-Upload ist damit nicht erfolgt.

## Paketinhalt und Herkunft

Das Repository ist derzeit privat. Ein Test-Build des vorherigen Stands mit
`uv build --out-dir /tmp/fem-143-baseline-dist` ergab ein Wheel mit 43 und
ein sdist mit 61 Einträgen. Beide enthielten **keine Lizenzdatei und keine
`License-Expression`**. Wheel und sdist enthielten keine Beispiel- oder
Testnetze. Die zwei `.msh`-Dateien unter `tests/fixtures` und
`examples/thermopiezoelectric` sind im Repository entstanden (Commits
`6ec4c9e` und `0c99c2b`); sie tragen keinen Drittquellenhinweis. Vor einer
externen Weitergabe ist die Rechtekette der Repository-Beiträge zu bestätigen.

Nach Aufnahme von `LICENSE`, `project.license = "MIT"` und
`project.license-files = ["LICENSE"]` ergab derselbe Build ein Wheel mit 44
und ein sdist mit 62 Einträgen. Beide enthalten den identischen Lizenztext
(SHA-256-Präfix `5ed34c5214b8`) und die Metadaten `License-Expression: MIT`
sowie `License-File: LICENSE`. Im Wheel liegt die Datei unter
`fem_toolbox-0.5.0.dist-info/licenses/LICENSE`, im sdist an der Wurzel.
Diese historische Prüfung verwendete Paketversion `0.5.0`; sie ändert frühere
veröffentlichte Archive nicht. Die finalen v1.0-Artefakte werden gemäß
[Release Notes](releases/v1.0.0.md) erneut geprüft.

## Direkt deklarierte Abhängigkeiten

Die folgenden Angaben stammen aus `pyproject.toml` und den lokal installierten
Distributionsmetadaten. Die unterstützten Versionen und der Update-Ablauf
stehen in [Doc/28](28_DEPENDENCY_COMPATIBILITY.md). Die Abhängigkeiten werden
beim Installieren bezogen;
ihre Dateien sind **nicht** Bestandteil des FEM-Toolbox-Wheels. Bei einer
Weitergabe eines gebündelten Installers sind die tatsächlich enthaltenen
Abhängigkeitsversionen und ihre Hinweise erneut zu inventarisieren.

| Paket | Rolle | Lokal geprüfte Version | Lizenzangabe der Distribution | Hinweis |
|---|---|---:|---|---|
| NumPy | Laufzeit | 2.4.6 | SPDX-Ausdruck `BSD-3-Clause AND 0BSD AND MIT AND Zlib AND CC0-1.0` | Distribution nennt mehrere Lizenzdateien; tatsächliches Wheel kann zusätzliche Komponenten enthalten |
| SciPy | Laufzeit | 1.18.1 | Copyright-/BSD-Text im `License`-Metadatenfeld | Binäre Bibliotheken und deren Notices je Ziel-Wheel prüfen |
| meshio | Laufzeit | 5.3.5 | MIT-Text, `LICENSE.txt` | Weitere transitive Abhängigkeiten je Installation prüfen |
| rich | transitiv über meshio | 15.0.0 | MIT-Metadaten und `LICENSE` | Lokal installiert, nicht Teil des Toolbox-Wheels |
| markdown-it-py | transitiv über rich | 4.2.0 | MIT-Klassifikator, `LICENSE` und `LICENSE.markdown-it` | Beide Hinweise bei gebündelter Weitergabe prüfen |
| mdurl | transitiv über markdown-it-py | 0.1.2 | MIT-Klassifikator | Kein `License-File`-Metadatenfeld; Notice bei Bündelung gesondert prüfen |
| Pygments | transitiv über rich | 2.21.0 | `BSD-2-Clause`, `AUTHORS` und `LICENSE` | Lokal installiert, nicht Teil des Toolbox-Wheels |
| pytest | Test/Entwicklung | 9.1.1 | MIT | Nicht Teil des FEM-Toolbox-Wheels |
| Ruff | Entwicklung | 0.16.8 | MIT | Nicht Teil des FEM-Toolbox-Wheels |
| mypy | Entwicklung | 2.3.1 | MIT | Nicht Teil des FEM-Toolbox-Wheels |

`setuptools>=77` ist das isolierte Build-Backend, keine deklarierte
Laufzeitabhängigkeit. Die vorliegende Bestandsaufnahme ist keine vollständige
Stückliste transitiver oder plattformspezifisch mitgelieferter Komponenten.
Für den finalen Release sind Wheel/sdist und die auf den Zielplattformen
installierten Abhängigkeiten erneut zu prüfen.

## Quellen für Paketmetadaten

Die [PyPA-Spezifikation](https://packaging.python.org/en/latest/specifications/pyproject-toml/#license)
beschreibt `project.license` als SPDX-Ausdruck und `project.license-files` als
Dateiliste. Ein Build-Backend muss passende Dateien in den Archiven ablegen;
das [Core-Metadata-Format](https://packaging.python.org/en/latest/specifications/core-metadata/#license-expression)
unterscheidet `License-Expression` und `License-File`. Die konkrete
Metadatenprüfung muss an **beiden** erzeugten Archiven erfolgen.
