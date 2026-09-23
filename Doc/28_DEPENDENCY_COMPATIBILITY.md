# Abhängigkeitsverträglichkeit und Updates für v1.0 (#144)

Stand: 23.09.2026. Offiziell geprüft wird der reguläre CPython-Build 3.11,
3.12 und 3.13 unter Linux. Die Paketmetadaten begrenzen die Installation auf
`>=3.11,<3.14`. macOS/ARM64 mit 3.13 wurde zusätzlich lokal geprüft; für
andere Plattformen und den experimentellen GIL-freien Python-Build besteht
keine eigene CI-Zusage. Das ist eine Testgrenze, keine Aussage, dass andere
Plattformen zwingend fehlschlagen.

## Installationsränder

Die CI-Jobfamilie `Compatibility` installiert je Python-Version in einer
frischen Umgebung zwei Kombinationen aus dem Paketquellstand: `oldest` mit
[`ci/oldest-constraints.txt`](../ci/oldest-constraints.txt) und `latest`
ohne Versionspins. `--only-binary=:all:` verlangt für externe Abhängigkeiten
verfügbare Wheels; es verhindert einen unbeabsichtigten lokalen NumPy- oder
SciPy-Quellbuild. `pip check`, ein Import **außerhalb** des Repository-Roots,
numerische Piezo-Benchmarks, das Gmsh-zu-VTU-Beispiel und Mesh-I/O-Tests sind
für beide Ränder Pflicht. `pip freeze` und die Paketmetadaten im Job-Log
halten die tatsächlich aufgelösten Versionen und `Requires-Dist` fest.

| CPython | Älteste geprüfte NumPy/SciPy/meshio-Kombination | Neueste Kombination |
|---|---|---|
| 3.11 | `1.26.0` / `1.12.0` / `5.3.5` | durch Resolver bei jedem CI-Lauf neu bestimmt und im Log ausgegeben |
| 3.12 | `1.26.0` / `1.12.0` / `5.3.5` | ebenso |
| 3.13 | `2.1.0` / `1.14.1` / `5.3.5` | ebenso |

Die Grenzen für 3.13 folgen aus den [NumPy-2.1-Release-Notizen](https://numpy.org/doc/2.2/release/2.1.0-notes.html)
und [SciPy-1.14.1-Release-Notizen](https://docs.scipy.org/doc/scipy-1.14.1/release/1.14.1-notes.html).
Die meshio-Untergrenze ist **gemessen**: Mit `meshio==5.3.0` und
`numpy==2.1.0` scheitern Import und alle Beispiele an `np.string_`, das in
NumPy 2 entfernt wurde. Die erste CI-Matrix zeigte den Fehler im Job
`Compatibility / Python 3.13 / oldest`; nach Wechsel auf `meshio==5.3.5`
liefen Wärme-, Piezo- und Gmsh/VTU-Smokes lokal durch. Ein unteres Limit in
`pyproject.toml` bedeutet nicht, dass beliebige selbst erzwungene Mischungen
kompatibel sind: Der Resolver muss zusätzlich die wechselseitigen
Versionsbedingungen der Pakete erfüllen.

## Paket- und Drittmetadaten

`pyproject.toml` deklariert die Python-/Laufzeitgrenzen und MIT mit einer im
Wheel/sdist enthaltenen `LICENSE`. `setuptools>=77` ist ausschließlich
Build-Abhängigkeit im isolierten Build-System. Die Laufzeitpakete sind NumPy,
SciPy und meshio; letztere zieht in der lokal geprüften Umgebung `rich`,
`markdown-it-py`, `mdurl` und `Pygments` nach. Das
[Lizenzinventar](27_DISTRIBUTION_INVENTORY.md) enthält Herkunft und Hinweise.
Die Toolbox-Wheels bündeln diese Fremdpakete nicht; ein späterer Installer mit
eingebetteten Wheels benötigt eine eigene vollständige Notice-Prüfung.

Für das Release-Issue [#28](https://github.com/maschin96/FEM-Toolbox/issues/28)
werden der grüne CI-Lauf, dessen sechs `Compatibility`-Jobs und die
Paketmetadaten des Zielcommits verlinkt. Das endgültige v1.0-Artefakt wird
erst im Release-Gate aus dem freizugebenden Commit gebaut und erneut geprüft.

## Sicherheitsupdates

`.github/dependabot.yml` prüft Python- und GitHub-Actions-Versionen
wöchentlich und eröffnet höchstens fünf Update-PRs je Ökosystem. Vor jedem
Release und nach einer relevanten Sicherheitsmeldung werden die für 3.11–3.13
aufgelösten Umgebungen mit einem aktuellen Advisory-Scanner (zum Beispiel
`pip-audit`) geprüft. Ein Fund wird mit Paket, Version, Advisory, betroffener
Funktion und verfügbarem Fix im **privaten** Release-Issue bewertet; Zugangsdaten
oder vertrauliche Exploit-Details werden nicht veröffentlicht. Kritische,
anwendbare Befunde blockieren den Release, bis Upgrade, Ersatz oder eine
begründete Ausnahme mit Auswirkung und Maßnahmen dokumentiert ist.

Bei dieser Prüfung waren Dependabot-**Sicherheitsalerts** für das private
Repository deaktiviert (API: HTTP 403); die wöchentliche Versionskonfiguration
ersetzt sie nicht. Ein lokaler `pip-audit`-Lauf am 23.09.2026 fand in der
Python-3.13-Mindestumgebung keine bekannten Befunde für die Laufzeitpakete;
das unveröffentlichte `fem-toolbox` wurde vom PyPI-basierten Scanner
übersprungen. Die ältere Entwicklungsumgebung enthielt `pip==25.2` mit
bekannten Advisories. `pip` gehört nicht zum Toolbox-Wheel; als Maßnahme
aktualisiert die CI das Installationswerkzeug vor dem Paketbezug auf
`pip>=26.2`. In der temporären Mindestumgebung wurde `pip==26.2.1` erneut
ohne bekannten Befund geprüft. Für die Python-3.11/3.12-Untergrenze wurden
NumPy 1.26.0, SciPy 1.12.0, meshio 5.3.5 und die lokal erfassten transitiven
Versionen zusätzlich als gepinnte Liste geprüft; auch dort meldete
`pip-audit` keinen bekannten Fund. Dieser zweite Scan löste keine
Abhängigkeiten neu auf und ersetzt daher keinen Scan der tatsächlich
installierten CI-Umgebung. Beide Ergebnisse sind Momentaufnahmen und kein
Ersatz für den erneuten Scan des Release-Kandidaten.

Bei einem Update wird zuerst die betroffene direkte oder transitive
Abhängigkeit aktualisiert. Danach laufen beide Ränder der CI-Matrix, die
gesamte Testsuite und der Artefaktcheck erneut. Ein oberes Limit wird nur
gesetzt, wenn eine konkrete neue Inkompatibilität reproduziert ist; ein
Sicherheitsfix allein ist kein Anlass für eine spekulative Obergrenze.
