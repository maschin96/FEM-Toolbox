# Benutzerhandbuch: Quellen und Kandidatenprüfung (#208)

Umfang: vorhandener M8-Funktionsstand und integrierte M9-Korrekturen,
Basiscommit `9daeffd`, Paketversion 1.3.0. Kandidat, keine M9-Release-Freigabe.
Dieser historische #208-Nachweis dokumentiert den damaligen vorgezogenen Stand.
Der Unsicherheitsworkflow #64 ist inzwischen im Kandidaten ergänzt; seine neuen
Nachweise stehen in Doc/67_UNCERTAINTY_WORKFLOW.md. Messvalidierung #63
wird weder durch synthetische Beispiele noch durch Surrogatfehler ersetzt.

## Quellen und Driftprüfung

`recipes.json` inventarisiert alle 63 ausführbaren Beispielinputs (Python,
TOML, Gmsh, CSV) mit SHA-256 und 22 vollständige Workfloweinträge. Die
byteidentischen Kopien unter `Doc/latex/examples/` machen das Quellenarchiv
unabhängig vom Checkout benutzbar. Alle Programme und deklarativen Modelle
sind im PDF vollständig eingebunden. Netzbytes stehen im Archiv bereit.

Die Kapitel verwenden kanonische Verträge: README/INSTALL, Doc/21/23/24,
Doc/33/34, Doc/38/39, Doc/40–55 sowie Doc/59–62. Schema-1-Dokumente sind
historische Teilverträge, ergänzt durch Schema 2/3 und Tet10; keine
implizite Migration wird behauptet. `coverage.json` enthält additive
`user_coverage`-Zuordnungen, ohne ursprüngliche Statuswerte zu überschreiben.

TDD: `test_user_handbook_sources.py` wurde vor Manifest/Kopien ausgeführt;
beide Tests scheiterten am fehlenden Rezeptmanifest (2 failed). Nach Anlegen
bestanden beide (2 passed). Sie sichern Inventar, Byteidentität, SHA-256,
Rezeptvollständigkeit und tatsächliche PDF-Einbindung. Fremde API-Änderungen
oder numerische Erwartungen wurden nicht eingeführt.

## Reproduzierbare Prüfung

Lokaler Lauf am 07.10.2026: neue virtuelle Umgebung unter `/tmp`, reguläre
Installation des gebauten Wheels, kein Editable-Install, kein Entwicklungs-
`PYTHONPATH`. Alle Python-Rezepte mit `-I`, `MPLBACKEND=Agg`; Importpfad
unter `site-packages` der neuen Umgebung ausdrücklich geprüft. Programm- und
Inputkopien in einem unabhängigen Arbeitsverzeichnis. Alle 22 Hauptbefehle
in `recipes.json` erfolgreich; zusätzlich zweimal Study-Resume, Piezo-Studie,
geladene Surrogatvorhersage, Runtime- und unabhängige Ergebnis-/Fehlerprüfung.
Jeder Befehl besitzt ein Budget von 90 s. Langsamster Hauptlauf war der erste
Matplotlib-Export (10.414 s einschließlich Fontinitialisierung); die übrigen
Hauptläufe lagen unter 1 s. Diese Zeiten sind Prüfbudgets, kein Leistungsbenchmark.

Artefakt: `fem_toolbox-1.3.0-py3-none-any.whl`, SHA-256
`bee582c33657938578516eb0ba6996e8b72b2a76432dca7c634dbce35e5c20ff`.
Build vom numerisch unveränderten M9-Basisstand `9daeffd`; nur die in diesem
PR dokumentierten Quellen/Tests wurden ergänzt. Der finale PDF-/CI-Build
identifiziert den tatsächlichen PR-Commit und alle Dokumentinput-Hashes.
Runtime: macOS 27.0.1 arm64, CPython 3.13.7; NumPy 2.5.3, SciPy 1.18.1,
meshio 5.3.5, Matplotlib 3.11.2, Pillow 12.3.0. `pip check` erfolgreich.
Die vollständig getrennten Wheel-/sdist-CI-Jobs unter Linux/Python 3.11–3.13
sind ein zweiter Runtime-Nachweis; ihr Resultat muss am aktuellen PR-Commit
geprüft werden. Eine zweite Person wird dadurch nicht behauptet.

| Workflow | Vorab festgelegte Erwartung | Beobachtung |
|---|---|---|
| Erstlauf Wärme | T=400−40x K, Abfluss 960 W; T-atol 3e−12 K, Bilanz 1e−10 W | Skript-Asserts bestanden; 320–400 K und 960 W |
| Thermopiezo | u_x=−9/440 m, phi=5/22 V, q_x=−20 W/m²; atol 1e−12 SI, rtol=0; Bilanzen 1e−10 SI | Bericht unabhängig gelesen und alle Größen/Bilanzen bestanden |
| Wärmezeitserie | Raster 0–10 s, endliche Temperaturen und zeitdiskrete Bilanz <1e−9 W | Endtemperaturen 296.5173832463 bis 296.5176819596 K; max. Bilanzdefekt 1.4734e−13 W |
| Studien-Resume | zwei Fälle, T_max=301/300.5 K (atol 1e−11 K); dritter Lauf keine neuen Versuche | beide succeeded, attempts jeweils 1; unveränderte Ergebnisse übernommen |
| Heat/Piezo-Lifecycle | analytische Felder, unabhängiger Direktlauf, Wiederaufnahme, Inputänderung und sichtbarer Fehlfall | beide `verification.json` erfolgreich, einschließlich Material-Negativfällen |
| Surrogat | separat geladene Q=5-Vorhersage [300.625 K, 5/6 W]; Einzelzielgrenzen 1e−10 K/1e−12 W | [300.625, 0.8333333333333333]; Skript prüft beide Held-out-FEM-Fälle und getrennte Toleranzen |
| Tet4-Moden | C=10 Pa, e=2 C/m², kappa=4 F/m, rho=6 kg/m³; omega=sqrt(55/3), u=sqrt(10), phi=sqrt(10)/2; rtol=1e−12 | separat gelöst: 4.281744192888376 rad/s, u=3.162277660168379, phi=1.58113883008419 |
| Tet10/Piezodynamik | unabhängige Feld-/Energie-/Moden-/diskrete Newmark-/komplexe Harmonikreferenzen mit den im eingebetteten Skript definierten SI-/rtol-Grenzen | alle neun Tet10-Programme und der Tet4-Dynamikworkflow erfolgreich |
| JPG | sechs Dateien, je 1200×900, feste gemeinsame Temperaturskala | alle Dimensionen separat per Pillow geprüft; Bildpaar im PDF enthalten |

Die Endtemperaturen des Aluminiumrezepts sind beobachtete numerische Werte,
keine vorher bekannten Messreferenzen und kein experimenteller Nachweis.
Patchtoleranzen sind bewusst an den jeweiligen kleinen analytischen Fall
gebunden; sie werden nicht zu einer universellen Anwendungstoleranz erklärt.

Nachgestellte Negativfälle: Leitfähigkeit −1, unbekanntes Modellschema 999,
entfernte elektrische Referenz, veränderte `values.csv` und Surrogatparameter
Q=11 außerhalb [0,10]. Konstruktion/Lesen/Vorhersage lehnten alle fünf Fälle
mit konkreter Diagnose ab. Tabellen-/Surrogatfehler sind `ValueError`;
Modellfehler sind `ModelSchemaError`. Kein Produktionsfehler wurde verdeckt
oder eine numerische Erwartung abgeschwächt.

Reproduktionsbefehle: Wheel bauen mit `uv build`, neue Umgebung mit
`python3.13 -m venv /tmp/fem-user-env`, Installation des Wheels inklusive
`[visualization]`, `pip check`. Ab `handbook/` alle `command`-Argumentlisten
und `continuation_commands` aus `recipes.json` in Reihenfolge ausführen;
vor `study-report` beide Resume-Aufrufe durchführen. `python` jeweils durch
den Interpreter der neuen Umgebung ersetzen. Ein neues Arbeitsverzeichnis
verwenden, da Exportziele nicht überschrieben werden. Die JPG-Darstellung
benötigt keine Änderung der Rohfelder. Die zwei mitgelieferten JPEGs sind
eigene, unretuschierte Ausgaben des Feldrezepts (Repository-MIT-Lizenz),
keine externen Fotos. Ihr visueller Ursprung und die gemeinsame Skala stehen
in der Bildunterschrift.

## Qualitäts- und PDF-Prüfungen

- Vollständige Suite: **773 passed, 83 skipped**, 85.40 s. Die übersprungenen
  Tests benötigen PETSc; die nativen real/complex-CI-Jobs sind separate Gates.
- `ruff check .`, `ruff format --check .` und `mypy src/fem_toolbox` bestanden
  (330 Python-Dateien formatiert, 67 typgeprüfte Kernmodule).
- Dokumentationsprüfungen nach Ergänzung der Matrix: **14 passed**.
- `git diff --check` und `build_handbook.py --check-matrix`: erfolgreich,
  **123** nachvollziehbar zugeordnete Quellen.
- `python3 scripts/build_handbook.py --self-test`: fünf PDFs und alle drei
  Negativfixtures erfolgreich; Quellenarchiv erneut entpackt und gebaut,
  Text-Roundtrip bestanden. Dokumentengine und vollständige Input-/Output-
  Hashes stehen im Buildmanifest.
- Benutzerband **69 Seiten**, Gesamthandbuch **104 Seiten**, keine Overfull-
  Meldungen in allen fünf Bänden. Die verbliebenen Theorie-Underfull-Hinweise
  sind sichtbar im Manifest; sie sind keine unterdrückten Inhaltsfehler.
- Alle 69 Benutzerseiten als Übersicht visuell geprüft; Tabellen,
  Vergleichsabbildung, lange Befehle, Codefortsetzung, Matrizen, Literatur und
  Register zusätzlich in lesbarer Auflösung. Gezielt Seiten 10, 14, 16, 19,
  38 und 63 geprüft. Kein abgeschnittener Inhalt; vollständige Quellen bleiben
  zum Kopieren im Archiv verfügbar. API- und V&V-Platzhalter sind weiterhin
  klar gekennzeichnet, kein Anspruch auf redaktionelle Gesamtabnahme #211.

Eigenprüfung durch den Autor. Unabhängige fachliche/redaktionelle Freigabe
ist offen; ein separat initialisierter Runtime ist keine zweite Person.
