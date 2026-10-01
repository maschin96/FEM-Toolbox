# Studienrezepte: Abnahme und Reproduktion (#61)

Stand: 01.10.2026. [Rezeptpaket](../examples/study_recipes/README.md),
[Wärme](../examples/study_recipes/heat/README.md),
[Piezo](../examples/study_recipes/piezo/README.md).

Die Abnahme verwendet zwei unabhängig **neu eingerichtete Umgebungen** auf
derselben Maschine und durch denselben Bearbeiter. Keine zweite Person und
kein unabhängiges fachliches Review werden behauptet. Das erfüllt die
Umgebungsalternative des Issues; Review vor Integration bleibt erforderlich.

## Identifizierter Kandidat und Befehle

Basiscommit des unveränderten numerischen Kerns:
`3b2571de22b985e8058ea34c7a267a9d42ee9d3d` (integrierter M8-Stand).
Rezeptskript und Model-/Netz-/Studieneingaben werden zusätzlich einzeln per
SHA-256 in [artifacts.json](evidence/issue61/artifacts.json) identifiziert.
Wheel/sdist stammen aus diesem Kernstand und den Issue-Rezepten; Dokumentation
wurde danach ergänzt. Version 1.2.0 bezeichnet hier einen unveröffentlichten
M8-Kandidaten, ausdrücklich nicht das stabile gleichnamige Release.

```sh
uv build --out-dir /tmp/fem61-artifacts
uv venv /tmp/fem61-wheel-env
uv pip install --python /tmp/fem61-wheel-env/bin/python /tmp/fem61-artifacts/fem_toolbox-1.2.0-py3-none-any.whl
uv venv /tmp/fem61-sdist-env
uv pip install --python /tmp/fem61-sdist-env/bin/python /tmp/fem61-artifacts/fem_toolbox-1.2.0.tar.gz
mkdir -p /tmp/fem61-reproduction
cp -R examples/study_recipes /tmp/fem61-reproduction/recipes
cd /tmp/fem61-reproduction
/tmp/fem61-wheel-env/bin/python -I recipes/run.py heat --output wheel-heat
/tmp/fem61-wheel-env/bin/python -I recipes/run.py piezo --output wheel-piezo
/tmp/fem61-sdist-env/bin/python -I recipes/run.py heat --output sdist-heat
/tmp/fem61-sdist-env/bin/python -I recipes/run.py piezo --output sdist-piezo
```

Alle vier Aufrufe Exit 0. Importpfade liegen in jeweiligem `site-packages`,
außerhalb des Checkouts; kein editable Install und kein PYTHONPATH.
CPython 3.13.5, macOS 27.0 arm64, NumPy 2.5.3, SciPy 1.18.1,
meshio 5.3.5. Die vollständige Plattformkennung steht in jedem JSON.

| Installation | Wärme (2 Fälle, 4 Knoten/1 Tet4) | Piezo (4 Fälle, 8 Knoten/6 Tet4) |
|---|---|---|
| Wheel | 0,0270 s | 0,0401 s |
| sdist | 0,0264 s | 0,0399 s |

Messgrenze: erster `Study.run(max_cases=1)`, Zwischenbericht, Fortsetzung und
unveränderte Wiederverwendung zusammen; keine Installation/Importzeit,
keine anschließenden direkten Referenz-, Fehler- oder Invalidierungsversuche.
Einzelmessung, keine Performancegarantie; Gesamtlaufbudget beider Rezepte
auf dieser Umgebung jeweils wenige Sekunden einschließlich aller Prüfungen.

## Numerische und operative Nachweise

Vorab begründete Toleranzen und analytische Sollwerte stehen direkt in beiden
Rezepten. Vollständige Soll-/Ist-Werte, separate direkte Model-Werte,
Bilanzwerte, Fall-IDs und Laufzeiten:

- [Wheel Wärme](evidence/issue61/wheel-heat.json)
- [Wheel Piezo](evidence/issue61/wheel-piezo.json)
- [sdist Wärme](evidence/issue61/sdist-heat.json)
- [sdist Piezo](evidence/issue61/sdist-piezo.json)

Jeder Lauf prüft alle Referenzkennwerte innerhalb seiner vorab festgelegten
absoluten/relativen Grenzen sowie die nodalen Felder gegen direkte Einzelmodelle.
`actual` enthält außerdem jede geprüfte globale Bilanz. Keine Anpassung der
Referenzen an gemessene Werte. Analytische Erwartungen sind unabhängig von
der Produktionsassembly; direkte Einzelmodelle prüfen zusätzlich Studienauflösung
und Bericht, verwenden aber bewusst denselben vorhandenen FEM-Solver.

Operativ: ein erfolgreicher Fall, übrige `pending` im Zwischenbericht;
Fortsetzung rechnet übrige Fälle, unveränderter Wiederholungskatalog identisch.
Ein zulässiger und ein absichtlich ungültiger Materialfall bleiben sichtbar.
Basismodelländerung hält IDs stabil, ändert Inputhashes, erhöht Versuchsanzahl
und rechnet alle betroffenen Fälle neu. Wiederherstellung rechnet die ursprünglichen
Eingaben erneut. Keine Behauptung eines getesteten harten Prozessabsturzes.

## Entwicklungsprüfung und Grenzen

`tests/regression/test_study_recipes.py` führt beide vollständigen Abläufe aus.
Red: beide Tests scheiterten am fehlenden Rezeptskript, nachdem die lokale
Importumgebung mittels `PYTHONPATH=src` korrigiert wurde. Der vorherige
`ModuleNotFoundError` war ein Umgebungsproblem und zählt nicht als Red-Nachweis.
Green: beide Tests bestanden; Formatierung danach ohne Verhaltensänderung.
Die sauberen Artefaktläufe oben benötigen keine solche Importkorrektur.

Es werden keine Solver-/Diskretisierungsoperatoren verändert. Daher keine
neuen allgemeinen Konvergenzstudien; bestehende wissenschaftliche Suite bleibt
maßgeblich. Diese Rezepte sind diskrete Wärme-/affine Piezo-Verifikation,
keine experimentelle Validierung. Plattformen Windows/Linux und Python 3.11/3.12
werden lokal nicht als reproduziert ausgewiesen; die reguläre CI prüft die
unterstützte Matrix. Die M8-Freigabe erfolgt separat über #62.

Lokale Gesamtprüfung: `PYTHONPATH=src .venv/bin/python -m pytest` —
696 bestanden, 83 optionale PETSc-Tests ausgelassen (PETSc lokal nicht vorhanden).
`python -m ruff check .`, `python -m ruff format --check .`,
`python -m mypy src/fem_toolbox` und `git diff --check` bestanden.
Die CI führt beide Rezepte zusätzlich aus übertragenen Wheel-/sdist-Artefakten
auf Python 3.11–3.13 ohne Checkout im Prüfjob aus.
