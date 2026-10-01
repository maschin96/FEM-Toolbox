# M8 / v1.3.0 — Freigabenachweis (#62)

Stand: 01.10.2026. **Freigabestand v1.3.0**.
Veröffentlichungsstatus, finaler Zielcommit und CI-Abschluss:
[GitHub Release v1.3.0](https://github.com/maschin96/FEM-Toolbox/releases/tag/v1.3.0)
und [Release-Issue #62](https://github.com/maschin96/FEM-Toolbox/issues/62).
Tag und GitHub Release werden erst nach Integration und erfolgreicher CI am
tatsächlichen main-Zielcommit erstellt. Dieser Bericht behauptet keine
experimentelle Validierung.

## Vollständiger Pflichtumfang

| Issue | Integrierter PR / Mergecommit | Vertrag und Verifikation |
|---|---|---|
| #59 Studien | [#212](https://github.com/maschin96/FEM-Toolbox/pull/212), `8b6da60` | [Manifest, Fingerprints, Resume](59_PARAMETER_STUDIES.md), `tests/unit/test_studies.py`, `tests/regression/test_parameter_study.py` |
| #60 Berichte | [#213](https://github.com/maschin96/FEM-Toolbox/pull/213), `db34a1d` | [Kennwerte, Herkunft, Fehlwerte und Qualitätsstatus](60_STUDY_REPORTS.md), `tests/unit/test_study_reports.py` |
| #205 Surrogate | [#214](https://github.com/maschin96/FEM-Toolbox/pull/214), `2962b04` | [Disjunkte Daten, Fehler, Persistenz und Extrapolation](61_SURROGATE_MODELS.md), `tests/unit/test_surrogates.py`, `tests/regression/test_surrogate_example.py` |
| #189 Feldbilder | [#215](https://github.com/maschin96/FEM-Toolbox/pull/215), `3b2571d` | [Tet4-Feldvertrag und JPEG](62_FIELD_VISUALIZATION.md), `tests/unit/test_visualization.py`, `tests/regression/test_field_visualization_example.py` |
| #61 Nutzerrezepte | [#216](https://github.com/maschin96/FEM-Toolbox/pull/216), `716c08a` | [Saubere Reproduktion und analytische Soll-/Ist-Werte](63_STUDY_RECIPE_ACCEPTANCE.md), `tests/regression/test_study_recipes.py` |

Alle fünf Issues sind geschlossen und ihre Beiträge in `codex/milestone-m8`
enthalten. Basis dieses Release-Issue-Branches: `716c08a39f5a099a01836c431cd9846946705c34`.
`origin/main` war beim Arbeitsbeginn bereits Vorfahr; `git merge origin/main`
meldete „Already up to date“. Vor Release-PR/Merge wird der Stand erneut geprüft.

## Kompatibilität, Grenzen und Lizenz

Additive neue APIs, keine entfernten öffentlichen Verträge und keine Änderung
numerischer Operatoren: SemVer **1.3.0** nach stabilem 1.2.0.
CPython 3.11–3.13, Model-Schemata 1–3, SI-/Vorzeichenkonventionen,
stationäre und dynamische Ergebnisformate sowie bisherige Tet4/Tet10-/PETSc-
Supportmatrix bleiben bestehen. Studien-/Tabellen-/Surrogatschema jeweils 1.
Versionsgebundene Studienfingerprints invalidieren frühere Ergebnisse bewusst.

Studien sind sequenziell, lokal und stationär, Model verwendet SciPy.
Feldbilder sind reelle stationäre Tet4-Außenflächen ohne komplexe Felder,
Tet10, Schnitt oder Zeitserie. Surrogate liefern empirische Testfehler und
keine Konfidenzintervalle/garantierte Fehlerschranken. Kein Optimierer,
keine Unsicherheits-/Adaptivitätsfunktion oder experimentelle Validierung.

MIT-Lizenz und Lizenzdatei werden in Wheel und sdist geprüft. Keine neuen
Core-Abhängigkeiten und keine vendorten Drittbibliotheken. Matplotlib ist
optional (BSD-kompatibel), NumPy/SciPy/meshio behalten ihre Lizenzen;
[Bestandsaufnahme](27_DISTRIBUTION_INVENTORY.md) und
[Adapterlizenz](62_FIELD_VISUALIZATION.md). Kein PyPI-Upload.

## Release-Verifikation

Versionstest Red: `tests/unit/test_package.py` erwartet 1.3.0 und scheiterte
an installierten 1.2.0-Metadaten. Green nach Änderung von `pyproject.toml`
und Neuinstallation: zwei Pakettests bestanden. Kein neues Solververhalten,
daher keine zusätzlichen künstlichen numerischen Red-Tests.

Die reguläre CI prüft vollständige Tests auf Python 3.11–3.13, Ruff,
Formatierung, mypy, älteste/aktuelle Abhängigkeiten, PETSc-real/complex sowie
sechs frische Wheel-/sdist-Installationen ohne Checkout im Prüfjob.
Alle bisherigen Beispiele, beide Studienrezepte, Wärme-Surrogat-Rundreise und
optionale JPEG-Bilder laufen aus diesen Artefakten. Die Core-Prüfung verlangt
explizit fehlendes Matplotlib/PETSc, anschließend wird das Extra separat installiert.

Die Rezeptprüfungen vergleichen nodale Felder und Kennwerte mit separaten
direkten Einzelmodellen sowie unabhängig hergeleiteten Tet4-/d₃₃-Erwartungen.
Vorab begründete Toleranzen bleiben unverändert. Der Ablauf prüft kontrollierten
Stopp zwischen Fällen, Wiederverwendung ohne neue Versuche, Materialfehler und
Basiseingabeinvalidierung. Kein neuer Nachweis eines harten Prozessabbruchs
während eines Solves. Bilanz-/Solverstatus und Netzkonvergenz werden getrennt.
Surrogattests prüfen zusätzlich Train-only-Skalierung, Rang, Extrapolation,
Artefaktkorruption und handberechnete Fehler; FEM-Parität ist keine Validierung.

Lokale Gesamtprüfung nach Versionswechsel: `PYTHONPATH=src .venv/bin/python -m pytest`
— **696 bestanden, 83 optionale PETSc-Tests ausgelassen**, 96,80 s.
Ruff, Formatprüfung, mypy (67 Module) und `git diff --check` bestanden.
PETSc fehlt lokal; native real/complex-Prüfung ist Teil der CI.

Neue Umgebungen: CPython 3.13.5, macOS 27.0 arm64, NumPy 2.5.3,
SciPy 1.18.1, meshio 5.3.5; Matplotlib 3.11.2/Pillow 12.3.0 erst nach
Core-Prüfung. Zwei neue Installationen aus dem 1.3.0-Wheel/sdist, außerhalb
des Checkouts, ohne editable Installation und mit `python -I`.
[Artefakthashes und Metadatenprüfung](evidence/issue62/artifacts.json).
Beide Rezeptpaare und beide Surrogatbeispiele Exit 0; pro Installation
sechs JPEGs, 1200×900 Pixel, Format JPEG geprüft. Temperaturbild zusätzlich
visuell auf Geometrie, Einheit und Colorbar geprüft.

```sh
uv build --out-dir /tmp/fem62-artifacts
uv venv /tmp/fem62-wheel-env
uv pip install --python /tmp/fem62-wheel-env/bin/python /tmp/fem62-artifacts/fem_toolbox-1.3.0-py3-none-any.whl
uv venv /tmp/fem62-sdist-env
uv pip install --python /tmp/fem62-sdist-env/bin/python /tmp/fem62-artifacts/fem_toolbox-1.3.0.tar.gz
mkdir -p /tmp/fem62-reproduction
cp -R examples/study_recipes examples/surrogate_heat examples/field_visualization /tmp/fem62-reproduction/
cd /tmp/fem62-reproduction
# Für beide Umgebungen; Ausgabeziele jeweils mit wheel-/sdist-Präfix:
/tmp/fem62-wheel-env/bin/python -I study_recipes/run.py heat --output wheel-heat
/tmp/fem62-wheel-env/bin/python -I study_recipes/run.py piezo --output wheel-piezo
/tmp/fem62-wheel-env/bin/python -I surrogate_heat/run.py --output wheel-surrogate
uv pip install --python /tmp/fem62-wheel-env/bin/python '/tmp/fem62-artifacts/fem_toolbox-1.3.0-py3-none-any.whl[visualization]'
/tmp/fem62-wheel-env/bin/python -I field_visualization/run.py --output wheel-images
```

Vollständige Soll-/Ist- und direkte Werte:
[Wheel Wärme](evidence/issue62/wheel-heat.json),
[Wheel Piezo](evidence/issue62/wheel-piezo.json),
[sdist Wärme](evidence/issue62/sdist-heat.json),
[sdist Piezo](evidence/issue62/sdist-piezo.json).
Surrogatberichte: [Wheel](evidence/issue62/wheel-surrogate.json),
[sdist](evidence/issue62/sdist-surrogate.json).
Temperatur-Testfehler 0 K, Leistungs-Maximalfehler 5,56e-16 W;
geladene Vorhersage bei Q=5 W/m³: 300,625 K und 5/6 W.
Diese Rundungswerte sind empirische Beobachtungen am Referenzfall.

Aktuelle CI-/Integrationslinks werden beim Abschluss im Release-Issue #62
festgehalten; das unten verlinkte Release enthält Zielcommit und gehashte
Artefakte aus der CI am tatsächlichen main-Stand.

## Review- und Integrationsstatus

Kein unabhängiger Reviewer verfügbar; keine unabhängige Freigabe behauptet.
Der verantwortliche Bearbeiter prüft Diff, Nachweise und aktuelle CI selbst.
Offene Review-Anmerkungen müssen vor Integration geklärt sein. Gemäß ADR-020
sind Rulesets kein zusätzliches Gate. Issue #62 bleibt bis Tag, geprüftem
main-Zielcommit und veröffentlichter Artefaktherkunft offen.
