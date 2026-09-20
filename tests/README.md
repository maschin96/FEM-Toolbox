# Teststruktur

Die Tests spiegeln die Verifikationsstufen aus
[`Doc/05_VERIFICATION_AND_VALIDATION.md`](../Doc/05_VERIFICATION_AND_VALIDATION.md):

- `unit/`: Algebra, Datenmodelle und kleine deterministische Funktionen
- `patch/`: Element- und Patch-Tests physikalischer Felder
- `convergence/`: räumliche und später zeitliche Konvergenzstudien
- `regression/`: fachlich begründete Referenzkennwerte und Fehlerregressionen

Tests werden aus dem Repository-Root mit `python -m pytest` ausgeführt.
