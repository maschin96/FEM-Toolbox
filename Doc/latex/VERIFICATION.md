# Abnahmeprotokoll der Buchstruktur (#206)

Prüfstand: 06.10.2026; Strukturkandidat auf stabilem M8-/v1.3.0-Stand.
Die exakte Version, Commit-ID, Quellzustand und SHA-256 der tatsächlich gebauten
Quellen werden in jedem Kandidatenmanifest erfasst. Dieses Protokoll ist kein
Releasebericht und keine unabhängige Freigabe.

## Ausgeführte Nachweise

- Red: sechs fatale Logdiagnosen, Layoutdiagnose und Matrixvalidierung scheiterten
  am noch nicht implementierten Verhalten (acht Tests); der Quellenmatrixtest
  meldete die noch fehlende Matrixdatei. Nach Implementierung neun grüne Tests.
- Zweiter Red-Schritt: ein ungültiges Ziel eines geplanten M9-Themas wurde noch
  nicht abgelehnt (1 fehlgeschlagen, 11 bestanden); nach Korrektur alle 12 grün.
- Vollständige Python-Suite: 771 bestanden, 83 optionale PETSc-Tests lokal
  ausgelassen. PETSc ist hier nicht installiert; CI prüft beide Scalar-Varianten.
- Ruff-Lint, Ruff-Format und mypy auf `src/fem_toolbox` sowie dem neuen
  Buildwerkzeug bestanden. Diff-/Quellenprüfung vor Commit.
- Container-PDF-Build: Gesamtdokument und alle vier Bände erfolgreich.
  Musterband: Gleichungen mit korrekt negativem äußeren Wärme-Randterm,
  Tabelle, selbst erstellte Abbildung, Python-Code, Querverweise und Zitat.
- Drei tatsächlich beschädigte LaTeX-Kopien scheitern für fehlende Datei,
  fehlende Referenz und fehlenden Literatureintrag. Diagnose-Logs werden als
  Artefakte erhalten; die positiven PDFs haben keine ungelösten Verweise/Zitate.
- PDF-Text/-Metadaten: Version, vollständiger Commit, Kandidatenstatus,
  deutsche Umlaute, Literatur, Sachregister und eingebettete Unicode-Schriften
  automatisiert geprüft. Lesezeichenstruktur ist vorhanden.
- Quellenarchiv in frischem Verzeichnis entpackt; alle fünf PDFs ohne zusätzliche
  Projektquellen gebaut. Der extrahierte PDF-Inhalt ist identisch zum positiven
  Repository-Build. Kein Byteidentitätsanspruch.
- Visuelle Eigenprüfung aller neun Musterbandseiten: Titel, Verzeichnisse,
  Symbol-/Abkürzungstabelle, Formeln, Grafik, Tabelle, Code, Bibliografie und
  Index lesbar, ohne abgeschnittene Inhalte. Die gemeinsamen Seitenvorlagen
  werden in allen Bänden verwendet. Layoutwarnungen bleiben im Manifest sichtbar.

## Grenzen und spätere Gates

Der Musterteil ist gegen Doc/02 auf Vorzeichen und SI-Einheiten geprüft.
Ein vollständiges mathematisches Fachreview und die Handbuchinhalte gehören zu
#207–#210. Keine neue Physik, keine Änderungen numerischer Toleranzen oder
Referenzdaten; Konvergenztests für neue Operatoren sind nicht anwendbar.

Die initiale Quellenmatrix ordnet den Bestand zu; sie behauptet keine bereits
abgeschlossene Übertragung in die Bände. #63–#65 bleiben geplante Inhalte.
CI-Nachweis am aktuellen PR-Commit und ein unabhängiges Review sind vor
Integration gesondert zu prüfen. Ein unabhängiger Reviewer steht bislang nicht
zur Verfügung; die visuelle Prüfung ist eine Eigenprüfung.
