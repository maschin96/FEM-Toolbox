# Versionierter Nachweiskatalog (#252)

14.10.2026 · Schema 1 · M15-Implementierungskandidat, keine neue Physikfreigabe.

`Doc/evidence/catalog/catalog.json` verbindet Fähigkeit, Analyse-/Elementumfang,
Kennwertkriterium, Commit, freigegebene Version, exakte Nachweisdateien,
Gültigkeitsgrenze und Reviewstatus. Die 17 Gruppen des veröffentlichten
v1.4.0-Handbuchs stehen neben dem M10-Wärmekandidaten, den ausschließlich im
Branch vorhandenen Quad-FV-Grundoperatoren und der weiterhin ungeprüften
experimentellen Piezo-Fallvalidierung. Ein geschlossenes Issue ist kein
Release-Nachweis. Die Quad-FV-Eintragung belegt keinen SIMPLE-Strömungssolver.

## Prüfung und Lesesicht

Aus jedem Arbeitsverzeichnis mit dem absoluten Pfad zum Skript:

```sh
python scripts/capability_catalog.py
python scripts/capability_catalog.py --format json
```

Die Standardausgabe ist eine Markdown-Tabelle mit Kennwertkriterien, Referenzen,
Umfang, Grenzen und Reviewstatus. Ihre Quelllinks zeigen auf den jeweiligen
unveränderlichen Commit. `--catalog` wählt eine alternative JSON-Datei;
relative Pfade beziehen sich auf die Repositorywurzel. Erfolgreiche Prüfung
liefert Exitcode 0; Schema-, Datei-, Hash- oder Provenienzfehler Exitcode 2 mit
konkretem Kontext. `missing_evidence` ist ein sichtbarer fachlicher Zustand,
keine Ausnahme und keine positive Prüfung.

## Evidenzvertrag

- `schema = 1`; unbekannte Formatversionen werden abgelehnt. Fähigkeiten besitzen
  eindeutige IDs und nichtleere Kriterien, Grenzen und Reviewerangaben.
- `released` verlangt eine vollständige Zuordnung zum archivierten Handbuch.
  Version und Commit müssen mit dem unveränderten Release-Manifest und der
  Release-Verifikation übereinstimmen; die damalige CI muss am selben Commit
  erfolgreich gewesen sein. Der Capability-Map-Hash stammt aus dem Manifest.
  Das archivierte Quellpaket muss zum veröffentlichten Artefakthash und seiner
  eigenen Paketmetadatenversion passen. Fehlende oder zusätzliche Release-Gruppen
  und gegenüber dem damaligen Handbuch erweiterte Umfänge werden abgelehnt.
- `branch_candidate` besitzt einen vollständigen Quellcommit und exakte
  Nachweisbytes, jedoch weder `released_version` noch eine Release-Zuordnung.
  Die Branch-Snapshots tragen ihren ursprünglichen Pfad; historische Dateien
  werden weder umgeschrieben noch durch heutige Messungen ersetzt.
- `planned` mit `missing_evidence` bildet fehlende Messdaten ausdrücklich ab.
  Schema 1 erlaubt ausschließlich `experimental_validation = not_performed`
  und `review_status = not_available`; es liegt kein unabhängiger Reviewbeleg vor.
  Eine zukünftige experimentelle Freigabe braucht einen separat beschlossenen
  Messdaten-/Holdout-Vertrag und eine neue Formatversion.
- Jede lokale Nachweisdatei muss innerhalb des Repositorys liegen und ihren
  SHA-256 erfüllen. Das erkennt Drift und falsche Zuordnungen; es ist keine
  kryptographische Signatur oder unabhängige fachliche Prüfung der Aussagen.
  Bei neuen Quellständen sind bisherige Release-Referenzen zu archivieren und
  neue Einträge/Versionen ausdrücklich zu reviewen, statt alte Hashes still zu
  aktualisieren. Änderungen des Katalogs erfordern einen geprüften PR.

Die Dateien `manifest.json`, `RELEASE_VERIFICATION.json` und das Quellpaket
wurden bytegenau aus dem GitHub-Release v1.4.0 übernommen; ihre Hashes entsprechen
den veröffentlichten Asset-Digests. Der historische Capability-Map-Snapshot stammt
aus dessen Commit und stimmt mit dem Manifest überein. Kandidatennachweise
stammen aus den im Katalog genannten Commits. Die Referenztests sind Nachweisquellen;
das Lesen des Katalogs führt sie nicht aus und erzeugt keine neue Numerikabnahme.

## Verifikation und Grenzen

Deterministische Prüfungen unterscheiden den freigegebenen Wärmeumfang, den
unveröffentlichten CFD-Operatorbranch, den M10-Kandidaten und den Messfall ohne
Daten. Negative Fälle prüfen falsche Version/Commit, fälschliche Promotion eines
Branches, fehlende/extern liegende Artefakte, falsche Hashes, erweiterten Umfang,
fehlende Release-Gruppen und CLI-Fehlerstatus. Die damaligen Release-/Handbuch-
Manifeste bleiben unverändert und benötigen keine Migration.

Keine neue FEM-Numerik, Messdatenvalidierung oder unabhängige Reviewerfreigabe.
Der operative Abnahmestand bleibt von den tatsächlichen Self-hosted-CI-Ergebnissen
des aktuellen Issue-PR abhängig; lokale Katalogtests ersetzen diese nicht.
