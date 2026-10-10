# M9: erneute Messdatensuche und entscheidbarer Releaseumfang (#63 / #66)

Stand: 10. Oktober 2026. **Recherchebefund und genehmigte Neuplanung; keine Freigabe.**
Der Eigentümer hat erneut bestätigt, dass keine eigenen Messdaten vorliegen,
und eine weitere Internetsuche sowie gegebenenfalls eine alternative Planung
beauftragt. Der Eigentümer hat anschließend **Weg B ausdrücklich genehmigt**:
„M9 wie vorgeschlagen neu planen“. Der
[Messdatenvertrag](65_VALIDATION_DATA_CONTRACT.md) gilt weiter für den offenen
experimentellen Folgeumfang; sein No-go blockiert nach der Neuplanung nicht
mehr die generische Software-/Dokumentationsfreigabe v1.4.0.

## 1. Ergebnis und Prüfgrenze

Die erneute Suche hat zusätzliche echte, öffentlich referenzierbare Messdaten
gefunden. Es wurde aber **noch kein Datensatz gefunden und angenommen, der alle
M9-Metadaten-, Unabhängigkeits- und Unsicherheitsgates erfüllt**. Datenzugang und
eine passende Datenlizenz sind notwendig, reichen dafür allein nicht aus.
Das Ergebnis behauptet weder weltweite Nichtverfügbarkeit noch die endgültige
Untauglichkeit sämtlicher Kandidaten.

Gesucht wurde nach quasistatischer Spannungs-/Verschiebungsantwort, d33,
statischer Kraft-/Ladungsantwort und als möglicher Umfangsalternative nach
Bimorph-/Quarz-Fällen. Quellen: Zenodo, Figshare, Mendeley, institutionelle
Archive, Primärpublikationen und ihre Datenverweise. Suchbegriffe umfassten
`piezoelectric quasistatic d33 displacement voltage experimental raw data`,
`piezoelectric static displacement dataset`, `piezoelectric bimorph experimental
dataset` und `quartz tuning fork static dataset`.

Vertieft geprüft: Repository-API-Metadaten für vier Datensätze, zwei kleine
CSV-Rohdateien, eine TXT-Rohdatei und die Blatt-/Textstruktur eines XLSX-Archivs.
Keine fremden Auswertungsskripte wurden ausgeführt, keine Kurven digitalisiert,
keine Materialwerte gefittet und keine Validierungsantworten mit FEM verglichen.
Die Sichtung ist eine Beschaffungsprüfung vor Fallannahme, keine abgeschlossene
statistische oder experimentelle Validierung.

## 2. Neue und nachgeprüfte Kandidaten

| Primärquelle | Nachgewiesen / Sichtung | Noch fehlend oder außerhalb des Vertrags |
| --- | --- | --- |
| [Green Piezo, Zenodo 7920056 v1](https://zenodo.org/records/7920056), zu [Monroe et al., 10.1002/admt.202300745](https://doi.org/10.1002/admt.202300745) | Die API bestätigt **CC BY 4.0**. Design-, Charakterisierungs- und Kraft-/Spannungsdateien werden beschrieben. | Vollständige unabhängige Materialtensoren, Unsicherheiten und quasistatischer Messschaltungsvertrag weiter nicht nachgewiesen. Das 776-MB-Archiv wurde nicht vollständig geladen/auditiert; Lizenzklärung ist keine Datenannahme. |
| [Light-Annealed Piezoelectric Films, Zenodo-Record 15526822](https://zenodo.org/records/15526822), [Cardoletti et al., Primärartikel](https://pmc.ncbi.nlm.nih.gov/articles/PMC12412497/) | **CC BY 4.0**, tatsächliche Feld-/Verschiebungsdatei `Figure_4d_540nm.txt`, Einheitenheader und Messfrequenz. Die Publikation beschreibt ein Glas/PZT-Cantilever mit interferometrischer Messung. | Die gesichtete Datei nennt **100 Hz**; das allein belegt keinen quasistatischen Zustand. Der Artikel leitet den effektiven Koeffizienten aus der gemessenen Verschiebung ab: diesen Wert wieder als Eingabe zu verwenden wäre keine unabhängige Prüfung derselben Antwort. Volltensor-, Unsicherheits- und unabhängiger Laufvertrag nicht nachgewiesen. |
| [Polaron-defect complexes, Zenodo 21454715 v2](https://zenodo.org/records/21454715), [Studienmitteilung der beteiligten Universität](https://mse-en.xjtu.edu.cn/info/1011/2134.htm) | **CC BY 4.0**. `Fig3-A displacement.csv` enthält räumliche Koordinaten und Verschiebung; `figS2-d33-U-f.csv` enthält Bias-, Amplituden- und Frequenzabschnitte mit Fehler-Spalten. | Die Universitätsquelle beschreibt eine durch DC-Bias aktivierte, elektrostriktive Kilohertz-Antwort. Die gesichteten CSVs liefern keinen vollständigen unabhängigen linearen quasistatischen Probenvertrag; Bedeutung der Fehler-Spalten bleibt unqualifiziert. Kein Zugriff auf den vollständigen Science-Artikel behauptet. |
| [Relaxor-Kristalle, Figshare 19346039 v1](https://doi.org/10.6084/m9.figshare.19346039.v1), [Yang et al., 2022](https://www.nature.com/articles/s41467-022-29962-6) | **CC BY 4.0**. Tatsächliches `Source Data.xlsx` enthält 14 Abbildungsblätter mit Material-/Polarisations-/Permittivitätscharakterisierung. | Die gesichtete Blatt-/Textstruktur und der Artikel weisen keine vollständige separate statische Proben-Verschiebungsreihe mit allen M9-Pflichtmetadaten nach. Materialcharakterisierung allein ist kein unabhängiger Gerätevergleich. Keine vollständige Auswertung sämtlicher Zellen behauptet. |
| [Quarz-Stimmgabel, Castellanos-Gomez et al., arXiv:1203.5767v1](https://arxiv.org/abs/1203.5767v1), [Volltext](https://arxiv.org/pdf/1203.5767) | Wissenschaftlich interessanter **linearer statischer** Alternativfall; Publikation beschreibt eine Kalibriermethode und Unsicherheit. | In der gesichteten Publikation kein qualifiziertes separat lizenziertes Rohdatenarchiv gefunden; vollständige Probengeometrie, Elektroden-/Kristallachsen- und unabhängiger Tensorvertrag nicht angenommen. Zudem Biegungs-/Kalibrierfall statt gewähltem monolithischem d33-Fall. |
| [Duval et al., 2007, Southampton](https://eprints.soton.ac.uk/263266/) | Primärnachweis von Verschiebungs- und Blockierkraftmessungen an PZT-Mikrocantilevern. | Repository-Volltext ist ausdrücklich auf Mitarbeiter eingeschränkt; zugängliche Rohdateien und ihre Nutzungsrechte wurden nicht nachgewiesen. Der Abstract genügt keiner Fallannahme. |
| [Gesing, Figshare 5766708](https://figshare.com/articles/dataset/Dataset_for_Figure_3/5766708) | Repositorybeschreibung nennt gemessene und simulierte Ladung pro Beschleunigung und **CC BY 4.0**. | Frequenzabhängiger Beschleunigungssensor; kein statischer d33-Verschiebungsfall. API-/Rohdateizugriff wurde in dieser Recherche nicht erfolgreich geprüft. |

Die [sechs älteren Kandidaten](65_VALIDATION_DATA_CONTRACT.md#5-geprüfte-öffentliche-kandidaten)
bleiben mit ihren damaligen Prüfgrenzen dokumentiert. Numerische oder DFT-Daten
werden auch bei offenem Download nicht als Messdaten angenommen.

## 3. Nachvollziehbare Beschaffungssichtung

[Die Evidenzliste](evidence/issue63/research-2026-10-10.json) hält abgerufene
URLs, Version/Lizenz, Dateigröße und SHA-256 der tatsächlich gesichteten Dateien
fest. Die kleinen, unveränderten Text-Rohdateien liegen im
[Sichtungsverzeichnis](evidence/issue63/README.md) mit Zuordnung und Attribution.
Das XLSX wird aus Platzgründen referenziert statt ins Repository kopiert.
API-Snapshots enthalten den relevanten Identitäts-/Lizenz-/Dateiausschnitt;
ihr Hash und der Hash der vollständigen abgerufenen Antwort sind getrennt.
Kein Snapshot wird als eingefrorener angenommener Validierungsbestand bezeichnet.

## 4. Drei tatsächlich unterschiedliche Wege

**A — Bestehenden M9-Umfang behalten.** Dokumentation weiterentwickeln,
v1.4.0 zurückhalten, fehlende Daten/Metadaten beschaffen. Bei Green Piezo,
Glas/PZT oder Quarz wäre die nächste Aktivität eine gezielte Klärung mit den
Autorinnen/Autoren. Eine Anfrage setzt einen separat autorisierten Kontakt
voraus; in dieser Recherche wurde niemand angeschrieben. Erfolg und Dauer
lassen sich nicht garantieren.

**B — M9 ausdrücklich neu planen (am 10.10.2026 genehmigt).**
v1.4.0 erhält den Titel **„Unsicherheitsworkflow und vollständige Dokumentation“**.
Die experimentelle Fallfreigabe wird als weiterhin offener Folgeumfang geführt.
Numerische Unsicherheitsfortpflanzung wird mit analytischen Referenzen und
synthetischen Lehrfällen verifiziert. Alle vier Fachbände dokumentieren diesen
tatsächlichen Stand; sie geben keinen bestandenen Messvergleich vor.
Eine solche Version ist fachlich sinnvoll und als kompatible Erweiterung
grundsätzlich mit einer Minor-Version vereinbar. Sie wäre jedoch eine **neue
Umfangsentscheidung**, nicht das Bestehen des bisherigen M9-Gates.

**C — Anderen realen Validierungsfall wählen.** Ein Quarz-/Bimorphfall oder
ein anderer unterstützter Einzelfeldfall könnte leichter beschaffbar sein.
Vor Umsetzung braucht er einen neuen vorab festgelegten Versuch, dieselben
Herkunfts-/Unsicherheits-/Unabhängigkeitsgates und ein positives Datengate.
Diese Recherche hat noch keinen solchen Ersatzfall qualifiziert; C ist daher
keine sofortige Release-Abkürzung.

## 5. Konkreter Änderungsumfang für B

Die folgende Liste ist der genehmigte Änderungsumfang. Umsetzung und Release
benötigen weiterhin die vollständige Abnahme; die Paketversion wird erst im
Releasepaket angepasst.

| Bestehendes Paket | Nach einer Zustimmung umzusetzen |
| --- | --- |
| M9 / #66 | Titel und Pflichtumfang auf Unsicherheitsworkflow plus vier vollständige Fachbände ändern. Expliziter Satz: keine experimentelle Validierungsfreigabe in v1.4.0. SemVer und Abnahme erneut prüfen. |
| #63 / #65 | Offen als experimenteller Folgeumfang erhalten; Recherchebefund/No-go dokumentieren. Nicht als erfolgreich schließen. Zuordnung zum späteren Horizont vor Umsetzung explizit entscheiden. |
| #64 | Abhängigkeit vom angenommenen Experiment für die **generische Software** lösen; explizite nutzerbegründete Intervalle/Verteilungen, begrenzter Korrelationsvertrag, reproduzierbares Sampling, Sensitivitäten, getrennte Fehlerbudgets und Provenienz implementieren. Kein automatisches Fitten oder erfundene Verteilungen für reale Materialien. |
| #64-Abnahme | Analytischen linearen Referenzfall, Erwartung/Varianz, gültige/ungültige Eingaben, feste Seeds und Stichprobenkonvergenz mit vorab gesetzten Kriterien prüfen. Ein synthetischer Piezo-Lehrfall wird entsprechend gekennzeichnet. |
| #208 / #209 | Bestehende Bände prüfen/integrieren und nach #64 um dessen reale API, Formate, Rezept und Driftprüfungen ergänzen. |
| #210 | Bestehende V&V konsolidieren, #64 detailliert nachweisen und dokumentiertes No-go der experimentellen Auswahl erklären. Vollständiger ehrlicher Befund statt erfundener Messvergleich. Die widersprüchlichen bisherigen Abhängigkeiten/Abnahmeformulierungen vorher konsistent ändern. |
| #211 | Alle Bände redaktionell abnehmen; Quellenmatrix, Querverweise, Rezepte, API und LaTeX-/Archiv-Rebuild prüfen, alle Seiten visuell kontrollieren. Keine Pflichtplatzhalter. |
| #66-Release | Alle neu vereinbarten Pflichtpakete integriert, relevante CI/Reviewbefunde geklärt, frische Wheel/sdist-Installationen geprüft, PDF-/Quellenartefakte mit Commit/Version/Prüfsummen erzeugt. Erst dann Release-PR nach main, Zielcommitprüfung und unveränderlicher annotierter Tag. |

Für B reicht weder dieser Vorschlag noch ein grüner Dokumentationsbuild zur
Freigabe: #64 und die vollständigen Bände fehlen weiterhin. Ein
Code-zu-Code-Vergleich kann einen zusätzlichen Verifikationsnachweis liefern,
aber keine unabhängigen Messdaten ersetzen. Herstellerkennlinien oder aus der
Antwort gefittete Koeffizienten erfüllen das bisherige Datengate ebenfalls nicht.

## 6. Prüfstatus

Eigene Quellen-/Einheiten-/Metadatenprüfung; kein unabhängiges fachliches Review.
Keine numerischen Produktionsänderungen. Quellenmatrixprüfung (133 Quellen),
12 Dokumentations-Gatetests, Ruff/Lint/Format (308 Dateien), mypy (67 Module),
lokale Markdown-Linkprüfung und Diffcheck bestanden. Sieben gespeicherte
Dateihashes nachgerechnet. Kein neuer PDF-Inhalt und kein TeX-/Buildwerkzeug
geändert; deshalb kein neuer PDF-/Layoutnachweis für diesen Recherchebeitrag.
Eine erneute numerische Gesamtsuite ist für diesen reinen Quellen-/Planungsdiff
nicht erforderlich und wird nicht als durchgeführt behauptet.
Entscheidung B ist genehmigt; Implementierung, Bände und Releasegate bleiben
bis zu ihren eigenen Nachweisen offen. Die drei gesichteten Fremd-Rohdateien
werden durch konkrete `.gitattributes`-Einträge byteidentisch bewahrt;
CRLF und originale nachgestellte Leerzeichen werden nicht bereinigt.
