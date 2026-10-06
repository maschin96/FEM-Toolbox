# M9: Validierungsfall, Messdatenvertrag und Datenauswahl (#63)

Stand: 6. Oktober 2026. **Vorbereitung; Datengate nicht erfüllt.**
Kanonische Grundlage: [V&V](05_VERIFICATION_AND_VALIDATION.md),
[Konventionen](02_MATHEMATICAL_MODEL.md),
[Piezo-Material](15_PIEZOELECTRIC_MATERIAL.md),
[Elektroden](20_ELECTRODES.md).

Der Eigentümer hat am 6. Oktober bestätigt, dass kein vorgesehener Datensatz
oder Messkontakt vorliegt. Die unten dokumentierte öffentliche Recherche
liefert bislang keinen für den gesamten Vertrag qualifizierten Datensatz.
Das ist keine Behauptung, dass geeignete Daten grundsätzlich nicht existieren.
Es wurden keine Messwerte erfunden, keine Kurven digitalisiert und keine
Materialparameter an Beobachtungen angepasst. #63, #64, #65 und das M9-Releasegate
bleiben offen. #208/#209 setzen die integrierte #64-Funktion ausdrücklich voraus.

## 1. Abgegrenzter Zielversuch

Vorgesehen ist ein kleiner, linearer quasistatischer Piezo-Körper im d33-Modus,
vorzugsweise ein homogenes, rechteckiges monolithisches Keramikelement.
Die Auswahl ist ein Beschaffungsziel, noch kein ausgewählter Probekörper.
Ein Multilayer-Stack ist nur zulässig, wenn einzelne Schichten, Polarisationen,
Elektroden und mechanische Kontakte vollständig modellierbar sind.

Primärvergleich: axiale Endflächenverschiebung gegenüber tatsächlich am
Probekörper anliegender Spannung bei dokumentierter mechanischer Lagerung.
Optional getrennte Fälle: blockierte Reaktionskraft bei Spannungsantrieb;
Leiterladung bei Kraftantrieb mit eindeutigem Kurzschluss/Leerlauf.
Die Endflächenverschiebung ist ein flächengewichtetes Mittel in festgelegter
Messrichtung relativ zur Referenzfläche, keine beliebige Maximalverschiebung.
Punkt-/Lasermessungen benötigen stattdessen ihren wirklichen Messpunkt,
Strahlrichtung und Referenzbezug; beide Messoperatoren dürfen nicht vermischt werden.

Der Anspruch gilt nur für den vorab festgelegten kleinen Spannungs-/Lastbereich,
Temperatur, Probenpopulation, Lagerung und elektrische Randbedingung.
Hysterese, zeitabhängiges Kriechen, Kontaktwechsel, Schaltlasten, Resonanz,
Depolarisation und nichtlineares Materialverhalten gehören nicht zum Modell.
Ein niedrigfrequenter Messlauf ist nicht allein durch seine Frequenz als
quasistatisch qualifiziert: Einschwingen, Haltezeit und Drift müssen belegt sein.
Ein bestandener Vergleich validiert diesen Fall, nicht die gesamte Toolbox.

## 2. Pflichtmetadaten und Herkunft

| Gruppe | Vor Freigabe erforderliche Angaben |
| --- | --- |
| Identität | DOI/Datensatzversion oder unveränderlicher Archivbezug, Autor/Institut, Messdatum, eindeutige Proben-/Versuchs-/Wiederholungs-IDs |
| Nutzungsrechte | Datenlizenz oder dokumentierte Erlaubnis einschließlich Verarbeitung, Weitergabe und abgeleiteter Artefakte; Publikationszugang allein genügt nicht |
| Originaldateien | Unveränderte Rohdaten, Dateiliste, Größen und SHA-256; Herkunftsbeleg und Beschreibung jeder Verarbeitung mit Quellversion |
| Geometrie | Dimensionen und Toleranzen, reale Schicht-/Elektrodendicken, Polarisation, Koordinatenursprung/Achsen, Messorte, aktive Elektrodenflächen |
| Material | Unabhängige C/e/κ-Daten oder vollständige s/d/κ-Daten, SI-Einheiten, Voigt-/Schubkonvention, Materialachsen, Temperatur und Messverfahren; Unsicherheiten und Korrelationen |
| Lagerung | Kontaktflächen, Klemmrichtung, Nachgiebigkeit/Preload und deren Unsicherheit; dokumentierte freie Oberflächen und Gegenkräfte |
| Elektrik | Gemessene Proben-Spannung statt bloßer Generatorvorgabe, Verstärker, Referenzelektrode, Leitung/Leerlauf/Kurzschluss, Eingangswiderstand und Messketteneinfluss |
| Messung | Sensor, Kalibrierung, Auflösung, Unsicherheit, Referenzbezug, Zeitsynchronisation, Filter, Haltezeiten, Wiederholungen und Umgebungsbedingungen |
| Datenwerte | Spaltennamen, SI-Einheiten, Vorzeichen, fehlende/ungültige Werte, Ausschlussgründe, Unsicherheitstyp (Standard/erweitert), Überdeckungsfaktor und Kovarianz |
| Unabhängigkeit | Herkunft von Material-/Kalibrierdaten, physisch getrennte Validierungsproben/-läufe, feste IDs der Partition und Sperrung vor Anpassung |

Auflösung ist keine Standardunsicherheit. Ein Hersteller-Nennwert ist ohne
Herkunft und Streuung kein vollständiger unabhängiger Materialnachweis.
Fehlende Unsicherheiten werden nicht mit null ersetzt; unbekannte Korrelation
ist nicht automatisch Unabhängigkeit. Die Eingabeverteilung für #64 benötigt
eine eigene sachliche Begründung, nicht nur bequeme Normalverteilungen.

## 3. Messwertformat und Beschaffung

Ein künftig qualifiziertes Archiv enthält unveränderte Originaldateien,
Lizenz/Herkunft, Metadaten, Materialquellen und eine abgeleitete SI-Tabelle.
Die Tabelle ist ein **geplanter Austauschvertrag**, noch kein implementierter
Produktionsparser oder bestehender Datensatz. Pro Zeile sind erforderlich:

`sample_id, run_id, repeat_id, time_s, voltage_V, force_N, displacement_m,
charge_C, temperature_K, observable, standard_uncertainty, partition, valid`

Nur tatsächlich gemessene Größen werden befüllt. Bei nicht gemessenen optionalen
Größen bleibt das Feld leer und der Messumfang wird im Metadatensatz erklärt.
Die Standardunsicherheit bezieht sich auf `observable` und trägt dessen
SI-Einheit. Mehrere Observablen benötigen getrennte Zeilen sowie eine
zugeordnete Kovarianzdatei; eine einzige Zahl ersetzt keine Kovarianzmatrix.
Generator- und Proben-Spannung erhalten getrennte Spalten, falls beide vorliegen.

Beschaffung: qualifizierte Datensatzversion festlegen, Lizenz prüfen,
Originale unverändert archivieren und Hashes dokumentieren; Einheiten,
Metadaten und Plausibilität prüfen; erst danach Datensatzannahme und Partition
freigeben. Fremde Auswertungsskripte werden vor Ausführung geprüft.
Bei großen Archiven genügt ein dauerhafter Beschaffungsbezug mit Dateihashes,
sofern die Lizenz die erforderliche Nutzung erlaubt. Es gibt bislang kein
angenommenes Archiv und daher keine erfundenen Downloadbefehle oder Hashes.

## 4. Vorab festzuschreibender Vergleich

Vor jeglicher Anpassung werden Proben-/Lauf-IDs, Observablen, Messoperatoren,
Vorzeichen, Auswertezeitfenster, Ausschlüsse, Materialquellen und Grenzen in
einer versionierten Annahmeakte eingefroren. Validierungsdaten bleiben
unangetastet; nachträgliche Auswahl eines günstigen linearen Kurvenstücks ist
keine unabhängige Validierung. Kalibrierung darf ausschließlich die vorab
zugewiesenen Kalibrierdaten verwenden. Bei fehlendem getrenntem Bestand wird
nicht kalibriert; dieselben Daten dürfen nicht beide Rollen erfüllen.

Für jeden reservierten Messpunkt lautet die Abweichung in der Zieleinheit
`delta = y_FEM - y_measured`. Eingabeunsicherheit, Messunsicherheit und
numerischer Raum-/Solverfehler werden getrennt berichtet. Für die
Unsicherheit der Differenz gilt mit nachgewiesener Kovarianz
`u_delta² = u_prediction² + u_measurement² - 2 cov(prediction, measurement)`.
Die Null-Kovarianz-Vereinfachung benötigt tatsächlich unabhängige Quellen.
Diskretisierungsfehler ist kein zufälliger Messfehler; er wird durch eine
separate Netzstudie mit fixen physikalischen Eingaben begrenzt.

Eine Akzeptanzakte muss zwei verschiedene Fragen vorab beantworten:

1. **Messkonsistenz:** Deckt das begründete Unsicherheitsintervall die
   Abweichung? Ein Faktor 2 darf nur mit dokumentiertem Überdeckungsmodell
   verwendet werden; er ist keine allgemeine 95-Prozent-Garantie.
2. **Anwendungstauglichkeit:** Liegt die Abweichung innerhalb einer vorab aus
   dem Verwendungszweck begründeten absoluten/relativen Toleranz? Relative
   Toleranzen allein sind bei Null oder kleinen Zielwerten ungeeignet.

Konkrete Spannungs-/Lastgrenzen, Messzeitfenster, Anforderungstoleranzen,
Überdeckungsmodell und maximaler numerischer Fehler sind **noch offen**:
Ohne Probekörper, Messkette und Nutzungsanforderung wären Zahlen unbegründet.
Sie werden vor Einsicht in die Validierungsantworten verbindlich festgelegt;
dieser Vorbereitungsstand behauptet keine fertig präregistrierte Abnahme.
Bei multiplen Messpunkten ist zusätzlich vorab festzulegen, ob punktweise
oder gemeinsame Überdeckung bewertet wird. Ein nach Ergebnis verbreitertes
Band oder nachträglich aufgewertetes Fehlerbudget ist nicht zulässig.

## 5. Geprüfte öffentliche Kandidaten

Recherche am 6. Oktober 2026: Datensatzbeschreibungen und Primärquellen,
keine vollständige Rohdaten-/Lizenzauditierung. Gründe unten sind eine
Eignungsbewertung für den abgegrenzten Fall, kein Qualitätsurteil über die Studie.

| Quelle | Beobachtung und Entscheidung |
| --- | --- |
| [EPFL Green Piezo, Zenodo v1, DOI 10.5281/zenodo.7920056](https://zenodo.org/records/7920056) | Gedruckte KNbO3-Geräte auf Papier, Kraftpulse/Spannungsantwort, Geometrie- und Charakterisierungsdateien. Weiterer Kandidat; vollständige unabhängige Tensoren, Messunsicherheiten und Abbildung der Messschaltung auf einen linearen quasistatischen Fall bisher nicht nachgewiesen. Lizenzfeld im gelesenen Text ohne auswertbaren Wert; daher keine Lizenzfreigabe behauptet. |
| [Serien-/Parallelschaltungen, DOI 10.5281/zenodo.13830840](https://zenodo.org/records/13830840) | Frequenz-/Lastwiderstandsabhängige Strom-/Spannungsdaten. In dieser Form dynamischer Schaltungsfall, kein statischer d33-Verschiebungsnachweis. |
| [Pz26-Stacks, Mendeley v1, DOI 10.17632/yfr68bw58c.1](https://data.mendeley.com/datasets/yfr68bw58c/1) | CC BY 4.0, Geometrie und Messgeräte beschrieben; frequenzgetriebene Platten, Glas/Kleber und nachgiebige Klemmung. Keine ausreichend bestimmte statische Lagerung und Unsicherheitsbeschreibung im gelesenen Eintrag; nicht angenommen. |
| [Preisach-Identifikation, Mendeley v9, DOI 10.17632/kzf668wggx.9](https://data.mendeley.com/datasets/kzf668wggx/9) | Hysterese-/Identifikationsdaten, Generator-Spannung vor Verstärkung, CC BY NC 3.0. Keine allgemeine Weitergabelizenz oder unabhängige lineare FEM-Validierung daraus abgeleitet. |
| [NPL DEPC-MPR 040, Lodeiro/Stewart/Cain 2005](https://eprintspublications.npl.co.uk/3372/) | Statische Kraft-/Verschiebungscharakterisierung, Messkette und Nachgiebigkeit; nichtlineares Verhalten/Hysterese erwähnt. Methodischer Beschaffungskandidat; Rohdatennutzung und vollständige unabhängige Materialmetadaten nicht nachgewiesen. |
| [Dryad, DOI 10.5061/dryad.n63m4](https://datadryad.org/dataset/doi%3A10.5061/dryad.n63m4) | First-principles-Materialberechnungen statt unabhängiger Probenmessungen; nicht als experimentelle Validierung zugelassen. |

Keine Kontakte wurden angeschrieben und keine eingeschränkten Daten weitergegeben.
Ein Kandidat kann nach zusätzlicher Rohdaten-/Metadatenprüfung geeignet werden;
eine ungeprüfte Teilmenge wird hier nicht als Messdatensatz angenommen.

## 6. Go/no-go und weitere Arbeit

| Gate | Zustand | Noch erforderlich |
| --- | --- | --- |
| Referenzierbare echte Messdaten und Nutzungsrechte | Offen | Qualifizierte Datensatzversion, Lizenzbeleg, Rohdateien/Hashes |
| Vollständiger linearer quasistatischer Probenvertrag | Offen | Materialtensoren/Unsicherheiten, Lagerung, Messoperator und elektrischer Messketteneinfluss |
| Vorab eingefrorene Grenzen | Offen | Anforderung, Amplituden-/Temperaturbereich, Auswertefenster und Abnahmeentscheidung vor Kalibrierung |
| Getrennte Kalibrier-/Validierungsdaten | Offen | Konkrete unabhängige IDs/Proben oder Verzicht auf Kalibrierung mit unabhängigem Materialbestand |

**No-go für den Abschluss von #63 und experimentelle M9-Freigabe.**
Der nächste fachliche Schritt ist die Beschaffung oder ergänzende Klärung eines
qualifizierten Datensatzes. Fehlende Daten werden nicht durch die vorhandenen
analytischen d33-Tests, Herstellerkennlinien oder synthetische Stichproben ersetzt.
Der vorbereitete Vertrag kann geprüft werden, schließt aber #63 nicht ab.

Prüfung dieses Beitrags: Inhalts-/Link-/Quellenmatrixkontrolle, Abgleich mit den
M5- und V&V-Konventionen. Keine neue ausführbare Funktion; TDD-Red-Test nicht
anwendbar. Unabhängiges wissenschaftliches Review liegt noch nicht vor.
