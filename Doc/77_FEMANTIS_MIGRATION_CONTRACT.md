# Femantis — Namensinventar und verbindlicher Migrationsvertrag

05.10.2026 · [Issue #235](https://github.com/maschin96/FEM-Toolbox/issues/235) · Paket A

[Umbenennungsplan](76_FEMANTIS_RENAMING_PLAN.md) ·
[Corporate Identity](75_FEMANTIS_CORPORATE_IDENTITY.md) ·
[ADR-053](07_RISKS_AND_DECISIONS.md#adr-053--femantis-namespace-und-kompatibilitätsvertrag-235)

## 1. Basis und exaktes Inventar

Basis ist der nach #234 integrierte main-Commit `fe6bdf3f765146f6a4fcef33c787647e83f81640`.
[Maschinenlesbares Inventar](brand/femantis-name-inventory.json) erfasst alle
372 getrackten Dateien dieses Commits: 267 Dateien mit Namensfundstellen,
1402 Textvorkommen und 67 bestehende Pythonmodule.
1 Binär-/nicht-UTF-8-Dateien sind mit Größe und SHA-256 separat gelistet.
Die gelieferte Bildtafel bleibt als binäre visuelle Markenreferenz erhalten.

Die frühere Schätzung von 259 Dateien in Doc/76 ist historisch und verwendet
andere Suchbedingungen. Dieses Inventar verwendet exakt den regulären Ausdruck
`fem[_ .-]?toolbox|femantis`, ohne Beachtung der Groß-/Kleinschreibung, auf Dateipfaden und
UTF-8-Textinhalten. Inhaltstreffer erhalten Zeile und Spalte (jeweils ab 1),
Originalschreibweise und Dateiprüfsumme. Dateipfadtreffer sind gesondert erfasst.
Neue Dokumente dieses PRs zählen noch nicht zur eingefrorenen Basis.

| Kategorie | Dateien |
|---|---:|
| documentation | 54 |
| examples-benchmarks | 30 |
| historical-evidence | 16 |
| implementation | 68 |
| packaging-ci | 5 |
| tests | 94 |

Reproduktion: `git ls-tree -r --name-only <basis_commit>` liefert die Pfade;
`git show <basis_commit>:<pfad>` liefert die jeweiligen Bytes. UTF-8-Dateien
mit dem obigen Python-regex durchsuchen; SHA-256 aus den unveränderten Bytes
bilden. NUL-haltige/nicht dekodierbare Dateien separat listen. Nicht gegen die
Arbeitskopie zählen, da sich Dokumentation und Code nach dem Basiscommit ändern.

Die Maßnahmen im JSON sind dateiweise Prüfaufträge, keine Freigabe für pauschales
Ersetzen: auch Implementierungsdateien enthalten persistierte Altkennungen.
Historische Dokumente enthalten aktuelle Hinweise neben datierter Evidenz;
diese Inhalte müssen bei Paket E abschnittsweise beurteilt werden.

## 2. Entscheidung: Namen, Version und Integration

Der Produktname ist **Femantis**, die gestaltete Wortmarke **FEMantis**.
Der kanonische Python-Namespace und Distributionsname lauten `femantis`.
GitHub-Ziel ist `maschin96/Femantis`. Die aktuelle Distribution 1.3.0 wird
nicht nachträglich umbenannt oder neu veröffentlicht.

Die Umbenennung wird der bereits geplanten **v2.0.0-Freigabe M14 / #113**
zugeordnet. Sie erhält dort zusätzliche verpflichtende Migrationsgates;
kein konkurrierendes zweites v2.0-Release. Bis die M14-Voraussetzungen gemäß
bestehender Roadmap erfüllt sind, entsteht kein `codex/milestone-m14`.
Die Codeumsetzung B–F beginnt erst auf dessen regelgerecht angelegtem
Integrationsbranch und geht über eigene Issue-PRs dorthin. Das bindet die
Migration an die bestehende Freigabeplanung und verschiebt keine unfertige
Milestone-Funktion nach main. Aktuelle Planungsdokumentation bleibt
milestone-unabhängige Repositorypflege.

Der Legacy-Import wird während **aller v2.x-Versionen** unterstützt. Entfernung
frühestens mit v3.0.0, nach explizitem Issue, Release-Notes-Ankündigung und
erneuter Kompatibilitätsentscheidung. Es gibt keinen kalenderbasierten
Ablauf und keinen automatischen Löschauftrag. Ein möglicher Legacy-
Distributionsadapter ist optional und darf keine konkurrierenden Paketdateien
besitzen; das neue Wheel enthält den Importadapter. Ohne publizierbare neue
Distribution ist das v2.0-Namensgate unerfüllt.

## 3. Kanonischer Code und Altimporte

Das Inventar listet jeden vorhandenen Modulpfad samt Zielpfad. Paket B muss
sämtliche dort gelisteten Altmodule erreichbar halten; daraus wird nicht
nachträglich eine stabile fachliche API für interne Module. Stabilität und
öffentliche Einstiegspunkte richten sich weiterhin nach Doc/29 und den
jeweiligen neueren Fachverträgen.

Eine einzige Implementierung liegt unter `src/femantis/`. Alt- und Neuimporte
referenzieren identische Modul-/Klassen-/Exceptionobjekte; kein kopierter Code
und kein doppelt ausgeführtes Modul. Der Adapter muss unter beiden
Importreihenfolgen und in frischen Prozessen funktionieren. Optionale PETSc-
und Visualisierungsabhängigkeiten werden nicht allein wegen Import des
Rootpakets geladen. Importwarnings dürfen optionale Imports nicht erzwingen.

`fem_toolbox` erzeugt eine `DeprecationWarning` mit Ziel `femantis`, einmal pro
Prozess beim Laden des Legacy-Rootpakets. Der Warnungstext benennt die
Unterstützung bis Ende v2.x. Tests prüfen Kategorie, Text und Auslöseort;
keine globale Warnungsunterdrückung. Typinformationen bleiben für beide
Importpfade nutzbar; Paket C prüft `py.typed` und mypy-Verbraucher.

B muss insbesondere Mesh/Materialobjekte, Exceptionidentität, `isinstance`,
Model-Läufe, gemischte Importreihenfolge und kanonische Serialisierungs-
identität prüfen. Bisherige Pickle-Nutzung ist kein dokumentierter
persistierter Datenvertrag; keine pauschale Zusage für fremde Pickles.
Tatsächlich vorhandene Serialisierungsverträge separat inventarisieren.

## 4. Persistenz bleibt kompatibel

Die Umbenennung allein erhöht keine Daten-Schemaversion und benennt keine
vorhandenen Keys um. Bereits veröffentlichte Fixtures, Referenzen und
Checksummen bleiben unverändert.

| Vertrag | Entscheidung für v2.0-Migration |
|---|---|
| VTU `fem_toolbox:*` einschließlich `set:`, `set-present:`, `external_id`, Tags | Reader/Writer behalten Kennungen; kein paralleles Präfix ohne neues Schema |
| Wärme-Serie `fem-toolbox-heat-series` | Kennung erhalten; Marke nur außerhalb des alten Identitätsvertrags ändern |
| Ergebnis-/Zeit-/Piezo-/Studien-/Surrogatfelder | Bestehende Schlüssel und Versionsregeln erhalten, einschließlich `toolbox_version` und `toolbox` |
| Historische `solver_name`-Werte und Provenienz | Alte Artefakte unverändert lesbar; neu erzeugte Solverkennung darf kanonisch sein, sofern das Schema sie als freien Wert erlaubt |
| Modell-/Studienhashes und Wiederaufnahme | Kein nachträgliches Umschreiben; Softwareversionswechsel kann bewusst neuen Hash erzeugen und muss dokumentiert werden |
| Alte Tags, Release Notes und V&V-Evidenz | Historische Namen und URLs erhalten; neue Hinweise ergänzen statt Evidenz ersetzen |

Eine spätere Änderung von Präfixen/Keys erfordert ein eigenes Schema-Issue,
Dual-Reader, Konfliktablehnung und unabhängige Roundtrip-Tests. Sie gehört
nicht zum reinen Namespacewechsel. Paket D prüft die vollständigen Varianten
von Ergebnisbündeln, Serien, Studien und Surrogaten anhand realer Altfixtures.

## 5. Öffentliche Namensbefunde am 05.10.2026

Es wurden ausschließlich lesende Abfragen durchgeführt. Kein Account,
Domainname oder Paket wurde angelegt/reserviert. Die HTTP-Befunde sind
Momentaufnahmen und müssen unmittelbar vor Veröffentlichung erneuert werden.

| Dienst / Abfrage | Befund | Aussage und Grenze |
|---|---|---|
| [PyPI JSON: femantis](https://pypi.org/pypi/femantis/json) | HTTP 404 | Kein öffentliches Projektobjekt; keine Garantie gegen Namenssperre, Verwechslung oder fehlende Publikationsrechte |
| [GitHub API: maschin96/Femantis](https://api.github.com/repos/maschin96/Femantis) | HTTP 404 mit authentifiziertem Aufruf | Kein für diese Sitzung sichtbares Zielrepo; Umbenennung erst nach erneuter Kollisionsprüfung |
| [Verisign RDAP: femantis.com](https://rdap.verisign.com/com/v1/domain/femantis.com) | HTTP 404 | Kein RDAP-Domainobjekt bei dieser Abfrage; keine Registrierungs-, Preis- oder Markenrechtsgarantie |
| Aktuelles GitHub-Repository | API meldet administrativen Zugriff, `has_pages=false`, keine Homepage | Sitzung kann Repository verwalten; daraus folgt keine Freigabe zum Umbenennen |

`femantis.com` ist eine untersuchte Option, keine beschlossene oder gekaufte
Domain. Andere TLDs wurden nicht geprüft. Markenregister-/Namensrechtsprüfung
und Rechte an der Bildtafel sind weiterhin offen; Verfügbarkeitsbefunde
belegen keine juristische Freigabe. Diese Nachweise blockieren externe
Veröffentlichung, nicht das dokumentierte lokale Architekturkonzept.

## 6. Externe Zuständigkeiten und Freigabegates

| Fläche | Zuständige Rolle | Nachweis vor Paket F |
|---|---|---|
| GitHub-Name, Remotes, Links, About | Repository-Eigentümer maschin96 | Umbenennung freigegeben, Ziel konfliktfrei, reale Clone-/Issue-/Release-Links getestet |
| PyPI/TestPyPI falls verwendet | Vom Eigentümer benannte Registry-Verantwortliche | Account-/Projektzugriff, vertrauenswürdige Publikationsanbindung und Zielnamenberechtigung nachgewiesen |
| Domain falls gewünscht | Vom Eigentümer benannte Domain-Verantwortliche | Domainwahl, Registrierung und DNS-Zugriff geklärt |
| Bild-/Markenrechte | Projekt-/Markenverantwortliche | Urheberschaft/Nutzungsrechte und Namensprüfung dokumentiert |
| Numerische und API-Freigabe | Integrationsverantwortliche plus unabhängiges Review | Parität, Altimporte, Artefakte und historische Daten geprüft |
| Release M14 | Verantwortliche von #113 | Zusätzliche Migrationsgates und vollständiger M14-Umfang abgenommen |

Registry, Domain, Webhooks, Secrets und externe Verbraucher sind nicht
vollständig aus dem öffentlichen Quellstand feststellbar. Es wurden keine
Secrets gelesen oder veröffentlicht. Interne Deployment-/Webhook-/Registry-
Anbindungen vor F mit dem Eigentümer inventarisieren; fehlende Angaben
explizit als offene Veröffentlichungsgates führen. Docker-Tags im vorhandenen
CI sind lokale Buildnamen, kein Beleg für eine veröffentlichte Registry.

## 7. Arbeitspaket-Gates

**A liefert:** eingefrorenes exaktes Inventar, überprüfte öffentliche
Namensbefunde, ADR und Release-/Kompatibilitätsentscheidungen. Die Entscheidung
steht im PR zur unabhängigen Prüfung; sie ist keine behauptete Reviewfreigabe.

**B startet nach Integration von A und Bereitstellung des M14-Branches:**
Red-Nachweis für fehlendes `femantis`, anschließend kanonischer Code und
identischer Legacyadapter gemäß Modulmatrix. Keine Paket- oder Datenmigration
als versteckter Nebenbeitrag. C/D/E/F folgen gemäß Doc/76, mit den hier
festgelegten Namen, Zeiträumen und Veröffentlichungsgates.

**F bleibt gesperrt**, solange Registry-/Markenrechte, Releasegates oder
Eigentümerfreigaben fehlen. A dokumentiert diese offenen externen Nachweise;
es behauptet sie nicht als erfüllt. Eine Domain ist keine Voraussetzung
für die Python-Distribution; sie ist nur bei tatsächlicher Nutzung ein Gate.
