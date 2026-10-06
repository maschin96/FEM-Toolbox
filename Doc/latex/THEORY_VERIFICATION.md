# Theorieband: fachliche Eigenprüfung und Herkunft (#207)

Stand: 6. Oktober 2026. Kandidat, keine unabhängige fachliche Freigabe.
Numerische Grundlage: M9-Integrationsstand `6a1fd36` nach #206/PR #240.
Dieser Beitrag ändert ausschließlich Dokumentation; Operatoren, Materialdaten,
Solver, tolerierte Fehler und historische Referenzwerte bleiben unverändert.

## Umfang und Zuordnung

Der separate `theory.tex` und der Gesamtband `main.tex` enthalten zehn Kapitel.
[Themenmatrix](theory-map.json) nennt ihre Gleichungslabels, kanonischen
Markdown-Quellen, Implementierung und unabhängigen bestehenden Prüfquellen.
[Quellenmatrix](coverage.json) bewahrt den initialen Herkunftsstatus und ergänzt
`theory_coverage` für ausgearbeitete Theorie. Das ist keine automatische Abnahme.

| Bereich | Inhalt und fachlicher Kontrollpunkt |
| --- | --- |
| Konventionen | SI, Engineering-Voigt, Arbeitskonjugation, orthogonale Materialachsen; Strain-/Stress-Charge-Konversion und absolute Permittivität |
| Geometrie | Affine Tet4/Tet10, Eck-/Kantenfolge, Jacobi-Transformation, Simplexmomente, Duffy-/Flächenquadratur und konsistente Lasten |
| Einzelfelder | Starke und schwache Formen von Wärme, Elektrostatik und Mechanik; auswärts gerichteter Wärmefluss, D-Fluss und Leiterladung unterscheiden |
| Thermomechanik | Einseitiger Ablauf, Eigenverzerrungs-RHS, freie/behinderte Expansion; Tet4-Mittelwertenergie von kontinuierlichem Energieintegral unterscheiden |
| Piezo/Thermopiezo | Enthalpievariation, symmetrisch-indefinite Blöcke, thermischer elektrischer Schalter; d33-Anregungen und Reziprozität |
| Algebra | Feld-/Elementassembly, Dirichlet-Elimination aus Originalsystem, Reaktionen, Nullräume und Gesamtladungs-Elektroden |
| Dynamik | Masse/Kapazität, P2-HRZ statt negativer Zeilensummen; Backward Euler mit bewegten Rändern, Energie-/Wärmebilanz |
| Piezo-Dynamik | Algebraische Potentialgleichung, positive Schur-Steifigkeit, Modal, Newmark, komplexe Harmonik und Resonanzdiagnostik |
| Solver/Indikatoren | Kongruente Skalierung, LU, PETSc COMM_SELF und FGMRES/FieldSplit; physikalische Residuen, Tet4-Wärmeindikator und begrenzter Kostenvergleich |
| Surrogate | Grad 1/2, Train-only-Skalierung, vollständige Basis, SVD-Least-Squares/Ridge ohne Intercept-Strafe, Rang und Holdout-/Boxgrenzen |

Die Handrechnungen verwenden bewusst einfache synthetische SI-Zahlen, keine
als reale Werkstoffmessungen ausgegebenen Materialdaten. Herleitungen werden
gegen `Doc/02_MATHEMATICAL_MODEL.md` und die jeweiligen Fachverträge geprüft.
Eigene Abbildung und Texte stehen unter der Repository-MIT-Lizenz. Die
Bibliografie enthält primäre externe Methodenquellen; externe Abbildungen,
Buchpassagen oder Messdaten werden nicht übernommen.

## Nachzurechnende Referenzen

| Beispiel | Unabhängig hergeleitete Erwartung |
| --- | --- |
| Einheits-Tet4, Wärme | V = 1/6 m³; T = x, k = 2 W/(m K): q = (-2,0,0) W/m² und TᵀKT = 1/3 W K |
| Einheits-Tet10, T = x² | Gradientennorm 2/15 W K; bloße Schwerpunktbewertung 1/12 W K ist kein Integral |
| Robin-Quader | k = 2, h = 4, T links = 400 K, Umgebung = 300 K: q = 400/3 W/m², T rechts = 1000/3 K |
| Kondensator | L = 0,002 m, A = 0,0012 m², κ = 2e-11 F/m: C = 12 pF, bei 100 V Q = 1,2 nC, U = 60 nJ |
| Elastischer Affinpatch | E = 100 Pa, ν = 0,25, εxx = 0,01, sonst null: σ = (1,2;0,4;0,4) Pa, U = 0,001 J |
| Behinderte Expansion | Eigenstrain xx = 0,02: σ = (-2,4;-0,8;-0,8) Pa, U = 0,004 J |
| Piezo-Affinpatch | Cxx = 10, e = 2, κ = 3, εxx = 0,2, φ′ = 0,3: σxx = 2,6 Pa, Dx = -0,5 C/m², H = 37/1200 J |
| Thermopiezo geklemmt | Cxx = 10, e = 2, κ = 4, αx = 0,01/K, ΔT = 10 K: bei χ = 1 φ(1) = -0,05 V, σxx = -1,1 Pa; bei χ = 0 φ = 0, σxx = -1 Pa |
| d33, Spannungsantrieb | L = 0,02 m, A = 1e-4 m², C = 1e9 Pa, e = 5 C/m², κ = 1e-8 F/m: bei 100 V u = -5e-7 m, Q = 1,75e-8 C, H = -8,75e-7 J |
| d33, Kraftantrieb | Bei 10 N kurzgeschlossen u = 2e-6 m, Q = -5e-8 C; offen u = 4e-6/7 m, φ = 2000/7 V |
| Zwei Federn | x0 = 2, Endlast 3: x = (2,5,8), Originalreaktion = (-3,0,0) |
| Piezo-Elektrode | e = 0,5, κ = 2, ε = 0,1, Länge 1: offen φ = 0,025 V; bei QL = 1,95 C φ = 1 V |
| Masse/Kapazität | Einheits-Tet4, ρ = 6: konsistente Diagonale 0,1 kg, Nebenwerte 0,05 kg, Gesamtmasse 1 kg; c = 2, ΔT = 10 K speichert 20 J |
| Backward Euler | λ = 1/s, θ0 = 10 K, Δt = 0,1 s: θ1 = 100/11 K gegenüber exakt 10 exp(-0,1) K |
| Ein-DOF-Piezo | M = 0,1 kg, Koffen = 17/9 N/m, ω²offen = 170/9 s⁻²; bei ω = 1 und F = 1 u = 90/161 m, φ = 60/161 V |
| Surrogat | y = 3 + 2x, xtrain = (-2,-1,0,1,2): skalierte Koeffizienten (3,4), Holdouts 0,5 und 1,5 liefern 4 und 6 |
| Ridge | ztrain = (-1,1), y = (1,5), λ = 2: Intercept 3, Steigung 1; Fehler (1,-2) liefern MAE 1,5, RMSE √(5/2), maxabs 2 |

## Prüfung und Grenzen

Die tatsächlichen Ausführungs- und PDF-Reviewergebnisse stehen unten. Für reine Dokumentation wird kein künstlicher Red-Test erzeugt.
Vorhandene Referenz-, Patch-, Bilanz- und Konvergenztests prüfen unveränderte
Produktionspfade mit ihren ursprünglichen fachlichen Toleranzen.

Historische Vier-Netze-Nachweise sind unter `Doc/25_THERMOPIEZO_CONVERGENCE.md` und
`Doc/55_TET10_ACCEPTANCE.md` an ihre jeweiligen Ursprungsstände gebunden.
Die Theorie ordnet diese Nachweise zu, erzeugt keine neuen Referenzdaten und
führt Kostenmessungen nicht als neue Performancebehauptung aus.
Native PETSc-Prüfungen benötigen die gesonderte CI-Laufzeit; lokales Skip ist
kein Erfolg. Eigenprüfung und erfolgreiche CI ersetzen kein unabhängiges
Review. Die weiteren Bände #208–#210, Gesamtabnahme #211 und Release #66 bleiben
außerhalb dieses abgeschlossenen Theorieumfangs offen.

## Ausgeführte Prüfungen

Alle Befehle laufen aus dem Repository-Root mit Entwicklungsumgebung und
`PYTHONPATH=src`; Build/Qualitätswerkzeuge sind bestehende unveränderte Werkzeuge.

- `python -m pytest`: **771 passed, 83 skipped**, 86,46 s, CPython 3.13.
  Enthalten sind die kleinen analytischen Referenzen sowie die vorhandenen
  Konvergenztests. Die 83 Skips betreffen die lokal fehlende native PETSc-Laufzeit;
  diese Prüfungen werden zusätzlich durch die separate PETSc-CI abgedeckt.
- `python -m ruff check .`: bestanden.
- `python -m ruff format --check .`: 306 Dateien bereits formatiert.
- `python -m mypy src/fem_toolbox`: bestanden, 67 Quelldateien.
- `git diff --check`: bestanden.
- `python3 scripts/build_handbook.py --check-matrix`: 122 Quellen zugeordnet.
- `python3 scripts/build_handbook.py --self-test`: fünf PDFs, drei korrekt
  abgelehnte Negativfälle (fehlende Datei, Referenz, Literatur); entpacktes
  Quellenarchiv baut alle fünf PDFs mit identischem extrahiertem Text.
- Separate Nachrechnung mit exakter rationaler Arithmetik bestätigt
  Wärmeenergien 1/3 und 2/15, Robinfluss 400/3, Kapazität 12 pF,
  mechanische Energien 1/1000 und 1/250, Piezo-Enthalpie 37/1200,
  Massendiagonale 1/10, Backward-Euler-Temperaturdifferenz 100/11 und
  harmonische Verschiebung 90/161 in ihren angegebenen SI-Einheiten.
  Das veröffentlichte Matrixlisting wurde unverändert ausgeführt und liefert
  37/1200 J bis 1e-15 J absolut. Das unabhängige Ridge-Beispiel wurde mit der
  Produktions-API nachgerechnet: Werte (2,3,4) bei z = (-1,0,1), absolute
  Toleranz 1e-13, keine Anpassung der bestehenden Testtoleranzen.
- Alle Kapitelpfade sowie Kapitel-, Vertrags-, Modul- und Testpfade der
  Themenmatrix wurden auf Existenz geprüft; Referenzen und Zitate durch den
  strikten LaTeX-Build geprüft. Die Basis bleibt Commit `6a1fd36`; das
  Buildmanifest nennt zusätzlich den tatsächlichen Dokumentationscommit und
  Eingabe-SHA-256, kein historischer Nachweis wird auf einen neuen Commit umetikettiert.

## Visuelles und fachliches Eigenreview

Alle Seiten des zunächst 39-seitigen Theoriebandes und die nach der
Stabilitätsherleitung geänderten Seiten des nun 40-seitigen Bandes wurden als gerenderte Seiten
kontrolliert: Titel/Status, Inhalts-, Abbildungs- und Tabellenverzeichnis,
Einheitentabelle, Gleichungsblöcke, Kopplungsdiagramm, lange Nachweistabelle,
Matrixlisting, Bibliografie und Sachregister. Keine abgeschnittenen Formeln,
Tabellen oder Pfade im Theorieband. Die 43-seitige Gesamtintegration wurde zusätzlich
an Titel, Verzeichnissen und den Übergängen zu den offenen Bänden geprüft.

Die Theorie enthält keine Overfull-Warnung. Sieben Underfull-Meldungen bleiben
für umbrechende Nachweistabellen sichtbar im Manifest/Log; die entsprechenden
Zellen sind vollständig lesbar und keine Inhalte fehlen. Jede Schrift wird
als eingebettet, unterteilt und mit Unicode-Abbildung geprüft; Textsuche,
Umlaute und mathematische Integral-/Gradientenzeichen funktionieren.
Die getrennten Anwendung/API/Nachweis-PDFs bleiben klar markierte Platzhalter.

Fachliche Eigenprüfung: Vorzeichen wurden aus partieller Integration und
Enthalpievariation neu abgeleitet; SI-Dimensionen, Randreaktionen, Leiterladung,
thermische Kopplungsrichtung, Schur-Vorzeichen, Anfangsdaten und komplexe
Phasoren gegen Verträge und Implementierung geprüft. Keine echte Abweichung
von Doc/02 gefunden, daher kein neuer Konventions-ADR. Insbesondere ist die
Tet4-Eigenstrain-Mittelwertenergie ausdrücklich abgegrenzt und die Tet10-
Schwerpunktprobe nicht als Energieintegration ausgegeben.

Es liegt **kein unabhängiges fachliches Review** vor. Der PR nennt diese
Prüflücke; CI und Eigenreview sind kein erfundener Reviewerentscheid.
