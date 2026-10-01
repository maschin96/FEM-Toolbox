# Piezo: Blocklänge und Kopplungsstärke

[Installation, Ausführung und Lebenszyklus](../README.md).
Aufruf: `python -I study_recipes/run.py piezo --output piezo-results` aus der
sauberen Umgebung außerhalb des Checkouts.

Die vorbereiteten Netze `short.msh`/`long.msh` bilden einen Block
0,01×0,01×L m mit L=0,02/0,04 m, A=1e-4 m² ab. Acht Knoten in
lexikografischer Reihenfolge `(x,y,z)` mit z als schnellster Koordinate,
sechs Tet4 über die Raumdiagonale, zwei untere und zwei obere Triangle3.
Region `specimen` enthält Zellzeilen 0–5; untere Knoten 0,2,4,6, obere
Knoten 1,3,5,7. Beide Netze haben identische Topologie. Reproduzierbares
Netzrezept: Einheitswürfelknoten aus `{0,1}³`, für jede der sechs
Permutationen der Achsen von (0,0,0) sukzessive eine Koordinate auf 1
setzen; die vier besuchten Knoten bilden ein Tet4. Danach Koordinaten
mit (0,01,0,01,L) skalieren. Gmsh-Import orientiert Tet4 positiv.

Synthetischer Benchmark aus [d₃₃-Nachweis](../../../Doc/18_PIEZO_D33_VERIFICATION.md):
C=diag(2e9,2e9,1e9,5e8,5e8,5e8) Pa, κ=diag(1e-8,1e-8,1e-8) F/m,
e hat nur e33=5 oder 10 C/m². Die zweite Variante ist eine explizite
synthetische Verdopplung, kein reales alternatives Keramikmaterial.
Materialachsen stimmen mit globalen x,y,z überein; Polarisationsachse +z,
Stress-Charge-Konvention und technische Schubverzerrungen wie im Projekt.
Die Studie variiert zuerst den Netzpfad, dann den gesamten 3×6-e-Tensor;
das kartesische Produkt liefert vier Fälle.

u_x=u_y=0 an allen Knoten, u_z=0 unten, obere/seitliche mechanische
Flächen lastfrei. Untere Elektrode `ground` bei 0 V, obere `driven`
bei 100 V; seitlicher elektrischer Fluss null. Keine Volumenladungen,
keine thermische Kopplung. Diagonales C vermeidet transversale Spannungen
in diesem bewusst eingeschränkten 1D-invarianten Modell.

Mit E_z=−V/L und σzz=0 folgt unabhängig:

```text
εzz = −e33 V/(C33 L)
u_top = −e33 V/C33
Q_top = A/L (κ33 + e33²/C33) V
H_mech = A/(2L) e33²/C33 V²
H_coupling = −A/L e33²/C33 V²
H_dielectric = −A/(2L) κ33 V²
H_total = −Q_top V/2
```

Q_top ist die **Elektrodenladung**, mit entgegengesetztem Vorzeichen
zur elektrischen Reaktion im Körper. Negative Enthalpie unter festem
Potential ist zulässig. Beiträge nicht als drei positive Energien addieren.

| L [m] | e33 [C/m²] | u_top [µm] | Q_top [nC] | H_mech [µJ] | H_coupling [µJ] | H_dielectric [µJ] |
|---|---|---|---|---|---|---|
| 0,02 | 5 | −0,5 | 17,5 | 0,625 | −1,25 | −0,25 |
| 0,02 | 10 | −1 | 55 | 2,5 | −5 | −0,25 |
| 0,04 | 5 | −0,5 | 8,75 | 0,3125 | −0,625 | −0,125 |
| 0,04 | 10 | −1 | 27,5 | 1,25 | −2,5 | −0,125 |

Bei doppelter Länge bleibt die freie Endverschiebung gleich, Ladung und
Enthalpiebeiträge halbieren sich. Doppeltes e verdoppelt u und vervierfacht
die mechanischen/Kopplungsbeiträge; der dielektrische Beitrag bleibt gleich.
Die Tabelle mittelt u_z nur über obere Knoten, liest die obere Elektrode
und summiert integrierte Energien ausschließlich über `specimen`.
Separate direkte `Model`-Läufe prüfen u/φ an allen Knoten, Ladung,
Energiebeiträge und globale Kraft-/Ladungsbilanz.

Vorab: `rtol=2e-11` wie im vorhandenen d₃₃-Benchmark; `atol=1e-16 m`
für u, `1e-18 C` für Ladung und `1e-16 J` für Energie. Für nullwertige
Bilanzen gelten absolut 1e-10 N und 1e-18 C. Die feldweise Skalierung
und kleine Sparse-LU-Systeme erwarten reine Rundungsabweichungen; Grenzen
liegen deutlich unter µm/nC/µJ-Referenzgrößen. Es werden keine Toleranzen
nachträglich erweitert. Budget: vier Fälle, je acht Knoten/sechs Tet4.
[Gemessene Laufzeiten](../../../Doc/63_STUDY_RECIPE_ACCEPTANCE.md).
Zusätzlich: vier direkte Modelle, zwei Fehlfallversuche und zweimal vier
Fälle für Spannungsänderung/Wiederherstellung.

Annahmen: quasistatisch, kleine Verzerrungen, lineare homogene Materialien,
perfekte Elektroden, keine Fransenfelder außerhalb des Körpers. Tet4 kann
die affinen Felder dieses Blocks exakt repräsentieren; das beweist keine
Genauigkeit für allgemeine Geometrien, Randlagen oder Materialgrenzflächen.
Keine allgemeine Netzkonvergenz und keine experimentelle Validierung.

[Fehler-/Resume-Anleitung](../README.md): negative definite Permittivität
führt zu einem sichtbaren Fehlfall; Spannungsänderung 100→50 V invalidiert
die ursprünglichen Ergebnisse (u/Q halbieren sich, Energien vierteln sich).
