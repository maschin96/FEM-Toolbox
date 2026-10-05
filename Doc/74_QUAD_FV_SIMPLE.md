# Stationäre laminare Quad-FV-Strömung (#230 / #226b)

Status: Implementiert; lokale Nachweise unten, Integration und unabhängiges
Review noch offen. Dieser Beitrag
bleibt auf M12; die stabile v1.3-API und die MAC-Referenz bleiben unverändert.
Keine experimentelle Validierung und keine RANS-/FSI-Freigabe.

## API und SI-Vertrag

`solve_quad_flow(FluidDomain, velocity_scale_m_s=..., length_scale_m=...)`
verwendet konstante positive Dichte und positive Zell-/Randviskositäten.
Optionale `force_N_m3` hat Form (2,) oder (ncell,2), `viscosity_Pa_s`
Form ()/(ncell,), `boundary_viscosity_Pa_s` Form ()/(nline2,).
Variable Viskositäten müssen zusammen geliefert werden; der vollständige
Tensor mu(grad U + grad Uᵀ) bleibt erhalten. Extrusionsbreite aus FluidDomain
multipliziert sämtliche Faceflächen und Zellvolumen.

U ist (ncell,2) in m/s, p (ncell,) in Pa und Massenfluss (nface,) in kg/s,
positiv nach außen zum Eigentümer. Arrays und Provenienz sind unveränderlich.
Periodische Partner speichern exakt denselben Fluss mit Gegenvorzeichen;
periodische Druckabfälle werden als äquivalente Volumenkraft ergänzt, das
Druckfeld enthält den periodischen Anteil mit Gauge. Geschlossene Domänen
haben genau eine vorgegebene Referenzzelle; ein Druckauslass schließt eine
zusätzliche Referenz aus. `pressure_reference_Pa` verschiebt allein die Gauge.

Wand: U=0; Einlass: vollständige vorgeschriebene Geschwindigkeit;
Druckauslass: vorgegebener Druck und ∂nU=0; Symmetrie: Un=0 und ∂nUt=0.
Der Wand-/Einlassdruck wird rekonstruiert, nicht mit künstlichem ∂np=0
festgelegt. Auslass-Rückströmung wird mit konkreter Face und kg/s abgelehnt.

## Diskretisierung und Iteration

Die konservativen Bausteine aus Doc/73 bilden sparse affine Operatoren.
Komponentenweise unrelaxierte integrierte Impulsdiagonalen a bilden V/a;
die Normalmobilität ist nᵀdiag(V/a)n. Die kompakte Druckkorrektur und die
Rhie-Chow-Korrektur verwenden dieselbe Mobilität und Facegeometrie.
Periodische Verschiebungen gehen in Nachbardistanzen und Rekonstruktionen ein.

Gewichtete Least Squares rekonstruieren kubische Polynome, verankert am
Zellwert. Zwei Face-Nachbarringe im Inneren und drei am physikalischen Rand
enthalten virtuelle periodische Bilder und explizite Dirichlet-/Neumann-Samples.
Daraus entstehen konsistente Zellgradienten und Facewerte einschließlich
Randableitungen und Drucktraces; rangdefiziente Stencils werden abgelehnt.
Eine isolierte quadratische oder kubische Randableitung bei linearer
Innenrekonstruktion erzeugte Übergangsfehler und erfüllte die Druckordnung
nicht. Der kubische Patchtest prüft Zellgradienten, Face-U, Randableitungen,
p=y³ und den Diffusionsoperator für U=y−y³ unabhängig; lineare Exaktheit
bleibt enthalten. Für die innere Spannungsableitung wird der Facegradient
aus beiden Polynomen ausgewertet und over-relaxed zum kompakten Zellunterschied
korrigiert. Dessen kubischer Term D³U[d,d,d]/24 wird zuvor abgezogen, damit
die innere Normalableitung zur kubischen Randableitung konsistent ist.
Konservativer gemeinsamer Facefluss und kompakter Anteil bleiben erhalten.

Die Picard-Schleife friert den konservativen Transport ein. Ein SIMPLE-
Vorkonditionierer löst Impulsblock und kompakten Druckkorrekturblock getrennt
mit sparse LU und korrigiert U um −G p′/a. GMRES beschleunigt diese Schritte;
es gibt keinen gekoppelten LU-Solve im produktiven Solver. Bis zu 20
Nachkorrekturen verwenden den gesamten nichtorthogonalen Oseen-Operator.
Innere Restbudgets werden aus den äußeren volumenbezogenen SI-Gates
abgeleitet. U und p werden anschließend separat unterrelaxiert.

`convection_scheme`: `donor-cell` (Standard), `linear-upwind` mit explizit
begrenzter Gradientenkorrektur oder `central` für Diffusionsdominanz.
Zentraltransport verlangt auf jedem Face Pe<2. Upwind ist konservativ;
die begrenzte Korrektur extrapoliert nicht über die Nachbarwerte hinaus.
Standardgrenzen: 500 äußere Iterationen, 500 SIMPLE-Schritte pro innerer
Nachkorrektur, 20 Nachkorrekturen und 50 000 U/p-Unbekannte.
`QuadFlowConvergenceError` bewahrt die vollständige äußere Residuenhistorie.

## Frische Rest- und Bilanzprüfung

Nach dem Update werden Transport, Diagonale, Druck-/Massenfluss und vollständige
Spannungen neu berechnet. Die direkte Face-Kontraktion wird zusätzlich gegen
den assemblierten nichtorthogonalen Spannungsoperator geprüft. Sämtliche
Kontinuitätszeilen einschließlich der bei Gauge-Elimination ersetzten Zeile
gehen in die finale Prüfung ein. Uref und Lref sind explizit positive SI-Werte,
Aref=Breite·Lref. Die unveränderten absoluten Gates aus Doc/71 lauten:

| Größe | Absolutes Gate |
| --- | --- |
| Jede lokale Massendivergenz | 1e−8 rho Uref/Lref |
| Jede lokale Impulskomponente, einschließlich Nichtorthogonalität | 1e−8 rho Uref²/Lref |
| Globale Masse | 1e−10 rho Uref Aref |
| Jede globale Impulskomponente | 1e−8 rho Uref² Aref |
| Diskrete Leistung | 1e−6 rho Uref³ Aref |

Die Leistungsbilanz verwendet eigenständig rekonstruierte Facearbeiten.
`gradient_dissipation_W` ist Σ2mu S:S V. Die diskrete viskose Arbeit
`viscous_dissipation_W` ist die gemeinsame innere Facekraft mal Geschwindigkeitssprung
plus äußere Arbeit (Uf−UP)·Fvisc. Ihre Differenz wird als
`viscous_discretization_power_W` ausgewiesen. Ebenso werden konservative
Transportarbeit minus Rand-Kinetiktransport und diskrete Druckarbeit minus
physikalischer Druckrandarbeit separat ausgegeben. Die geprüfte diskrete Bilanz
ist Körperarbeit + viskose Randarbeit + Druckrandarbeit + Druckinterpolationsarbeit
− Gradientendissipation − viskose Diskretisierungsarbeit − Rand-Kinetiktransport
− numerische Transportarbeit. Diese zusätzlichen Arbeiten sind Diskretisierungs-
befunde, keine nachträglich angepassten physikalischen Quellen. Das kleine
diskrete Leistungsresiduum allein bestätigt keine genaue Kontinuumsdissipation.

## Nachweise und Grenzen

`tests/unit/test_quad_flow.py` prüft Kanalfelder, Referenzverschiebung,
20°-Scherung mit Symmetrie, affinen hydrostatischen Druck, unabhängigen
gekoppelten LU-Vergleich, MAC-Parität, Randpatches, Rückströmung und Fehlerfälle.
Die vier Netze in `tests/convergence/test_quad_flow_convergence.py` verwenden die
unveränderten Doc/70/71-Budgets für Poiseuille, nichtlineare MMS und variable mu.
Die MMS verwendet wörtlich psi=sin²(pi x)sin²(pi y); die ältere MAC-MMS
verwendet doppelte Amplitude. Die unabhängigen analytischen Quellen passen
zur jeweiligen Amplitude; keine Fehlergrenze wird deshalb verändert.

Der [Dateikanal](../examples/quad_flow/README.md) läuft aus installiertem
Wheel/sdist. Input- und Geometriehash, Schemes, SI-Skalen, Relaxation,
Ressourcenlimits, NumPy/SciPy-Versionen und vollständige Rohresiduen werden
in der Provenienz gespeichert. SA ist #87; FSI, bewegte Netze, Transienz und
3D sind nicht Bestandteil dieses Beitrags.

## Lokaler Verifikationsstand (05.10.2026)

Python 3.13.5, NumPy/SciPy-Versionen werden je Lauf in Provenienz erfasst.
Die unveränderten Fehlerbudgets für Zentral- und begrenzte Upwind-MMS wurden
auf 16/32/64/128 bestanden:

| Verfahren / Fehler | 16 | 32 | 64 | 128 |
| --- | ---: | ---: | ---: | ---: |
| Zentral, RMS-U (m/s) | 0,0121328 | 0,00309771 | 0,000773520 | 0,000193219 |
| Zentral, RMS-p (Pa) | 0,0308871 | 0,00314491 | 0,000428110 | 0,0000908579 |
| Linear-Upwind, RMS-U (m/s) | 0,0113543 | 0,00300362 | 0,000761719 | 0,000191743 |
| Linear-Upwind, RMS-p (Pa) | 0,0683501 | 0,0161901 | 0,00405735 | 0,00102213 |

Alle drei Netzverdopplungen überschreiten die jeweiligen Mindestordnungen.
Poiseuille auf 8/16/32/64 erreicht feinsten maximalen U-Fehler 1,17e−13 m/s
und RMS-Druckgradientfehler 1,94e−12 Pa/m; diese Zahlen sind durch
Polynomexaktheit und numerische Lösung begrenzt, keine allgemeine
Genauigkeitsbehauptung für andere Felder.

Der uniforme Variable-mu-Kanal ist ebenfalls polynomexakt. Mit Standard-
Relaxation lag der RMS-U-Fehler auf allen vier Netzen bei etwa 3,4e−10 m/s,
unterhalb des Iterationsfehlers; ein Quotient solcher Fehler misst keine
Diskretisierungsordnung. Die ursprünglichen vier uniformen Netze bleiben als
strenger Exaktheitsnachweis (<1e−9 auf **jedem** Netz) erhalten. Zusätzlich
prüfen vier glatt verformte Netze desselben physikalischen Einheitsquadrats
mit y=eta+0,08 sin(2pi xi)sin(pi eta) dieselbe Lösung, mu und Quelle mit
unverändertem feinsten Fehlerbudget <1e−3 und Ordnung >1,8. Zentraltransport
ist bei dieser Diffusionsdominanz zulässig; alpha_U=alpha_p=1 wird explizit
protokolliert. Die RMS-U-Fehler der verformten Netze betragen
9,37e−5 / 2,38e−5 / 5,96e−6 / 1,49e−6 m/s (Ordnungen 1,98 / 2,00 / 2,00).
Mit alpha=1 liegen die uniformen Fehler zwischen 1,46e−15 und 6,47e−14 m/s.
Der unabhängige volle-Spannungs-Patch auf orthogonalen und
20°-gescherten Netzen prüft -d[(1+y)(1−2y)]/dy=1+4y bis 1e−10 N/m³.
Diese Behandlung des exakten Referenzfalls ist vor Integration zu reviewen;
keine Änderung der numerischen Fehlergrenzen in Doc/70/71.

TDD: fehlende Solver-API zunächst ImportError; anschließend Kanalfelder grün.
Quadratische/kubische Rand-, Zell-/Face-Rekonstruktion und kubischer
Diffusionspatch reproduzierten echte Diskretisierungsfehler vor Korrektur.
Fehlende Geometrie-/Residuenprovenienz zunächst KeyError, anschließend grün.
Dateilauf zunächst fehlend, nach Implementierung samt Überschreibschutz grün.

Wheel und sdist wurden mit isoliertem Build gebaut, getrennt in frische
Umgebungen installiert und mit `python -I` außerhalb des Checkouts ausgeführt.
Beide Dateiläufe liefern identische Reports (NumPy 2.5.3/SciPy 1.18.1), globale
Masse −3,97e−12 kg/s, x-Impuls 6,41e−11 N und Leistung −5,54e−11 W.
83 optionale native PETSc-Prüfungen sind lokal mangels PETSc übersprungen;
CI enthält deren gesonderte Real-/Complex-Jobs. CI und unabhängiges Review
sind vor Integration am aktuellen PR-Commit separat zu prüfen.

## Sichtbare und getrennte CI-Abnahme

Der erste gemeinsame Testjob überschritt auf allen drei GitHub-Python-Versionen
das 30-Minuten-Budget während der Quad-FV-Abnahme. Deshalb laufen die fünf
`cfd_convergence`-Prüfungen in verpflichtenden separaten Jobs: Poiseuille,
zentrale MMS, begrenztes Upwind und beide Variable-mu-Netzfamilien, jeweils
unter Python 3.11/3.12/3.13. Jeder Job zeigt den Beginn und Fehlerwerte jeder
Netzstufe, bleibt auf 20 Minuten begrenzt und speichert seinen JUnit-Nachweis.
BLAS-/OpenMP-Threads sind in den Testjobs auf eins begrenzt, um konkurrierende
Rechenthreads bei den vielen kleinen Krylov-Vektoroperationen zu vermeiden.
Der übrige Testjob hat wieder ein 15-Minuten-Budget. Diese Aufteilung lässt
keine Abnahme aus; sämtliche CFD-Jobs müssen vor Integration grün sein.

```sh
python -m pytest -m "not cfd_convergence" --durations=10
python -m pytest -m cfd_convergence -vv -s -o junit_family=legacy
```

`python -m pytest` ohne Filter führt weiterhin sämtliche Prüfungen aus.
