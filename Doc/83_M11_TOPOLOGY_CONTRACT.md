# M11: Struktur-Topologieoptimierung (#76)

19.10.2026 · Vertrag vor Implementierung/Nachweisläufen; keine Freigabe.

Stationäre kleine Verzerrungen auf Tet4, lineare positive Materialtensoren,
entwurfsunabhängige Totlasten, homogene Verschiebungsvorgaben. Alle Lastfälle
teilen Netz, Tensoren und Dichte. Ziel C=sum(w_l f_l^T u_l) [J], w_l>0,
ist die doppelte gewichtete Formänderungsenergie. Keine thermischen Lasten,
Tet10, bewegten Lagerungen, dichteabhängigen Lasten oder gekoppelte Optimierung.

x_e in [0,1]; feste Voll-/Leerzellen besitzen x=1/0. Volumen v=sum(V rho)/sum(V)
bezieht sich auf das gesamte Gebiet einschließlich fester Zellen.
Unmögliche Budgets werden vor dem Lauf diagnostiziert.
Dichtefilter H_ej=max(0,r-|c_e-c_j|) V_j, P_ej=H_ej/sum(H_ek), r>0 [m].
Randzeilen werden normalisiert; feste Eingangsvariablen bleiben vorgeschrieben.
Nach Filterung werden physische feste Zeilen exakt auf 1/0 überschrieben;
ihre Ableitungszeilen und feste Variablenspalten sind null.
K_e=a_e K_e^0, a=epsilon+(1-epsilon)rho^p, Standard p=3, epsilon=1e-3;
p>=1, 0<epsilon<1. rho=0 besitzt numerische Reststeifigkeit, keine entfernte
Geometrie. Kein Heaviside-Schritt innerhalb der Optimierung. Der Filterradius
regularisiert; er garantiert keine harte Mindestwandstärke. Zwischenmaterial:
volumengewichteter Anteil 0.1<rho<0.9 und sum(4rho(1-rho)V)/sum(V).
Ungefilterte SIMP und nicht konsistente Sensitivitätsfilter sind nicht Standard.

dC/drho_e=-p(1-epsilon)rho^(p-1) 2U_e/a_e [J], U gewichtete aktuelle
Zellenergie; Rückführung über P^T. dv/dx=P^T V/sum(V) [1]. Feste Zeilen/Spalten
werden in der Ableitung eliminiert. Der bestehende Mechaniksolver rechnet
jeden Zustand. SciPy SLSQP verwendet vollständige Gradienten, Bounds und
lineare Volumenrestriktion; Zielskalierung mit C_start. Zulässiger konstanter
Start über monotone Volumeninversion. Explizite Starts müssen zulässig sein.
Status: converged, iteration_limit, evaluation_limit, solver_failed.
Konvergenz ist lokal, keine globale Optimalitätsgarantie oder Fertigungsfreigabe.

## Vorab festgelegte Nachweise und Gates

- Vollmaterialparität: relative Felder/C <=2e-10, absolut u<=1e-12 m.
- Zentrale Gradienten h=1e-3,1e-4,1e-5; relative Normabweichung <=2e-4 bei
  h=1e-5, absolute Volumenableitung <=2e-9. Zwei Lastfälle, feste Gebiete,
  Randzellen, ungleiche Volumina und nahe leere/volle Entwürfe.
- Zustandsresiduum <=1e-9; relative Arbeit/Energieabweichung und normierte
  Kräftebilanz <=1e-8. Grenzen erlauben Rundung/Differenzierung kleiner
  konditionierter Systeme; keine nachträgliche Grenzweitung.
- Referenzquader 2x1x1 m, E=1000 Pa, nu=0.3; x=0 fest, getrennte
  Resultierende -y/-z je 1 N bei x=2, Gewichte 0.5/0.5. Eine Zellschicht an
  Last/Klemme voll. v=0.6; Netze (4,2,2),(6,3,3),(8,4,4), r=0.3/0.6 m.
  C-Reduktion >=1%, v<=0.6+1e-8; reproduzierbares C relativ <=1e-10.
  Maximal 100 Iterationen/300 Auswertungen, skaliertes SLSQP-ftol=1e-7.
  Budgetstopp ist kein Konvergenznachweis. Zwischenmaterial, Kosten und
  Netz-/Filtereinfluss werden explizit ausgewiesen; keine behauptete Netzrate.
- Geometrie: Schwellen 0.25/0.35/0.45; selektierte vorhandene Tet4 und exakt
  rekonstruierte Randflächen, ohne Glättung/CAD. Flächenverbundenheit,
  Lastflächenerhalt und ausreichende Lagerung sind Pflicht. Vollmaterial-Reanalyse
  mit unveränderten Lasten: C_geo/C_density<=2, Volumen<=0.9 Ausgangsvolumen;
  mindestens eine Schwelle je Referenzlauf muss bestehen. Andere Schwellen
  bleiben No-go. Dies ersetzt nicht die ursprüngliche 0.6-Volumenrestriktion.
  Von-Mises-Spannung und kritische Zellen werden berichtet; ideale Klemm- und
  Lastsingularitäten erlauben keine Festigkeitsfreigabe.

## API, Provenienz und Release

Separater Topologiebaustein delegiert an vorhandene Mechanik. Additiver
Model-Snapshot statt zweitem Material-/Lastparser. TOML referenziert bestehende
Model-Dateien; direkte Python-API ebenfalls möglich. Neues Ergebnis-Schema 1:
Netz, Material, Lastfälle, feste Gebiete, Dichten, Parameter, Historie,
Status, Eingabehashes und RunProvenance; VTU für Dichten und Geometrie/Spannung.
Bestehende Persistenz bleibt unverändert. Ungültige/nicht endliche Eingaben,
Nichtnull-Lagerwerte und unmögliche Volumen liefern ValidationError;
Zustandsfehler und Budgets werden eindeutig protokolliert. Geometrie-No-go separat.

#83 zuletzt: #76–#82, unabhängige Gradientenreferenzen, Geometrie-Reanalyse,
Artefaktinstallation, aktuelle CI und Dokumentations-/Versions-/Reviewnachweise.
Offene Gates verhindern Release/Tag. M11 ist die Themenkennung; die im Backlog
vorläufige v1.9.0 ist keine automatisch festgelegte Veröffentlichungsnummer.

## Implementierungsbefund zur Variablenskalierung (19.10.2026)

Der erste Matrixlauf ohne Koordinatenskalierung erreichte auf feineren Netzen
innerhalb von 100 Iterationen nur Budgetstopps. Diese Läufe sind No-go und
werden nicht umetikettiert. SLSQP verwendet deshalb z=x/sqrt(n_free): die
quadratische Schrittnorm entspricht der mittleren statt der aufsummierten
Zellvariation. Gradienten, lineare Restriktion und Schranken werden vollständig
zurücktransformiert; Ziel C/C_start, ftol, Budgets und Gates bleiben unverändert.
Dies ist eine Konditionierung des gleichen Problems, keine Grenzweitung.
Nichtkonvexe lokale Minima können sich durch die Suchkoordinaten unterscheiden.

## Expliziter zweiter Geometrieweg vor seiner Verifikation

Reine Schwellengeometrien bleiben mit ihrem ursprünglichen No-go sichtbar.
Als gesonderter, optionaler Weg wird nach Schwellenbildung ein deterministischer
Flächenpfad ergänzt: Dijkstra über den ursprünglichen Tet4-Flächengraphen,
Zielzellkosten V*(1-rho+1e-3), stabile Zellindexreihenfolge, Verbindung der
Komponente mit kleinstem Index zur billigsten weiteren Komponente bis zur
Zusammenhängigkeit. Feste Leergebiete werden nie aufgefüllt. Alle hinzugefügten
Zellen werden ausgewiesen. Schwelle, Material, Lasten, Reanalyse und Grenzen
(C-Verhältnis<=2, Gesamtvolumen<=0.9) bleiben unverändert. Das ist eine
explizite Geometrieänderung, keine behauptete Reanalyse des reinen Schwellenfeldes.
Die Mindestabnahme erlaubt einen erfolgreichen der beiden dokumentierten Wege;
fehlgeschlagene reine Schwellen werden weder gelöscht noch als bestanden gezählt.
