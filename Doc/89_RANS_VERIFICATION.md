# M12 — getrennte RANS-Code- und Modellprüfung (#88)

Stand 20.10.2026, Entwicklungskandidat 1.7.0.dev1. Dies ist keine M12-Freigabe.
Der [Vertrag](71_M12_RANS_FSI_CONTRACT.md) und der vorab fixierte
[Prüfplan](88_SA_NOFT2_IMPLEMENTATION.md) bleiben maßgeblich.
Eigene Messwerte mit Quellhashes stehen im
[Nachweisarchiv](evidence/m12/v1.7.0.dev1-rans-verification.json).

## Codeverifikation

Die unabhängig differenzierte SA-MMS erreicht auf 16/32/64/128 Zellen je
Richtung relative RMS-Fehler 9.50617e-4, 2.48245e-4, 6.36585e-5 und
1.61304e-5; die beobachteten Ordnungen sind 1.937, 1.963 und 1.981.
Alle Zustände sind nichtnegativ, die frischen SA-Residuen erfüllen das
absolute Gate. Der virtuelle MMS-Wandabstand ist kein Wandmodellnachweis.

Der kontinuierliche Halbkanal-BVP integriert die Modellgleichungen mit eigener
Kollokation und ohne Produktions-CFD-Importe. Zweite Initialisierung und
strengere Toleranz bestätigen denselben nichttrivialen Ast. Die numerische
Lösung benutzt diesen BVP ausschließlich zur Fehlerauswertung.

| Zellen/Halbkanal | relativer RMS-U-Fehler | Bulk U [m/s] | maximales y+ |
| --- | --- | --- | --- |
| 32 | 0.00347335 | 18.6450185 | 0.498003 |
| 64 | 0.00187375 | 18.6119767 | 0.240251 |
| 128 | 0.000212285 | 18.5807321 | 0.118014 |
| 256 | 0.0000881530 | 18.5784062 | 0.0584885 |

Der feinste Fehler beträgt 0.00882 % gegenüber dem Gate 0.5 %.
Wandkraftfehler liegen unter 1.8e-13 relativ (Gate 0.1 %).
Viskosität +/-5 %, Wandverteilung beta=5 statt 4.5 und Wiederholung mit
SA-Relaxation 0.08 bestehen die unveränderten Feld-/Bilanzgates.
Die Wandverteilungsänderung verändert den Bulk um ungefähr 0.006 %.
Eine geometrisch ähnliche SI-Variante mit H=0.1 m, u_tau=0.01 m/s,
Uref=0.2 m/s und Lref=0.1 m erreicht auf 32 Zellen/Halbkanal
0.336 % relativen RMS-Fehler und y+=0.498. Ein unabhängiger
BVP-Skalierungstest prüft außerdem bei festem Volumenstrom und konstantem nu:
H'=gamma H, G'=G/gamma³, U'=U/gamma, a'=a.
Das ist eine Vorbereitung der FSI-Referenz, kein Kopplungsnachweis.

## Unabhängiger Modellvergleich

Die Originalquelle ist [Moser/Kim/Mansour (1999)](https://turbulence.oden.utexas.edu/MKM_1999.html),
[Profil chan590.means](https://turbulence.oden.utexas.edu/data/MKM/chan590/profiles/chan590.means).
Der Header nennt Re_tau=587.19 und Normierung U_tau,h. Der Rundname 590
wird nicht als exakter Reynoldsparameter übernommen. Herkunft, Downloadtag,
Dateigröße, Spalten und SHA-256 stehen im [Quellenmanifest](evidence/m12/reference_sources.json).
Originaldaten werden wegen fehlender expliziter Weitergabelizenz weder im
Repository noch in Paketen, Handbucharchiven oder Release-Assets verbreitet.
Die dedizierte CI lädt sie temporär und verlangt den unveränderten Hash.
Ein lokaler Lauf ohne M12_DNS_PROFILE überspringt diesen Test sichtbar.

| Zellen/Halbkanal | RMS-U+-Abweichung | Bulkabweichung |
| --- | --- | --- |
| 128 | 1.07751 % | 0.39493 % |
| 256 | 1.07465 % | 0.40740 % |

Beide Grenzen von 5 % sind erfüllt. Die letzte Profiländerung beträgt
0.01628 %, unter 1 % und unter einem Fünftel des Modellbudgets.
Der DNS-Bulk ist 18.6544039 U_tau; das Modell ist nicht an diese Daten angepasst.
Interpolation benutzt ausschließlich numerische Geschwindigkeiten, tatsächliche
Wandschubspannung und Wand-/Symmetriebedingungen. Der Zentrumwert wird aus
drei benachbarten numerischen Reihen symmetrisch extrapoliert.

Ein erster Vergleich scheiterte durch einen Auswertungsfehler: letzte
Gleitkommastellen gleicher y-Reihen erzeugten künstliche Zusatzreihen und
einen schlecht konditionierten Interpolanten. Ein parabolischer Regressionstest
mit Koordinatenrauschen reproduzierte dies. Reihen werden nun allein innerhalb
128 eps mal geometrischer Längenskala volumenbezogen zusammengefasst;
physisch verschiedene Wandlagen bleiben getrennt. Daten, Solver, Modell und
Abnahmegrenzen wurden nicht geändert. Fünf Offline-Quellen-/Interpolationstests
und der erneute vollständige Modellvergleich bestehen.

## Aussagegrenzen und Reproduktion

Dieser DNS-Vergleich beurteilt das Modell für den glatten voll entwickelten
Kanal bei Re_tau=587.19. Er ist keine experimentelle Validierung und deckt
weder Transition noch Rauheit, starke Ablösung oder allgemeine 3D-Strömung ab.
Die eigenständige BVP-Parität prüft den Code, nicht die physikalische Modellgüte.
FSI, Netzbewegung und deren unabhängige Referenzen bleiben eigene Gates.
Unabhängiges menschliches Fachreview liegt nicht vor.

Gezielte Läufe: tests/convergence/test_sa_transport_convergence.py,
test_rans_channel.py und test_rans_dns.py, jeweils mit explizitem
pytest-Marker cfd_convergence. Der DNS-Lauf benötigt die gehashte lokale
Quelldatei über M12_DNS_PROFILE. Unit-Tests laden keine Netzressourcen.
Die dedizierte CI führt den Quellenvergleich auf allen unterstützten
Python-Versionen zusätzlich zu den bisherigen CFD-Netzfolgen aus.
