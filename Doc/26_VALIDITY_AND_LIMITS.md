# Gültigkeitsbereich und numerische Nutzungsgrenzen für v1.0 (#142)

Die stationären Grenzen bleiben in v1.1.0 gültig. Für die additiven transienten
und dynamischen Modi gilt [Dokument 35](35_M6_RELEASE_VERIFICATION.md).
Die folgende Tabelle beschreibt ausdrücklich den stationären v1.0-Umfang.

Status: Gültigkeitsbereich der **stabilen stationären Forschungs-API v1.0.0**.
Vor einem Lauf sind Modellannahmen, Randbedingungen und die unten genannten
Prüfungen gegen den konkreten Anwendungsfall abzugleichen. Die
[Release Notes](releases/v1.0.0.md) benennen den Freigabenachweis.

## Entscheidungsübersicht

| Frage | Innerhalb des dokumentierten Bereichs | Außerhalb oder noch ohne Nachweis |
|---|---|---|
| Physik | Stationäre Wärme `T`, lineare statische Elastizität `u`, Elektrostatik `φ`, sequenziell `T → u`, monolithisch quasistatisch `u ↔ φ`, sequenziell `T → (u, φ)` | Transienten, Dynamik, harmonische Felder, Rückwirkung von `u`/`φ` auf `T`, nichtlineare Stoffe oder große Verformungen |
| Raum und Netz | 3D, gerade lineare Tet4, konforme Triangle3-Randflächen; ein gemeinsames Netz für gekoppelte Felder | Tet10/Hex8, gekrümmte oder nichtkonforme Grenzflächen, Kontakt und automatische Vernetzung |
| Materialien | Stückweise konstante lineare isotrope oder anisotrope Tensoren je Tet4; positive zulässige Stoffwerte; orthonormale Materialachsen | Ortsabhängigkeit innerhalb eines Elements, Hysterese, Plastizität, Bruch und experimentell unbestätigte Werkstoffparameter als Vorhersagegrundlage |
| Grenzflächen | Perfekte Bindung und stetige Primärfelder auf konformem Netz | Grenzflächenwiderstand, Delamination, Kontakt, Sprünge der Primärfelder |
| Thermische Kopplung | Temperatur zuerst, dann thermische Eigenverzerrung in der Mechanik; in `T → (u, φ)` elektrischer thermischer Beitrag nur mit `electric_thermal_strain=true` | Voll reziproke oder zeitabhängige thermische Kopplung |
| Einheiten und Vorzeichen | Intern SI; Längenimport mit `length_scale_to_m`; Temperatur K oder explizit umgerechnete `degC`; Stress-Charge und Engineering-Voigt `xx,yy,zz,yz,xz,xy` | Implizit angenommene Gmsh-Längeneinheit, unmarkierte relative Permittivität oder andere Tensor-/Polarisationskonvention |
| Solver | SciPy Sparse Direct/SuperLU; bei Piezo feldweise Blockskalierung voreingestellt; feste oder schwebende Equipotential-Elektroden | PETSc, verteiltes Rechnen und eine garantierte Skalierung auf große Modelle |

Das [Modellschema](21_MODEL_SCHEMA.md) enthält die genauen Eingabefelder und
Einheiten. [Mathematisches Modell](02_MATHEMATICAL_MODEL.md),
[Materialregionen](19_MATERIAL_REGIONS.md),
[Elektroden](20_ELECTRODES.md) und
[Thermopiezoelektrizität](22_THERMOPIEZOELECTRIC.md) definieren die zugehörigen
Konventionen. Die Modellprüfung erkennt unter anderem fehlende Referenzen,
Starrkörpermoden, widersprüchliche Regionen und Randbedingungen
([Preflight](23_MODEL_PREFLIGHT.md)); sie beweist keine Eignung des Modells.

## Numerische Grenzen und konkrete Prüfung

| Grenze | Mögliche Auswirkung | Erkennbares Signal | Empfohlene Prüfung |
|---|---|---|---|
| Lineare Tet4 bei Biegung oder starken Gradienten | Zu steife Antwort, langsame Konvergenz | Verschiebungen, Energien oder Reaktionen ändern sich bei Verfeinerung deutlich | Mindestens drei systematisch verfeinerte, formähnliche Netze vergleichen; für Biegung gegen eine unabhängige Referenz oder ein höherwertiges Element prüfen |
| Nahe Inkompressibilität | Volumetrisches Locking und zu steife Verschiebung | Bei isotropem Material wird `ν ≥ 0,49` abgelehnt; darunter kann Verfeinerung dennoch unzureichend sein | Abstand zu `ν=0,5` und Netzkonvergenz dokumentieren; für diesen Bereich eine geeignete gemischte Formulierung verwenden |
| Flache oder entartete Tet4 | Große Gradientenfehler oder instabile Matrix | Geometrieprüfung lehnt Entartung ab; Qualitätsindikatoren verschlechtern sich vor der Grenze | Qualitätsbericht nach [Tet4-Formulierung](03_TET4_FORMULATION.md) prüfen und Netz verbessern; Bestehen der Mindestgrenze allein genügt nicht |
| Punktlasten, Ecken und scharfe Material-/Randwechsel | Lokale Singularität oder nicht konvergenter Spitzenwert | Peak wächst bei Netzverfeinerung, während integrale Größen stabilisieren | Flächen- statt Punktlast verwenden, wenn physikalisch passend; Mittel-/Integralwerte und Abstand zur Singularität berichten |
| Rohe Zellspannung und Zellflüsse | Sprünge über Elementgrenzen; lokale Maxima sind netzabhängig | Nachbarzellen haben unterschiedliche Werte; nodale Darstellung glättet sichtbar | Rohwerte auf mehreren Netzen und Bilanzen vergleichen; Mittelung nur als Visualisierung kennzeichnen |
| Schlechte Stoffskalierung oder unzureichende Lagerung/Gauge | Singuläres oder schlecht konditioniertes System, unzuverlässige Felder | Preflight-/LU-Fehler, große feldweise Residuen oder Bilanzfehler | Lagerung, Potentialreferenz und Einheiten prüfen; Blockskalierung verwenden; Residuen **und** physikalische Bilanzen auswerten |
| Nicht validierte Stoffdaten oder ideale Grenzfläche | Rechnerisch konvergentes, physikalisch falsches Ergebnis | Sensitivität gegenüber Stoffwerten oder Grenzflächenannahme | Materialherkunft und Unsicherheiten erfassen; für reale Vorhersagen experimentell validieren |

Die Akzeptanzgrenzen einzelner Verifikationstests gelten für **deren**
synthetische Fälle, nicht pauschal für andere Geometrien, Lasten oder
Materialien. Für einen neuen Fall sind Zielgröße, Fehlermaß und erforderliche
Netzkonvergenz vor der Auswertung festzulegen. Ein kleines algebraisches
Residuum allein ist kein Genauigkeitsnachweis.

## Beleglage und Ressourcen

Patch-, Bilanz- und analytische Tests der Einzelfelder stehen in
[Elastizität](10_LINEAR_ELASTICITY.md),
[Wärme](11_STEADY_HEAT.md) und
[Elektrostatik](12_ELECTROSTATICS.md), die Piezo-Referenzen in
[Doc/17](17_PIEZOELECTRIC_SOLVER.md) und
[Doc/18](18_PIEZO_D33_VERIFICATION.md). Der gekoppelte
[Vier-Netze-Nachweis](25_THERMOPIEZO_CONVERGENCE.md) prüft `n=2,4,6,8`
(48–3072 Tet4) mit einer unabhängig integrierten geschichteten 1D-Lösung
auf einem 3D-Netz. Für die dortigen glatten Felder liegen die gemessenen
Konvergenzordnungen innerhalb der vorab gesetzten Fenster `L²: 1,75–2,25`,
`H¹: 0,75–1,25`; Wärme-, Kraft-, Moment- und Ladungsbilanz bleiben unter
`10⁻⁸` in der jeweiligen SI-Einheit. Der Nachweis deckt insbesondere keine
gekrümmten Grenzflächen, allgemeine exakte Quergradienten oder reale
Materialdaten ab.

Zusätzliche **Ressourcenbeobachtung** am 22.09.2026 auf macOS 27.0,
ARM64, 10 logischen CPUs, Python 3.13.7, NumPy 2.4.6, SciPy 1.18.1:

| Reguläres Netz | Tet4 | Knoten | `Model.solve()` [s] | Spitzen-RSS des Prozesses [MiB] |
|---:|---:|---:|---:|---:|
| `n=2` | 48 | 27 | 0,028 | 68,8 |
| `n=8` | 3072 | 729 | 1,519 | 211,2 |

Reproduktion aus dem Repository-Root mit installierten Entwicklungsabhängigkeiten;
je Netz ein frischer Prozess. Der `Model`-Aufbau erfolgt vor der Zeitmessung,
ist aber im Spitzen-RSS enthalten. `ru_maxrss` wird auf macOS in Byte
geliefert; angegeben ist Byte / 2²⁰. Messbefehl (für `n=2` oder `8`):

```sh
PYTHONPATH=src .venv/bin/python -c 'import resource,sys,time; from tests.conftest import cube_mesh; from tests.convergence.test_thermopiezo_model_convergence import _model; n=int(sys.argv[1]); model=_model(cube_mesh.__wrapped__(),n,voltage=0.7); start=time.perf_counter(); model.solve(); print(f"n={n} solve_s={time.perf_counter()-start:.3f} peak_rss_mib={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1048576:.1f}")' 2
```

Dies sind Einzelmessungen eines synthetischen Falls, keine Laufzeitprognose,
RAM-Obergrenze oder allgemeine Freigabe für 3072 Elemente. Sparse-LU-Fill-in
hängt von Netz, DOFs, Material und Randbedingungen ab. Für größere Fälle
Speicher und Laufzeit auf der Zielhardware separat messen.

## Nachweisstatus vor einer realen Anwendung

**Code-Verifikation** prüft die Umsetzung der Gleichungen;
**Lösungsverifikation** prüft die diskrete Genauigkeit des konkreten Falls.
Beide sind von **Validierung** gegen Messdaten zu unterscheiden
([V&V-Konzept](05_VERIFICATION_AND_VALIDATION.md)). Die genannten Tests und
synthetischen Referenzen stützen Code- und Lösungsverifikation innerhalb ihres
Testbereichs. Eine experimentelle Validierung realer Werkstoffe oder Bauteile
ist für v1.0 nicht beansprucht.

Vor einer physikalischen Schlussfolgerung: Modellannahmen und SI-Eingaben
prüfen, Netzqualität und Randbedingungen kontrollieren, mehrere Netze für die
Zielgrößen vergleichen, feldweise Residuen sowie Wärme-/Kraft-/Ladungsbilanzen
auswerten und bei realen Anwendungen Material- und Messunsicherheit erfassen.
