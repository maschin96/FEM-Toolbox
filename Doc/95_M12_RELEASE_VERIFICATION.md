# M12 / v1.7.0 — Implementierungs- und Freigabenachweis (#94)

19.11.2026 · Freigabevorbereitung; tatsächlicher Abschluss in
[Release-Issue #94](https://github.com/maschin96/FEM-Toolbox/issues/94).
Diese Datei behauptet keinen erfolgten Tag und keine abgeschlossene
Veröffentlichung.

## Umfang und Versionsentscheidung

M12 liefert ebene stationäre inkompressible isotherme Strömung auf dem
Quad-FV/SIMPLE-Pfad, die SA-noft2-Schließung mit aufgelöster Wand und eine
stationäre bidirektionale Dirichlet–Neumann-FSI mit spannweiteninvarianter
extrudierter Tet4-Kleinverformungsstruktur. Der 2D-MAC-Kern bleibt Referenz.

Version: **1.7.0**. M12 folgt auf v1.6.0; alle M12-Einstiegspunkte
(`fem_toolbox.cfd`, `fem_toolbox.coupling`) sind additiv, bestehende Solver,
Model-Schema 1–3 und Persistenzformate bleiben unverändert. FluidDomain-,
FluidMesh-, FlowProblem- und `fsi_schema_version` 1 sind eigenständige neue
Formate. Das begründet einen kompatiblen Minor-Schritt. Die M12-Module
bleiben als Forschungs-API klassifiziert (`api-map.json`, Doc/29) und sind
keine stabile SemVer-Oberfläche; versioniert sind ihre Dateischemata. Das vorläufige
Planungsziel v1.11.0 der Planungswelle 3 beruhte auf einer Releasefolge
M8–M11 mit eigenen Versionen; nach der tatsächlichen Folge (M11 = v1.6.0)
wäre jede höhere Nummer ohne inkompatible Änderung nicht begründet.

## Integrierte Pflichtissues

| Issue | Inhalt | Integration | Nachweis |
| --- | --- | --- | --- |
| #84 | RANS-/FSI-Vertrag | [PR #227](https://github.com/maschin96/FEM-Toolbox/pull/227) | [Doc/71](71_M12_RANS_FSI_CONTRACT.md) |
| #85 | Quad4-Fluidnetz, Wandqualität, BC-Preflight | [PR #228](https://github.com/maschin96/FEM-Toolbox/pull/228) | [Doc/72](72_FLUID_MESH_ADAPTER.md) |
| #86 | stationärer MAC-Referenzkern | [PR #223](https://github.com/maschin96/FEM-Toolbox/pull/223) | [Doc/70](70_STATIONARY_FLOW_CORE.md) |
| #226 / #229 | Quad-FV-Grundoperatoren | [PR #231](https://github.com/maschin96/FEM-Toolbox/pull/231) | [Doc/73](73_QUAD_FV_OPERATORS.md) |
| #226 / #230 | Quad-FV-SIMPLE und Strömungsabnahme | [PR #232](https://github.com/maschin96/FEM-Toolbox/pull/232) | [Doc/74](74_QUAD_FV_SIMPLE.md) |
| #87 | SA-noft2 und gemeinsame RANS-Iteration | [PR #289](https://github.com/maschin96/FEM-Toolbox/pull/289) | [Doc/88](88_SA_NOFT2_IMPLEMENTATION.md) |
| #88 | RANS-Codeverifikation, DNS-Modellvergleich | [PR #290](https://github.com/maschin96/FEM-Toolbox/pull/290) | [Doc/89](89_RANS_VERIFICATION.md) |
| #89 | konservative Lastübertragung | [PR #291](https://github.com/maschin96/FEM-Toolbox/pull/291) | [Doc/90](90_CONSERVATIVE_FSI_INTERFACE.md) |
| #90 | Netzbewegung und Geometrieprüfung | [PR #292](https://github.com/maschin96/FEM-Toolbox/pull/292) | [Doc/91](91_STATIONARY_MESH_MOTION.md) |
| #91 | stationäre bidirektionale Kopplung | [PR #294](https://github.com/maschin96/FEM-Toolbox/pull/294) | [Doc/92](92_STATIONARY_FSI_COUPLING.md) |
| #92 | FSI-Benchmark und Bilanzen | [PR #295](https://github.com/maschin96/FEM-Toolbox/pull/295) | [Doc/93](93_STATIONARY_FSI_VERIFICATION.md) |
| #93 | deklarativer Workflow und Provenienz | [PR #296](https://github.com/maschin96/FEM-Toolbox/pull/296) | [Doc/94](94_FSI_WORKFLOW.md) |

[PR #288](https://github.com/maschin96/FEM-Toolbox/pull/288) und dieser
Release-PR bringen M12 auf den stabilen main-Stand (v1.6.0 und #287).

## Gates (vorab festgelegt in Doc/71)

- **CFD-Codeverifikation:** Vier-Netz-Poiseuille, nichtlineare MMS und
  Variable-mu-Kanal mit unveränderten Budgets und MAC-Parität (Doc/74);
  SA-MMS 16–128 Zellen, beobachtete Ordnung 1.94–1.98 (Doc/89).
- **RANS-Modell:** unabhängiger SA-Halbkanal-BVP, feinster RMS-U-Fehler
  0.0088 % (Gate 0.5 %); DNS Moser/Kim/Mansour Re_τ = 587.19: RMS-U⁺
  1.07 %, Bulk 0.41 % (Gates je 5 %). Netz-/Wandsensitivität: Viskosität
  ±5 %, Wandverteilung β = 5, SA-Relaxation 0.08, y⁺ ≤ 0.5.
- **Lastübertragung:** Transferaudit und virtuelle Arbeit ≤ 1e-10 relativ,
  Strukturknotenkräfte in N, Kleinverformung ≤ 0.02 (Doc/90).
- **Netzbewegung:** harmonische Verschiebung, Nichtorthogonalität ≤ 60°,
  Skewness ≤ 0.25, Flächenverhältnis ≥ 0.5, Diagnose statt Mittelung (Doc/91).
- **Bidirektionale FSI:** geführter Deckel, je vier Netzpaare. Laminar
  Ordnung 2.0, feinster Weg −0.048 %, bidirektionaler Effekt 1.7471 % gegen
  1.7479 % analytisch; SA feinster Weg −0.119 %, Effekt 2.8553 % gegen
  2.8585 %, y⁺ ≤ 0.5 (Doc/93). Massen-, Kraft-/Reaktions-, Arbeits- und
  Schnittstellenbilanzen pro Netzpaar geprüft.
- **Workflow:** beide Beispiele aus Dateien, Wiederholung exakt,
  Abweichung zu #92 ≤ 2e-9 m, CI-Ausführung aus Wheel und sdist (Doc/94).

## Finale lokale Prüfung von Version 1.7.0

Aus einem sauberen `git archive` des Release-Commits, frische Umgebung,
Python 3.13, macOS arm64:

- Suite ohne die CI-Konvergenzjobs: **1205 bestanden, 83 optionale Tests
  übersprungen**. Drei Baseline-Worker benötigen `git rev-parse` und laufen
  nur im Repository (CI); zwei Handbuch-Abdeckungsfehler für die neuen
  Releasedokumente wurden behoben.
- Ruff, Formatprüfung, striktes mypy über 95 Quelldateien und API-Vertragsprüfung bestanden.
- Wheel und sdist getrennt frisch installiert: isolierter Import, Version
  1.7.0, MIT, NumPy 2.5.5/SciPy 1.18.1; beide `examples/fsi_lid`-Fälle
  konvergieren aus der Installation.
- Konvergenzstudien (#88, #92) und DNS-Vergleich laufen als eigene CI-Jobs.
- [Lokale Belege mit Quellhashes](evidence/m12/v1.7.0-local-verification.json).
  Vorabprüfungen, nicht die später aus main erzeugten Veröffentlichungsartefakte.

## Bekannte Grenzen

Nur ebene 2D-Strömung, stationär, inkompressibel, isotherm; SA-noft2 mit
aufgelöster Wand (y⁺ ≤ 1), keine Wandfunktionen. FSI nur mit expliziter
Breite, spannweiteninvarianter Spur und kleinen Verformungen; konstante
Relaxation. Keine allgemeine 3D-Fluid-, URANS-, LES-/DNS-, dynamische FSI-,
kompressible oder Mehrphasenfreigabe. Der SA-Vergleich ist ein
Modellvergleich mit DNS, keine experimentelle Validierung; die DNS-Daten
werden nicht verbreitet.

## Reviewstatus und Integrationsgate

Unabhängiges menschliches Fachreview liegt für M12 nicht vor und wird nicht
behauptet. #230 nennt ein unabhängiges Review vor Integration als
Akzeptanzkriterium; die technischen Gates von #226/#229/#230 sind belegt
(Doc/73/74, CI), das Review fehlt. Gemäß Doc/09 wird das hier benannt und
nicht als Freigabe gewertet. Der Abschluss von #226/#230 mit diesem Vermerk
und der Tag v1.7.0 erfordern eine ausdrückliche Eigentümerentscheidung.

Vor dem main-Merge: aktuelle Release-PR-CI vollständig grün. Danach main-CI
am exakten Zielcommit prüfen; erst dann annotierter Tag v1.7.0 und
Veröffentlichung. #94 wird zuletzt geschlossen; Commit-, Lauf-, Tag- und
Artefaktlinks werden dort gebunden.
