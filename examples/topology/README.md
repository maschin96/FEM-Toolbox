# Struktur-Topologieoptimierung (M11 / v1.6.0)

Nach Paketinstallation: `python -I examples/topology/run.py --output results/topology`.
Die vorbereitete TOML-Studie rechnet zwei Lastfälle eines 2x1x1-m-Quaders,
minimiert gefilterte SIMP-Nachgiebigkeit bei 60% Dichtevolumen und rechnet
Schwellengeometrien mit Vollmaterial erneut. Neue Ausgabeordner sind Pflicht.

`python -I examples/topology/verify.py --output results/topology-verification`
prüft drei Netze und zwei Radien, Wiederholung, Einzellastvariante und
Gradienten; `--quick` ist ausschließlich ein Smoke-Test.
`prepare.py --output <neuer Ordner> --divisions 3 --radius 0.6` erzeugt andere
vorbereitete Netze/Modelle deterministisch.

[API, Formate, Einheiten, Fehler und Grenzen](../../Doc/84_M11_TOPOLOGY_WORKFLOW.md).
[Vorab festgelegte Referenzen/Toleranzen](../../Doc/83_M11_TOPOLOGY_CONTRACT.md).
[Verifikation und Freigabestatus](../../Doc/85_M11_VERIFICATION.md).
Keine globale Optimalität, garantierte Mindestwandstärke oder Fertigungsfreigabe.
