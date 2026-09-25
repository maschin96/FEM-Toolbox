# Tet10-/Triangle6-Topologie (#35a)

M7 erweitert das Netzmodell additiv um eine quadratische Simplex-Familie.
Dieser Teil implementiert Topologie, Metadaten und I/O; quadratische
Formfunktionen und Quadratur ergänzt [#35b](41_TET10_GEOMETRY.md).
Operatoren und die Model-Integration folgen in #35c/d.

## Datenvertrag

`Mesh(points=..., tetra10=..., triangle6=...)` verwendet nullbasierte Indizes.
Tet10-Zeilen haben zehn Einträge: vier Ecken 0,1,2,3, dann Kantenknoten
01,12,20,03,13,23. Triangle6 verwendet drei Ecken 0,1,2, dann 01,12,20.
Dies entspricht der internen VTK-/meshio-Reihenfolge. Die abweichende rohe
Gmsh-Tet10-Reihenfolge wird genau einmal durch
[meshio 5.3.5](https://github.com/nschloe/meshio/blob/v5.3.5/src/meshio/gmsh/common.py)
konvertiert. Eine handgeschriebene Gmsh-4.1-Datei prüft diese Grenze unabhängig
vom Exportadapter.

Die neuen Arrays und optionalen `tetra10_ids`/`triangle6_ids` sind wie bisher
kopiert und unveränderlich. `cells`, `facets`, `cell_ids`, `facet_ids`,
`cell_type`, `facet_type`, `n_cells` und `n_facets` adressieren die aktive Familie.
`cell_tags`/`volume_sets` beziehen sich auf deren Volumenzeilen,
`facet_tags`/`facet_sets` auf deren Oberflächenzeilen. Überlappende Gruppen
bleiben erhalten. Materialzuweisungen verwenden dieselben Volumenzeilen.

Tet4-Konstruktoren bleiben kompatibel; nicht verwendete Konnektivitäten haben
null Zeilen mit der jeweiligen Breite. Ein Netz darf nur eine Familie enthalten,
auch bei räumlich getrennten Teilnetzen. Gemischte Tet4/Tet10- und
Triangle3/Triangle6-Netze werden ausdrücklich abgelehnt. Gemeinsame Kanten müssen
denselben Kantenknoten besitzen; ein Kantenknoten darf weder Ecke noch Knoten
einer anderen Kante sein. Doppelte Eck-Simplizes sind ungültig.

## Geometrie und Randtopologie

`Mesh` speichert Topologie; der Import und die Randtopologie prüfen zusätzlich
mittige Kantenknoten. Die Abweichung darf höchstens 1e-10 der jeweiligen
Kantenlänge betragen. Relative Koordinaten vermeiden eine Aufweitung durch
Translation. Gekrümmte Geometrie wird abgelehnt. Volumen, Flächen und Normalen
verwenden danach die bestehende affine Eckgeometrie; dies führt keine
quadratische Integrationsregel ein.

`inspect_boundary_topology` gibt bei Tet10 alle sechs orientierten Flächenknoten
zurück und prüft deren Zuordnung zum eindeutigen Besitzer. Die Inzidenzschlüssel
bleiben die drei sortierten Eckindizes. Innere, fremde, nichtkonforme und
nichtmannigfaltige Randflächen sind ungültig; vollständige Randabdeckung bleibt
optional. Eine negative Gmsh-Volumenorientierung wird mit allen Kantenknoten
konsistent umgeordnet.

## Import und Export

`read_gmsh` unterstützt Tet10/Triangle6 inklusive SI-Skalierung und physischen
Gruppen. Originale Gmsh-Entitätsnummern stellt meshio weiterhin nicht als
Toolbox-IDs bereit. `write_vtu` schreibt aktive Zelltypen, Ergebnisarrays und
Toolbox-Metadaten. `read_vtu` liest das Netz ohne Orientierungsänderung zurück:
optionale externe IDs, Tags sowie überlappende und leere Auswahlsätze werden
für vorhandene Zelltypen erhalten. Präsenzmasken unterscheiden fehlende Werte
von numerischem Null. Metadaten eines Zelltyps ohne Zeilen sind nicht übertragbar.

Ergebnisarrays liest man weiterhin über meshio; `read_vtu` liefert ausschließlich
`Mesh`. Der Präfix `fem_toolbox:` ist für Metadaten reserviert. Fremde VTU-Dateien
ohne diese Metadaten liefern ein Netz ohne entsprechende IDs oder Gruppen.

## Verifikation und Grenzen

Unit-Tests prüfen unveränderliche Speicherung, zwei benachbarte Tet10, gemeinsame
Kanten, orientierte Sechsknotenflächen, Materialregionen, Gmsh 2.2/4.1,
Orientierungswechsel, Einheiten und VTU-Rundreisen einschließlich Tet4.
Ungültige Mischfamilien, gekrümmte Kanten und falsche Randzuordnungen werden
abgelehnt. Die bisherige Tet4-Suite bleibt der Kompatibilitätsnachweis.

Tet10 ist noch kein berechenbares Model. `Model` lehnt quadratische Netze früh
ab; bestehende Tet4-Solver erhalten keine still reduzierte Konnektivität.
Patch-, Energie- und Konvergenznachweise quadratischer Operatoren gehören zu
#35b/c und sind durch diese Topologietests nicht erbracht.
