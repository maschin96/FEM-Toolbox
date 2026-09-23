# Kleine thermische Parameterbibliothek

Stand: 23.09.2026; Katalogversion 1; Teil von #29. Die Bibliothek arbeitet offline.

| Schlüssel | Material/Zustand | Dichte kg/m³ | c J/(kg K) | k W/(m K) | Quelle |
|---|---|---:|---:|---:|---|
| `aluminium_6060_t6_20c` | EN AW-6060, T6, 20 °C | 2700 | 901 | 200 | [Hydro, Rev.2 01/2019, S.2](https://www.hydro.com/globalassets/08-about-hydro/hydro-worldwide/germany/extrusion-germany/alloy-data-sheets/hydro-en-aw-6060.pdf) |
| `copper_cda_20c` | generisches Referenzkupfer, 20 °C | 8920 | 386 | 394 | [CDA, The Copper Advantage, Tabelle 3, S.9](https://www.copper.org/publications/pub_list/pdf/a1360.pdf) |
| `stainless_14301_rt` | EN 1.4301 / 304, Raumtemperatur | 7900 | 500 | 15 | [Outokumpu, Steel Grades, Properties and Global Standards, Tabelle 13, S.13](https://www.outokumpu.com/-/media/from-sharepoint/documents/product/outokumpu-steel-grades-properties-global-standards.pdf) |

Hydro nennt typische Kennwerte bei 20 °C und den T6-Zustand für die
Wärmeleitfähigkeit. Der am Recherchetag indizierte Herstellertext enthält
2,70 g/cm³, 901 J/(kg K), 200 W/(m K). Der direkte PDF-Abruf ergab 404;
diese eingeschränkte Quellenverfügbarkeit ist auch im Katalogeintrag vermerkt.
Die Dichte wird mit Faktor 1000 in kg/m³ umgerechnet.

CDA Tabelle 3 nennt in metrischer Darstellung 8,92 g/cm³, 0,386 J/(g °C)
und 3,94 W·cm/(cm² °C). Umrechnungsfaktoren sind 1000, 1000 und 100.
Der Eintrag ist ausdrücklich kein vollständiger, chargenzertifizierter
C11000-Datensatz. Outokumpu nennt für 1.4301 bei RT 7,9 kg/dm³,
500 J/(kg °C), 15 W/(m °C); nur die Dichte benötigt Faktor 1000.
RT wird in den Metadaten als 293,15 K repräsentiert.

Alle drei Einträge sind **typische konstante Referenzwerte**. Die Quellen
begründen kein allgemeines Temperaturintervall, keine Unsicherheitsverteilung
und keine Temperaturabhängigkeit. Insbesondere werden die Werte nicht für
Phasenwechsel oder große Temperaturbereiche extrapoliert. Der Solver hält die
Parameter konstant; die Eignung für den konkreten Temperaturbereich bleibt
am Materialdatenblatt zu prüfen. Die Bibliothek enthält thermische Parameter,
keinen frei erfundenen piezoelektrischen Volltensor.

```python
from fem_toolbox.materials.catalog import material_catalog, thermal_material

for entry in material_catalog():
    print(entry.key, entry.source_url)
material = thermal_material("aluminium_6060_t6_20c")
# material ist direkt ein ThermalMaterial, einschließlich density/specific_heat.
```

Eigene Daten lassen sich mit
`ThermalMaterial.isotropic(k, density=rho, specific_heat=c)` oder mit einem
anisotropen 3×3-Leitfähigkeitstensor im Konstruktor verwenden. Dichte und
spezifische Wärmekapazität müssen gemeinsam positiv und endlich vorliegen.
Stationäre Konstruktoren ohne Speicherparameter bleiben unverändert zulässig.
`assign_region_materials` kann diese Materialien weiterhin Regionen zuordnen;
seine Provenienz enthält Speicherparameter und Quellenidentität.
