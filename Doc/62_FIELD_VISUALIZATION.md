# Stationäre Tet4-Feldbilder (#189)

Die optionale öffentliche API liegt in `fem_toolbox.post.visualization`:
`plot_field(mesh, values, *, association, name, unit, nodes=None,
component=None, limits=None, cmap="viridis", view=(25, -60), figsize=(8, 6))`
liefert `FieldPlot` mit `figure`, `axes`, `artist`, `faces`, `face_values` und
`close()`. `save_jpeg(plot, path, *, dpi=150, overwrite=False)` speichert separat.
Figure/Axes können im eigenen Skript weiter gestaltet werden. Der Aufrufer
schließt den Plot in `finally`; es gibt kein implizites `show()`.

## Installation und Umfang

```sh
python -m pip install '.[visualization]'
python examples/field_visualization/run.py --output field_images
```

Bei Installation aus einem Paketarchiv das Extra `visualization` an dessen
Pfad anhängen. Matplotlib >=3.9,<4 ist optional; Solver und der Import des
Adapters benötigen es nicht. Bei fehlender Installation nennt `plot_field`
den Installationsbefehl. Die CI prüft Kernartefakte ohne Matplotlib und
Rendering auf Python 3.11–3.13 mit Extra. Lokal geprüft: Python 3.13 und
Matplotlib 3.9.4 und 3.11.2. Matplotlib verwendet eine BSD-kompatible Lizenz;
[Original-Lizenz](https://matplotlib.org/stable/users/project/license.html).
Rendering verwendet den expliziten
[Agg-Canvas](https://matplotlib.org/stable/api/backend_agg_api.html), ohne globale
Backend- oder Stiländerungen. Auf Systemen ohne beschreibbaren Benutzer-Cache
kann `MPLCONFIGDIR` auf ein beschreibbares Verzeichnis gesetzt werden.

| Eingabe | Auswahl | Zuordnung |
|---|---|---|
| T [K], phi [V], andere reelle Skalare | `component=None` | Knoten oder Zellen |
| u [m], q [W/m²], E [V/m], D [C/m²] | x, y, z, magnitude | Knoten oder Zellen |
| Spannung [Pa] | xx, yy, zz, yz, xz, xy | Zellen, vorhandene Voigt-Reihenfolge |

Nur reelle stationäre Tet4-Daten auf unverformter Geometrie werden zugesagt.
Tet10, komplexe Felder, Zeit-/Frequenzserien, Schnitte und Volumenrendering sind
nicht unterstützt. Der generische Arrayadapter kennt keine Ergebnisart;
der Aufrufer verantwortet Feldname, Einheit und stationäre Datenherkunft.
Der numerische Kern erhält keine neue Abhängigkeit oder Interpolation.

## Feldvertrag und Darstellung

`association="node"`: Zeilen folgen `nodes`, standardmäßig den sortierten
aktiven Volumenknoten. Jeder aktive Knoten muss genau einmal vorkommen;
inaktive Knoten und Export-Padding werden abgelehnt. Solver liefern die
Zuordnung als `node_ids`. Bei `association="cell"` folgen Zeilen den Tet4-Zellen;
`nodes` ist dann unzulässig. Skalare haben Form `(n,)`, Vektoren `(n,3)` und
Voigt-Felder `(n,6)`. Werte müssen reell und endlich sein.

Außenflächen werden aus der Volumentopologie gewonnen, auch bei fehlenden
Randdreiecken; gemeinsame innere Facetten bleiben unsichtbar. Eingetragene
Randflächen unterliegen dem bestehenden Topologie-Preflight.
Nodale Komponenten oder Beträge werden zuerst je Knoten ausgewertet und dann
arithmetisch über die drei Flächenknoten gemittelt. Dies ist eine ausdrücklich
benannte Darstellungsreduktion, kein exakter kontinuierlicher Feldverlauf.
Zellfarben stammen unmittelbar vom Eigentümerelement, ohne Materialglättung.

Achsen zeigen Meter, Geometrieproportionen bleiben erhalten. `view` enthält
Elevation und Azimut in Grad; `figsize` Breite/Höhe in Zoll. Colorbar und Titel
nennen Feld, Komponente und Reduktion, die Colorbar zusätzlich die Einheit.
`limits=(min,max)` ist endlich und streng aufsteigend. Ohne Grenzen werden die
Flächenfarben verwendet; konstante Felder erhalten eine symmetrische Erweiterung
von `max(1,abs(value))*1e-6`. Außerhalb expliziter Grenzen liegende Farben werden
abgeschnitten und durch Colorbar-Endmarkierungen sichtbar. Gleiche Grenzen und
Colormap erlauben Vergleiche. Rohwerte bleiben in `face_values` zugänglich.

JPEG verlangt `.jpg`/`.jpeg`, ein vorhandenes Zielverzeichnis und positive
endliche DPI. Standardmäßig erfolgt exklusive Dateierzeugung, vorhandene Dateien
lösen `FileExistsError` aus. `overwrite=True` erlaubt Ersatz. Schreibfehler
werden als OS-Fehler mit Pfad weitergegeben. Pixelgröße folgt Zollgröße × DPI;
Bytes sind zwischen Font-/Matplotlib-Versionen nicht reproduzierbar zugesagt.

## Nachweis und Beispiel

[Ausführbarer Model-Workflow](../examples/field_visualization/README.md) zeigt
Temperatur, Verschiebungsbetrag, zellbasierten Wärmestrom, Potential und Spannung
sowie zwei unabhängige Modellläufe mit gleicher Temperaturskala.
Tests prüfen sechs Außenflächen zweier benachbarter Tetraeder, getrennte
Zellfarben 2/8, lineares skalares Feld und Vektorbetrag `(3,4,0) → 5`, sämtliche
Komponenten, aktive Knoten, Fehlerfälle, konstante Farben, Figure-Lebensdauer
und JPEG-Format/-Abmessungen über Pillow. Zuordnungen sind exakt; analytische
Mittel/Beträge verwenden `atol=1e-14, rtol=0` für Größenordnung eins.
Die Bilder des Solverbeispiels wurden lokal auf Geometrie, Labels, Colorbar
und Farbzuordnung angesehen. Das belegt den Adapter, keine Netzkonvergenz
oder experimentelle Validierung. Numerische Operatoren werden nicht verändert;
neue Solver-Konvergenztests sind daher nicht erforderlich.
