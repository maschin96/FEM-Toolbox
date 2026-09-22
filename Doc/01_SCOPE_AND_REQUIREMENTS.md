# Umfang und Anforderungen

Status: Entwurf

## 1. Ziel

Entwickelt wird eine nachvollziehbare Python-FEM-Toolbox für lineare gekoppelte
3D-Feldprobleme. Der wissenschaftliche Kern ist die konsistente gemeinsame
Behandlung von Verschiebung, Temperatur und elektrischem Potential auf einem
Tet4-Netz. Direkter und inverser piezoelektrischer Effekt müssen aus derselben
konstitutiven Beschreibung folgen.

Die erste Version ist ein Referenzcode für Forschung, Lehre und die Prüfung
neuer Formulierungen. Sie ist zunächst kein allgemeiner industrieller
Produktionssolver.

## 2. Modellannahmen der ersten Version

- dreidimensionaler euklidischer Raum
- infinitesimale Verzerrungen und kleine Verschiebungen
- lineare Materialien, stückweise konstant je Volumenelement
- lineare Formfunktionen auf geraden Vierknoten-Tetraedern
- quasistatische Mechanik und Elektrostatik
- stationäre Wärmeleitung im ersten Ausbau
- SI-Einheiten intern
- perfekte Grenzflächenbindung mit stetigen Primärvariablen
- keine freien Volumenladungen im Piezo-Basistest, sofern nicht explizit gesetzt

## 3. Funktionale Anforderungen

### 3.1 Netz und Gebiete

- Knotenkoordinaten und Tet4-Konnektivität speichern
- orientierungsunabhängige Volumenberechnung; negative Orientierung erkennen und
  kontrolliert korrigieren
- degenerierte Elemente mit konfigurierbarer Toleranz ablehnen
- Volumen- und Oberflächengruppen über stabile ganzzahlige Tags adressieren
- Gmsh-Netze importieren und VTU-Ergebnisse exportieren

### 3.2 Materialien

- isotrope und vollständig anisotrope Elastizität `C^E`
- anisotrope Wärmeleitfähigkeit `k`
- anisotrope Permittivität `kappa^S`
- piezoelektrische Kopplungsmatrix `e`
- thermischer Ausdehnungsvektor `alpha`
- Dichte `rho` und Wärmekapazität `c_p` für spätere Transienten
- Materialachsen durch explizite orthonormale Orientierung
- Validierung von Form, Symmetrie, Endlichkeit und physikalisch zulässigen
  Eigenwerten

### 3.3 Randbedingungen und Lasten

- komponentenweise Dirichlet-Bedingungen für `u`, `T`, `phi`
- Volumenkraft, Oberflächenlast und Punktlast
- Wärmequelle, Wärmestrom und Konvektion
- freie Ladungsdichte und Oberflächenladung
- geerdete und vorgegebene Elektroden
- Equipotential-Elektrode mit entweder vorgegebenem Potential oder vorgegebener
  Gesamtladung als spätere Erweiterung

### 3.4 Analysen

1. lineare statische Elastizität
2. stationäre Wärmeleitung
3. Elektrostatik
4. stationäre Thermoelastizität als sequenzielle Kopplung `T -> u`
5. monolithische quasistatische Piezoelektrizität `u <-> phi`
6. kombiniertes stationäres Modell `T -> (u, phi)`

Zeitabhängige Wärmeleitung, elastische Dynamik, harmonische Piezoanalyse und
reziproke Thermoelastizität folgen erst nach diesen sechs Analysen.

## 4. Qualitätsanforderungen

- reproduzierbare Ergebnisse für identische Eingaben
- klare Fehlermeldungen mit Gebiet, Element oder Freiheitsgrad
- keine stillen Einheitenkonversionen
- getrennte Element-, Assemblierungs-, Randbedingungs- und Solvertests
- Residuum, Randreaktionen und globale Bilanzen als Standardausgabe
- vollständige Provenienz: Eingabe, Toolbox-Version, Solverparameter und
  Konvergenzinformationen im Ergebnis
- Vektorisierung erst nach verifizierter skalarer Referenzimplementierung

## 5. Nichtziele der ersten Version

- geometrische oder materielle Nichtlinearität
- Kontakt, Bruch, Plastizität und Hysterese
- elektromagnetische Wellenausbreitung
- Fluid-Struktur-Kopplung
- adaptive Netzverfeinerung
- automatische CAD-Vernetzung im Kernpaket
- GPU- oder verteiltes Rechnen
- Tet10-, Hex8- oder gemischte Elemente

Diese Punkte dürfen die Architektur beeinflussen, aber den ersten verifizierten
Solver nicht verzögern.

## 6. Abnahmekriterien für den stationären Ausbau bis v1.0

Der stationäre Zielumfang bis v1.0 ist erreicht, wenn:

1. alle drei Einzelfeldsolver Patch- und Bilanztests bestehen;
2. der piezoelektrische Blockoperator bis auf die dokumentierte
   Sattelpunkt-Signatur symmetrisch ist;
3. direkter und inverser Piezoeffekt an mindestens einem Benchmark quantitativ
   gegen eine unabhängige Referenz übereinstimmen;
4. eine Netzverfeinerungsstudie die erwartete asymptotische Konvergenz zeigt;
5. ein vollständiges Beispiel von Gmsh-Datei bis VTU-Ergebnis reproduzierbar
   läuft;
6. alle Annahmen, Einheiten und Materialkonventionen im Ergebnisbericht stehen.

Die ursprüngliche Bezeichnung „Version 0.1“ meinte diesen fachlichen Zielumfang.
Im eingeführten Release-Zyklus ist v0.1.0 dagegen M0 (Projektfundament).
Die Kriterien werden schrittweise in M2–M5 erfüllt. Langfristige Erweiterungen
und ihre Nichtziele beschreibt die [Produktvision](13_PRODUCT_VISION.md).
