# Verifikation und Validierung

Status: Entwurf

## 1. Terminologie

- **Code-Verifikation:** Lösen wir die implementierten Gleichungen richtig?
- **Lösungsverifikation:** Ist die diskrete Lösung hinreichend konvergiert?
- **Validierung:** Beschreiben die Gleichungen die beobachtete Realität für den
  vorgesehenen Anwendungsfall ausreichend?

Diese Stufen werden nicht vermischt. Übereinstimmung mit einem anderen Solver
allein ist keine vollständige Code-Verifikation.

## 2. Testpyramide

### 2.1 Algebraische Unit-Tests

- Tet4-Volumen und Formfunktionsgradienten
- `B`- und `G`-Operatoren
- isotrope Steifigkeitsmatrix aus `E, nu`
- Tensorrotationen für `C`, `e`, `kappa`, `k`, `alpha`
- lokale Matrixsymmetrie und Kopplungsreziprozität
- Flächengeometrie und nach außen orientierte Normalen
- DOF-Zuordnung und symmetrieerhaltende Dirichlet-Elimination

### 2.2 Einzelelement- und Patch-Tests

| Physik | Test | Erwartung |
|---|---|---|
| Mechanik | Starrkörpertranslation/-rotation | Nullspannung, Nullenergie |
| Mechanik | affine Verschiebung | konstante exakte Dehnung/Spannung |
| Wärme | lineares Temperaturfeld | konstanter exakter Wärmestrom |
| Elektrostatik | lineares Potential | konstantes exaktes E- und D-Feld |
| Piezo | affine `u` und `phi` | konstante gekoppelte Felder, exaktes Residuum |

Patch-Tests werden auch auf verzerrten, unterschiedlich orientierten
Tetraedern und nach Knotenumordnung durchgeführt.

### 2.3 Analytische Benchmarks

- Zugstab und dreidimensionaler uniaxialer Spannungszustand
- Patch unter hydrostatischer Last
- stationäre Wärmeleitung im Quader
- Plattenkondensator mit analytischer Kapazität
- frei expandierender und vollständig behinderter Thermokörper
- piezoelektrischer Block im `d_33`-Modus:
  - prescribed voltage -> Dehnung/Verschiebung
  - prescribed pressure/force -> Ladung/Spannung

Für den piezoelektrischen Benchmark werden Kurzschluss (`phi=0` an beiden
Elektroden) und Leerlauf (eine Referenzelektrode, andere Elektrode mit
vorgegebener Gesamtladung null) unterschieden. Punktweise `D dot n = 0` ist kein
Ersatz für eine leitende Equipotential-Elektrode.

Der implementierte [d₃₃-Nachweis](18_PIEZO_D33_VERIFICATION.md) verifiziert
direkte und inverse Wirkung sowie einen 1D-invarianten Leerlaufgrenzfall.
Die leitende Floating Electrode ist im M5-Solver implementiert. Der
[Vier-Netze-Nachweis des vollständigen stationären Modells](25_THERMOPIEZO_CONVERGENCE.md)
prüft sie zusätzlich bei Materialgrenze und thermischer Anregung.

### 2.4 Manufactured Solutions und Konvergenz

Glatte vorgegebene Felder `u*(x)`, `T*(x)`, `phi*(x)` werden in die starken
Gleichungen eingesetzt; daraus entstehen konsistente Quellen und
Randbedingungen. Mindestens vier systematisch verfeinerte Netze werden geprüft.

Erwartet werden für hinreichend glatte Lösungen und reguläre Netze ungefähr:

- `L2`-Fehler der linearen Primärfelder: Ordnung 2,
- Energie-/`H1`-Fehler: Ordnung 1.

Die beobachtete Ordnung wird aus drei aufeinanderfolgenden Netzen berechnet.
Asymptotische Ordnung, nicht nur kleiner absoluter Fehler, ist das Kriterium.

### 2.5 Globale Bilanzen

Jeder gelöste Fall berichtet:

- Summe äußerer Kräfte und Lagerreaktionen,
- zugeführte und abgeführte Wärmeleistung,
- freie und induzierte Gesamtladung,
- algebraisches Residuum je Feld,
- mechanische, elektrische und gekoppelte Energiebeiträge.

Eine Lösung mit kleinem linearem Residuum, aber verletzter physikalischer Bilanz
gilt als fehlgeschlagen.

## 3. Unabhängige Referenzen

Mindestens ein gekoppelter Fall wird gegen zwei unabhängige Referenzen geprüft,
beispielsweise eine analytische 1D-Reduktion und einen etablierten FEM-Code.
Verglichen werden nicht nur Maximalwerte, sondern Feldnormen, Reaktionen,
Gesamtladung und Energien. Materialkonventionen und elektrische Randbedingungen
müssen identisch transformiert sein.

## 4. Regression

Regressionstests speichern kleine, menschenprüfbare Referenzgrößen:

- wenige ausgewählte DOFs,
- Bilanzen und Energien,
- Matrixdimension, Nichtnullstruktur und Symmetriefehler,
- Solverstatus und normierte Residuen.

Große binäre Golden Files werden vermieden. Eine fachlich begründete Änderung
der Referenzwerte benötigt einen dokumentierten Entscheid.

## 5. Validierung mit Experimenten

Erst nach Code- und Lösungsverifikation folgt ein einfacher Piezoexperiment-
Datensatz, etwa ein PZT-Block oder -Aktuator mit bekannter Geometrie,
Polarisationsrichtung, Elektrodenkonfiguration und Temperatur. Zu erfassen sind
Unsicherheiten in Materialdaten, Klemmung, Elektrodengeometrie, Messkette und
Geometrie. Das Validierungsergebnis ist eine Abweichung mit Unsicherheitsband,
kein binäres „richtig/falsch“.

## 6. Freigabegates

Ein Meilenstein wird nur freigegeben, wenn:

1. Unit- und Patch-Tests grün sind;
2. Bilanzfehler unter dokumentierter Toleranz liegen;
3. mindestens ein Benchmark unabhängig bestätigt ist;
4. keine unerklärte Verschlechterung der Konvergenzordnung vorliegt;
5. die Beispielrechnung aus frischer Umgebung reproduzierbar ist.

## M12: getrennte Code-, Modell- und Kopplungsnachweise

[Der RANS-/FSI-Vertrag #84](71_M12_RANS_FSI_CONTRACT.md) setzt vor Umsetzung
separate Gates für die neue FV-Druckkopplung, variable Viskosität, SA-Terme,
Wandabstände, nicht passende Schnittstellen und stationäre Netzbewegung.
Der unabhängige 1D-SA-Codevergleich prüft Modellimplementation; der DNS-Kanal
prüft den Modellfehler und ist keine experimentelle Validierung. Laminarer
und SA-Deckeltest sind beide für RANS-FSI erforderlich. Neue Gates werden
erst nach tatsächlich dokumentierter Prüfung als bestanden bezeichnet.
