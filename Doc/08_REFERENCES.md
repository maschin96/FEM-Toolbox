# Quellen und Referenzen

Status: Entwurf

## Fachbücher

1. O. C. Zienkiewicz, R. L. Taylor, J. Z. Zhu: *The Finite Element Method: Its
   Basis and Fundamentals*. Elsevier.
2. T. J. R. Hughes: *The Finite Element Method: Linear Static and Dynamic Finite
   Element Analysis*. Dover.
3. J. N. Reddy: *An Introduction to the Finite Element Method*. McGraw-Hill.
4. A. Benjeddou: *Advances in Piezoelectric Finite Element Modeling of Adaptive
   Structural Elements: A Survey*. Computers & Structures 76 (2000).
5. IEEE Standard on Piezoelectricity, ANSI/IEEE Std 176-1987. Relevant für
   Materialkonventionen; Zugriff und Lizenz müssen projektseitig geklärt werden.

## Verifikation

6. K. Salari, P. Knupp: *Code Verification by the Method of Manufactured
   Solutions*, SAND2000-1444, Sandia National Laboratories.
   https://www.osti.gov/servlets/purl/759450-wLI4Ux/

## Primäre Softwaredokumentation

7. Gmsh Reference Manual — Elementtypen, Knotenanordnung und Python-API:
   https://gmsh.info/doc/texinfo/
8. SciPy `sparse.linalg` — direkte und iterative Sparse-Solver:
   https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html
9. PETSc KSP/PC — Block- und Schur-Komplement-Verfahren:
   https://petsc.org/main/manual/ksp/
10. meshio — Netzformat-Ein-/Ausgabe:
    https://github.com/nschloe/meshio

## Vergleichbare offene Implementierungen

11. SfePy: Beispiel zur linearen Piezoelastizität:
    https://sfepy.org/doc/examples/multi_physics-piezo_elasticity.html
12. SfePy: Beispiel zur Piezoelastodynamik:
    https://sfepy.org/doc/examples/multi_physics-piezo_elastodynamic.html

## Theorie-Quervergleich

13. COMSOL Structural Mechanics Theory, Piezoelectricity — Stress-Charge- und
    Strain-Charge-Form sowie Materialkonventionen:
    https://doc.comsol.com/6.3/doc/com.comsol.help.sme/sme_ug_theory.06.043.html
14. COMSOL Thermoelasticity Theory — thermodynamisch reziproke Kopplung:
    https://doc.comsol.com/6.4/doc/com.comsol.help.mems/mems_ug_sme.08.7.html
15. Abaqus Theory Guide, triangular/tetrahedral/wedge elements — Eigenschaften
    linearer Tetraeder:
    https://docs.software.vt.edu/abaqusv2025/English/SIMACAETHERefMap/simathe-c-tritetwedge.htm

## Nutzung der Quellen

Externe Dokumentation dient zur Gegenprüfung, nicht als Spezifikation des
eigenen Codes. Vor jeder numerischen Referenzrechnung werden Voigt-Reihenfolge,
Schubkonvention, elektrische Vorzeichen, Materialachsen, Permittivitätsdefinition
und Elektrodenrandbedingungen explizit abgeglichen.
