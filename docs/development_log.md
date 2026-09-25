# Development log

How the project evolved, reconstructed from the working notes and the files in
this repository. Dates are 2026.

## March: literature and planning
* Literature review of five papers ([literature_review.md](literature_review.md)),
  research gaps and objectives, and a first presentation draft.
* Wing parameters chosen: **NACA 2412**, semi-span **1.2 m**, root chord
  **0.25 m**, tip chord **0.20 m** (taper 0.8), no sweep at the leading edge, no
  twist.
* Airfoil coordinates (Selig format) imported into Fusion 360 with a small Python
  script ([`tools/fusion_import_airfoil.py`](../tools/fusion_import_airfoil.py)).
  Fixes needed along the way: a UTF-8 byte-order mark in the `.dat` file, and
  `ObjectCollection` instead of a Python list for the spline points.
* Loft between root and tip sections, after closing the trailing edge so the
  profiles are closed regions.

## March–April: internal structure in Fusion 360
* Front spar at 25 % of the root chord (x = 62.5 mm), 20 × 5 mm. Rear spar at
  70 % chord, 17 × 5 mm.
* Spars extruded "To Object" along the span. "Join" failed until the spar
  profiles were fully inside the section.
* The rear spar stuck out of the thinner tip section, so it was **tapered** to
  follow the 70 % chord line.
* Ten ribs at 109 mm pitch (span / 11). Rib profiles were made with Intersect
  because the wing tapers. They were first joined into one solid, then kept as
  separate bodies.

## 19–20 April: first Ansys attempt (solid model), "failed attempt"
* The all-solid STEP (`cad/archive/v1_solid_wing_2026-04-19.step`) imported into
  Workbench: 41k nodes, 35 contact regions (105k contact elements), 80 Pa.
* The solver stopped with **"small pivot"** errors at rib 5 (UX): parts were not
  properly connected, and thin solids had one element through the thickness.
* A merged single-solid variant solved: 0.55 mm and 1.49 MPa at 100 Pa. The
  decision was then to rebuild as a **shell model**.

## 30 April – 9 May: shell model
* The skin was converted to a surface in Fusion (`new geo.step`, then
  `wing_final.step`). The root and tip closing faces became separate surface
  bodies (Body44, Body45).
* Meshed with SHELL181 (skin and caps) and SOLID186 (ribs and spars), plus 16
  bonded contacts. Mesh warnings were accepted.
* XFLR5: wing built, VLM2 analysis at 20 m/s for α = 4–8°. The spanwise load
  could not be exported from the UI version in use. Coordinate-based pressure
  functions in Mechanical also failed ("Y is undefined"), so a uniform
  **800 Pa** was used instead.
* Pivot error on Body45 (a free ROTX DOF): weak springs were switched on.
* Iterations on shell thickness gave tip deflections of 0.012 → 0.055 → 0.21 mm,
  then an unstable run, then a 11.4 m "rigid-body" result, and finally the saved
  configuration: **8 mm skin, 0.1 mm end caps → 7.14 mm, 154.9 MPa**.
* Modal analysis (6 modes): 13.9 to 87.7 Hz.

## 10 May: CFRP and fatigue
* The static and modal systems were duplicated with "CFRP": E = 70 GPa,
  ν = 0.3, ρ = 1600 kg/m³ (from paper 4). Results: 7.137 mm, 157.9 MPa,
  frequencies × 1.30, mass −41 %.
* Fatigue on the aluminium static system: S-N curve (300 MPa at 10³ cycles down
  to 90 MPa at 10⁷), UTS 310 MPa, zero-based loading, Goodman. Minimum life
  **3.33 × 10⁶ cycles**.
* Conclusions and a flutter discussion were written, and the project was added
  to the CV.

## September: repository, audit and re-analysis
* Everything organised into this repository.
* The solver decks and result files were audited
  ([model_audit.md](model_audit.md)). The root boundary condition and the load
  definition turned out to drive the reported stress, deflection and
  frequencies.
* Independent re-analysis with a verified shell FE code, a clamped root, a
  derived design load, sizing with buckling, a composite laminate study and
  fatigue ([validation/README.md](../validation/README.md)).
