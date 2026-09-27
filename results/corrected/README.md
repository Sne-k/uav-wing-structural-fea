# Corrected Ansys results

Solved in Ansys Workbench 2026 R1 Student from `ansys/corrected/wing_fea_corrected.wbpj`,
a copy of the original project built by
[`validation/ansys_rerun/build_corrected_workbench.py`](../../validation/ansys_rerun/build_corrected_workbench.py).

Changes in every system:
* the 800 Pa pressure is suppressed;
* two APDL **Commands** objects are added: `root_clamp.inp` clamps the whole root
  section, and `design_load.inp` applies the design ultimate load (static systems
  only; 147.15 N per semi-span on the upper skin). Both are in
  [`validation/ansys_rerun/snippets/`](../../validation/ansys_rerun/snippets);
* weak springs are off;
* the modal systems find 12 modes.

| System | Result | Value |
|---|---|---|
| Aluminium static | max total deformation | **0.814 mm** |
| | max equivalent stress | **5.43 MPa** (local peak in rib 2; skin ≈ 3.4 MPa at the root) |
| | fatigue life (minimum) | 10⁷ cycles everywhere (end of the S-N curve: no fatigue damage) |
| "CFRP" static (isotropic, as in the project) | deformation / stress | 0.816 mm / 5.45 MPa |
| Aluminium modal | modes 1–12 | 22.77 flap 1 · 36.98 local · 64.05 local · 68.12 local · **123.23 flap 2** · **141.29 in-plane 1** · 157.74 local · 176.19 local · 211.11 local · **221.35 torsion 1** · 223.16 local · 266.20 local Hz |
| "CFRP" modal | modes 1–12 | every frequency is 1.2976 × the aluminium value (a density-only change) |

The local modes are vibrations of the 0.1 mm tip cap and of front-spar segments
between ribs. The spars are joined to the skin only by contact, not to the ribs.

Files:
* `*_total_deformation.png`, `*_equivalent_stress.png`, `al_static_fatigue_life.png`
* `*_modeNN_<f>Hz.png`
* `values_*.txt`: the values read from Mechanical

The modal "maximum deformation" values are mass-normalised and have no physical
size.
