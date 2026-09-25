# Model audit: verification and known issues

This page checks the saved Ansys model against:

* what the solver actually received (`ds.dat`) and produced (`solve.out`, `file.rst`);
* hand calculations;
* an independent, benchmark-verified shell FE model of the same wing with a
  correctly clamped root ([`validation/`](../validation/README.md)).

Every number can be reproduced:

```bash
python tools/inspect_ds_dat.py ansys/wing_fea_files/dp0/SYS-2/MECH/ds.dat   # what the solver received
python validation/read_ansys_results.py                                      # Ansys results (pip install ansys-mapdl-reader)
python tools/beam_check.py                                                   # beam theory
python validation/run_validation.py                                          # independent re-analysis
```

## 1. What the solver received

| Item | Value in the solver deck |
|---|---|
| Skin (`Wing_body`) | SHELL181, **8 mm**, mid-surface on the outer mould line, 12 870 elements, all triangles |
| Root cap (`Body44`) and tip cap (`Body45`) | SHELL181, 0.1 mm, 62 and 37 elements |
| Ribs 1 to 10 | SOLID186, about 2 mm thick, 41 to 59 elements each, at y = 0.109 to 1.091 m |
| Spars | SOLID186; spar 1 straight at x = 62–67 mm (25 % of root chord); spar 2 tapered, on the 70 % chord line |
| Connections | 16 bonded contacts (listed in 3.1); no shared topology, so every body is meshed separately |
| Fixed Support | **87 nodes, all on `Body44` (the root cap)** |
| Load | 800 Pa by components along +Z on **every skin facet (upper and lower surfaces)**; total reaction **434.6 N** |
| Solver options | linear static, weak springs on; modal: 6 modes |
| Mass | 12.419 kg (aluminium), 7.359 kg ("CFRP"); the 8 mm skin alone is about 11.95 kg |

## 2. Cross-check: three independent answers

| Quantity | Ansys (as saved) | Beam theory, clamped root | Shell FE re-analysis, clamped root (converged) |
|---|---|---|---|
| Tip deflection | **7.14 mm** | 2.86 mm | **2.70 mm** |
| Peak von Mises | **154.9 MPa** | 10.2 MPa (root bending) | **10.7 MPa** at the web/skin root corner; 8.3 MPa elsewhere |
| 1st flapwise bending | 13.89 Hz | 22.5 Hz (Rayleigh) | **22.50 Hz** |
| 1st in-plane (chordwise) bending | 23.86 Hz | – | **141.4 Hz** |
| 1st torsion | not in the first 6 modes | – | **214.0 Hz** |
| Total load / reaction | 434.6 N | 441 N | 440.9 N |
| Mass | 12.42 kg | 11.95 kg (skin only) | 12.85 kg (webs are full depth) |

Beam theory and the independent shell model agree with each other to within
5–6 % on deflection and stress, and to within 0.1 % on frequency. The Ansys
model is 2.6 times more flexible, its peak stress is 15 times higher, and its
in-plane bending frequency is 6 times lower. The reasons follow.

## 3. Known issues, ordered by impact

### 3.1 The wing is supported only through the spar end faces (critical)

The contact map from `ds.dat`:

| Contact regions | Joins |
|---|---|
| 1–10 | skin to rib 1 … rib 10 |
| 11, 12 | skin to spar 1, skin to spar 2 |
| 13, 14 | **root cap (`Body44`) to spar 1 (13 nodes), and to spar 2 (8 nodes)** |
| 15, 16 | tip cap (`Body45`) to spar 1, and to spar 2 |

The Fixed Support is on the root cap, and **no contact joins the skin to the
root cap**. The skin's root edge is therefore free. The entire root bending
moment (about 255 N·m) goes from the skin into the spars, and then through two
small contact patches (8 and 13 nodes) on the spar end faces into the support.

The result file confirms it (`validation/results/ansys_rst_readback.txt`). The
maximum von Mises, **154.9 MPa, is at the root end of spar 2** (y = 0), on the
8-node contact patch. Along that spar it falls to 91 MPa within 20 mm, 22 MPa at
50–100 mm, and 6 MPa beyond 200 mm. It is a local artefact of the connection,
not a wing stress.

Consequences:

* The reported "maximum stress near the root" and the factor of safety (1.61)
  describe the connection patch.
* The whole wing rotates about that soft root joint, which gives the 7.14 mm tip
  deflection instead of about 2.7 mm. It also gives the low 13.9 Hz bending and
  23.9 Hz in-plane frequencies instead of 22.5 Hz and 141 Hz.
* The fatigue minimum life (3.33 × 10⁶ cycles) sits on the same patch. It checks
  out exactly: zero-based loading, σa = σm = 77.5 MPa, Goodman with Su = 310 MPa
  gives 103 MPa, and on the S-N curve that is 3.3 × 10⁶ cycles. So the fatigue
  result describes the patch, not the wing.
* The "thickness optimisation" iterations were steered toward target
  deflection and stress ranges rather than design requirements. Their numbers
  can't be used for sizing.

**Fix:** clamp the complete root section (skin edge and both spar root faces),
or share topology between the root cap and skin. Then switch off weak springs
and confirm the solution still converges. The re-analysis does exactly this.

### 3.2 The load is doubled and not physically derived (high)

* The 800 Pa is defined **by components (+Z)** on the **whole** skin, so the
  lower surface is pushed upward too. The total is 434.6 N per semi-span, equal to
  1600 Pa over the 0.27 m² planform: twice the intended load.
* At the stated 20 m/s cruise (q = 245 Pa), that would need C_L ≈ 6.6. For the
  stated 5 kg UAV, one semi-span carries about 25 N in level flight, so the model
  applies about 17 g.
* XFLR5 (VLM2, α = 4–8°) was run, but the spanwise lift was never extracted. The
  800 Pa was a placeholder after the coordinate-based pressure expression failed.

**Fix:** a stated design case (mass, load factor, safety factor). Take the
spanwise lift from lifting-line theory, XFLR5 or a Schrenk distribution, and
apply it to the upper skin. The re-analysis uses a 5 kg MTOW, n_limit = 4 and a
1.5 safety factor, with lifting-line lift (Part B).

### 3.3 The "CFRP" case only changes the density (high)

Both materials use E = 70 GPa. Only ν (0.30 vs 0.33) and ρ (1600 vs 2700) differ.
So the static results are identical (7.137 vs 7.143 mm, 157.9 vs 154.9 MPa), and
every frequency rises by exactly √(2700/1600) = 1.299, a pure mass effect. The
re-analysis reproduces this: 1.298.

A carbon laminate has no yield point and fails ply by ply, so von Mises stress,
a 500 MPa "yield" and a factor of safety are not meaningful for it. These
results support **"same stiffness, 41 % lighter"**. They do **not** support
"stiffer", "less deformation" or "better flutter resistance".

**Fix:** orthotropic plies and a real layup, evaluated with Tsai-Wu or Hashin.
The re-analysis uses T300/5208 plies with classical lamination theory (Part D).

### 3.4 Modal results were mis-labelled (medium)

Participation factors from `solve.out` (aluminium; CFRP scales by 1.299):

| Mode | f [Hz] | Dominant effective mass | What it actually is |
|---|---|---|---|
| 1 | 13.89 | Z 67 %, ROTX 99 % | 1st **flapwise** bending |
| 2 | 23.86 | X 72 %, ROTZ 94 % | 1st **in-plane** bending (6 times too low, see 3.1) |
| 3 | 36.98 | ≈ 0 in all directions | local mode (normalised amplitude 149) |
| 4 | 63.93 | ≈ 0 | local mode |
| 5 | 68.12 | ≈ 0 | local mode at the tip (normalised amplitude 341) |
| 6 | 87.75 | Z 18 % | 2nd **flapwise** bending (87.75 / 13.89 = 6.3; textbook cantilever 6.27) |

No torsion mode appears in the first six, so the modal results say nothing about
bending–torsion coupling or flutter. With the root clamped, the first eight
modes are: 22.5 (flap 1), 123.4 (flap 2), 141.4 (in-plane 1), 214.0
(**torsion 1**), 317.5, 567.8, 574.9 and 700.7 Hz.

### 3.5 Sizing is far from a UAV wing (medium)

Under the design ultimate load (6 g on a 5 kg aircraft), the saved design reaches
only 3.4 MPa (reserve factor 82 on yield) and deflects 0.8 mm. A semi-span weighs
12.8 kg, so the two wings weigh about five times the whole aircraft. Part C of
the re-analysis sizes the skin with strength and buckling checks.

### 3.6 Mesh and model quality (medium)

* The skin mesh is 100 % degenerate triangular SHELL181. Ansys recommends
  triangles only as filler because they are overly stiff; use a quad-dominant
  mesh.
* There is no mesh-convergence study. (The re-analysis includes one: tip
  deflection changes 0.07 % between the last two meshes.)
* An 8 mm skin with mid-surface offset sits 4 mm outside the aerodynamic surface
  and overlaps the solid spars and ribs by 4 mm.
* Ribs are bonded only to the skin, not to the spars.
* Weak springs hide constraint problems. The April solid model failed with a
  pivot error at rib 5 for the same kind of reason.
* Aluminium has no yield strength defined in Engineering Data, so Mechanical
  cannot compute a safety factor. The FoS of 1.61 was computed by hand from the
  artefact stress.

### 3.7 Documentation consistency (low)

* The airfoil is **NACA 2412**. The root mesh nodes match it: 12.0 % thick at
  30 % chord, 2 % camber at 40 % chord. Some earlier write-ups say
  "NACA 63(1)-412".
* The elliptical-load expression tried in Mechanical assumed a 0.35 m semi-span
  instead of 1.2 m. The XFLR5 file isn't in the project folder, so it can't be
  checked that the aerodynamic and structural models match.
* The "5 kg UAV" mass is used but never justified.

## 4. What is solid

* The CAD and FE pipeline works end to end: airfoil import script, tapered loft,
  ribs, a tapered rear spar, STEP export, a shell-plus-solid idealisation, and
  linked static, modal and fatigue systems.
* Within the saved model, the numbers are internally consistent. The reaction
  equals the applied load; the CFRP frequency ratio is exactly 1.299; the fatigue
  life follows exactly from Goodman and the S-N curve; and mode 6 / mode 1 = 6.3,
  like a textbook cantilever.
* The fatigue method is sound: zero-based loading, Goodman and an S-N curve with
  a UTS entered. It just needs the correct stress field.
