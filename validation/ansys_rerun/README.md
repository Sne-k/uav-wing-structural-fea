# Re-solving the saved decks in Ansys MAPDL

`rerun_in_mapdl.py` works as follows:

1. It takes the solver decks that Workbench wrote
   (`ansys/wing_fea_files/dp0/SYS-2|SYS-3/MECH/ds.dat`).
2. It writes modified copies into a temporary folder.
3. It runs **Ansys MAPDL 2026 R1 (Student)** in batch mode.
4. It reads back the results that the /POST1 commands appended to each copy
   write out.

The student's mesh, contacts and materials are used unchanged. Only the
corrections listed below are inserted just before `solve`; the original files
are never modified.

```bash
python validation/ansys_rerun/rerun_in_mapdl.py   # ~5-10 min; set AWP_ROOT261 if Ansys is elsewhere
```

| Case | Change to the saved deck |
|---|---|
| R0 / M0 | none (static / modal deck exactly as saved) |
| R1 / M1 | every node in the root plane (y = 0) clamped: skin edge, spar root faces and root cap; modal with 12 modes |
| R2 | as R1, and the 800 Pa pressure replaced by the design ultimate load: lifting-line lift for a 5 kg UAV at 6 g (147.15 N per semi-span), centre of pressure at 25 % chord, as nodal forces on the upper-skin nodes |

## Static ([static_results.csv](static_results.csv))

| Case | Max total deformation | Tip deflection | Max von Mises (where) | Skin max | Ribs max | Spars max | Load from reactions |
|---|---|---|---|---|---|---|---|
| R0 as saved | **7.143 mm** | 7.124 mm | **154.91 MPa** (rear-spar root, y = 0) | 52.5 MPa (y = 0.009 m) | 31.0 MPa (rib 1) | 154.9 MPa (y = 0) | 434.4 N |
| R1 root clamped | **2.698 mm** | 2.689 mm | 17.7 MPa (inside rib 2, local) | **10.60 MPa** (root) | 17.7 MPa (rib 2) | 7.4 MPa | 434.5 N |
| R2 root clamped, design load | **0.818 mm** | 0.812 mm | 5.5 MPa (rib 2) | **3.45 MPa** (root) | 5.5 MPa (rib 2) | 2.2 MPa | 147.1 N |

R0 reproduces the Workbench results exactly. Clamping the root is the only
change between R0 and R1, and it removes the 154.9 MPa artefact and 62 % of the
deflection.

| | Ansys, root clamped (R1/R2) | wingfe shell model | Beam theory |
|---|---|---|---|
| Tip deflection, 800 Pa | 2.698 mm | 2.698 mm | 2.86 mm |
| Skin root stress, 800 Pa | 10.6 MPa | 10.7 MPa | 10.2 MPa |
| Tip deflection, design load | 0.818 mm | 0.814 mm | – |
| Skin root stress, design load | 3.45 MPa | 3.36 MPa | – |

## Modal ([modal_results.csv](modal_results.csv))

**M0, as saved.** The six frequencies reproduce exactly (13.892 … 87.749 Hz).
The located peaks show what the local modes are:
- **Mode 3 (37.0 Hz):** at the tip leading edge.
- **Mode 5 (68.1 Hz):** at the tip trailing edge.
- **Mode 4 (63.9 Hz):** on the front spar between ribs 2 and 3 (x = 0.065 m,
  y = 0.31 m).

**M1, root clamped, 12 modes.**

| Mode type | Ansys, root clamped | wingfe shell model |
|---|---|---|
| Flapwise bending 1 | 22.77 Hz | 22.50 Hz |
| Flapwise bending 2 | 123.23 Hz | 123.42 Hz |
| In-plane bending 1 | 141.29 Hz | 141.37 Hz |
| Torsion 1 | 221.35 Hz | 213.96 Hz |
| Local modes | 36.98, 64.05, 68.12, 157.7, 176.2, 211.1, 223.2, 266.2 Hz | none below 700 Hz |

The local modes in the Ansys model sit at the 0.1 mm tip cap, or on the front
spar between ribs. The spars are joined to the skin only by bonded contact, and
not to the ribs, so short spar segments can vibrate on their own. The wingfe
model has shared nodes and a 2 mm tip rib, so these modes do not appear in it.

## Workbench check

`ansys/wing_fea.wbpj` (renamed from `complete.wbpj` together with its `_files`
folder) opens in Workbench 2026 R1 Student. This was checked in batch mode on a
copy. The systems listed were: "dummy", "static and fatigue aluminium",
"Modal aluminum", "CFRP static", "CFRP modal".

## Corrected Workbench project

`build_corrected_workbench.py` applies the same corrections inside a copy of
the Workbench project and re-solves it; see
[`results/corrected/`](../../results/corrected/README.md).

It does this with two mesh-independent APDL Commands snippets in
[`snippets/`](snippets): the root clamp, and the design load on the upper skin.
Workbench gives, for aluminium:
* 0.814 mm and 5.43 MPa under the design load (MAPDL re-run above: 0.818 mm and
  5.45 MPa);
* 22.77, 123.23, 141.29 and 221.35 Hz for the main modes, identical to M1.
