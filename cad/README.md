# CAD (Autodesk Fusion 360)

| File | What it is |
|---|---|
| `wing_final.f3d` / `wing_final.step` | Geometry used by the final Ansys project. Skin as a surface body (`Wing_body`), root and tip closing faces (`Body44`, `Body45`), 10 solid ribs, 2 solid spars. |
| `naca2412.dat` | NACA 2412 coordinates in Selig format (generated from the NACA 4-digit equations) for [`tools/fusion_import_airfoil.py`](../tools/fusion_import_airfoil.py). |
| `archive/v1_solid_wing_2026-04-19.step` | First all-solid model (wing, spars, ribs), used in the failed April solve. |
| `archive/v2_skin_ribs_spar_2026-04-30.step` | Intermediate: skin surface, ribs, one spar. |
| `archive/v5_variant_unused_2026-05-09.*` | A variant with extra surface bodies and no spars. It was tried, then reverted. |

## Geometry

Values below were measured from the Ansys mesh, so they are what was actually
analysed.

| Item | Value |
|---|---|
| Airfoil | NACA 2412 (12 % thick at 30 % chord, 2 % camber at 40 % chord), finite TE 0.26 % c |
| Semi-span | 1.200 m |
| Root / tip chord | 0.250 m / 0.200 m (taper 0.8), leading edge straight at x = 0 |
| Sweep, twist, dihedral | 0 |
| Front spar | straight, x = 62.5–67.5 mm, 20 mm deep (25 % of root chord, 31 % at the tip) |
| Rear spar | on the 70 % chord line (tapered), 5 mm wide, 17 mm deep at the root |
| Ribs | 10, about 2 mm thick, at y = k × 1.2/11 (109 mm pitch) |

Axes: x chordwise (aft), y spanwise (root → tip), z up.
