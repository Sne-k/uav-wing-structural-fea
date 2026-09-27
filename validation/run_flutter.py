"""Preliminary flutter and divergence estimate (typical section at 75 % semi-span).

The section properties come from the shell FE model:
  * mass per span, CG and polar inertia: lumped nodal masses in the rib bay around y = 0.9 m
  * elastic (flexural) axis: chord position where a tip load produces no twist
  * omega_h, omega_alpha: first flapwise-bending and first torsion modes of the FE model
Aerodynamics: Theodorsen, V-g method (see wingfe/flutter.py).

Run:  python validation/run_flutter.py   (after run_validation.py)
Writes validation/results/F_flutter.csv
"""
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import run_validation as rv  # noqa: E402
from wingfe.aero import RHO_SL  # noqa: E402
from wingfe.flutter import divergence_speed, flutter_speed  # noqa: E402
from wingfe.materials import T300_5208, laminate_section  # noqa: E402
from wingfe.wing import SPAN, Mesh, chord  # noqa: E402

Y_TS = 0.75 * SPAN
V_DIVE_ASSUMED = 40.0   # m/s, assumed design dive speed of a ~20 m/s cruise UAV


def section_properties(res):
    m, info = res["model"], res["info"]
    pitch = SPAN / 11
    lo, hi = Y_TS - pitch / 2, Y_TS + pitch / 2
    # Element masses and centroids in the rib bay around Y_TS. The bay edges fall on
    # node rows, so selecting elements by centroid (not nodes by position) avoids a
    # boundary row being counted twice or not at all.
    me, ce = [], []
    for (conn, grp), (dofs, T, data) in zip(m.elements, m.edata):
        c = m.nodes[conn].mean(axis=0)
        if lo <= c[1] < hi:
            me.append(m.sections[grp].rho_t * data["area"])
            ce.append(c)
    me, ce = np.array(me), np.array(ce)
    mass_span = me.sum() / pitch
    x, z = ce[:, 0], ce[:, 2]
    xcg = (me * x).sum() / me.sum()
    zcg = (me * z).sum() / me.sum()

    # flexural axis from two unit tip loads (at the leading and trailing edge nodes of the tip)
    tip = info["tip"]
    le, te = tip[np.argmin(m.nodes[tip, 0])], tip[np.argmax(m.nodes[tip, 0])]
    j = int(np.argmin(np.abs(info["y_stations"] - Y_TS)))
    s_le, s_te = info["grid_u"][j][0], info["grid_u"][j][-1]
    twist = []
    for n in (le, te):
        F = np.zeros(m.ndof)
        F[6 * n + 2] = 1.0
        uu = m.solve_static(F)
        twist.append((uu[6 * s_te + 2] - uu[6 * s_le + 2]) / (m.nodes[s_te, 0] - m.nodes[s_le, 0]))
    x_le, x_te = m.nodes[le, 0], m.nodes[te, 0]
    x_fa_tip = x_le + (x_te - x_le) * twist[0] / (twist[0] - twist[1])
    x_ea = x_fa_tip / chord(SPAN) * chord(Y_TS)          # same chord fraction at the typical section
    I_ea = (me * ((x - x_ea)**2 + (z - zcg)**2)).sum() / pitch
    c = chord(Y_TS)
    b = c / 2
    return dict(b=b, mass_per_span=mass_span, x_cg_frac=xcg / c, x_ea_frac=x_ea / c,
                a=(x_ea - b) / b, x_a=(xcg - x_ea) / b, r2=I_ea / (mass_span * b * b),
                mu=mass_span / (math.pi * RHO_SL * b * b))


def run_design(name, sections):
    L = rv.DESIGN["semi_span_lift_ult_N"]
    r = rv.solve_case(sections, lambda m: rv.load_upper_skin(m, rv.design_pressure(L)), Mesh(), modes=8)
    f, lab = r["freqs"], r["mode_labels"]
    fh = next((x for x, l in zip(f, lab) if l.startswith("flapwise")), None)
    fa = next((x for x, l in zip(f, lab) if l.startswith("torsion")), None)
    if fh is None or fa is None:
        missing = "flapwise bending" if fh is None else "torsion"
        print(f"{name}: no {missing} mode among the first {len(f)} modes; skipped")
        return [name] + ["n/a"] * 13
    p = section_properties(r)
    wa = 2 * math.pi * fa
    sigma = fh / fa
    fl = flutter_speed(p["mu"], sigma, p["r2"], p["a"], p["x_a"])
    dv = divergence_speed(p["mu"], p["r2"], p["a"])
    VF = fl[0] * p["b"] * wa if fl else float("inf")
    VD = dv * p["b"] * wa if dv else float("inf")
    return [name, f"{fh:.1f}", f"{fa:.1f}", f"{p['mass_per_span']:.3f}", f"{p['mu']:.1f}",
            f"{p['x_ea_frac']:.3f}", f"{p['x_cg_frac']:.3f}", f"{p['a']:.3f}", f"{p['x_a']:.3f}",
            f"{p['r2']:.3f}", f"{sigma:.3f}", f"{VF:.0f}", f"{VD:.0f}", f"{VF / V_DIVE_ASSUMED:.1f}"]


def main():
    summ = json.load(open(os.path.join(rv.OUT, "summary.json")))
    ts = summ["chosen_aluminium_skin_mm"]
    designs = {
        "original (8 mm Al skin, 5 mm webs, 2 mm ribs)": rv.al_sections(0.008, 0.005, 0.002),
        f"sized Al ({ts} mm skin, 1.5 mm webs, 1.0 mm ribs)": rv.al_sections(ts * 1e-3, 1.5e-3, 1.0e-3),
        "CFRP QI [45/-45/0/90]s skin (1.0 mm)": {
            "skin": laminate_section(T300_5208, [45, -45, 0, 90, 90, 0, -45, 45]),
            "web": laminate_section(T300_5208, [45, -45, -45, 45]),
            "rib": laminate_section(T300_5208, [0, 90, 90, 0])},
    }
    rows = []
    for name, secs in designs.items():
        row = run_design(name, secs)
        rows.append(row)
        print(f"{name}: f_h {row[1]} Hz, f_alpha {row[2]} Hz, mu {row[4]}, EA {row[5]} c, CG {row[6]} c, "
              f"V_flutter {row[11]} m/s, V_divergence {row[12]} m/s ({row[13]} x V_dive {V_DIVE_ASSUMED:.0f} m/s)")
    rv.write_csv("F_flutter.csv", rows,
                 ["design", "f_bending_Hz", "f_torsion_Hz", "mass_per_span_kg_m", "mu", "EA_x_over_c",
                  "CG_x_over_c", "a", "x_alpha", "r_alpha^2", "sigma", "V_flutter_m_s", "V_divergence_m_s",
                  "V_flutter_over_V_dive"])


if __name__ == "__main__":
    main()
