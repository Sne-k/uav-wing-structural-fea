"""Re-analysis of the project wing with a corrected set-up.

Part A  Replication: original geometry, thicknesses and load (800 Pa along +z on
        the whole skin), but with the root section clamped. Mesh convergence,
        comparison with the Ansys run and with beam theory, modal classification.
Part B  Design load case: lift from lifting-line theory for a 5 kg UAV at the
        ultimate load factor, applied to the upper skin.
Part C  Aluminium skin sizing sweep with linear buckling.
Part D  CFRP laminates (T300/5208, CLT, Tsai-Wu) versus aluminium.
Part E  Fatigue at the critical location of the sized aluminium design.

Run:  python validation/run_validation.py        (about 10-15 minutes)
Results: validation/results/*.csv, *.json, *.png
"""
import csv
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from wingfe.aero import lifting_line, schrenk  # noqa: E402
from wingfe.materials import (AL6061_T6, T300_5208, Isotropic, engineering_constants,  # noqa: E402
                              isotropic_section, laminate_section)
from wingfe.model import modal_participation  # noqa: E402
from wingfe.wing import (SPAN, Layout, Mesh, build_wing, chord, isotropic_stresses,  # noqa: E402
                         laminate_strength, load_component_z_all_skin, load_upper_skin)

OUT = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUT, exist_ok=True)
G0 = 9.81

# ----------------------------------------------------------------- design inputs
DESIGN = dict(
    mtow_kg=5.0,             # all-up mass assumed during the project
    n_limit=4.0,             # limit manoeuvre load factor (CS-23 normal category ~3.8)
    safety_factor=1.5,       # ultimate = 1.5 x limit
    cl_max=1.2,
    note="Inertia relief from the wing's own mass is ignored (conservative).",
)
DESIGN["n_ult"] = DESIGN["n_limit"] * DESIGN["safety_factor"]
W = DESIGN["mtow_kg"] * G0
S_REF = 2 * (0.25 + 0.20) / 2 * SPAN
DESIGN["semi_span_lift_limit_N"] = DESIGN["n_limit"] * W / 2
DESIGN["semi_span_lift_ult_N"] = DESIGN["n_ult"] * W / 2
DESIGN["stall_speed_ms"] = math.sqrt(2 * W / (1.225 * S_REF * DESIGN["cl_max"]))
DESIGN["manoeuvre_speed_ms"] = DESIGN["stall_speed_ms"] * math.sqrt(DESIGN["n_limit"])

ANSYS = dict(tip_deflection_mm=7.143, max_vm_MPa=154.91, mass_kg=12.419,
             freqs_Hz=[13.892, 23.863, 36.981, 63.935, 68.116, 87.749],
             cfrp_freqs_Hz=[18.053, 30.997, 47.779, 83.042, 88.201, 113.86])
BEAM = dict(tip_deflection_mm=2.86, root_stress_MPa=10.2, f1_Hz=22.5)

LOG = []


def log(msg):
    print(msg, flush=True)
    LOG.append(msg)


def write_csv(name, rows, header):
    with open(os.path.join(OUT, name), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


# ----------------------------------------------------------------- helpers
def al_sections(skin, web, rib, mat=AL6061_T6):
    return {"skin": isotropic_section(mat, skin), "web": isotropic_section(mat, web),
            "rib": isotropic_section(mat, rib)}


def solve_case(sections, F_fn, mesh=Mesh(), modes=0, buckling=False, layout=Layout()):
    t0 = time.time()
    m, info = build_wing(mesh, layout, sections)
    m.build()
    m.fix_nodes(info["root"])
    F = F_fn(m)
    u = m.solve_static(F)
    res = dict(model=m, info=info, u=u, F=F, dof=m.ndof, mass=m.mass())
    tip = info["tip"]
    res["tip_w_mm"] = u[6 * tip + 2].mean() * 1e3
    lt = tip[np.argmin(m.nodes[tip, 0])]
    tt = tip[np.argmax(m.nodes[tip, 0])]
    res["tip_twist_deg"] = math.degrees((u[6 * tt + 2] - u[6 * lt + 2]) /
                                        (m.nodes[tt, 0] - m.nodes[lt, 0]))
    res["reaction_z"] = -m.reactions[6 * info["root"] + 2].sum()
    res["applied_z"] = F[2::6].sum()
    if modes:
        f, phi = m.solve_modes(modes)
        res["freqs"], res["phi"] = f, phi
        res["participation"] = modal_participation(m, phi)
        res["mode_labels"] = classify_modes(m, info, phi, res["participation"])
    if buckling:
        lam, bphi = m.solve_buckling(u, k=4)
        res["buckling"], res["buckling_phi"] = lam, bphi
    res["time_s"] = time.time() - t0
    return res


def classify_modes(m, info, phi, P):
    """Label modes from the motion of the tip section (translation, twist, fore-aft)."""
    tip = info["tip"]
    le, te = tip[np.argmin(m.nodes[tip, 0])], tip[np.argmax(m.nodes[tip, 0])]
    half_c = 0.5 * (m.nodes[te, 0] - m.nodes[le, 0])
    labels, count = [], {}
    for j in range(phi.shape[1]):
        p = phi[:, j]
        w = p[6 * tip + 2].mean()
        ux = p[6 * tip + 0].mean()
        twist = (p[6 * te + 2] - p[6 * le + 2]) / (2 * half_c)
        tip_motion = max(abs(w), abs(ux), abs(twist) * half_c)
        peak = np.abs(np.stack([p[0::6], p[1::6], p[2::6]])).max()
        if tip_motion < 0.1 * peak:
            kind = "local / panel"
        elif abs(ux) >= max(abs(w), abs(twist) * half_c):
            kind = "in-plane (chordwise) bending"
        elif abs(twist) * half_c > abs(w):
            kind = "torsion"
        else:
            kind = "flapwise bending"
        count[kind] = count.get(kind, 0) + 1
        labels.append(f"{kind} #{count[kind]}")
    return labels


def stress_summary(res, groups=None, away=0.05):
    g, a = isotropic_stresses(res["model"], res["u"], groups)
    far = a[:, 1] > away
    k = a[:, 3].argmax()
    kf = np.where(far)[0][a[far, 3].argmax()]
    return dict(max_vm_MPa=a[k, 3] / 1e6, max_vm_at=(g[k], *np.round(a[k, :3], 4)),
                max_vm_far_MPa=a[kf, 3] / 1e6, g=g, a=a)


def design_pressure(total_lift):
    ys, shape, CLa, e = lifting_line(chord, SPAN)

    def dp(x, y):
        c = chord(y)
        xi = min(max(x / c, 0.0), 1.0)
        l = total_lift * np.interp(y, ys, shape)          # N/m
        return l / c * 3 * (1 - xi)**2                     # centre of pressure at 25 % chord
    return dp


# ================================================================= PART A
def part_a():
    log("\n=== PART A: original geometry and load, root clamped ===")
    secs = al_sections(0.008, 0.005, 0.002)
    rows = []
    conv = {}
    for fct in (0.5, 1.0, 1.5, 2.0):
        mesh = Mesh().scaled(fct)
        r = solve_case(secs, lambda m: load_component_z_all_skin(m, 800.0), mesh, modes=8)
        s = stress_summary(r)
        fl = r["freqs"]
        lab = r["mode_labels"]
        f_tor = next(f for f, l in zip(fl, lab) if l.startswith("torsion"))
        f_inp = next(f for f, l in zip(fl, lab) if l.startswith("in-plane"))
        rows.append([fct, r["dof"], f"{r['tip_w_mm']:.4f}", f"{s['max_vm_MPa']:.3f}",
                     f"{s['max_vm_far_MPa']:.3f}", f"{fl[0]:.3f}", f"{f_inp:.2f}", f"{f_tor:.2f}",
                     f"{r['mass']:.3f}", f"{r['time_s']:.0f}"])
        log(f"mesh x{fct}: dof {r['dof']}, tip w {r['tip_w_mm']:.3f} mm, max vM {s['max_vm_MPa']:.2f} MPa "
            f"(y>0.05 m: {s['max_vm_far_MPa']:.2f}), f1 {fl[0]:.2f} Hz, in-plane {f_inp:.1f} Hz, "
            f"torsion {f_tor:.1f} Hz, mass {r['mass']:.3f} kg, reaction {r['reaction_z']:.1f} N")
        conv[fct] = (r, s)
    write_csv("A_mesh_convergence.csv", rows,
              ["mesh_factor", "dof", "tip_w_mm", "max_vm_MPa", "max_vm_y>0.05_MPa", "f1_Hz",
               "f_inplane_Hz", "f_torsion_Hz", "mass_kg", "time_s"])
    r, s = conv[2.0]
    modes = [[i + 1, f"{f:.2f}", lab, *[f"{p:.3f}" for p in pp]]
             for i, (f, lab, pp) in enumerate(zip(r["freqs"], r["mode_labels"], r["participation"]))]
    write_csv("A_modes_aluminium.csv", modes, ["mode", "f_Hz", "type", "X", "Y", "Z", "RX", "RY", "RZ"])
    for row in modes:
        log(f"  mode {row[0]}: {row[1]} Hz  {row[2]}")

    # CFRP as defined in the project: isotropic, E = 70 GPa, nu = 0.3, rho = 1600
    cf = Isotropic(E=70e9, nu=0.3, rho=1600.0, name="CFRP (isotropic, project)")
    rc = solve_case(al_sections(0.008, 0.005, 0.002, cf),
                    lambda m: load_component_z_all_skin(m, 800.0), Mesh().scaled(2.0), modes=8)
    sc = stress_summary(rc)
    log(f"project 'CFRP': tip w {rc['tip_w_mm']:.3f} mm, max vM {sc['max_vm_MPa']:.2f} MPa, "
        f"f1 {rc['freqs'][0]:.2f} Hz (ratio {rc['freqs'][0] / r['freqs'][0]:.3f}), mass {rc['mass']:.3f} kg")
    return conv, rc, sc


# ================================================================= PART B
def part_b():
    log("\n=== PART B: design load case ===")
    for k, v in DESIGN.items():
        log(f"  {k}: {v if isinstance(v, str) else round(v, 3)}")
    ys, sh_llt, CLa, e = lifting_line(chord, SPAN)
    ys2, sh_s = schrenk(chord, SPAN)
    log(f"  lifting line: AR {(2 * SPAN)**2 / S_REF:.2f}, span efficiency e = {e:.3f}, CL_alpha = {CLa:.3f} /rad")
    rows = [[f"{y:.3f}", f"{a * DESIGN['semi_span_lift_ult_N']:.3f}", f"{b * DESIGN['semi_span_lift_ult_N']:.3f}"]
            for y, a, b in zip(ys, sh_llt, np.interp(ys, ys2, sh_s))]
    write_csv("B_spanwise_lift_ultimate.csv", rows, ["y_m", "lift_llt_N_per_m", "lift_schrenk_N_per_m"])
    secs = al_sections(0.008, 0.005, 0.002)
    L = DESIGN["semi_span_lift_ult_N"]
    r = solve_case(secs, lambda m: load_upper_skin(m, design_pressure(L)), Mesh().scaled(1.5))
    s = stress_summary(r)
    log(f"  original 8 mm design at ULTIMATE design load ({L:.1f} N/semi-span, applied {r['applied_z']:.1f} N): "
        f"tip w {r['tip_w_mm']:.3f} mm, max vM {s['max_vm_MPa']:.2f} MPa, "
        f"reserve factor on yield {AL6061_T6.Sy / 1e6 / s['max_vm_MPa']:.0f}, semi-span mass {r['mass']:.2f} kg "
        f"({2 * r['mass'] / DESIGN['mtow_kg'] * 100:.0f} % of the 5 kg MTOW for both wings)")
    return ys, sh_llt, ys2, sh_s, r, s


# ================================================================= PART C
def part_c():
    log("\n=== PART C: aluminium skin sizing (webs 1.5 mm, ribs 1.0 mm) ===")
    L_ult = DESIGN["semi_span_lift_ult_N"]
    L_lim = DESIGN["semi_span_lift_limit_N"]
    rows, results = [], {}
    for ts in (0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5):
        secs = al_sections(ts * 1e-3, 1.5e-3, 1.0e-3)
        r = solve_case(secs, lambda m: load_upper_skin(m, design_pressure(L_ult)), Mesh().scaled(2.0),
                       modes=8, buckling=True)
        s = stress_summary(r)
        lam = r["buckling"][0]
        f = r["freqs"]
        lab = r["mode_labels"]
        f_tor = next((x for x, l in zip(f, lab) if l.startswith("torsion")), float("nan"))
        tip_lim = r["tip_w_mm"] * L_lim / L_ult
        rf_yield = AL6061_T6.Sy / 1e6 / s["max_vm_MPa"]
        ok = lam >= 1.0 and rf_yield >= 1.0
        rows.append([ts, f"{r['mass']:.3f}", f"{tip_lim:.2f}", f"{s['max_vm_MPa']:.1f}", f"{rf_yield:.2f}",
                     f"{lam:.3f}", f"{f[0]:.2f}", f"{f_tor:.1f}", "yes" if ok else "no"])
        log(f"  skin {ts} mm: mass {r['mass']:.3f} kg, tip w @limit {tip_lim:.1f} mm, max vM @ult "
            f"{s['max_vm_MPa']:.1f} MPa (RF {rf_yield:.2f}), buckling factor {lam:.3f} x ultimate, "
            f"f1 {f[0]:.1f} Hz, torsion {f_tor:.1f} Hz -> {'OK' if ok else 'fails'}")
        results[ts] = (r, s)
    write_csv("C_aluminium_sizing.csv", rows,
              ["skin_mm", "semi_span_mass_kg", "tip_w_limit_mm", "max_vm_ult_MPa", "RF_yield",
               "buckling_factor_vs_ultimate", "f1_Hz", "f_torsion_Hz", "meets_strength_and_buckling"])
    return results


# ================================================================= PART D
def part_d():
    log("\n=== PART D: CFRP laminates (T300/5208, ply 0.125 mm) ===")
    L_ult = DESIGN["semi_span_lift_ult_N"]
    L_lim = DESIGN["semi_span_lift_limit_N"]
    web = laminate_section(T300_5208, [45, -45, -45, 45])            # 0.5 mm +/-45 shear webs
    rib = laminate_section(T300_5208, [0, 90, 90, 0])                 # 0.5 mm ribs
    skins = {
        "QI [45/-45/0/90]s (1.0 mm)": [45, -45, 0, 90, 90, 0, -45, 45],
        "[45/-45/0/0]s (1.0 mm)": [45, -45, 0, 0, 0, 0, -45, 45],
        "[45/-45/0]s (0.75 mm)": [45, -45, 0, 0, -45, 45],
    }
    rows, results = [], {}
    for name, lay in skins.items():
        sk = laminate_section(T300_5208, lay)
        ec = engineering_constants(sk)
        secs = {"skin": sk, "web": web, "rib": rib}
        r = solve_case(secs, lambda m: load_upper_skin(m, design_pressure(L_ult)), Mesh().scaled(2.0),
                       modes=8, buckling=True)
        sr, where = laminate_strength(r["model"], r["u"], groups=("skin_upper", "skin_lower", "web1", "web2"))
        lam = r["buckling"][0]
        f = r["freqs"]
        lab = r["mode_labels"]
        f_tor = next((x for x, l in zip(f, lab) if l.startswith("torsion")), float("nan"))
        tip_lim = r["tip_w_mm"] * L_lim / L_ult
        ok = lam >= 1.0 and sr >= 1.0
        rows.append([name, f"{ec['Ex'] / 1e9:.1f}", f"{ec['Gxy'] / 1e9:.1f}", f"{r['mass']:.3f}",
                     f"{tip_lim:.2f}", f"{sr:.2f}", f"{lam:.3f}", f"{f[0]:.2f}", f"{f_tor:.1f}",
                     "yes" if ok else "no"])
        log(f"  {name}: Ex {ec['Ex'] / 1e9:.1f} GPa, Gxy {ec['Gxy'] / 1e9:.1f} GPa, mass {r['mass']:.3f} kg, "
            f"tip w @limit {tip_lim:.1f} mm, Tsai-Wu strength ratio @ult {sr:.2f} "
            f"({where['group']}, {where['angle']} deg ply), buckling {lam:.3f} x ult, f1 {f[0]:.1f} Hz, "
            f"torsion {f_tor:.1f} Hz -> {'OK' if ok else 'fails'}")
        results[name] = r
    write_csv("D_cfrp_laminates.csv", rows,
              ["skin_layup", "Ex_GPa", "Gxy_GPa", "semi_span_mass_kg", "tip_w_limit_mm",
               "tsai_wu_strength_ratio_ult", "buckling_factor_vs_ultimate", "f1_Hz", "f_torsion_Hz",
               "meets_strength_and_buckling"])
    return results


# ================================================================= PART E
SN = [(1e3, 300e6), (1e4, 250e6), (1e5, 180e6), (1e6, 120e6), (1e7, 90e6)]   # project S-N curve


def sn_life(sa_eq):
    cyc, s = zip(*SN)
    if sa_eq >= s[0]:
        return cyc[0]
    if sa_eq <= s[-1]:
        return float("inf")
    return 10**np.interp(-math.log10(sa_eq), [-math.log10(v) for v in s], [math.log10(c) for c in cyc])


def part_e(al_results, chosen):
    log(f"\n=== PART E: fatigue of the sized aluminium design (skin {chosen} mm) ===")
    r, s = al_results[chosen]
    L_ult = DESIGN["semi_span_lift_ult_N"]
    g, a = s["g"], s["a"]
    lower = np.char.startswith(g.astype(str), "skin_lower")
    per_g = 1.0 / (DESIGN["n_ult"])                    # stress per 1 g relative to the ultimate case
    sig_ult_tension = a[lower, 4].max()                # spanwise membrane tension, lower skin
    sig_ult_vm = a[:, 3].max()
    rows = []
    for label, base, amp in (("0 -> limit (4 g) manoeuvre cycles", DESIGN["n_limit"] / 2, DESIGN["n_limit"] / 2),
                             ("1 g +/- 0.5 g gust cycles", 1.0, 0.5),
                             ("1 g +/- 1.0 g gust cycles", 1.0, 1.0)):
        for comp, sig in (("lower-skin spanwise tension", sig_ult_tension), ("von Mises (project method)", sig_ult_vm)):
            sm = sig * per_g * base
            sa = sig * per_g * amp
            sa_eq = sa / (1 - sm / AL6061_T6.Su)       # Goodman
            N = sn_life(sa_eq)
            rows.append([label, comp, f"{sm / 1e6:.2f}", f"{sa / 1e6:.2f}", f"{sa_eq / 1e6:.2f}",
                         "> 1e7 (below the lowest S-N point)" if N == float("inf") else f"{N:.3g}"])
            log(f"  {label}, {comp}: mean {sm / 1e6:.1f} MPa, amplitude {sa / 1e6:.1f} MPa, Goodman "
                f"{sa_eq / 1e6:.1f} MPa -> life {'> 1e7' if N == float('inf') else f'{N:.3g}'} cycles")
    write_csv("E_fatigue.csv", rows, ["cycle", "stress_component", "mean_MPa", "amplitude_MPa",
                                      "goodman_equivalent_MPa", "life_cycles"])


# ================================================================= plots
def plots(conv, ys, sh_llt, ys2, sh_s, al_results, cf_results, chosen):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1 mesh convergence
    fct = sorted(conv)
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    ax[0].plot([conv[f][0]["dof"] for f in fct], [conv[f][0]["tip_w_mm"] for f in fct], "o-")
    ax[0].axhline(BEAM["tip_deflection_mm"], ls="--", c="gray", label="beam theory")
    ax[0].set(xlabel="DOF", ylabel="tip deflection [mm]", title="Tip deflection")
    ax[1].plot([conv[f][0]["dof"] for f in fct], [conv[f][1]["max_vm_MPa"] for f in fct], "o-", label="max")
    ax[1].plot([conv[f][0]["dof"] for f in fct], [conv[f][1]["max_vm_far_MPa"] for f in fct], "s-", label="y > 0.05 m")
    ax[1].axhline(BEAM["root_stress_MPa"], ls="--", c="gray", label="beam theory")
    ax[1].set(xlabel="DOF", ylabel="von Mises [MPa]", title="Peak stress")
    ax[2].plot([conv[f][0]["dof"] for f in fct], [conv[f][0]["freqs"][0] for f in fct], "o-")
    ax[2].axhline(BEAM["f1_Hz"], ls="--", c="gray", label="Rayleigh")
    ax[2].set(xlabel="DOF", ylabel="f1 [Hz]", title="1st bending frequency")
    for a in ax:
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
    fig.suptitle("Part A - mesh convergence (original design, 800 Pa, clamped root)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "A_mesh_convergence.png"), dpi=150)
    plt.close(fig)

    # 2 spanwise deflection and skin stress, corrected model vs beam theory
    r, s = conv[2.0]
    m, u, info = r["model"], r["u"], r["info"]
    top = [row[-1] for row in info["webs"][1]]
    y = m.nodes[top, 1]
    w = u[6 * np.array(top) + 2] * 1e3
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    ax[0].plot(y, w, label="wingfe shell model (root clamped)")
    ax[0].plot([0, SPAN], [0, ANSYS["tip_deflection_mm"]], "r--", label="Ansys tip value 7.14 mm")
    ax[0].set(xlabel="span y [m]", ylabel="deflection w [mm]", title="Deflection along the front spar")
    g, a = s["g"], s["a"]
    for name, mk in (("skin_upper", "min"), ("skin_lower", "max")):
        sel = np.char.startswith(g.astype(str), name)
        bins = np.linspace(0, SPAN, 25)
        idx = np.digitize(a[sel, 1], bins)
        vals = [getattr(np, mk)(a[sel, 4][idx == k]) / 1e6 for k in range(1, len(bins))]
        ax[1].plot(0.5 * (bins[1:] + bins[:-1]), vals, label=f"{name.replace('_', ' ')} spanwise stress")
    ax[1].axhline(0, c="k", lw=0.5)
    ax[1].set(xlabel="span y [m]", ylabel="stress [MPa]", title="Skin bending stress")
    for a_ in ax:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "A_deflection_and_skin_stress.png"), dpi=150)
    plt.close(fig)

    # 3 load distributions
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    L = DESIGN["semi_span_lift_ult_N"]
    ax.plot(ys, sh_llt * L, label="lifting line (used)")
    ax.plot(ys2, sh_s * L, "--", label="Schrenk")
    uni = np.array([chord(v) for v in ys]) * 800 * 2.02
    ax.plot(ys, uni, ":", label="project: 800 Pa on whole skin")
    ax.set(xlabel="span y [m]", ylabel="lift per span [N/m]",
           title=f"Spanwise load, ultimate {DESIGN['n_ult']:.0f} g for 5 kg UAV")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "B_spanwise_load.png"), dpi=150)
    plt.close(fig)

    # 4 sizing trade
    ts = sorted(al_results)
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    ax[0].plot(ts, [al_results[t][0]["mass"] for t in ts], "o-", label="aluminium")
    ax[1].plot(ts, [al_results[t][1]["max_vm_MPa"] for t in ts], "o-")
    ax[1].axhline(AL6061_T6.Sy / 1e6, ls="--", c="r", label="yield 276 MPa")
    ax[2].plot(ts, [al_results[t][0]["buckling"][0] for t in ts], "o-")
    ax[2].axhline(1.0, ls="--", c="r", label="buckling at ultimate")
    for name, rr in cf_results.items():
        ax[0].plot([float(name.split("(")[1].split()[0])], [rr["mass"]], "s", label=name)
    ax[0].set(xlabel="skin thickness [mm]", ylabel="semi-span mass [kg]", title="Mass")
    ax[1].set(xlabel="skin thickness [mm]", ylabel="max von Mises @ ultimate [MPa]", title="Strength")
    ax[2].set(xlabel="skin thickness [mm]", ylabel="buckling load factor", title="Skin buckling")
    for a_ in ax:
        a_.grid(alpha=0.3)
        a_.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "C_aluminium_sizing.png"), dpi=150)
    plt.close(fig)

    # 5 mode shapes (original geometry, clamped) and first buckling mode of the chosen design
    def surface_plot(ax, m, disp, title, scale):
        for conn, grp in m.elements:
            if not grp.startswith("skin_upper"):
                continue
            X = m.nodes[conn] + scale * disp[conn]
            ax.plot_trisurf(X[:, 0], X[:, 1], X[:, 2], color="tab:blue", alpha=0.6, linewidth=0)
        ax.set_title(title, fontsize=9)
        ax.set_box_aspect((0.25, 1.2, 0.25))
        ax.set_axis_off()

    phi = r["phi"]
    fig = plt.figure(figsize=(12, 4))
    for k in range(4):
        d = np.stack([phi[0::6, k], phi[1::6, k], phi[2::6, k]], axis=1)
        d /= np.abs(d).max()
        ax = fig.add_subplot(1, 4, k + 1, projection="3d")
        surface_plot(ax, m, d, f"mode {k + 1}: {r['freqs'][k]:.1f} Hz\n{r['mode_labels'][k]}", 0.08)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "A_mode_shapes.png"), dpi=150)
    plt.close(fig)

    rb, sb = al_results[chosen]
    bphi = rb["buckling_phi"][:, 0]
    d = np.stack([bphi[0::6], bphi[1::6], bphi[2::6]], axis=1)
    d /= np.abs(d).max()
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection="3d")
    surface_plot(ax, rb["model"], d, f"1st buckling mode, {chosen} mm Al skin\n"
                 f"load factor {rb['buckling'][0]:.2f} x ultimate", 0.03)
    ax.view_init(elev=35, azim=-60)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "C_buckling_mode.png"), dpi=150)
    plt.close(fig)


def main():
    t0 = time.time()
    conv, rc, sc = part_a()
    ys, sh_llt, ys2, sh_s, rb, sb = part_b()
    al = part_c()
    ok = [t for t in sorted(al) if al[t][0]["buckling"][0] >= 1.0 and
          AL6061_T6.Sy / 1e6 / al[t][1]["max_vm_MPa"] >= 1.0]
    chosen = ok[0] if ok else max(al)
    log(f"\nlightest aluminium skin meeting yield and buckling at ultimate: {chosen} mm")
    cf = part_d()
    part_e(al, chosen)
    plots(conv, ys, sh_llt, ys2, sh_s, al, cf, chosen)

    rA, sA = conv[2.0]
    summary = dict(
        design=DESIGN,
        ansys_original=ANSYS,
        beam_theory=BEAM,
        corrected_original=dict(tip_w_mm=rA["tip_w_mm"], max_vm_MPa=sA["max_vm_MPa"],
                                max_vm_far_MPa=sA["max_vm_far_MPa"], mass_kg=rA["mass"],
                                reaction_z_N=rA["reaction_z"],
                                modes=[dict(f=float(f), type=l) for f, l in zip(rA["freqs"], rA["mode_labels"])]),
        corrected_project_cfrp=dict(tip_w_mm=rc["tip_w_mm"], max_vm_MPa=sc["max_vm_MPa"],
                                    mass_kg=rc["mass"], f1=float(rc["freqs"][0])),
        design_load_on_original=dict(tip_w_mm=rb["tip_w_mm"], max_vm_MPa=sb["max_vm_MPa"], mass_kg=rb["mass"]),
        chosen_aluminium_skin_mm=chosen,
        runtime_min=(time.time() - t0) / 60,
    )
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    with open(os.path.join(OUT, "run_log.txt"), "w") as f:
        f.write("\n".join(LOG))
    log(f"\ndone in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
