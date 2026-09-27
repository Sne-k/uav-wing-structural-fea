"""Re-solve the saved Workbench decks in Ansys MAPDL, as saved and with the corrections.

Cases (all use the student's own mesh, contacts and materials from ansys/):
  R0  static deck exactly as saved                  -> must reproduce 7.14 mm / 154.9 MPa
  R1  static, whole root section clamped (y = 0)     -> compare with the wingfe re-analysis
  R2  as R1, but the 800 Pa is replaced by the design ultimate load
      (lifting-line lift for a 5 kg UAV at 6 g, on the upper skin only)
  M0  modal deck exactly as saved                   -> must reproduce 13.9 ... 87.7 Hz
  M1  modal, whole root section clamped, 12 modes

The original decks are never modified: each case writes a modified copy into a
temporary working folder, runs MAPDL in batch mode, and reads back a small
results file written by appended /POST1 commands.

Run:  python validation/ansys_rerun/rerun_in_mapdl.py     (needs Ansys 2026 R1, ~5-10 min)
Writes validation/ansys_rerun/static_results.csv and modal_results.csv
"""
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DP0 = os.path.join(REPO, "ansys", "wing_fea_files", "dp0")
sys.path.insert(0, os.path.join(REPO, "validation"))
from wingfe.aero import lifting_line  # noqa: E402
from wingfe.wing import SPAN, chord  # noqa: E402

AWP = os.environ.get("AWP_ROOT261", r"D:\apps\ANSYS\ANSYS Inc\ANSYS Student\v261")
MAPDL = os.path.join(AWP, "ansys", "bin", "winx64", "MAPDL.exe")
L_ULT = 6.0 * 5.0 * 9.81 / 2          # ultimate semi-span lift [N] (same design case as run_validation.py)

CLAMP = """! --- correction: clamp the whole root section (skin edge, spar root faces, root cap)
nsel,s,loc,y,-1e-6,1e-6
d,all,all
allsel,all
"""

STATIC_POST = """
/post1
set,last
allsel,all
nsort,u,sum
*get,umax,sort,0,max
nle=node(0,1.2,0)
nte=node(0.2,1.2,0)
*get,uzle,node,nle,u,z
*get,uzte,node,nte,u,z
shell,top
nsort,s,eqv
*get,sm_t,sort,0,max
*get,nm_t,sort,0,imax
shell,bot
nsort,s,eqv
*get,sm_b,sort,0,max
*get,nm_b,sort,0,imax
nm=nm_t
*if,sm_b,gt,sm_t,then
nm=nm_b
*endif
*get,ym,node,nm,loc,y
nsel,s,loc,y,0.05,1.3
shell,top
nsort,s,eqv
*get,sf_t,sort,0,max
shell,bot
nsort,s,eqv
*get,sf_b,sort,0,max
shell,top
allsel,all
! total vertical reaction: sum RF over all nodes of the root plane
! (Workbench does not write element nodal loads, so FSUM cannot be used)
nsel,s,loc,y,-1e-6,1e-6
*get,nn,node,0,count
rfz=0
nd=0
*do,k,1,nn
nd=ndnext(nd)
*get,rr,node,nd,rf,fz
rfz=rfz+rr
*enddo
allsel,all
%GROUPS%
*cfopen,static_results,txt
*vwrite,umax,uzle,uzte,sm_t,sm_b,ym,sf_t,sf_b,rfz
(9E16.7)
*vwrite,g1,g1y,g2,g2y,g3,g3y,g4,g4y
(8E16.7)
*cfclos
finish
"""

GROUP_TEMPLATE = """{select}
nsle,s
shell,top
nsort,s,eqv
*get,{g}t,sort,0,max
*get,{g}tn,sort,0,imax
shell,bot
nsort,s,eqv
*get,{g}b,sort,0,max
*get,{g}bn,sort,0,imax
{g}={g}t
{g}n={g}tn
*if,{g}b,gt,{g}t,then
{g}={g}b
{g}n={g}bn
*endif
*get,{g}y,node,{g}n,loc,y
shell,top
allsel,all
"""
GROUPS = (("g1", "skin", ("Wing_body",)), ("g2", "end caps", ("Body44", "Body45")),
          ("g3", "ribs", ("Rib1",) + tuple(f"rib{k}" for k in range(2, 11))), ("g4", "spars", ("spar1", "spar2")))


def group_post(body_types):
    """APDL that finds the max von Mises (top/bottom) and its y position per group of bodies."""
    out = []
    for g, _, names in GROUPS:
        types = [body_types[n] for n in names if n in body_types]
        sel = "\n".join(("esel,s" if i == 0 else "esel,a") + f",type,,{t}" for i, t in enumerate(types))
        out.append(GROUP_TEMPLATE.format(select=sel, g=g))
    return "".join(out)

MODAL_POST = """
/post1
allsel,all
nle=node(0,1.2,0)
nte=node(0.2,1.2,0)
*get,nset,active,0,set,nset
*cfopen,modal_results,txt
*do,i,1,nset
set,1,i
*get,fr,active,0,set,freq
*get,a,node,nle,u,z
*get,b,node,nte,u,z
*get,c,node,nle,u,x
nsort,u,sum
*get,um,sort,0,max
*get,umn,sort,0,imax
*get,ux_,node,umn,loc,x
*get,uy_,node,umn,loc,y
*get,uz_,node,umn,loc,z
*vwrite,i,fr,a,b,c,um,ux_,uy_,uz_
(F5.0,8E16.7)
*enddo
*cfclos
finish
"""


def naca2412_camber(xi):
    m, p = 0.02, 0.4
    return np.where(xi < p, m / p**2 * (2 * p * xi - xi**2), m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * xi - xi**2))


def read_deck(path):
    lines = open(path, errors="ignore").read().splitlines()
    nodes, facets, surf_types, cur = {}, [], set(), None
    body_types, body = {}, None
    i = 0
    while i < len(lines):
        s = lines[i]
        m = re.match(r"/com,\*+ Elements for Body \d+ '(.+)' \*+", s)
        if m:
            body = m.group(1).split("|")[-1]
        m = re.match(r"et,(\d+),(\d+)", s)
        if m and body:
            body_types[body] = int(m.group(1))
            body = None
        if s.lower().startswith("nblock"):
            i += 2
            while not lines[i].strip().startswith("-1"):
                p = lines[i].split()
                nodes[int(p[0])] = np.array([float(v) for v in p[1:4]])
                i += 1
        m = re.match(r"et,(\d+),154", s)
        if m:
            surf_types.add(int(m.group(1)))
        m = re.search(r"TYPE,(\d+)", s)
        if m:
            cur = int(m.group(1))
        if s.startswith("eblock") and "COMPACT" in s.upper():
            i += 2
            while not lines[i].strip().startswith("-1"):
                ids = [int(x) for x in lines[i].split()][1:]
                if cur in surf_types:
                    facets.append(list(dict.fromkeys(ids[:4])))
                i += 1
        i += 1
    return lines, nodes, facets, surf_types, body_types


def design_forces(nodes, facets):
    """Nodal forces (+z) on upper-skin nodes: lifting-line lift, centre of pressure at 25 % chord."""
    ys, shape, _, _ = lifting_line(chord, SPAN)
    F = {}
    for f in facets:
        P = np.array([nodes[n] for n in f[:3]])
        x, y, z = P.mean(axis=0)
        c = chord(y)
        xi = min(max(x / c, 0.0), 1.0)
        if z <= naca2412_camber(xi) * c:
            continue                                   # lower surface
        if len(f) == 4:
            P4 = np.array([nodes[n] for n in f])
            area_z = 0.5 * abs(np.cross(P4[2] - P4[0], P4[3] - P4[1])[2])
        else:
            area_z = 0.5 * abs(np.cross(P[1] - P[0], P[2] - P[0])[2])
        dp = L_ULT * np.interp(y, ys, shape) / c * 3 * (1 - xi)**2
        for n in f:
            F[n] = F.get(n, 0.0) + dp * area_z / len(f)
    total = sum(F.values())
    scale = L_ULT / total                              # remove the small faceting error
    return {n: v * scale for n, v in F.items()}, total


def make_deck(src, extra_before_solve="", modes=None, post=""):
    lines = open(src, errors="ignore").read().splitlines()
    idx = [i for i, s in enumerate(lines) if s.strip().lower() == "solve"]
    assert len(idx) == 1, f"expected one 'solve' in {src}"
    if modes:
        lines = [re.sub(r"^modopt,lanb,\d+", f"modopt,lanb,{modes}", s) for s in lines]
    out = lines[:idx[0]] + extra_before_solve.splitlines() + lines[idx[0]:] + post.splitlines()
    return "\n".join(out) + "\n"


def run(case, deck_text, work, extra_files=None):
    d = os.path.join(work, case)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "input.dat"), "w") as f:
        f.write(deck_text)
    for name, text in (extra_files or {}).items():
        with open(os.path.join(d, name), "w") as f:
            f.write(text)
    r = subprocess.run([MAPDL, "-b", "-np", "2", "-dir", d, "-j", "file", "-i", "input.dat", "-o", "solve.out"],
                       cwd=d, capture_output=True, text=True)
    res = os.path.join(d, "static_results.txt" if case.startswith("R") else "modal_results.txt")
    if not os.path.exists(res):
        tail = open(os.path.join(d, "solve.out"), errors="ignore").read()[-3000:] if os.path.exists(
            os.path.join(d, "solve.out")) else r.stdout[-3000:] + r.stderr[-3000:]
        raise RuntimeError(f"{case}: no results file (exit {r.returncode})\n{tail}")
    return [[float(v) for v in ln.split()] for ln in open(res) if ln.strip()]


def classify(fr, uz_le, uz_te, ux, umax):
    half_c = 0.1
    w = 0.5 * (uz_le + uz_te)
    twist = (uz_te - uz_le) / (2 * half_c)
    tip = max(abs(w), abs(ux), abs(twist) * half_c)
    if tip < 0.1 * umax:
        return "local / panel"
    if abs(ux) >= max(abs(w), abs(twist) * half_c):
        return "in-plane (chordwise) bending"
    if abs(twist) * half_c > abs(w):
        return "torsion"
    return "flapwise bending"


def main():
    if not os.path.exists(MAPDL):
        sys.exit(f"MAPDL not found at {MAPDL} (set AWP_ROOT261)")
    static_src = os.path.join(DP0, "SYS-2", "MECH", "ds.dat")
    modal_src = os.path.join(DP0, "SYS-3", "MECH", "ds.dat")
    _, nodes, facets, surf_types, body_types = read_deck(static_src)
    post = STATIC_POST.replace("%GROUPS%", group_post(body_types))
    F, raw_total = design_forces(nodes, facets)
    forces = "\n".join(f"f,{n},fz,{v:.8e}" for n, v in sorted(F.items())) + "\n"
    remove_pressure = "".join(f"esel,s,type,,{t}\nsfedele,all,all,pres\n" for t in sorted(surf_types)) + "allsel,all\n"
    print(f"design load: {len(F)} upper-skin nodes, lift {L_ULT:.2f} N (faceting error before scaling "
          f"{(raw_total / L_ULT - 1) * 100:+.2f} %)")

    work = tempfile.mkdtemp(prefix="wing_mapdl_")
    print("working folder:", work)
    static_rows, modal_rows = [], []
    for case, text, extra in (
            ("R0_as_saved", make_deck(static_src, post=post), None),
            ("R1_root_clamped", make_deck(static_src, CLAMP, post=post), None),
            ("R2_clamped_design_load", make_deck(static_src, CLAMP + remove_pressure + "/input,design_forces,inp\n",
                                                  post=post), {"design_forces.inp": forces})):
        v, grp = run(case, text, work, extra)[:2]
        umax, uzle, uzte, smt, smb, ym, sft, sfb, rfz = v
        row = [case, f"{umax * 1e3:.3f}", f"{0.5 * (uzle + uzte) * 1e3:.3f}", f"{max(smt, smb) / 1e6:.2f}",
               f"{ym:.3f}", f"{max(sft, sfb) / 1e6:.2f}", f"{-rfz:.2f}"]
        for k in range(len(GROUPS)):
            row += [f"{grp[2 * k] / 1e6:.2f}", f"{grp[2 * k + 1]:.3f}"]
        static_rows.append(row)
        print(f"{case}: max total deformation {row[1]} mm, tip w {row[2]} mm, max von Mises {row[3]} MPa "
              f"at y = {row[4]} m, beyond y = 0.05 m {row[5]} MPa, applied load (from reactions) {row[6]} N")
        print("   per part: " + ", ".join(f"{GROUPS[k][1]} {row[7 + 2 * k]} MPa at y = {row[8 + 2 * k]} m"
                                          for k in range(len(GROUPS))))
    for case, text in (("M0_as_saved", make_deck(modal_src, post=MODAL_POST)),
                       ("M1_root_clamped", make_deck(modal_src, CLAMP, modes=12, post=MODAL_POST))):
        for i, fr, a, b, c, um, px, py, pz in run(case, text, work):
            label = classify(fr, a, b, c, um)
            modal_rows.append([case, int(i), f"{fr:.3f}", label, f"{px:.3f}", f"{py:.3f}", f"{pz:.3f}"])
            print(f"{case} mode {int(i)}: {fr:.3f} Hz  {label:<30} peak at x = {px:.3f}, y = {py:.3f}, z = {pz:.3f} m")

    with open(os.path.join(HERE, "static_results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "max_total_deformation_mm", "tip_w_mm", "max_von_mises_MPa", "y_of_max_m",
                    "max_von_mises_beyond_y0.05_MPa", "applied_Fz_from_reactions_N"] +
                   [f"{s}_{g[1].replace(' ', '_')}" for g in GROUPS for s in ("max_vm_MPa", "y_m")])
        w.writerows(static_rows)
    with open(os.path.join(HERE, "modal_results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "mode", "f_Hz", "type", "peak_x_m", "peak_y_m", "peak_z_m"])
        w.writerows(modal_rows)
    shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
