"""Read the original Ansys result files (file.rst) without an Ansys installation.

Prints, for the aluminium and "CFRP" static runs: the total reaction force, the
maximum nodal von Mises stress per body, and how the von Mises stress decays
along the spars and skin from the root. It also prints the natural frequencies
of the modal runs. Output goes to validation/results/ansys_rst_readback.txt.

Requires:  pip install ansys-mapdl-reader
Note: the nodal displacement records of these 2026 R1 files are stored in a
compressed format that ansys-mapdl-reader does not decode (it returns zeros or
huge values). Stresses, reactions and frequencies read correctly and match the
Workbench screenshots exactly.
"""
import os
import re
import warnings

import numpy as np

warnings.filterwarnings("ignore")
from ansys.mapdl.reader import read_binary  # noqa: E402

HERE = os.path.dirname(__file__)
DP0 = os.path.join(HERE, "..", "ansys", "wing_fea_files", "dp0")
OUT = os.path.join(HERE, "results", "ansys_rst_readback.txt")


def bodies(ds_path):
    """Element type number -> body name, and the node set of each body (from ds.dat)."""
    lines = open(ds_path, errors="ignore").read().splitlines()
    names, nodes, cur, name, i = {}, {}, None, None, 0
    while i < len(lines):
        s = lines[i]
        m = re.match(r"/com,\*+ Elements for Body \d+ '(.+)' \*+", s)
        if m:
            name = m.group(1).split("|")[-1]
        m = re.match(r"et,(\d+),(\d+)", s)
        if m and name:
            names[int(m.group(1))] = name
            name = None
        m = re.search(r"TYPE,(\d+)", s)
        if m:
            cur = int(m.group(1))
        if s.startswith("eblock") and "COMPACT" in s.upper():
            i += 2
            while not lines[i].strip().startswith("-1"):
                ids = [int(x) for x in lines[i].split()][1:]
                nodes.setdefault(cur, set()).update(n for n in ids if n > 0)
                i += 1
        i += 1
    return {names[t]: np.array(sorted(nodes[t])) for t in names}


def von_mises(s):
    sx, sy, sz, sxy, syz, sxz = s.T
    return np.sqrt(0.5 * ((sx - sy)**2 + (sy - sz)**2 + (sz - sx)**2) + 3 * (sxy**2 + syz**2 + sxz**2))


def static_report(sysdir, label, out):
    r = read_binary(os.path.join(DP0, sysdir, "MECH", "file.rst"))
    bn = bodies(os.path.join(DP0, sysdir, "MECH", "ds.dat"))
    xyz = r.mesh.nodes
    idx = {n: i for i, n in enumerate(r.mesh.nnum)}
    sn, s = r.nodal_stress(0)
    vm = von_mises(s)
    pos = {n: i for i, n in enumerate(sn)}
    out.append(f"\n== {label} ({sysdir}) ==")
    rf, _, dof = r.nodal_reaction_forces(0)
    tot = {d: rf[dof == d].sum() for d in np.unique(dof)}
    out.append("total reaction by DOF (1 UX, 2 UY, 3 UZ): " +
               ", ".join(f"{int(k)}: {v:.3f} N" for k, v in tot.items() if int(k) <= 3))
    out.append(f"{'body':<10} {'max vM [MPa]':>12}  at y [m]")
    for b, nodes in bn.items():
        ii = [pos[n] for n in nodes if n in pos]
        v = vm[ii]
        if np.all(np.isnan(v)) or np.nanmax(v) == 0:
            continue
        j = int(np.nanargmax(v))
        out.append(f"{b:<10} {np.nanmax(v) / 1e6:12.2f}  {xyz[idx[nodes[j]]][1]:.3f}")
    bands = ((0, 0.005), (0.005, 0.02), (0.02, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.6),
             (0.6, 0.9), (0.9, 1.2))
    out.append("max vM [MPa] in spanwise bands from the root:")
    out.append(f"{'':<10} " + " ".join(f"{a:.3f}-{b:.2f}" for a, b in bands))
    for b in ("spar1", "spar2", "Wing_body"):
        nodes = bn[b]
        Y = np.array([xyz[idx[n]][1] for n in nodes])
        V = np.array([vm[pos[n]] if n in pos else np.nan for n in nodes])
        out.append(f"{b:<10} " + " ".join(f"{np.nanmax(V[(Y >= a) & (Y < c)]) / 1e6:10.1f}" for a, c in bands))
    # Only the invariant von Mises value is reported: shell stress components in
    # the result file may be stored in element axes rather than global axes.


def main():
    out = []
    static_report("SYS-2", "static, aluminium", out)
    static_report("SYS-4", "static, CFRP (isotropic)", out)
    for sysdir, label in (("SYS-3", "modal, aluminium"), ("SYS-5", "modal, CFRP")):
        r = read_binary(os.path.join(DP0, sysdir, "MECH", "file.rst"))
        out.append(f"\n== {label} ({sysdir}) ==  frequencies [Hz]: " +
                   ", ".join(f"{f:.3f}" for f in r.time_values))
    text = "\n".join(out)
    print(text)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
