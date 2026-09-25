"""Summarise an Ansys Mechanical solver input deck (ds.dat).

Prints, for every body in the model: element type, element count, node count
and bounding box; which bodies carry the Fixed Support nodes; which bodies the
pressure (SURF154) elements sit on; the loaded area and the total force; and
the shell section thicknesses.

Usage:
    python tools/inspect_ds_dat.py ansys/wing_fea_files/dp0/SYS-2/MECH/ds.dat

Only the Python standard library is needed.
"""
import math
import re
import sys
from collections import defaultdict


def parse(path):
    lines = open(path, errors="ignore").read().splitlines()
    nodes, bodies, sections = {}, [], []
    elem_nodes = defaultdict(set)
    elem_count = defaultdict(int)
    surf_elems, fixed, sfe = [], [], {}
    body_name = elem_type = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.lower().startswith("nblock"):
            i += 2
            while not lines[i].strip().startswith("-1"):
                p = lines[i].split()
                nodes[int(p[0])] = tuple(float(v) for v in p[1:4])
                i += 1
        m = re.match(r"/com,\*+ Elements for Body \d+ '(.+)' \*+", line)
        if m:
            body_name = m.group(1).split("|")[-1]
        m = re.match(r"et,(\d+),(\d+)", line)
        if m and body_name:
            bodies.append((body_name, int(m.group(1)), int(m.group(2))))
            body_name = None  # later et commands (contact, SURF154) are not bodies
        m = re.search(r"TYPE,(\d+)", line)
        if m:
            elem_type = int(m.group(1))
        if line.startswith("secdata"):
            sections.append(float(line.split(",")[1]))
        if line.startswith("eblock") and "COMPACT" in line.upper():
            i += 2
            while not lines[i].strip().startswith("-1"):
                ids = [int(x) for x in lines[i].split()][1:]
                elem_count[elem_type] += 1
                elem_nodes[elem_type].update(n for n in ids if n > 0)
                if elem_type == 48:
                    surf_elems.append(ids)
                i += 1
        if line.startswith("CMBLOCK,_FIXEDSU"):
            vals, j = [], i + 2
            while j < len(lines) and re.match(r"^\s*-?\d+(\s+-?\d+)*\s*$", lines[j]):
                vals += [int(t) for t in lines[j].split()]
                j += 1
            k = 0
            while k < len(vals):  # "a -b" means the range a..b
                if k + 1 < len(vals) and vals[k + 1] < 0:
                    fixed += range(vals[k], -vals[k + 1] + 1)
                    k += 2
                else:
                    fixed.append(vals[k])
                    k += 1
        m = re.match(r"sfe,all,(\d),pres,1,([-0-9.eE+]+)", line)
        if m:
            sfe[int(m.group(1))] = float(m.group(2))
        i += 1
    return nodes, bodies, sections, elem_nodes, elem_count, surf_elems, fixed, sfe


def facet_area(pts):
    def tri(a, b, c):
        u = [b[k] - a[k] for k in range(3)]
        v = [c[k] - a[k] for k in range(3)]
        cr = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
        return 0.5 * math.sqrt(sum(x * x for x in cr))
    if len(pts) == 3:
        return tri(*pts)
    return tri(pts[0], pts[1], pts[2]) + tri(pts[0], pts[2], pts[3])


def main(path):
    nodes, bodies, sections, elem_nodes, elem_count, surf_elems, fixed, sfe = parse(path)
    print(f"{path}\n{len(nodes)} nodes")
    xs, ys, zs = zip(*nodes.values())
    print(f"model bounding box  X {min(xs):.4f}..{max(xs):.4f}  Y {min(ys):.4f}..{max(ys):.4f}  Z {min(zs):.4f}..{max(zs):.4f} m\n")
    print(f"{'body':<10} {'ANSYS type':<10} {'elements':>8} {'nodes':>7}   bounding box [m]")
    for name, et, etype in bodies:
        ns = elem_nodes.get(et)
        if not ns:
            continue
        pts = [nodes[n] for n in ns]
        lo = [min(p[k] for p in pts) for k in range(3)]
        hi = [max(p[k] for p in pts) for k in range(3)]
        print(f"{name:<10} {etype:<10} {elem_count[et]:>8} {len(ns):>7}   "
              f"X[{lo[0]:.3f},{hi[0]:.3f}] Y[{lo[1]:.3f},{hi[1]:.3f}] Z[{lo[2]:.3f},{hi[2]:.3f}]")
    print(f"\nshell section thicknesses (secdata): {sections} m")
    fixed_set = set(fixed)
    print(f"\nFixed Support: {len(fixed_set)} nodes")
    for name, et, _ in bodies:
        shared = len(fixed_set & elem_nodes.get(et, set()))
        if shared:
            print(f"  on {name}: {shared} nodes")
    if surf_elems:
        area = 0.0
        for ids in surf_elems:
            corners = list(dict.fromkeys(ids[:4]))
            area += facet_area([nodes[n] for n in corners])
        loaded = {n for ids in surf_elems for n in ids if n > 0}
        on = [name for name, et, _ in bodies if loaded & elem_nodes.get(et, set())]
        comp = ", ".join(f"face {k}: {v:g} Pa" for k, v in sorted(sfe.items()))
        print(f"\nPressure (SURF154): {len(surf_elems)} facets on {on}, area {area:.4f} m^2 ({comp})")
        p = max(abs(v) for v in sfe.values()) if sfe else 0.0
        print(f"  resultant if applied to full area: {p * area:.1f} N")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
