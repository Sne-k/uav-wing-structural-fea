"""Benchmarks for the wingfe shell element against closed-form / published results.

Run:  python validation/benchmarks.py
Writes validation/results/benchmarks.csv and prints a table.
"""
import csv
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from wingfe.materials import Isotropic, Ply, isotropic_section, laminate_section  # noqa: E402
from wingfe.model import ShellModel  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "results")


def grid(nx, ny, f):
    """Structured quad grid; f(i/nx, j/ny) -> xyz."""
    ids = np.arange((nx + 1) * (ny + 1)).reshape(nx + 1, ny + 1)
    nodes = np.array([f(i / nx, j / ny) for i in range(nx + 1) for j in range(ny + 1)])
    elems = [[ids[i, j], ids[i + 1, j], ids[i + 1, j + 1], ids[i, j + 1]]
             for i in range(nx) for j in range(ny)]
    return nodes, elems, ids


def consistent_edge_load(ids_edge, nodes, ndof_dir, total, ndof):
    """Distribute a total force uniformly along a straight node line (trapezoidal)."""
    F = np.zeros(ndof)
    n = len(ids_edge)
    w = np.ones(n)
    w[0] = w[-1] = 0.5
    w /= w.sum()
    for k, nid in enumerate(ids_edge):
        F[6 * nid + ndof_dir] += total * w[k]
    return F


def cantilever_plate(nx=20, ny=2):
    L, b, t = 1.0, 0.1, 0.01
    mat = Isotropic(E=70e9, nu=0.0, rho=2700.0)
    nodes, elems, ids = grid(nx, ny, lambda s, r: (s * L, r * b, 0.0))
    m = ShellModel(nodes, [(e, "p") for e in elems], {"p": isotropic_section(mat, t)},
                   {"p": np.array([1.0, 0, 0])}).build()
    m.fix_nodes(ids[0, :])
    P = 100.0
    F = consistent_edge_load(ids[-1, :], nodes, 2, P, m.ndof)
    u = m.solve_static(F)
    I = b * t**3 / 12
    exact = P * L**3 / (3 * mat.E * I) + P * L / (5 / 6 * mat.G * b * t)
    fe = u[6 * ids[-1, :] + 2].mean()
    f, _ = m.solve_modes(3)
    f_exact = 1.87510**2 / (2 * math.pi * L**2) * math.sqrt(mat.E * I / (mat.rho * b * t))
    return [("cantilever plate: tip deflection [m]", fe, exact),
            ("cantilever plate: 1st bending frequency [Hz]", f[0], f_exact),
            ("cantilever plate: rigid-body check (max |K r| / max K_ii)", m.rigid_body_check().max(), 0.0)]


def inplane_beam(nx=10, ny=1):
    L, h, t = 1.0, 0.1, 0.01
    mat = Isotropic(E=70e9, nu=0.3, rho=2700.0)
    nodes, elems, ids = grid(nx, ny, lambda s, r: (s * L, r * h, 0.0))
    m = ShellModel(nodes, [(e, "w") for e in elems], {"w": isotropic_section(mat, t)},
                   {"w": np.array([1.0, 0, 0])}).build()
    m.fix_nodes(ids[0, :])
    P = 1000.0
    F = consistent_edge_load(ids[-1, :], nodes, 1, P, m.ndof)
    u = m.solve_static(F)
    I = t * h**3 / 12
    exact = P * L**3 / (3 * mat.E * I) + P * L / (5 / 6 * mat.G * h * t)
    fe = u[6 * ids[-1, :] + 1].mean()
    return [(f"in-plane (web) cantilever, {nx}x{ny} mesh: tip deflection [m]", fe, exact)]


def ss_plate_buckling(n=16):
    a = 1.0
    t = 0.01
    mat = Isotropic(E=70e9, nu=0.3, rho=2700.0)
    nodes, elems, ids = grid(n, n, lambda s, r: (s * a, r * a, 0.0))
    m = ShellModel(nodes, [(e, "p") for e in elems], {"p": isotropic_section(mat, t)},
                   {"p": np.array([1.0, 0, 0])}).build()
    edge = set(ids[0, :]) | set(ids[-1, :]) | set(ids[:, 0]) | set(ids[:, -1])
    fixed = {6 * k + 2 for k in edge}                     # w = 0 on all edges (simple support)
    fixed |= {6 * k + 0 for k in ids[0, :]}               # u = 0 on x = 0
    fixed.add(6 * ids[0, n // 2] + 1)                      # v = 0 at one point
    fixed |= {6 * k + 5 for k in range(len(nodes))}       # drilling rotations (flat plate)
    m.fixed = np.array(sorted(fixed))
    free = np.ones(m.ndof, bool)
    free[m.fixed] = False
    m.free = np.where(free)[0]
    Nx = 1.0  # N/m compression
    F = consistent_edge_load(ids[-1, :], nodes, 0, -Nx * a, m.ndof)
    u = m.solve_static(F)
    lam, _ = m.solve_buckling(u, k=3)
    D = mat.E * t**3 / (12 * (1 - mat.nu**2))
    exact = 4 * math.pi**2 * D / a**2
    return [("simply supported square plate: buckling Nx_cr [N/m]", lam[0] * Nx, exact)]


def box_torsion(nb=10, nh=6, nl=40):
    """Closed rectangular tube, clamped at one end, torque at the other (Bredt)."""
    bw, hh, L, t = 0.2, 0.1, 1.0, 0.002
    mat = Isotropic(E=70e9, nu=0.33, rho=2700.0)
    per = []  # perimeter points (x, z) counter-clockwise
    for i in range(nb):
        per.append((-bw / 2 + bw * i / nb, -hh / 2))
    for i in range(nh):
        per.append((bw / 2, -hh / 2 + hh * i / nh))
    for i in range(nb):
        per.append((bw / 2 - bw * i / nb, hh / 2))
    for i in range(nh):
        per.append((-bw / 2, hh / 2 - hh * i / nh))
    npp = len(per)
    nodes = np.array([(x, L * j / nl, z) for j in range(nl + 1) for (x, z) in per])
    elems = []
    for j in range(nl):
        for i in range(npp):
            a0, a1 = j * npp + i, j * npp + (i + 1) % npp
            elems.append([a0, a1, a1 + npp, a0 + npp])
    m = ShellModel(nodes, [(e, "s") for e in elems], {"s": isotropic_section(mat, t)},
                   {"s": np.array([0.0, 1.0, 0])}).build()
    m.fix_nodes(range(npp))
    T = 100.0
    Aenc = bw * hh
    q = T / (2 * Aenc)  # shear flow
    F = np.zeros(m.ndof)
    tip0 = nl * npp
    for i in range(npp):  # each perimeter segment carries q * length along its direction
        p0, p1 = np.array(per[i]), np.array(per[(i + 1) % npp])
        seg = p1 - p0
        f = q * seg  # force vector (x, z) along the segment, rotating about +y
        for nid in (tip0 + i, tip0 + (i + 1) % npp):
            F[6 * nid + 0] += -0.5 * f[0]
            F[6 * nid + 2] += -0.5 * f[1]
    u = m.solve_static(F)
    # twist from displacement of the four corners at the tip
    tip = nodes[tip0:tip0 + npp]
    ux = u[6 * np.arange(tip0, tip0 + npp)]
    uz = u[6 * np.arange(tip0, tip0 + npp) + 2]
    theta = np.mean((tip[:, 2] * ux - tip[:, 0] * uz) / (tip[:, 0]**2 + tip[:, 2]**2))
    J = 4 * Aenc**2 * t / (2 * (bw + hh))
    exact = T * L / (mat.G * J)
    return [("closed box, torque at tip: twist [rad] vs Bredt-Batho", abs(theta), exact)]


def scordelis_lo(n=12):
    R, L, phi0, t = 25.0, 50.0, math.radians(40.0), 0.25
    mat = Isotropic(E=4.32e8, nu=0.0, rho=1.0)
    nodes, elems, ids = grid(n, n, lambda s, r: (R * math.sin(s * phi0), r * L / 2, R * math.cos(s * phi0)))
    m = ShellModel(nodes, [(e, "r") for e in elems], {"r": isotropic_section(mat, t)},
                   {"r": np.array([0.0, 1.0, 0])}).build()
    fixed = set()
    for k in ids[:, -1]:          # diaphragm at y = L/2: ux = uz = 0
        fixed |= {6 * k + 0, 6 * k + 2}
    for k in ids[:, 0]:           # symmetry y = 0: uy = rotx = rotz = 0
        fixed |= {6 * k + 1, 6 * k + 3, 6 * k + 5}
    for k in ids[0, :]:           # symmetry x = 0 (crown): ux = roty = rotz = 0
        fixed |= {6 * k + 0, 6 * k + 4, 6 * k + 5}
    m.fixed = np.array(sorted(fixed))
    free = np.ones(m.ndof, bool)
    free[m.fixed] = False
    m.free = np.where(free)[0]
    F = np.zeros(m.ndof)
    for (conn, _), (dofs, T, data) in zip(m.elements, m.edata):
        area = data["area"]
        for nid in conn:
            F[6 * nid + 2] -= 90.0 * area / 4
    u = m.solve_static(F)
    fe = -u[6 * ids[-1, 0] + 2]
    return [(f"Scordelis-Lo roof ({n}x{n}): midside free-edge deflection", fe, 0.3024)]


def ply_orientation():
    """[0]8 and [90]8 plate strips must reproduce E1 and E2 bending stiffness."""
    ply = Ply(E1=181e9, E2=10.3e9, G12=7.17e9, nu12=0.28, G13=7.17e9, G23=3.5e9, rho=1600,
              Xt=1, Xc=1, Yt=1, Yc=1, S12=1, t=0.125e-3)
    out = []
    for ang, E in ((0, ply.E1), (90, ply.E2)):
        L, b = 0.2, 0.02
        sec = laminate_section(ply, [ang] * 8)
        nodes, elems, ids = grid(20, 2, lambda s, r: (s * L, r * b, 0.0))
        m = ShellModel(nodes, [(e, "p") for e in elems], {"p": sec}, {"p": np.array([1.0, 0, 0])}).build()
        m.fix_nodes(ids[0, :])
        F = consistent_edge_load(ids[-1, :], nodes, 2, 1.0, m.ndof)
        u = m.solve_static(F)
        h = 8 * ply.t
        # plate-strip: use D11^-1 compliance (free anticlastic curvature)
        d = np.linalg.inv(sec.D)
        exact = 1.0 * L**3 / (3 * b / d[0, 0])
        out.append((f"[{ang}]8 laminate strip: tip deflection [m] (E ~ {E / 1e9:.0f} GPa)",
                    u[6 * ids[-1, :] + 2].mean(), exact))
    return out


def lifting_line_elliptic():
    """Elliptic planform: lifting line must give e = 1 and CL_a = a0 / (1 + a0 / (pi AR))."""
    from wingfe.aero import lifting_line
    s, c0 = 1.2, 0.25
    _, _, CLa, e = lifting_line(lambda y: c0 * math.sqrt(max(0.0, 1 - (y / s)**2)), s)
    AR = (2 * s)**2 / (math.pi * s * c0 / 2)
    a0 = 2 * math.pi
    return [("lifting line, elliptic wing: CL_alpha [1/rad]", CLa, a0 / (1 + a0 / (math.pi * AR))),
            ("lifting line, elliptic wing: span efficiency e", e, 1.0)]


def theodorsen_function():
    """Hankel-function form used in wingfe vs the independent modified-Bessel form K1/(K0+K1)."""
    from scipy.special import kv
    from wingfe.flutter import theodorsen
    rows = []
    for k in (0.1, 0.5, 1.0):
        ref = kv(1, 1j * k) / (kv(0, 1j * k) + kv(1, 1j * k))
        c = theodorsen(k)
        rows.append((f"Theodorsen C(k = {k}): real part", c.real, ref.real))
        rows.append((f"Theodorsen C(k = {k}): imaginary part", c.imag, ref.imag))
    return rows


def flutter_trend():
    """Typical section: the flutter speed must rise as the CG moves forward (mass balancing)."""
    from wingfe.flutter import flutter_speed
    speeds = [flutter_speed(20, 0.4, 6 / 25, -0.2, xa)[0] for xa in (0.2, 0.1, 0.0)]
    return [("typical section, x_a = 0.1: flutter speed U_F/(b w_alpha) (no reference)", speeds[1], 0.0),
            ("typical section: flutter speed rises as the CG moves forward (1 = yes)",
             float(speeds[0] < speeds[1] < speeds[2]), 1.0)]


def wing_rigid_body():
    """Assembled project wing: rigid-body motions must produce no forces."""
    from wingfe.wing import Layout, Mesh, build_wing
    mat = Isotropic(E=70e9, nu=0.33, rho=2700.0)
    secs = {"skin": isotropic_section(mat, 0.008), "web": isotropic_section(mat, 0.005),
            "rib": isotropic_section(mat, 0.002)}
    m, _ = build_wing(Mesh(), Layout(), secs)
    m.build(want_mass=False)
    return [("assembled wing: rigid-body check (max |K r| / max K_ii)", m.rigid_body_check().max(), 0.0)]


def main():
    rows = []
    for fn in (cantilever_plate, inplane_beam, lambda: inplane_beam(10, 2), ss_plate_buckling,
               box_torsion, scordelis_lo, lambda: scordelis_lo(20), ply_orientation,
               lifting_line_elliptic, theodorsen_function, flutter_trend, wing_rigid_body):
        rows += fn()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "benchmarks.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case", "wingfe", "reference", "error_%"])
        for name, fe, ref in rows:
            err = (fe - ref) / ref * 100 if ref else float("nan")
            w.writerow([name, f"{fe:.6g}", f"{ref:.6g}", f"{err:.2f}"])
            print(f"{name:<75} {fe:12.5g} {ref:12.5g} {err:8.2f} %" if ref else f"{name:<75} {fe:12.3g}")


if __name__ == "__main__":
    main()
