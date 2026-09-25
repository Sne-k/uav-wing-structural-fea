"""Parametric shell model of the project wing: skin, two spar webs and ribs.

Geometry (read back from the Ansys mesh, see docs/model_audit.md):
    NACA 2412, semi-span 1.2 m, chord 0.25 m (root) -> 0.20 m (tip), straight
    unswept leading edge at x = 0, no twist or dihedral
    spar 1: straight web at x = 0.0625 m (25 % of the root chord)
    spar 2: web on the 70 % chord line (tapered)
    10 ribs at y = k * 1.2 / 11 (k = 1..10), plus an optional tip rib

All parts share nodes along their junctions (no contact), and the whole root
section (skin edge + both webs) is clamped.
Axes: x chordwise (aft), y spanwise, z up. Units: SI.
"""
from dataclasses import dataclass, field

import numpy as np

from .model import ShellModel

SPAN, C_ROOT, C_TIP = 1.2, 0.25, 0.20
SPAR1_X, SPAR2_FRAC = 0.0625, 0.70


def chord(y):
    return C_ROOT + (C_TIP - C_ROOT) * y / SPAN


def naca2412(xi, closed_te=True):
    """Upper and lower surface z/c at chord fraction xi (vertical-thickness form)."""
    m, p, t = 0.02, 0.4, 0.12
    a4 = -0.1036 if closed_te else -0.1015
    yt = 5 * t * (0.2969 * np.sqrt(xi) - 0.1260 * xi - 0.3516 * xi**2 + 0.2843 * xi**3 + a4 * xi**4)
    yc = np.where(xi < p, m / p**2 * (2 * p * xi - xi**2), m / (1 - p)**2 * ((1 - 2 * p) + 2 * p * xi - xi**2))
    return yc + yt, yc - yt


@dataclass
class Mesh:
    n_le: int = 6      # elements per surface, leading edge -> spar 1
    n_mid: int = 10    # spar 1 -> spar 2
    n_te: int = 5      # spar 2 -> trailing edge
    n_web: int = 4     # through the depth of webs and ribs
    n_bay: int = 4     # spanwise elements per rib bay

    def scaled(self, f):
        return Mesh(*(max(1, int(round(v * f))) for v in
                      (self.n_le, self.n_mid, self.n_te, self.n_web, self.n_bay)))


@dataclass
class Layout:
    tip_rib: bool = True
    rib_stations: list = field(default_factory=lambda: [k * SPAN / 11 for k in range(1, 11)])


class _Registry:
    def __init__(self):
        self.key = {}
        self.xyz = []

    def add(self, p):
        k = tuple(np.round(p, 9))
        if k not in self.key:
            self.key[k] = len(self.xyz)
            self.xyz.append(np.array(p, float))
        return self.key[k]


def _chord_fracs(y, mesh):
    xi1 = SPAR1_X / chord(y)
    le = xi1 * (1 - np.cos(0.5 * np.pi * np.arange(mesh.n_le + 1) / mesh.n_le))
    mid = np.linspace(xi1, SPAR2_FRAC, mesh.n_mid + 1)
    te = np.linspace(SPAR2_FRAC, 1.0, mesh.n_te + 1)
    return np.concatenate([le, mid[1:], te[1:]]), mesh.n_le, mesh.n_le + mesh.n_mid


def build_wing(mesh=Mesh(), layout=Layout(), sections=None):
    """Return (ShellModel, info) — sections: {'skin', 'web', 'rib'} -> Section."""
    reg = _Registry()
    y_st = [0.0]
    ribs = sorted(layout.rib_stations)
    for y0, y1 in zip([0.0] + ribs, ribs + [SPAN]):
        y_st += list(np.linspace(y0, y1, mesh.n_bay + 1)[1:])
    y_st = np.array(y_st)

    def surf(xi, y, level):
        """Point at chord fraction xi, height fraction level (0 lower skin .. 1 upper skin)."""
        c = chord(y)
        zu, zl = naca2412(xi)
        return np.array([xi * c, y, (zl + level * (zu - zl)) * c])

    elems = []
    upper_ids, web_ids = [], {1: [], 2: []}
    grid_u, grid_l, webs = [], [], {1: [], 2: []}
    for y in y_st:
        xs, i1, i2 = _chord_fracs(y, mesh)
        grid_u.append([reg.add(surf(x, y, 1.0)) for x in xs])
        grid_l.append([reg.add(surf(x, y, 0.0)) for x in xs])
        for s, i in ((1, i1), (2, i2)):
            webs[s].append([reg.add(surf(xs[i], y, k / mesh.n_web)) for k in range(mesh.n_web + 1)])

    nx = len(grid_u[0])
    for j in range(len(y_st) - 1):
        for k in range(nx - 1):
            elems.append(([grid_u[j][k + 1], grid_u[j][k], grid_u[j + 1][k], grid_u[j + 1][k + 1]], "skin_upper"))
            elems.append(([grid_l[j][k], grid_l[j][k + 1], grid_l[j + 1][k + 1], grid_l[j + 1][k]], "skin_lower"))
        for s in (1, 2):
            for k in range(mesh.n_web):
                elems.append(([webs[s][j][k], webs[s][j + 1][k], webs[s][j + 1][k + 1], webs[s][j][k + 1]], f"web{s}"))

    rib_y = list(ribs) + ([SPAN] if layout.tip_rib else [])
    for y in rib_y:
        xs, _, _ = _chord_fracs(y, mesh)
        col = [[reg.add(surf(x, y, k / mesh.n_web)) for k in range(mesh.n_web + 1)] for x in xs]
        for i in range(len(xs) - 1):
            for k in range(mesh.n_web):
                conn = [col[i][k], col[i + 1][k], col[i + 1][k + 1], col[i][k + 1]]
                if len(set(conn)) >= 3:
                    elems.append((conn, "rib"))

    nodes = np.array(reg.xyz)
    sec = {"skin_upper": sections["skin"], "skin_lower": sections["skin"],
           "web1": sections["web"], "web2": sections["web"], "rib": sections["rib"]}
    ref = {"skin_upper": np.array([0.0, 1, 0]), "skin_lower": np.array([0.0, 1, 0]),
           "web1": np.array([0.0, 1, 0]), "web2": np.array([0.0, 1, 0]), "rib": np.array([1.0, 0, 0])}
    model = ShellModel(nodes, elems, sec, ref)
    root = np.where(np.abs(nodes[:, 1]) < 1e-9)[0]
    tip = np.where(np.abs(nodes[:, 1] - SPAN) < 1e-9)[0]
    info = dict(y_stations=y_st, root=root, tip=tip, grid_u=grid_u, grid_l=grid_l, webs=webs)
    return model, info


# ---------------------------------------------------------------- loads
def load_component_z_all_skin(model, p):
    """Ansys-style load: pressure p [Pa] along +z on the full area of every skin facet."""
    F = np.zeros(model.ndof)
    for (conn, grp), (dofs, T, data) in zip(model.elements, model.edata):
        if not grp.startswith("skin"):
            continue
        for (xi, eta), (*_, w) in zip(_GP, data["gp"]):
            n = _N(xi, eta)
            for a, nid in enumerate(conn):
                F[6 * nid + 2] += p * w * n[a]
    return F


def load_upper_skin(model, dp_func):
    """Vertical aerodynamic load on the upper skin: dp_func(x, y) [Pa per planform area]."""
    F = np.zeros(model.ndof)
    for (conn, grp), (dofs, T, data) in zip(model.elements, model.edata):
        if grp != "skin_upper":
            continue
        X = model.nodes[conn]
        nz = abs(data["R"][2, 2])  # projection factor of the facet onto the xy plane
        for (xi, eta), (*_, w) in zip(_GP, data["gp"]):
            n = _N(xi, eta)
            x, y, _ = n @ X
            q = dp_func(x, y) * w * nz
            for a, nid in enumerate(conn):
                F[6 * nid + 2] += q * n[a]
    return F


_g = 1 / np.sqrt(3)
_GP = [(-_g, -_g), (_g, -_g), (_g, _g), (-_g, _g)]


def _N(xi, eta):
    return 0.25 * (1 + xi * np.array([-1, 1, 1, -1])) * (1 + eta * np.array([-1, -1, 1, 1]))


# ---------------------------------------------------------------- post-processing
def isotropic_stresses(model, u, groups=None):
    """Von Mises stress at top/bottom surface of every Gauss point.

    Returns arrays (group, y, x, vm_max, sig_span_membrane) per Gauss point.
    """
    out = []
    for (conn, grp), (dofs, T, data), e in zip(model.elements, model.edata, model.element_results(u)):
        if groups and grp not in groups:
            continue
        sec = model.sections[grp]
        t = np.sqrt(12 * sec.D[0, 0] / sec.A[0, 0])
        X = model.nodes[conn]
        for (xi, eta), ek in zip(_GP, e):
            N = sec.A @ ek[:3] + sec.B @ ek[3:]
            M = sec.B @ ek[:3] + sec.D @ ek[3:]
            vm = 0.0
            for s in (1, -1):
                sx, sy, txy = N / t + s * 6 * M / t**2
                vm = max(vm, np.sqrt(sx * sx - sx * sy + sy * sy + 3 * txy * txy))
            x, y, z = _N(xi, eta) @ X
            out.append((grp, x, y, z, vm, N[0] / t))
    g = np.array([o[0] for o in out])
    a = np.array([o[1:] for o in out], float)
    return g, a  # a columns: x, y, z, von Mises, local-x membrane stress (spanwise for skin/webs)


def laminate_strength(model, u, groups=("skin_upper", "skin_lower")):
    """Minimum Tsai-Wu strength ratio over all plies and Gauss points of laminate groups."""
    from .materials import ply_stresses, tsai_wu_strength_ratio
    worst = (np.inf, None)
    for (conn, grp), e in zip(model.elements, model.element_results(u)):
        if grp not in groups or model.sections[grp].plies is None:
            continue
        X = model.nodes[conn]
        for (xi, eta), ek in zip(_GP, e):
            for ang, z, s1, s2, t12, ply in ply_stresses(model.sections[grp], ek):
                r = tsai_wu_strength_ratio(ply, s1, s2, t12)
                if r < worst[0]:
                    worst = (r, dict(group=grp, angle=ang, z=z, xyz=_N(xi, eta) @ X,
                                     s1=s1, s2=s2, t12=t12))
    return worst
