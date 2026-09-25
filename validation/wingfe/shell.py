"""Flat-facet 4-node shell element (6 DOF per node).

* membrane:  bilinear Q4 plus Wilson/Taylor incompatible modes (QM6), so thin
             webs in in-plane bending do not lock
* bending:   Reissner-Mindlin plate with MITC4 assumed transverse shear
             (Bathe & Dvorkin 1985), so thin plates do not shear-lock
* drilling:  small penalty tying the drilling rotation to the membrane rotation
             (theta_z - omega), which leaves rigid-body rotations energy-free
* laminates: the section is given as A, B, D (3x3) and As (2x2) matrices, so an
             isotropic plate and a composite layup use the same code path

A degenerate quad (two nodes repeated) is accepted and behaves as a triangle.

Local DOF order per node: u, v, w, theta_x, theta_y, theta_z.
Sign convention: u = z*theta_y, v = -z*theta_x (right-hand rotations).
"""
import numpy as np

_G = 1.0 / np.sqrt(3.0)
GAUSS = [(-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G)]
XI_N = np.array([-1.0, 1.0, 1.0, -1.0])
ETA_N = np.array([-1.0, -1.0, 1.0, 1.0])


def shape(xi, eta):
    n = 0.25 * (1 + xi * XI_N) * (1 + eta * ETA_N)
    dxi = 0.25 * XI_N * (1 + eta * ETA_N)
    deta = 0.25 * ETA_N * (1 + xi * XI_N)
    return n, dxi, deta


def local_frame(X, ref_dir):
    """Rotation matrix (rows = local axes) and in-plane nodal coordinates.

    Local x is the projection of ref_dir onto the element plane (this is also
    the material 0-degree direction), local z is the facet normal.
    """
    e3 = np.cross(X[2] - X[0], X[3] - X[1])
    e3 /= np.linalg.norm(e3)
    r = ref_dir - (ref_dir @ e3) * e3
    if np.linalg.norm(r) < 0.2 * np.linalg.norm(ref_dir):
        r = X[1] - X[0]
        r = r - (r @ e3) * e3
    e1 = r / np.linalg.norm(r)
    e2 = np.cross(e3, e1)
    R = np.vstack([e1, e2, e3])
    c = X.mean(axis=0)
    xy = (X - c) @ R[:2].T
    warp = (X - c) @ e3
    return R, xy, warp


def _jac(xy, dxi, deta):
    J = np.array([[dxi @ xy[:, 0], dxi @ xy[:, 1]],
                  [deta @ xy[:, 0], deta @ xy[:, 1]]])
    return J, np.linalg.det(J)


class Section:
    """Shell section: stiffness matrices, mass per area and rotary inertia per area."""

    def __init__(self, A, B, D, As, rho_t, rho_i, name="", plies=None):
        self.A, self.B, self.D, self.As = map(np.asarray, (A, B, D, As))
        self.rho_t, self.rho_i = rho_t, rho_i
        self.name = name
        self.plies = plies  # for laminates: list of (angle_deg, z_bot, z_top, ply)

    @property
    def g_drill(self):
        return 1e-3 * self.A[2, 2]


def _membrane_rows(dNx, dNy):
    Bm = np.zeros((3, 24))
    for i in range(4):
        Bm[0, 6 * i] = dNx[i]
        Bm[1, 6 * i + 1] = dNy[i]
        Bm[2, 6 * i] = dNy[i]
        Bm[2, 6 * i + 1] = dNx[i]
    return Bm


def _bending_rows(dNx, dNy):
    Bb = np.zeros((3, 24))
    for i in range(4):
        Bb[0, 6 * i + 4] = dNx[i]
        Bb[1, 6 * i + 3] = -dNy[i]
        Bb[2, 6 * i + 4] = dNy[i]
        Bb[2, 6 * i + 3] = -dNx[i]
    return Bb


def _mitc4_shear(xy, xi, eta):
    """Assumed transverse shear strain rows (2x24) at (xi, eta)."""
    def covariant(tp_xi, tp_eta, which):
        n, dxi, deta = shape(tp_xi, tp_eta)
        d = dxi if which == "xi" else deta
        dx = d @ xy[:, 0]
        dy = d @ xy[:, 1]
        row = np.zeros(24)
        for i in range(4):
            row[6 * i + 2] = d[i]
            row[6 * i + 4] = n[i] * dx
            row[6 * i + 3] = -n[i] * dy
        return row

    gA = covariant(0.0, 1.0, "xi")
    gC = covariant(0.0, -1.0, "xi")
    gB = covariant(-1.0, 0.0, "eta")
    gD = covariant(1.0, 0.0, "eta")
    g_xi = 0.5 * (1 + eta) * gA + 0.5 * (1 - eta) * gC
    g_eta = 0.5 * (1 + xi) * gD + 0.5 * (1 - xi) * gB
    _, dxi, deta = shape(xi, eta)
    J, _ = _jac(xy, dxi, deta)
    return np.linalg.solve(J, np.vstack([g_xi, g_eta]))


def _incompatible_rows(xy, xi, eta, J0inv, detJ0, detJ):
    """Membrane strain rows for the 4 incompatible-mode amplitudes (3x4)."""
    dP = np.array([[-2 * xi, 0.0], [0.0, -2 * eta]])  # rows: P1, P2 ; cols: d/dxi, d/deta
    f = detJ0 / detJ
    dPx = np.zeros(2)
    dPy = np.zeros(2)
    for k in range(2):
        g = f * J0inv @ dP[k]
        dPx[k], dPy[k] = g
    Ba = np.zeros((3, 4))
    Ba[0, 0], Ba[0, 1] = dPx[0], dPx[1]           # eps_x  <- a1, a2 (u)
    Ba[1, 2], Ba[1, 3] = dPy[0], dPy[1]           # eps_y  <- a3, a4 (v)
    Ba[2, 0], Ba[2, 1] = dPy[0], dPy[1]           # gamma_xy
    Ba[2, 2], Ba[2, 3] = dPx[0], dPx[1]
    return Ba


def element_matrices(X, sec, ref_dir, want_mass=True):
    """Return (K_global 24x24, M_global 24x24, T, data) for one element.

    data keeps what is needed for stress recovery and geometric stiffness.
    """
    R, xy, warp = local_frame(X, ref_dir)
    C = np.block([[sec.A, sec.B], [sec.B, sec.D]])
    _, dxi0, deta0 = shape(0.0, 0.0)
    J0, detJ0 = _jac(xy, dxi0, deta0)
    J0inv = np.linalg.inv(J0)

    Kuu = np.zeros((24, 24))
    Kua = np.zeros((24, 4))
    Kaa = np.zeros((4, 4))
    M = np.zeros((24, 24))
    gp_data = []
    area = 0.0
    for xi, eta in GAUSS:
        n, dxi, deta = shape(xi, eta)
        J, detJ = _jac(xy, dxi, deta)
        g = np.linalg.solve(J, np.vstack([dxi, deta]))
        dNx, dNy = g
        Bm = _membrane_rows(dNx, dNy)
        Bb = _bending_rows(dNx, dNy)
        Ba = _incompatible_rows(xy, xi, eta, J0inv, detJ0, detJ)
        Bs = _mitc4_shear(xy, xi, eta)
        E6 = np.vstack([Bm, Bb])
        Ea = np.vstack([Ba, np.zeros((3, 4))])
        w = detJ  # Gauss weights are 1 for 2x2
        area += w
        Kuu += w * (E6.T @ C @ E6 + Bs.T @ sec.As @ Bs)
        Kua += w * (E6.T @ C @ Ea)
        Kaa += w * (Ea.T @ C @ Ea)
        # drilling penalty: (theta_z - omega)^2, omega = (v,x - u,y)/2
        dr = np.zeros(24)
        for i in range(4):
            dr[6 * i + 5] = n[i]
            dr[6 * i] += 0.5 * dNy[i]
            dr[6 * i + 1] -= 0.5 * dNx[i]
        Kuu += w * sec.g_drill * np.outer(dr, dr)
        if want_mass:
            for i in range(4):
                for j in range(4):
                    mij = w * n[i] * n[j]
                    for d in range(3):
                        M[6 * i + d, 6 * j + d] += sec.rho_t * mij
                    for d in range(3, 6):
                        M[6 * i + d, 6 * j + d] += sec.rho_i * mij
        gp_data.append((Bm, Bb, Ba, dNx, dNy, w))

    Kaa_inv = np.linalg.inv(Kaa)
    K = Kuu - Kua @ Kaa_inv @ Kua.T
    T = np.zeros((24, 24))
    for b in range(8):
        T[3 * b:3 * b + 3, 3 * b:3 * b + 3] = R
    data = dict(R=R, xy=xy, warp=np.abs(warp).max(), area=area, gp=gp_data,
                recover=Kaa_inv @ Kua.T, C=C)
    return T.T @ K @ T, T.T @ M @ T, T, data


def generalized_strains(data, T, u_glob):
    """Mid-surface strains and curvatures at the 4 Gauss points (local axes)."""
    ul = T @ u_glob
    a = -data["recover"] @ ul
    out = []
    for Bm, Bb, Ba, *_ in data["gp"]:
        out.append(np.concatenate([Bm @ ul + Ba @ a, Bb @ ul]))
    return np.array(out)  # (4, 6): eps_x, eps_y, gam_xy, k_x, k_y, k_xy


def geometric_stiffness(data, T, N_gp):
    """Geometric (initial-stress) stiffness in global axes.

    N_gp: (4, 3) membrane resultants Nx, Ny, Nxy at the Gauss points.
    The in-plane stress acts on the gradients of all three translations.
    """
    Kg = np.zeros((24, 24))
    for (Bm, Bb, Ba, dNx, dNy, w), (nx, ny, nxy) in zip(data["gp"], N_gp):
        S = np.array([[nx, nxy], [nxy, ny]])
        G = np.vstack([dNx, dNy])
        k4 = w * G.T @ S @ G
        for d in range(3):
            idx = [6 * i + d for i in range(4)]
            Kg[np.ix_(idx, idx)] += k4
    return T.T @ Kg @ T
