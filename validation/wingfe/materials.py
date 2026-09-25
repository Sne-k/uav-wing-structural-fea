"""Section properties: isotropic plates and classical lamination theory (CLT).

Ply angles are measured from the element local x axis, which the mesh sets to
the spanwise direction for skins and webs.
"""
from dataclasses import dataclass

import numpy as np

from .shell import Section

SHEAR_CORR = 5.0 / 6.0


@dataclass
class Isotropic:
    E: float
    nu: float
    rho: float
    Sy: float = np.nan      # yield strength [Pa]
    Su: float = np.nan      # ultimate strength [Pa]
    name: str = ""

    @property
    def G(self):
        return self.E / (2 * (1 + self.nu))


@dataclass
class Ply:
    """Unidirectional ply (material axes 1 = fibre, 2 = transverse)."""
    E1: float
    E2: float
    G12: float
    nu12: float
    G13: float
    G23: float
    rho: float
    Xt: float
    Xc: float
    Yt: float
    Yc: float
    S12: float
    t: float
    name: str = ""

    def Q(self):
        nu21 = self.nu12 * self.E2 / self.E1
        d = 1 - self.nu12 * nu21
        return np.array([[self.E1 / d, self.nu12 * self.E2 / d, 0],
                         [self.nu12 * self.E2 / d, self.E2 / d, 0],
                         [0, 0, self.G12]])


# T300/5208 carbon/epoxy (Kaw, Mechanics of Composite Materials, 2nd ed., Table 2.1)
T300_5208 = Ply(E1=181e9, E2=10.3e9, G12=7.17e9, nu12=0.28, G13=7.17e9, G23=3.5e9,
                rho=1600.0, Xt=1500e6, Xc=1500e6, Yt=40e6, Yc=246e6, S12=68e6,
                t=0.125e-3, name="T300/5208")

AL6061_T6 = Isotropic(E=70e9, nu=0.33, rho=2700.0, Sy=276e6, Su=310e6, name="Al 6061-T6")


def isotropic_section(mat, t, name=""):
    D0 = mat.E / (1 - mat.nu**2)
    Qm = D0 * np.array([[1, mat.nu, 0], [mat.nu, 1, 0], [0, 0, (1 - mat.nu) / 2]])
    A = Qm * t
    D = Qm * t**3 / 12
    As = SHEAR_CORR * mat.G * t * np.eye(2)
    return Section(A, np.zeros((3, 3)), D, As, mat.rho * t, mat.rho * t**3 / 12,
                   name=name or f"{mat.name} {t * 1e3:g} mm")


def _T_sigma(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c * c, s * s, 2 * s * c],
                     [s * s, c * c, -2 * s * c],
                     [-s * c, s * c, c * c - s * s]])


def qbar(ply, theta_deg):
    th = np.radians(theta_deg)
    T = _T_sigma(th)
    Rr = np.diag([1, 1, 2])
    return np.linalg.inv(T) @ ply.Q() @ Rr @ T @ np.linalg.inv(Rr)


def laminate_section(ply, angles_deg, name=""):
    """CLT section for a stack of identical plies listed bottom to top."""
    n = len(angles_deg)
    h = n * ply.t
    z = np.linspace(-h / 2, h / 2, n + 1)
    A = np.zeros((3, 3))
    B = np.zeros((3, 3))
    D = np.zeros((3, 3))
    As = np.zeros((2, 2))
    plies = []
    for k, ang in enumerate(angles_deg):
        Qb = qbar(ply, ang)
        A += Qb * (z[k + 1] - z[k])
        B += Qb * (z[k + 1]**2 - z[k]**2) / 2
        D += Qb * (z[k + 1]**3 - z[k]**3) / 3
        c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
        Qs = np.array([[ply.G13 * c * c + ply.G23 * s * s, (ply.G13 - ply.G23) * c * s],
                       [(ply.G13 - ply.G23) * c * s, ply.G13 * s * s + ply.G23 * c * c]])
        As += SHEAR_CORR * Qs * (z[k + 1] - z[k])
        plies.append((ang, z[k], z[k + 1], ply))
    rho_t = ply.rho * h
    rho_i = ply.rho * h**3 / 12
    return Section(A, B, D, As, rho_t, rho_i,
                   name=name or f"{ply.name} [{'/'.join(str(a) for a in angles_deg)}]",
                   plies=plies)


def engineering_constants(sec):
    """Effective in-plane moduli of a laminate from the A matrix."""
    h = sum(p[2] - p[1] for p in sec.plies) if sec.plies else None
    a = np.linalg.inv(sec.A)
    return dict(Ex=1 / (a[0, 0] * h), Ey=1 / (a[1, 1] * h), Gxy=1 / (a[2, 2] * h),
                nuxy=-a[0, 1] / a[0, 0], h=h)


def tsai_wu_index(ply, s1, s2, t12):
    F1 = 1 / ply.Xt - 1 / ply.Xc
    F2 = 1 / ply.Yt - 1 / ply.Yc
    F11 = 1 / (ply.Xt * ply.Xc)
    F22 = 1 / (ply.Yt * ply.Yc)
    F66 = 1 / ply.S12**2
    F12 = -0.5 * np.sqrt(F11 * F22)
    return F1 * s1 + F2 * s2 + F11 * s1**2 + F22 * s2**2 + F66 * t12**2 + 2 * F12 * s1 * s2


def tsai_wu_strength_ratio(ply, s1, s2, t12):
    """Load multiplier R at which the Tsai-Wu index reaches 1 (R > 1 is safe)."""
    F1 = 1 / ply.Xt - 1 / ply.Xc
    F2 = 1 / ply.Yt - 1 / ply.Yc
    F11 = 1 / (ply.Xt * ply.Xc)
    F22 = 1 / (ply.Yt * ply.Yc)
    F66 = 1 / ply.S12**2
    F12 = -0.5 * np.sqrt(F11 * F22)
    a = F11 * s1**2 + F22 * s2**2 + F66 * t12**2 + 2 * F12 * s1 * s2
    b = F1 * s1 + F2 * s2
    if a <= 0:
        return np.inf if b <= 0 else 1 / b
    return (-b + np.sqrt(b * b + 4 * a)) / (2 * a)


def ply_stresses(sec, eps_kappa):
    """Ply stresses in material axes at the top and bottom of every ply.

    eps_kappa: (6,) mid-surface strains and curvatures in the section (local) axes.
    Returns a list of (angle, z, s1, s2, t12, ply).
    """
    e0, k = eps_kappa[:3], eps_kappa[3:]
    out = []
    for ang, zb, zt, ply in sec.plies:
        Qb = qbar(ply, ang)
        T = _T_sigma(np.radians(ang))
        for z in (zb, zt):
            sig_xy = Qb @ (e0 + z * k)
            s1, s2, t12 = T @ sig_xy
            out.append((ang, z, s1, s2, t12, ply))
    return out
