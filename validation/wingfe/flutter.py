"""Typical-section (2-DOF bending-torsion) flutter with Theodorsen aerodynamics, V-g method.

Notation (Hodges & Pierce, Introduction to Structural Dynamics and Aeroelasticity):
    b      semi-chord
    a      elastic-axis position aft of mid-chord, in semi-chords
    x_a    CG position aft of the elastic axis, in semi-chords
    r2     (radius of gyration about the elastic axis / b)^2
    mu     m / (pi rho b^2)
    sigma  omega_h / omega_alpha
h is positive down and alpha is positive nose up.
"""
import numpy as np
from scipy.special import hankel2


def theodorsen(k):
    H0, H1 = hankel2(0, k), hankel2(1, k)
    return H1 / (H1 + 1j * H0)


def vg_analysis(mu, sigma, r2, a, x_a, ks=None):
    """Return list of (k, [(V/(b w_a), w/w_a, g) for both branches])."""
    if ks is None:
        ks = np.concatenate([np.linspace(3.0, 0.2, 400), np.geomspace(0.2, 0.005, 400)[1:]])
    out = []
    for k in ks:
        C = theodorsen(k)
        lh = 1 - 2j * C / k
        la = -(a + 1j / k + 2 * C / k**2 + 2j * C * (0.5 - a) / k)
        mh = -a + 2j * (a + 0.5) * C / k
        ma = (1 / 8 + a * a - 1j * (0.5 - a) / k + 2 * (a + 0.5) * C / k**2
              + 2j * (a + 0.5) * (0.5 - a) * C / k)
        A = np.array([[mu + lh, mu * x_a + la], [mu * x_a + mh, mu * r2 + ma]])
        B = np.array([[mu * sigma**2, 0], [0, mu * r2]])
        Z = np.linalg.eigvals(np.linalg.solve(B, A))      # Z = (w_a / w)^2 (1 + i g)
        branches = []
        for z in Z:
            if z.real <= 0:
                branches.append((np.nan, np.nan, np.nan))
                continue
            w_ratio = 1 / np.sqrt(z.real)
            g = z.imag / z.real
            branches.append((w_ratio / k, w_ratio, g))
        out.append((k, sorted(branches, key=lambda t: (np.nan_to_num(t[1], nan=9e9)))))
    return out


def flutter_speed(mu, sigma, r2, a, x_a):
    """Lowest reduced speed V/(b w_alpha) where a branch's damping g crosses zero upward."""
    res = vg_analysis(mu, sigma, r2, a, x_a)
    best = None
    for br in range(2):
        prev = None
        for k, bs in res:
            V, w, g = bs[br]
            if np.isnan(g):
                prev = None
                continue
            if prev is not None and prev[2] < 0 <= g:
                f = -prev[2] / (g - prev[2])
                Vf = prev[0] + f * (V - prev[0])
                wf = prev[1] + f * (w - prev[1])
                if best is None or Vf < best[0]:
                    best = (Vf, wf)
            prev = (V, w, g)
    return best  # (V_F / (b w_a), w_F / w_a) or None


def divergence_speed(mu, r2, a):
    """V_D / (b w_alpha) for the typical section (None if the EA is at/ahead of the quarter chord)."""
    if 1 + 2 * a <= 0:
        return None
    return np.sqrt(mu * r2 / (1 + 2 * a))
