"""Prandtl lifting-line theory (Glauert's monoplane equation) for a straight wing.

Used to turn a total design lift into a spanwise lift distribution, in place of
the XFLR5 output that could not be exported during the project.
"""
import numpy as np

RHO_SL = 1.225
_trapz = getattr(np, "trapezoid", None) or np.trapz   # renamed in NumPy 2.0


def lifting_line(chord_fn, semi_span, a0=2 * np.pi, n_terms=40):
    """Symmetric loading on an untwisted wing.

    chord_fn(y) for 0 <= y <= semi_span. Returns (y, l_shape, CL_alpha, e) where
    l_shape is the spanwise lift per unit span normalised so that its integral
    over the semi-span is 1 (shape is independent of alpha for an untwisted wing).
    """
    b = 2 * semi_span
    # collocation over the half span, theta in (0, pi/2]
    th = np.pi / 2 * (np.arange(1, n_terms + 1)) / n_terms      # excludes the tip (theta = 0)
    n = 2 * np.arange(n_terms) + 1                                # odd terms only
    y = semi_span * np.cos(th)
    c = np.array([chord_fn(abs(v)) for v in y])
    Amat = (np.sin(np.outer(th, n)) * (4 * b / (a0 * c))[:, None]
            + np.sin(np.outer(th, n)) * n[None, :] / np.sin(th)[:, None])
    An = np.linalg.solve(Amat, np.ones(n_terms))                 # per unit (alpha - alpha_0)
    S = 2 * _trapz([chord_fn(v) for v in np.linspace(0, semi_span, 401)],
                         np.linspace(0, semi_span, 401))
    AR = b * b / S
    CLa = np.pi * AR * An[0]
    delta = np.sum(n[1:] * (An[1:] / An[0])**2)
    ys = np.linspace(0, semi_span, 241)
    ths = np.arccos(np.clip(ys / semi_span, -1, 1))
    gam = np.sin(np.outer(ths, n)) @ An                          # ~ circulation
    shape = gam / _trapz(gam, ys)
    return ys, shape, CLa, 1 / (1 + delta)


def schrenk(chord_fn, semi_span, n=241):
    ys = np.linspace(0, semi_span, n)
    S_half = _trapz([chord_fn(v) for v in ys], ys)
    ell = 4 * S_half / (np.pi * semi_span) * np.sqrt(np.clip(1 - (ys / semi_span)**2, 0, None))
    trap = np.array([chord_fn(v) for v in ys])
    shape = 0.5 * (ell + trap)
    return ys, shape / _trapz(shape, ys)
