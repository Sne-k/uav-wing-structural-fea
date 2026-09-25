"""Hand-calculation check of the wing FEA using Euler-Bernoulli beam theory.

The wing is idealised as a tapered cantilever, clamped at the root:
  * NACA 2412 section, chord 0.25 m (root) to 0.20 m (tip), semi-span 1.2 m
  * thin-walled aluminium skin on the outer mould line (same idealisation as the
    SHELL181 skin in the FE model) plus the two rectangular spars
  * 800 Pa acting vertically on the whole skin area (upper + lower surfaces),
    which is how the pressure is applied in the Ansys model

It prints the root bending moment, root bending stress, tip deflection and a
Rayleigh estimate of the first flapwise frequency. These are the values the
FE model should approach if the root were properly clamped.

Usage:
    python tools/beam_check.py            # 8 mm skin (as modelled)
    python tools/beam_check.py 1.0        # any other skin thickness in mm
"""
import math
import sys

E = 70e9            # Pa, aluminium (the "CFRP" case in the model uses the same E)
SPAN = 1.2          # m
ROOT_CHORD, TIP_CHORD = 0.25, 0.20
PRESSURE = 800.0    # Pa
FE_MASS = {"aluminium": 12.419, "CFRP (isotropic)": 7.3593}  # total mass reported by the solver [kg]


def naca2412(x, m=0.02, p=0.4, t=0.12):
    yt = 5 * t * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)
    yc = m / p**2 * (2 * p * x - x * x) if x < p else m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x * x)
    return yc + yt, yc - yt


def section(chord, skin):
    """Return (I about the chordwise neutral axis, wetted perimeter, max fibre distance)."""
    n = 400
    xs = [0.5 * (1 - math.cos(math.pi * k / n)) for k in range(n + 1)]
    pieces = []  # (area, centroid z, own I)
    perimeter = 0.0
    for surf in (0, 1):
        pts = [(x * chord, naca2412(x)[surf] * chord) for x in xs]
        for a, b in zip(pts[:-1], pts[1:]):
            ds = math.hypot(b[0] - a[0], b[1] - a[1])
            perimeter += ds
            pieces.append((skin * ds, 0.5 * (a[1] + b[1]), 0.0))
    for x_c, h_root in ((0.25, 0.020), (0.70, 0.017)):   # spars: 5 mm wide, height scales with chord
        h = h_root * chord / ROOT_CHORD
        zu, zl = naca2412(x_c)
        pieces.append((0.005 * h, 0.5 * (zu + zl) * chord, 0.005 * h**3 / 12))
    area = sum(a for a, _, _ in pieces)
    z_na = sum(a * z for a, z, _ in pieces) / area
    inertia = sum(a * (z - z_na) ** 2 + i0 for a, z, i0 in pieces)
    c_max = max(max(abs(naca2412(x)[0] * chord - z_na), abs(naca2412(x)[1] * chord - z_na)) for x in xs)
    return inertia, perimeter, c_max + skin / 2


def run(skin):
    n = 240
    dy = SPAN / n
    ys = [SPAN * i / n for i in range(n + 1)]
    secs = [section(ROOT_CHORD + (TIP_CHORD - ROOT_CHORD) * y / SPAN, skin) for y in ys]
    w = [PRESSURE * s[1] for s in secs]                 # N/m
    shear, moment = [0.0] * (n + 1), [0.0] * (n + 1)
    for i in range(n - 1, -1, -1):                      # integrate from the free tip
        shear[i] = shear[i + 1] + 0.5 * (w[i] + w[i + 1]) * dy
        moment[i] = moment[i + 1] + 0.5 * (shear[i] + shear[i + 1]) * dy
    curv = [moment[i] / (E * secs[i][0]) for i in range(n + 1)]
    slope, defl = [0.0] * (n + 1), [0.0] * (n + 1)
    for i in range(1, n + 1):                           # clamped root: w = w' = 0
        slope[i] = slope[i - 1] + 0.5 * (curv[i - 1] + curv[i]) * dy
        defl[i] = defl[i - 1] + 0.5 * (slope[i - 1] + slope[i]) * dy

    inertia, perimeter, c_max = secs[0]
    print(f"skin thickness            {skin * 1e3:.1f} mm")
    print(f"total load on semi-span   {shear[0]:.0f} N")
    print(f"root bending moment       {moment[0]:.0f} N m")
    print(f"root bending stress       {moment[0] * c_max / inertia / 1e6:.1f} MPa")
    print(f"tip deflection            {defl[-1] * 1e3:.2f} mm")
    skin_mass = sum(s[1] for s in secs) * dy * skin * 2700
    print(f"aluminium skin mass       {skin_mass:.2f} kg")
    if abs(skin - 0.008) < 1e-9:   # Rayleigh quotient with the static deflected shape
        num = sum(E * secs[i][0] * curv[i] ** 2 for i in range(n + 1)) * dy
        per_sum = sum(s[1] for s in secs) * dy
        for label, total in FE_MASS.items():
            m = [total * s[1] / per_sum for s in secs]      # FE mass, distributed like the skin
            den = sum(m[i] * defl[i] ** 2 for i in range(n + 1)) * dy
            print(f"1st flap frequency ({label}, Rayleigh)  {math.sqrt(num / den) / (2 * math.pi):.1f} Hz")


if __name__ == "__main__":
    run(float(sys.argv[1]) / 1e3 if len(sys.argv) > 1 else 0.008)
