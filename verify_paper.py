#!/usr/bin/env python3
"""
verify_paper.py -- high-precision check of "Discreteness as the late-time limit of a smooth flow"

Every identity, inequality and number stated in the paper is re-checked here in
arbitrary-precision arithmetic (mpmath, 50 digits by default) on random positions and
random times of EITHER sign.  The maps are re-implemented from the formulas in the paper,
independently of discreteness.py, so the two files check each other (see --no-code).

Usage
-----
  python verify_paper.py                    # 50 digits, 150 samples per check
  python verify_paper.py --dps 80 --samples 400 --seed 7
  python verify_paper.py --no-code          # skip the cross-check against discreteness.py

Requirements: Python >= 3.9, mpmath.  The optional cross-check also needs numpy and a
discreteness.py in the same directory.  Exit status is 0 if every check passes, 1 otherwise.

What is checked (paper section / result)
----------------------------------------
  2.2   cell form == global closed form (either sign of t)
  2.4   D^2 + ((1-k) sin 2 pi x)^2 = 2E ;  D, E > 0
  R1    R_0 = id ;  dR/dt = -sin(2 pi R) ;  sin(2 pi R) = 2k sin(2 pi x)/E
  R2    dR/dx = 2k/E ;  R_t(x+n) = R_t(x) + n
  R3    slopes e^{-2 pi t} at integers, e^{+2 pi t} at half-integers ; both fixed
  R4    R_{s+t} = R_s o R_t ;  R_{-t} o R_t = id                         (either sign)
  R5    gap to the rounding map stays 1/2 at every finite t (non-uniform limit)
  R6    d/dt V(R_t) = -sin^2(2 pi R_t)
  R7    tan(pi R_t) = e^{-2 pi t} tan(pi x)
  R8    monotone approach to the integers; closed cells mapped into themselves   (t >= 0, and
        separately for ALL real s <= t of either sign, as Result 8 is now stated)
  R9    e^{-2 pi t}|x-y| <= |R_t(x)-R_t(y)| <= e^{2 pi t}|x-y|                    (t >= 0)
  5     Moebius form (z+r)/(1+rz), r = tanh(pi t); multiplier at +1; velocity addition;
        hyperbolic distance 2 artanh r = 2 pi t; Poisson-kernel push-forward
  3/6   resolved fraction 2 R_{-t}(delta) = (2/pi) arctan(tan(pi delta)/k), checked against
        a root-find of the forward flow, plus its small-k asymptotics
  6     resolution time: the three quoted pairs, the asymptotic estimate, and the exact
        t* = ln[1/(tan(pi delta) tan(pi eps))] / (2 pi); the leading error (pi/6)(eps^2+delta^2)
        of the estimate; the unresolved radius eps_c(t) = (1/pi) arctan(k/tan(pi delta)), i.e.
        t*(eps_c(t), delta) = t and 1 - (resolved fraction) = 2 eps_c(t)
  code  (optional) float64 flow() and resolved_fraction() in discreteness.py vs this reference
"""
from __future__ import annotations

import argparse
import random
import sys

from mpmath import (atan, atanh, cos, diff, exp, fabs, findroot, log, mp, mpc, mpf, nint, nstr, pi,
                    sin, tan, tanh)

TP = None  # 2*pi; set in main() once the working precision is known


# --------------------------------------------------------------------------- #
# The maps, written from the formulas in the paper                            #
# --------------------------------------------------------------------------- #
def kfac(t):
    return exp(-TP * t)


def D_(x, t):
    c = cos(TP * x)
    return (1 + c) + kfac(t) * (1 - c)


def E_(x, t):
    c = cos(TP * x)
    return (1 + c) + kfac(t) ** 2 * (1 - c)


def R_closed(x, t):
    """Section 2.1: x - (1/pi) arctan((1-k) sin(2 pi x) / D)."""
    k = kfac(t)
    return x - atan((1 - k) * sin(TP * x) / D_(x, t)) / pi


def R_cell(x, t):
    """Section 2.2: n + (1/pi) arctan(k tan(pi (x-n))), valid off the half-integers."""
    n = nint(x)
    return n + atan(kfac(t) * tan(pi * (x - n))) / pi


def V(x):
    return -cos(TP * x) / TP


# --------------------------------------------------------------------------- #
# Bookkeeping                                                                 #
# --------------------------------------------------------------------------- #
class Report:
    def __init__(self):
        self.rows = []  # (name, value, threshold, kind)

    def residual(self, name, value, tol):
        self.rows.append((name, value, tol, "res"))

    def count(self, name, bad, total):
        self.rows.append((name, bad, 0, f"cnt:{total}"))

    def flag(self, name, ok):
        self.rows.append((name, 0 if ok else 1, 0, "flag"))

    def show(self):
        w = max(len(r[0]) for r in self.rows)
        print(f"{'check':<{w}}  {'result':>14}   status")
        print("-" * (w + 28))
        allok = True
        for name, val, tol, kind in self.rows:
            if kind == "res":
                ok, res = val < tol, nstr(val, 3)
            elif kind.startswith("cnt"):
                ok, res = val == 0, f"{int(val)}/{kind.split(':')[1]} viol."
            else:
                ok, res = val == 0, "ok" if val == 0 else "FAILED"
            allok &= ok
            print(f"{name:<{w}}  {res:>14}   {'PASS' if ok else 'FAIL'}")
        print("-" * (w + 28))
        print("ALL CHECKS PASSED" if allok else "SOME CHECKS FAILED")
        return allok


class MaxTracker(dict):
    """Keeps the largest |value| seen per check name (insertion-ordered)."""

    def update_max(self, name, v):
        v = fabs(v)
        if name not in self or v > self[name]:
            self[name] = v


# --------------------------------------------------------------------------- #
# Checks                                                                      #
# --------------------------------------------------------------------------- #
def check_identities(rep, rng, N, tol):
    m = MaxTracker()
    for _ in range(N):
        x, t, s = mpf(rng.uniform(-3, 3)), mpf(rng.uniform(-1.2, 1.2)), mpf(rng.uniform(-1.2, 1.2))
        k = kfac(t)
        y = R_closed(x, t)
        th = TP * y

        m.update_max("2.2  cell form == closed form (t of either sign)", R_closed(x, t) - R_cell(x, t))
        m.update_max("2.4  D^2 + ((1-k) sin)^2 = 2E", D_(x, t) ** 2 + ((1 - k) * sin(TP * x)) ** 2 - 2 * E_(x, t))
        m.update_max("2.4  D and E strictly positive (violation size)",
                     max(mpf(0), -D_(x, t)) + max(mpf(0), -E_(x, t)))
        m.update_max("R1   R_0 = id", R_closed(x, 0) - x)
        m.update_max("R1   dR/dt = -sin(2 pi R)", diff(lambda tt: R_closed(x, tt), t) + sin(th))
        m.update_max("R1   sin(2 pi R) = 2k sin(2 pi x)/E", sin(th) - 2 * k * sin(TP * x) / E_(x, t))
        m.update_max("R2   dR/dx = 2k/E", diff(lambda xx: R_closed(xx, t), x) - 2 * k / E_(x, t))
        m.update_max("R2   shift equivariance R(x+2) = R(x)+2", R_closed(x + 2, t) - R_closed(x, t) - 2)
        m.update_max("R4   group law R_{s+t} = R_s o R_t (either sign)", R_closed(x, s + t) - R_closed(y, s))
        m.update_max("R4   inverse R_{-t} o R_t = id", R_closed(y, -t) - x)
        m.update_max("R6   dV/dt = -sin^2(2 pi R)", diff(lambda tt: V(R_closed(x, tt)), t) + sin(th) ** 2)
        if fabs(cos(pi * x)) > mpf("1e-3"):          # tan has a pole at half-integers
            m.update_max("R7   tan(pi R) = k tan(pi x)", tan(pi * y) - k * tan(pi * x))

        # Section 5: geometry
        z, r, rs = exp(mpc(0, 1) * TP * x), tanh(pi * t), tanh(pi * s)
        m.update_max("5    Moebius form (z+r)/(1+rz), r = tanh(pi t)",
                     exp(mpc(0, 1) * th) - (z + r) / (1 + r * z))
        m.update_max("5    multiplier at +1: (1-r)/(1+r) = e^{-2 pi t}", (1 - r) / (1 + r) - k)
        m.update_max("5    velocity addition for r_{s+t}", tanh(pi * (s + t)) - (rs + r) / (1 + rs * r))
        m.update_max("5    hyperbolic distance 2 artanh r = 2 pi t", 2 * atanh(r) - TP * t)
        m.update_max("5    Poisson kernel: E/(2k) = P_r(2 pi R)",
                     E_(x, t) / (2 * k) - (1 - r * r) / (1 - 2 * r * cos(th) + r * r))
    for name, v in m.items():
        rep.residual(name, v, tol)


def check_fixed_points(rep, rng, tol):
    m = MaxTracker()
    for _ in range(40):
        t = mpf(rng.uniform(-1.2, 1.2))
        n = rng.randint(-50, 50)
        h = mpf(n) + mpf(1) / 2
        m.update_max("R3   integers fixed", R_closed(mpf(n), t) - n)
        m.update_max("R3   half-integers fixed", R_closed(h, t) - h)
        m.update_max("R3   slope at integer = e^{-2 pi t}", diff(lambda xx: R_closed(xx, t), mpf(n)) - kfac(t))
        m.update_max("R3   slope at half-integer = e^{+2 pi t}", diff(lambda xx: R_closed(xx, t), h) - 1 / kfac(t))
    for name, v in m.items():
        rep.residual(name, v, tol)


def check_inequalities(rep, rng, N):
    slack = mpf(10) ** (-(mp.dps - 8))               # allow for rounding at the working precision only
    v8 = v8c = v9 = 0
    for _ in range(N):
        x = mpf(rng.uniform(-0.5, 0.5))
        s = mpf(rng.uniform(0, 1.2))
        t = s + mpf(rng.uniform(0, 1.2))
        if fabs(R_cell(x, t)) > fabs(R_cell(x, s)) + slack:
            v8 += 1
        if fabs(R_cell(x, t)) > mpf(1) / 2 + slack:
            v8c += 1
        y = x + mpf(rng.uniform(-0.3, 0.3))
        u = mpf(rng.uniform(0, 1.5))
        if x != y:
            ratio = fabs(R_closed(x, u) - R_closed(y, u)) / fabs(x - y)
            k = kfac(u)
            if not (k * (1 - slack) <= ratio <= (1 / k) * (1 + slack)):
                v9 += 1
    rep.count("R8   |R_t - n| <= |R_s - n| for 0 <= s <= t", v8, N)
    rep.count("R8   closed cell mapped into itself", v8c, N)
    rep.count("R9   e^{-2 pi t} <= |dR|/|dx| <= e^{2 pi t}", v9, N)


def check_monotone_either_sign(rep, rng, N):
    """Result 8 as stated in the README: for ALL real s <= t (either sign) the distance to the
    nearest integer is non-increasing in t, and every closed cell is mapped into itself at every time.
    Uses its own random stream (see main) so that adding this check leaves the samples of all other
    checks unchanged.  Cells n = -5..5 are used, not only n = 0."""
    slack = mpf(10) ** (-(mp.dps - 8))
    v8 = v8c = 0
    for _ in range(N):
        n = rng.randint(-5, 5)
        x = mpf(n) + mpf(rng.uniform(-0.5, 0.5))
        s = mpf(rng.uniform(-1.2, 1.2))
        t = s + mpf(rng.uniform(0, 2.4))
        ds, dtt = fabs(R_cell(x, s) - n), fabs(R_cell(x, t) - n)
        if dtt > ds + slack:
            v8 += 1
        if ds > mpf(1) / 2 + slack or dtt > mpf(1) / 2 + slack:
            v8c += 1
    rep.count("R8   |R_t - n| <= |R_s - n| for ALL real s <= t (either sign)", v8, N)
    rep.count("R8   closed cell mapped into itself at times of either sign", v8c, N)


def check_nonuniform_limit(rep):
    """Result 5 remark: at every finite t there are starts whose distance to round(x) is still ~1/2."""
    worst = mpf(0)
    for t in (mpf(1), mpf(3), mpf(6)):
        eps = mpf("1e-40")                           # start 1e-40 below the half-integer 1/2; round(x) = 0
        gap = fabs(R_cell(mpf(1) / 2 - eps, t))
        worst = max(worst, fabs(gap - mpf(1) / 2))
    rep.residual("R5   gap to rounding is 1/2 at finite t (eps=1e-40)", worst, mpf("1e-10"))


def check_resolved_fraction(rep, rng, tol):
    """The measure of {x in a unit cell : |R_t(x) - n| <= delta} equals 2 R_{-t}(delta)
    = (2/pi) arctan(tan(pi delta)/k).  Independent check: find, by root-finding on the FORWARD
    flow, the start that lands exactly at delta."""
    m, asym = MaxTracker(), MaxTracker()
    rtol = mpf(10) ** (-(mp.dps - 5))
    for _ in range(25):
        t = mpf(rng.uniform(0, 3))
        dl = mpf(rng.uniform(0.005, 0.45))
        formula = (2 / pi) * atan(tan(pi * dl) / kfac(t))
        edge = findroot(lambda xx: R_cell(xx, t) - dl, formula / 2, tol=rtol, maxsteps=200)
        m.update_max("3/6  resolved fraction = 2 x (start landing at delta)", formula - 2 * edge)
        m.update_max("3/6  resolved fraction = 2 R_{-t}(delta)", formula - 2 * R_cell(dl, -t))
    for t, dl in ((mpf(4), mpf("0.01")), (mpf(6), mpf("0.001"))):
        k = kfac(t)
        unresolved = 1 - (2 / pi) * atan(tan(pi * dl) / k)
        asym.update_max("3/6  unresolved ~ 2k/(pi^2 delta): relative error", unresolved / (2 * k / (pi ** 2 * dl)) - 1)
    for name, v in m.items():
        rep.residual(name, v, tol)
    for name, v in asym.items():
        rep.residual(name, v, mpf("1e-3"))


def check_resolution_time(rep):
    quoted = [  # (eps, delta, exact, approx) exactly as printed in Section 6 (six decimals)
        (mpf("1e-3"), mpf("1e-3"), "1.834427", "1.834428"),
        (mpf("1e-6"), mpf("1e-4"), "3.300299", "3.300299"),
        (mpf("1e-9"), mpf("1e-6"), "5.132638", "5.132638"),
    ]
    extra = [(mpf("0.01"), mpf("0.1"), None, None), (mpf("0.2"), mpf("0.2"), None, None)]
    approx = lambda e, d: log(1 / (pi ** 2 * e * d)) / TP
    exact_closed = lambda e, d: log(1 / (tan(pi * d) * tan(pi * e))) / TP

    print("Section 6 -- resolution time t*(eps, delta); the last two rows are outside the paper's regime")
    print(f"  {'eps':>8} {'delta':>8} | {'exact (root-find)':>18} {'exact (closed)':>16} {'paper approx':>14} | {'approx - exact':>14}")
    worst, digits_ok = mpf(0), True
    for e, d, q_exact, q_approx in quoted + extra:
        x0 = mpf(1) / 2 - e
        te = findroot(lambda t: R_cell(x0, t) - d, approx(e, d), tol=mpf(10) ** (-(mp.dps - 8)), maxsteps=200)
        tc, ta = exact_closed(e, d), approx(e, d)
        worst = max(worst, fabs(te - tc))
        print(f"  {nstr(e, 2):>8} {nstr(d, 2):>8} | {nstr(te, 12):>18} {nstr(tc, 12):>16} {nstr(ta, 12):>14} | {nstr(ta - te, 3):>14}")
        if q_exact is not None:
            digits_ok &= f"{float(te):.6f}" == q_exact and f"{float(ta):.6f}" == q_approx
    print()
    rep.residual("6    exact t* (closed form) vs root-find", worst, mpf(10) ** (-(mp.dps - 10)))
    rep.flag("6    the three quoted (exact, approx) pairs, 6 decimals", digits_ok)
    rep.flag("6    t* > 0 iff delta + eps < 1/2",
             exact_closed(mpf("0.2"), mpf("0.2")) > 0 and exact_closed(mpf("0.3"), mpf("0.3")) < 0)


def check_resolution_links(rep, rng, tol):
    """Section 6 claims beyond the table: (i) the leading error (pi/6)(eps^2 + delta^2) of the
    small-scale estimate; (ii) eps_c(t) = (1/pi) arctan(k / tan(pi delta)) is the radius of the
    still-unresolved neighbourhood of a half-integer, so t*(eps_c(t), delta) = t, and the two ends
    of each cell give the unresolved fraction 1 - 2 R_{-t}(delta) = 2 eps_c(t)."""
    m = MaxTracker()
    for _ in range(25):
        t, dl = mpf(rng.uniform(0, 3)), mpf(rng.uniform(0.01, 0.45))
        k = kfac(t)
        ec = atan(k / tan(pi * dl)) / pi
        m.update_max("6    t*(eps_c(t), delta) = t", log(1 / (tan(pi * dl) * tan(pi * ec))) / TP - t)
        m.update_max("6    unresolved fraction = 2 eps_c(t)", 1 - (2 / pi) * atan(tan(pi * dl) / k) - 2 * ec)
    for name, v in m.items():
        rep.residual(name, v, tol)
    worst = mpf(0)                                   # relative deviation from the stated leading term
    for e, d in ((mpf("1e-3"), mpf("1e-3")), (mpf("1e-4"), mpf("1e-3")), (mpf("1e-5"), mpf("1e-4"))):
        approx = log(1 / (pi ** 2 * e * d)) / TP
        exact = log(1 / (tan(pi * d) * tan(pi * e))) / TP
        worst = max(worst, fabs((approx - exact) / (pi / 6 * (e * e + d * d)) - 1))
    rep.residual("6    estimate - exact = (pi/6)(eps^2+delta^2), relative", worst, mpf("1e-4"))


def check_code(rep, rng):
    """Optional: float64 discreteness.py vs the high-precision reference above."""
    try:
        import numpy as np
        import discreteness as dz
    except Exception as exc:  # noqa: BLE001
        print(f"(cross-check against discreteness.py skipped: {exc.__class__.__name__}: {exc})\n")
        return
    N = 400
    xs = [rng.uniform(-3, 3) for _ in range(N)]
    ts = [rng.uniform(-1.0, 6.0) for _ in range(N)]
    got = dz.flow(np.array(ts), np.array(xs))
    abs_err = rel_err = 0.0
    for x, t, g in zip(xs, ts, got):
        xm, tm = mpf(x), mpf(t)                      # the exact binary values the code received
        ref = R_cell(xm, tm)
        d = fabs(ref - nint(xm))                     # distance to the nearest integer
        abs_err = max(abs_err, float(fabs(mpf(float(g)) - ref)))
        # Relative accuracy of the distance to the integer.  Only meaningful in the cell n = 0:
        # flow() returns the absolute position n + d, and float64 cannot hold a distance d below
        # ~1e-16 next to an integer n != 0 (so the absolute error above is the right measure there).
        # Also skip starts within 1e-3 of a half-integer, where rounding of the float input dominates.
        if nint(xm) == 0 and d > mpf("1e-300") and fabs(mpf(1) / 2 - fabs(xm)) > mpf("1e-3"):
            rel_err = max(rel_err, float(fabs(mpf(float(g)) - ref) / d))
    rep.residual("code flow(): max abs error vs reference (t in [-1,6])", mpf(abs_err), mpf("1e-12"))
    rep.residual("code flow(): max relative error to the integer (cell n=0)", mpf(rel_err), mpf("1e-10"))
    if hasattr(dz, "resolved_fraction"):
        worst = max(float(fabs(mpf(dz.resolved_fraction(t)) - (2 / pi) * atan(tan(pi * mpf("0.02")) / kfac(mpf(t)))))
                    for t in (0.0, 0.1, 0.4, 1.0, 1.5, 3.0))
        rep.residual("code resolved_fraction(): max abs error", mpf(worst), mpf("1e-12"))
    if hasattr(dz, "resolution_time"):
        pairs = [(1e-3, 1e-3), (1e-6, 1e-4), (1e-9, 1e-6), (1e-2, 1e-1), (0.2, 0.2)]
        worst = max(float(fabs(mpf(dz.resolution_time(e, d))
                               - log(1 / (tan(pi * mpf(d)) * tan(pi * mpf(e)))) / TP)) for e, d in pairs)
        rep.residual("code resolution_time(): max abs error", mpf(worst), mpf("1e-12"))


# --------------------------------------------------------------------------- #
def main(argv=None):
    global TP
    ap = argparse.ArgumentParser(description="High-precision verification of the paper's results.")
    ap.add_argument("--dps", type=int, default=50, help="working precision in decimal digits (default 50)")
    ap.add_argument("--samples", type=int, default=150, help="random samples per check (default 150)")
    ap.add_argument("--seed", type=int, default=1, help="random seed (default 1)")
    ap.add_argument("--no-code", action="store_true", help="skip the cross-check against discreteness.py")
    args = ap.parse_args(argv)
    if args.dps < 40:
        ap.error("--dps must be at least 40 (the pass thresholds assume it)")

    mp.dps = args.dps
    TP = 2 * pi
    rng = random.Random(args.seed)

    print(f"working precision {mp.dps} digits, {args.samples} samples per check, seed {args.seed}; "
          f"exact identities must agree to < 1e-30\n")
    rep = Report()
    check_identities(rep, rng, args.samples, mpf(10) ** (-30))
    check_fixed_points(rep, rng, mpf(10) ** (-30))
    check_inequalities(rep, rng, max(args.samples, 1000))
    check_monotone_either_sign(rep, random.Random(args.seed + 1000), max(args.samples, 1000))
    check_nonuniform_limit(rep)
    check_resolved_fraction(rep, rng, mpf(10) ** (-30))
    check_resolution_time(rep)
    check_resolution_links(rep, rng, mpf(10) ** (-30))
    if not args.no_code:
        check_code(rep, rng)
    return 0 if rep.show() else 1


if __name__ == "__main__":
    sys.exit(main())
