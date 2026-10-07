#!/usr/bin/env python3
"""
oim_exact_flow_study.py
=======================
Does the exact flow of  dx/dt = -sin(2 pi x)  (README: "Discreteness as the
late-time limit of a smooth flow") survive noise and non-factorable couplings,
and is there any practical advantage over the usual approaches?

The testbed is an oscillator Ising machine (OIM, Wang & Roychowdhury 2019):

    dphi_i/dt = -K * sum_j J_ij sin(phi_i - phi_j)  -  Ks(t) * sin(2 phi_i)  [+ sigma * noise]

with J = -A (max-cut on a random sparse graph).  With x = phi/pi and
t_paper = S/pi, S = integral of Ks dt, the injection term is exactly the README flow.

Experiments (everything printed is measured, nothing is hard-coded)
-------------------------------------------------------------------
  0  Exactness checks: README map == RK4/DOP853; group law merges sub-steps;
     NEW (not in README): constant detuning keeps an exact closed form.
  1  COUPLING.  How fast does "use the exact injection flow, ignore the coupling"
     degrade with rho = K*lambda_max/Ks?  And does the README's closed-form
     derivative (2k/E) give a clean first-order correction?
  2  NOISE.  (a) the noise-blurred basin boundary (exact committor)
             (b) the floor under the unresolved fraction (exact von Mises law)
             (c) the noise-capped resolution time
  3  BENEFIT?  Exact-injection Strang splitting vs Euler / RK4 / LSODA / BDF (baseline = best explicit scheme):
     (a) cost to reach a fixed fidelity vs injection stiffness r = Ks_max/(K*lambda)
     (b) wall-clock vs graph size at strong injection
     (c) what happens to the benefit once noise is present
  4  Verdict computed from the numbers above, plus a figure.

Usage
-----
    python oim_exact_flow_study.py              # full run (about 1.5-2 min; timings vary by machine)
    python oim_exact_flow_study.py --quick      # smaller samples (about 30 s)
    python oim_exact_flow_study.py --no-plot    # skip the figure
    python oim_exact_flow_study.py --seed 7 --outdir results

Requirements: Python >= 3.9, numpy, scipy; matplotlib for the figure.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from scipy.integrate import quad, solve_ivp
from scipy.special import i0e, i1e
from scipy.stats import norm

TWO_PI = 2.0 * math.pi
PI = math.pi


# =========================================================================== #
# 0.  The exact injection flow (README) and its relatives                     #
# =========================================================================== #
def R_paper(t, x):
    """README closed form (Sec. 2.1):  R_t(x) = x - atan((1-k) sin 2pi x / D)/pi."""
    t = np.asarray(t, float)
    x = np.asarray(x, float)
    k = np.exp(-TWO_PI * t)
    c, s = np.cos(TWO_PI * x), np.sin(TWO_PI * x)
    D = (1 + c) + k * (1 - c)
    return x - np.arctan((1 - k) * s / D) / PI


def R_cell(t, x):
    """README cell form (Sec. 2.2), written with atan2 (stable for large t)."""
    k = np.exp(-TWO_PI * np.asarray(t, float))
    x = np.asarray(x, float)
    n = np.round(x)
    y = x - n
    return n + np.arctan2(k * np.sin(PI * y), np.cos(PI * y)) / PI


def inject(phi, S):
    """Exact flow of  dphi/dt = -Ks(t) sin(2 phi)  over a window with S = int Ks dt.

    Identical to the README map under x = phi/pi, t = S/pi (k = e^{-2 pi t} = e^{-2 S}).
    """
    k = np.exp(-2.0 * np.asarray(S, float))
    n = np.round(phi / PI)
    y = phi - n * PI
    return n * PI + np.arctan2(k * np.sin(y), np.cos(y))


def dinject(phi, S):
    """d inject / d phi = k / (cos^2 y + k^2 sin^2 y)   (= 2k/E of README Result 2)."""
    k = np.exp(-2.0 * np.asarray(S, float))
    y = phi - np.round(phi / PI) * PI
    return k / (np.cos(y) ** 2 + k ** 2 * np.sin(y) ** 2)


def inject_detuned(phi, S, d):
    """NEW (derived here, checked numerically below).  Exact flow of

            dphi/dS = d - sin(2 phi),     |d| < 1     (d = detuning / Ks)

    The tangent of the shifted angle obeys a Riccati equation, which is linear in
    homogeneous coordinates (p, q) = (sin psi, cos psi):  p' = -w p,  q' = -2d p + w q,
    w = sqrt(1-d^2).  The fixed points move to phi_s = asin(d)/2 (stable) and
    pi/2 - phi_s (unstable); the cell boundaries move with them.
    """
    w = math.sqrt(1.0 - d * d)
    phi_s = 0.5 * math.asin(d)
    psi = phi - phi_s
    psi_u = PI / 2 - 2 * phi_s
    m = np.floor((psi - psi_u) / PI) + 1
    p0, q0 = np.sin(psi), np.cos(psi)
    p = p0 * math.exp(-w * S)
    q = (q0 - d / w * p0) * math.exp(w * S) + d / w * p0 * math.exp(-w * S)
    al = np.arctan2(p, q)
    centre = psi_u + (m - 0.5) * PI
    return al + TWO_PI * np.round((centre - al) / TWO_PI) + phi_s


def experiment0_exactness(rng):
    print("\n[0] Exactness checks")
    x = rng.uniform(-3, 3, 400)
    t = rng.uniform(0, 2, 400)
    e_cell = np.max(np.abs(R_cell(t, x) - R_paper(t, x)))
    e_map = np.max(np.abs(PI * R_cell(t, x) - inject(PI * x, PI * t)))
    ks = 1.3
    ref = np.array([solve_ivp(lambda s, p: -ks * np.sin(2 * p), (0, ti / ks * PI), [PI * xi],
                              method="DOP853", rtol=1e-13, atol=1e-14).y[0, -1]
                    for xi, ti in zip(x[:100], t[:100])])
    e_ode = np.max(np.abs(inject(PI * x[:100], PI * t[:100]) - ref))
    s1, s2 = rng.uniform(0, 1, 2)
    e_grp = np.max(np.abs(inject(inject(PI * x, s1), s2) - inject(PI * x, s1 + s2)))
    print(f"    README cell form vs README global closed form      {e_cell:9.2e}")
    print(f"    inject(phi,S) vs README map (x=phi/pi, t=S/pi)     {e_map:9.2e}")
    print(f"    inject vs DOP853 integration of the ODE            {e_ode:9.2e}")
    print(f"    group law: inject(S2) o inject(S1) = inject(S1+S2) {e_grp:9.2e}")
    worst_det = 0.0
    for d in (0.0, 0.3, -0.6, 0.95):
        for S in (0.1, 0.7, 3.0):
            x0 = rng.uniform(-6, 6, 60)
            ex = np.array([solve_ivp(lambda s, p: d - np.sin(2 * p), (0, S), [a], method="DOP853",
                                     rtol=1e-13, atol=1e-14).y[0, -1] for a in x0])
            worst_det = max(worst_det, np.max(np.abs(inject_detuned(x0, S, d) - ex)))
    print(f"    NEW: detuned closed form vs DOP853 (d up to 0.95)  {worst_det:9.2e}")
    return dict(ok=max(e_cell, e_map, e_ode, e_grp, worst_det) < 1e-9, worst_detuned=worst_det)


# =========================================================================== #
# 1.  The OIM testbed                                                         #
# =========================================================================== #
def make_graph(N, rng, deg=3):
    """Random sparse undirected graph (about `deg` neighbours per node)."""
    rows, cols = [], []
    for i in range(N):
        for j in rng.choice(N, deg // 2 + 1, replace=False):
            if i != j:
                rows += [i, j]
                cols += [j, i]
    A = sps.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(N, N)).tocsr()
    A.data[:] = 1.0
    return A


class OIM:
    """dphi/dt = -K sum_j J_ij sin(phi_i-phi_j) - Ks(t) sin(2 phi_i),  J = -A (max-cut)."""

    def __init__(self, A, K=1.0):
        self.A, self.J, self.K = A.tocsr(), -A.tocsr(), K
        self.N = A.shape[0]
        self.lam = float(abs(spla.eigsh(self.J, k=1, which="LM", return_eigenvectors=False)[0]))
        self.edges = A.nnz / 2

    def coupling(self, phi):
        s, c = np.sin(phi), np.cos(phi)
        return -self.K * (s * (self.J @ c) - c * (self.J @ s))

    def cut(self, phi):
        s = np.where(np.cos(phi) > 0, 1.0, -1.0)
        return 0.5 * (self.edges - 0.5 * s @ (self.A @ s))

    def jac(self, phi, Ks):
        Jc = self.J.tocoo()
        c = np.cos(phi[Jc.row] - phi[Jc.col])
        off = sps.coo_matrix((self.K * Jc.data * c, (Jc.row, Jc.col)), shape=(self.N, self.N))
        diag = -np.bincount(Jc.row, weights=self.K * Jc.data * c, minlength=self.N) - 2 * Ks * np.cos(2 * phi)
        return (off + sps.diags(diag)).tocsc()


class Ramp:
    """Ks(t): linear ramp 0 -> Ks_max over [0, T], then constant.  G(t) = int_0^t Ks."""

    def __init__(self, Ks_max, T=4.0, t_end=6.0):
        self.Kmax, self.T, self.t_end = Ks_max, T, t_end

    def Ks(self, t):
        return self.Kmax * min(1.0, t / self.T)

    def G(self, t):
        return 0.5 * self.Kmax * t * t / self.T if t <= self.T else 0.5 * self.Kmax * self.T + self.Kmax * (t - self.T)

    def S(self, a, b):
        return self.G(b) - self.G(a)


# --- fixed-step integrators; each returns (phi_final, number_of_coupling_evaluations) --- #
def euler(m, rp, p0, h):
    p, n = p0.copy(), int(round(rp.t_end / h))
    for i in range(n):
        p = p + h * (m.coupling(p) - rp.Ks(i * h) * np.sin(2 * p))
    return p, n


def rk4(m, rp, p0, h):
    f = lambda t, q: m.coupling(q) - rp.Ks(t) * np.sin(2 * q)
    p, n = p0.copy(), int(round(rp.t_end / h))
    for i in range(n):
        t = i * h
        k1 = f(t, p)
        k2 = f(t + h / 2, p + h / 2 * k1)
        k3 = f(t + h / 2, p + h / 2 * k2)
        k4 = f(t + h, p + h * k3)
        p = p + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return p, 4 * n


def strang_exact(m, rp, p0, h):
    """Strang splitting:  A(h/2) B(h) A(h/2),  A = EXACT injection flow, B = Heun on the coupling.

    The group law lets consecutive half-steps merge, so A is called once per step."""
    n = int(round(rp.t_end / h))
    p = inject(p0, rp.S(0.0, h / 2))
    for i in range(n):
        f1 = m.coupling(p)
        f2 = m.coupling(p + h * f1)
        p = p + 0.5 * h * (f1 + f2)
        a, b = (i * h + h / 2, (i + 1) * h + h / 2) if i < n - 1 else (i * h + h / 2, (i + 1) * h)
        p = inject(p, rp.S(a, b))
    return p, 2 * n


def reference(m, rp, p0):
    return solve_ivp(lambda t, q: m.coupling(q) - rp.Ks(t) * np.sin(2 * q), (0, rp.t_end), p0,
                     method="DOP853", rtol=1e-10, atol=1e-11).y[:, -1]


def lsoda(m, rp, p0, rtol=1e-3):
    s = solve_ivp(lambda t, q: m.coupling(q) - rp.Ks(t) * np.sin(2 * q), (0, rp.t_end), p0,
                  method="LSODA", rtol=rtol, atol=rtol * 1e-3,
                  jac=lambda t, q: m.jac(q, rp.Ks(t)).toarray())
    return s.y[:, -1], s.nfev, s.njev


def bdf_sparse(m, rp, p0, rtol=1e-3):
    pat = (abs(m.J) + sps.identity(m.N)).tocsr()
    s = solve_ivp(lambda t, q: m.coupling(q) - rp.Ks(t) * np.sin(2 * q), (0, rp.t_end), p0,
                  method="BDF", rtol=rtol, atol=rtol * 1e-3,
                  jac=lambda t, q: m.jac(q, rp.Ks(t)), jac_sparsity=pat)
    return s.y[:, -1]


def spins(p):
    return np.cos(p) > 0


# =========================================================================== #
# 2.  Experiment 1: coupling                                                  #
# =========================================================================== #
def experiment1_coupling(rng, quick):
    """Ignore the coupling (zeroth order) vs add the closed-form first-order correction.

    Variation of constants:  phi(t) ~ R(t,phi0) + int_0^t R'_{t-s}(R_s(phi0)) F(R_s(phi0)) ds,
    where F is the coupling force along the zeroth-order trajectories and R' = 2k/E (README).
    """
    print("\n[1] Coupling: error of the exact-injection description vs rho = K*lambda_max/Ks")
    N, M = 80, (16 if quick else 40)
    m = OIM(make_graph(N, rng))
    Ks, tf = 1.0, 0.7
    ics = [rng.uniform(-3 * PI, 3 * PI, N) for _ in range(M)]
    rhos = [0.01, 0.02, 0.05, 0.1, 0.2, 0.4]
    ss = np.linspace(0, tf, 161)
    w = np.full(ss.size, ss[1] - ss[0])
    w[[0, -1]] *= 0.5
    out = dict(rho=rhos, e0=[], e1=[], flip=[])
    print(f"    (N={N}, lambda_max={m.lam:.2f}, Ks={Ks}, t={tf}, {M} random initial phase vectors)")
    print(f"    {'rho':>6s} {'basin flips':>12s} {'err 0th order':>14s} {'err 1st order':>14s} {'0th/rho':>9s} {'1st/rho^2':>10s}")
    for rho in rhos:
        K = rho * Ks / m.lam
        m.K = K
        e0, e1, flips = [], [], []
        for p in ics:
            ex = solve_ivp(lambda t, q: m.coupling(q) - Ks * np.sin(2 * q), (0, tf), p, method="DOP853",
                           rtol=1e-12, atol=1e-13).y[:, -1]
            z0 = inject(p, Ks * tf)
            n0, nf = np.round(p / PI), np.round(ex / PI)
            corr = np.zeros(N)
            for s, wi in zip(ss, w):
                zs = inject(p, Ks * s)
                corr += wi * dinject(zs, Ks * (tf - s)) * m.coupling(zs)
            y0 = p - n0 * PI
            keep = (np.abs(y0) < PI / 2 - 0.35) & (n0 == nf)    # away from the saddle layer, no basin flip
            flips.append(np.mean(n0 != nf))
            e0.append(np.abs(ex - z0)[keep])
            e1.append(np.abs(ex - (z0 + corr))[keep])
        e0, e1 = np.median(np.concatenate(e0)), np.median(np.concatenate(e1))
        out["e0"].append(e0)
        out["e1"].append(e1)
        out["flip"].append(float(np.mean(flips)))
        print(f"    {rho:6.3f} {np.mean(flips):12.4f} {e0:14.3e} {e1:14.3e} {e0 / rho:9.3f} {e1 / rho ** 2:10.4f}")
    lr = np.log(rhos)
    out["slope0"] = np.polyfit(lr, np.log(out["e0"]), 1)[0]
    out["slope1"] = np.polyfit(lr, np.log(out["e1"]), 1)[0]
    out["gain"] = out["e0"][3] / out["e1"][3]
    print(f"    fitted exponents:  err0 ~ rho^{out['slope0']:.2f},  err1 ~ rho^{out['slope1']:.2f}")
    print("    (basin flips are rare events -- a handful of oscillators per run -- so only their size, not a scaling law, is reported)")
    return out


# =========================================================================== #
# 3.  Experiment 2: noise  (single oscillator, README units: dx = -sin(2 pi x) dt + sigma dW)
# =========================================================================== #
def noisy_flow(x0, sigma, T, h, rng, record=None):
    """Strang splitting with the README map as the exact drift flow.  Half-steps merge (group law)."""
    x = R_cell(h / 2, x0)
    n = int(round(T / h))
    root = sigma * math.sqrt(h)
    for i in range(n):
        x = x + root * rng.standard_normal(x.shape)
        x = R_cell(h if i < n - 1 else h / 2, x)
        if record is not None:
            record(i + 1, x)
    return x


def committor(x, sigma):
    """Exact P(reach integer 1 before integer 0 | start x in (0,1)) for dx=-sin(2 pi x)dt+sigma dW.

    1-D diffusion theory: scale density s'(y) = exp(2V(y)/sigma^2) = exp(-cos(2 pi y)/(pi sigma^2))."""
    f = lambda y: math.exp(-(math.cos(TWO_PI * y) + 1.0) / (PI * sigma ** 2))
    return quad(f, 0, x, limit=200)[0] / quad(f, 0, 1, limit=200)[0]


def t_star(eps, delta):
    """README Sec. 6 exact resolution time."""
    return np.log(1.0 / (np.tan(PI * delta) * np.tan(PI * eps))) / TWO_PI


def experiment2_noise(rng, quick):
    out = {}
    print("\n[2] Noise  (units of the README: a = 1, so the saddle expands as e^{+2 pi t})")
    print("    Linearising at the half-integer (README Result 3: slope e^{2 pi t}) predicts a blur width s = sigma/sqrt(4 pi).")

    # (a) noisy basin boundary
    npt = 8000 if quick else 20000
    ms = np.array([-2, -1, -0.5, 0, 0.5, 1, 2.0])
    print("  (a) which integer is reached from x = 1/2 + m*s ?   simulated vs exact committor")
    rows = []
    for sigma in (0.10, 0.25):
        s = sigma / math.sqrt(4 * PI)
        x0 = np.repeat(0.5 + ms * s, npt)
        xf = noisy_flow(x0, sigma, 2.5, 4e-3, rng)
        p_sim = (np.round(xf) == 1).reshape(len(ms), npt).mean(axis=1)
        p_ex = np.array([committor(0.5 + mm * s, sigma) for mm in ms])
        p_g = norm.cdf(ms)
        z = (p_sim - p_ex) / np.sqrt(np.maximum(p_ex * (1 - p_ex), 1e-4) / npt)
        rows.append((sigma, s, p_sim, p_ex, p_g, z))
        print(f"      sigma={sigma}: s={s:.4f}   max |sim-exact| = {np.max(np.abs(p_sim - p_ex)):.4f}   "
              f"max |z| = {np.max(np.abs(z)):.2f}   max |gauss-exact| = {np.max(np.abs(p_g - p_ex)):.4f}")
    out["basin"] = dict(ms=ms, rows=rows, max_z=max(np.max(np.abs(r[5])) for r in rows), npt=npt)

    # (b) unresolved fraction and its floor
    sigma, delta = 0.10, 0.10
    M = 60000 if quick else 200000
    kappa = 1.0 / (PI * sigma ** 2)
    num = quad(lambda y: math.exp(kappa * (math.cos(TWO_PI * y) - 1)), -delta, delta)[0]
    den = quad(lambda y: math.exp(kappa * (math.cos(TWO_PI * y) - 1)), -0.5, 0.5)[0]
    floor = 1 - num / den
    checkpoints = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.5]
    x0 = rng.uniform(0, 1, M)
    h = 2e-3        # splitting samples the post-flow state: width biased by ~ -pi*h*100%, i.e. ~ -8% in a 3.5-sigma tail at h=2e-3
    seen = {}

    def rec(i, x):
        tt = round(i * h, 6)
        if any(abs(tt - c) < 1e-9 for c in checkpoints):
            seen[min(checkpoints, key=lambda c: abs(c - tt))] = float(np.mean(np.abs(x - np.round(x)) > delta))

    noisy_flow(x0, sigma, max(checkpoints), h, rng, record=rec)
    det = {c: (2 / PI) * math.atan(math.exp(-TWO_PI * c) / math.tan(PI * delta)) for c in checkpoints}
    print(f"  (b) unresolved fraction P(|x - round(x)| > delta), sigma={sigma}, delta={delta}, uniform start")
    print(f"      exact stationary floor (von Mises law)   = {floor:.3e}")
    print(f"      {'t':>5s} {'noisy sim':>11s} {'README formula':>15s} {'README + floor':>15s}")
    for c in checkpoints:
        print(f"      {c:5.1f} {seen[c]:11.4e} {det[c]:15.4e} {det[c] + floor:15.4e}")
    se = math.sqrt(floor * (1 - floor) / M)
    out["floor"] = dict(t=checkpoints, sim=[seen[c] for c in checkpoints], det=[det[c] for c in checkpoints],
                        floor=floor, z_late=(seen[checkpoints[-1]] - floor) / se,
                        max_excess=max(seen[c] / (det[c] + floor) - 1 for c in checkpoints), M=M)
    print(f"      late-time sim vs floor: z = {out['floor']['z_late']:+.2f};  "
          f"largest excess over 'README + floor' in the crossover region: {100 * out['floor']['max_excess']:.0f}%")

    # (c) resolution time
    s = sigma / math.sqrt(4 * PI)
    eps_list = [1e-6, 1e-2, 5e-2]
    npt = 6000 if quick else 20000
    x0 = np.concatenate([np.full(npt, 0.5 - e) for e in eps_list])
    first = np.full(x0.size, np.nan)
    h = 4e-3

    def rec2(i, x):
        hit = (np.abs(x - np.round(x)) < delta) & np.isnan(first)
        first[hit] = i * h

    noisy_flow(x0, sigma, 3.0, h, rng, record=rec2)
    print(f"  (c) time to get within delta={delta} of the destination integer, start at 1/2 - eps")
    print("      prediction = README t* with eps -> |eps + eta|, eta ~ N(0, s^2);  quantiles 10/50/90%")
    res = []
    for j, e in enumerate(eps_list):
        q_sim = np.nanpercentile(first[j * npt:(j + 1) * npt], [10, 50, 90])
        eta = s * rng.standard_normal(400000)
        q_pred = np.percentile(np.maximum(t_star(np.abs(e + eta), delta), 0), [10, 50, 90])
        det_t = t_star(e, delta)
        res.append((e, q_sim, q_pred, det_t, first[j * npt:(j + 1) * npt].copy()))
        print(f"      eps={e:7.0e}: sim {np.round(q_sim, 3)}   pred {np.round(q_pred, 3)}   noise-free t* = {det_t:.3f}")
    out["restime"] = dict(rows=res, max_rel_err=max(np.max(np.abs(r[1] - r[2]) / r[1]) for r in res))
    return out


# =========================================================================== #
# 4.  Experiment 3: is there a benefit?                                       #
# =========================================================================== #
def fidelity(m, rp, ics, refs, fn, h):
    full = frac = 0.0
    ev = 0
    cut = []
    for p, r in zip(ics, refs):
        with np.errstate(all="ignore"):
            q, ev = fn(m, rp, p, h)
        if not np.all(np.isfinite(q)):
            cut.append(np.nan)
            continue
        ok = spins(q) == spins(r)
        full += ok.all()
        frac += ok.mean()
        cut.append(m.cut(q))
    M = len(ics)
    return full / M, frac / M, ev, float(np.nanmean(cut))


def min_cost(m, rp, ics, refs, fn, grid, target=0.95):
    """Coarse-to-fine scan; return (h, evals) for the coarsest h that reaches `target`
    full-configuration agreement with the reference, confirmed at the next finer h too."""
    prev = None
    for i, h in enumerate(grid):
        full, _, ev, _ = fidelity(m, rp, ics, refs, fn, h)
        if full >= target:
            if i + 1 < len(grid) and fidelity(m, rp, ics, refs, fn, grid[i + 1])[0] < target:
                continue
            return h, ev
    return None, None


def experiment3a_cost_vs_stiffness(rng, quick):
    print("\n[3a] Cost to reach >= 95% identical final spin configurations (vs a 1e-10 reference)")
    N, M = 100, (10 if quick else 20)
    m = OIM(make_graph(N, rng))
    rs = [1, 10, 100] if quick else [1, 3, 10, 30, 100]
    grids = dict(euler=[0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001],
                 rk4=[0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002],
                 strang=[0.4, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005])
    fns = dict(euler=euler, rk4=rk4, strang=strang_exact)
    print(f"    N={N}, lambda_max={m.lam:.2f}, {M} random initial conditions; cost = coupling evaluations (2 sparse matvecs each)")
    print(f"    {'r=Ks/(K lam)':>12s} | {'Euler':>16s} | {'RK4':>16s} | {'exact-inj. Strang':>18s} | {'LSODA (nfev / njev)':>20s}")
    out = dict(r=rs, euler=[], rk4=[], strang=[], lsoda=[], h={})
    for r in rs:
        rp = Ramp(r * m.K * m.lam)
        ics = [rng.uniform(-PI, PI, N) for _ in range(M)]
        refs = [reference(m, rp, p) for p in ics]
        cell = {}
        for name in ("euler", "rk4", "strang"):
            h, ev = min_cost(m, rp, ics, refs, fns[name], grids[name])
            cell[name] = (h, ev)
            out[name].append(ev)
            out["h"][(name, r)] = h
        nf, nj, ok = [], [], 0
        for p, rf in zip(ics, refs):
            q, a, b = lsoda(m, rp, p)
            nf.append(a)
            nj.append(b)
            ok += (spins(q) == spins(rf)).all()
        out["lsoda"].append((np.mean(nf), np.mean(nj), ok / M))
        fmt = lambda c: f"h={c[0]:<6g}{c[1]:>6d}" if c[0] else "      > grid     "
        print(f"    {r:>12g} | {fmt(cell['euler']):>16s} | {fmt(cell['rk4']):>16s} | {fmt(cell['strang']):>18s} | "
              f"{np.mean(nf):7.0f} / {np.mean(nj):4.0f}  ({100 * ok / M:.0f}%)")
    return out


def experiment3b_wallclock(rng, quick, h_rk4, h_euler, h_split=0.05, r=100):
    print(f"\n[3b] Wall-clock per run at strong injection (r = {r}), growing graph size")
    print(f"    splitting h={h_split} (a conservative choice; 3a says coarser works), Euler h={h_euler}, RK4 h={h_rk4} (steps found in 3a);")
    print("    BDF/LSODA rtol=1e-3 with analytic Jacobians")
    sizes = [(100, 5), (400, 4)] if quick else [(100, 8), (400, 6), (1500, 6)]
    out = []
    print(f"    {'N':>5s} | {'method':<22s} {'ms/run':>9s} {'speedup':>8s} {'per-spin':>9s} {'identical':>10s} {'cut/ref':>8s}")
    for N, M in sizes:
        m = OIM(make_graph(N, rng))
        rp = Ramp(r * m.K * m.lam)
        ics = [rng.uniform(-PI, PI, N) for _ in range(M)]
        refs = [reference(m, rp, p) for p in ics]
        cut_ref = np.mean([m.cut(q) for q in refs])
        methods = [("exact-inj. Strang", lambda p: strang_exact(m, rp, p, h_split)[0]),
                   ("Euler", lambda p: euler(m, rp, p, h_euler)[0]),
                   ("RK4", lambda p: rk4(m, rp, p, h_rk4)[0]),
                   ("BDF (sparse Jacobian)", lambda p: bdf_sparse(m, rp, p))]
        if N <= 400:
            methods.append(("LSODA (dense Jacobian)", lambda p: lsoda(m, rp, p)[0]))
        rec = {}
        for name, fn in methods:
            t0 = time.perf_counter()
            res = [fn(p) for p in ics]
            dt = (time.perf_counter() - t0) / M * 1e3
            frac = np.mean([np.mean(spins(q) == spins(r_)) for q, r_ in zip(res, refs)])
            full = np.mean([(spins(q) == spins(r_)).all() for q, r_ in zip(res, refs)])
            cutr = np.mean([m.cut(q) for q in res]) / cut_ref
            rec[name] = (dt, frac, full, cutr)
        base = rec["exact-inj. Strang"][0]
        for name, (dt, frac, full, cutr) in rec.items():
            sp = "1.0x" if name.startswith("exact") else f"{dt / base:.1f}x slower"
            note = "   <- unfaithful: timing not comparable" if frac < 0.999 else ""
            print(f"    {N:>5d} | {name:<22s} {dt:9.1f} {sp:>8s} {frac:9.4f} {100 * full:9.0f}% {cutr:8.4f}{note}")
        out.append((N, rec))
    return out


def experiment3c_noise_bias(rng, quick):
    """Single oscillator dphi = -Ks sin 2phi dt + sigma dW.  Exact stationary law p ~ exp((Ks/sigma^2) cos 2phi)."""
    print("\n[3c] With noise: stationary-law bias of Euler-Maruyama vs exact-injection splitting")
    Ks, sigma, M = 1.0, 0.3, (30000 if quick else 100000)
    kappa = Ks / sigma ** 2
    D_exact = 1 - i1e(kappa) / i0e(kappa)               # <1 - cos 2phi>, exact
    ksh = [0.02, 0.05, 0.1, 0.25, 0.5, 1.0]
    out = dict(ksh=ksh, em=[], st=[])
    print(f"    exact <1-cos 2phi> = {D_exact:.5f}   (kappa={kappa:.1f});  relative error of the simulated value:")
    print(f"    {'Ks*h':>6s} {'Euler-Maruyama':>16s} {'exact-inj. splitting':>22s}")
    for v in ksh:
        h = v / Ks
        n = int(max(400, 12 / v))
        p = np.zeros(M)
        q = np.zeros(M)
        with np.errstate(all="ignore"):
            for _ in range(n):
                p = p - h * Ks * np.sin(2 * p) + sigma * math.sqrt(h) * rng.standard_normal(M)
        q = inject(q, Ks * h / 2)
        for i in range(n):
            q = q + sigma * math.sqrt(h) * rng.standard_normal(M)
            q = inject(q, Ks * h if i < n - 1 else Ks * h / 2)
        em = np.mean(1 - np.cos(2 * p)) / D_exact - 1
        st = np.mean(1 - np.cos(2 * q)) / D_exact - 1
        out["em"].append(em)
        out["st"].append(st)
        print(f"    {v:6.2f} {100 * em:+15.1f}% {100 * st:+21.1f}%")
    return out


# =========================================================================== #
# 5.  Verdict and figure                                                      #
# =========================================================================== #
def verdict(e0, e1, e2, e3a, e3b, e3c):
    L = []
    ok = lambda b: "NOT COMPARABLE" if b is None else ("CONFIRMED     " if b else "NOT CONFIRMED ")
    L.append(("Exactness (README flow, group law, detuned closed form)", e0["ok"], f"worst detuned-flow error {e0['worst_detuned']:.1e}"))
    L.append(("Coupling: zeroth-order error is O(rho)", abs(e1["slope0"] - 1) < 0.15, f"fitted exponent {e1['slope0']:.2f}"))
    L.append(("Coupling: closed-form first-order fix is O(rho^2)", abs(e1["slope1"] - 2) < 0.2,
              f"exponent {e1['slope1']:.2f}; {e1['gain']:.0f}x smaller error at rho=0.1"))
    L.append(("Noise: basin blur = exact committor, width sigma/sqrt(4 pi)", e2["basin"]["max_z"] < 4.0,
              f"max |z| = {e2['basin']['max_z']:.2f} over {2 * len(e2['basin']['ms'])} points"))
    L.append(("Noise: unresolved fraction saturates at the von Mises floor", abs(e2["floor"]["z_late"]) < 4.0,
              f"late-time z = {e2['floor']['z_late']:+.2f}"))
    L.append(("Noise: resolution time = README t* with eps -> |eps+eta|", e2["restime"]["max_rel_err"] < 0.10,
              f"max quantile error {100 * e2['restime']['max_rel_err']:.1f}%"))
    # benefit: baseline = cheapest of the explicit fixed-step schemes (Euler and RK4 share the same stability limit)
    rs = e3a["r"]
    ratios = []
    for r, ev_e, ev_r, ev_s in zip(rs, e3a["euler"], e3a["rk4"], e3a["strang"]):
        best = min([v for v in (ev_e, ev_r) if v], default=None)
        ratios.append((r, best / ev_s if (best and ev_s) else None))
    win = next((r for i, (r, _) in enumerate(ratios) if all(y is not None and y >= 2 for _, y in ratios[i:])), None)
    top = ratios[-1]
    L.append(("Benefit: >=2x fewer coupling evals than the BEST explicit scheme", win is not None,
              f"{top[1]:.0f}x at r={top[0]}; sustained from r={win if win else 'never'}  (ratios: "
              + ", ".join(f"r={r:g}:{y:.1f}x" for r, y in ratios if y) + ")"))
    lo = ratios[0]
    L.append(("Benefit also at weak injection (r ~ 1)?", lo[1] is not None and lo[1] > 1.0, f"best-explicit/splitting = {lo[1]:.2f} at r={lo[0]}"))
    FAITHFUL = 0.999       # per-spin agreement with the reference needed for a timing comparison to be meaningful
    N, rec = e3b[-1]
    base = rec["exact-inj. Strang"][0]
    L.append((f"Benefit: wall-clock vs RK4 / BDF at N={N}, r=100",
              rec["RK4"][0] / base > 3 and rec["BDF (sparse Jacobian)"][0] / base > 1,
              f"{rec['RK4'][0] / base:.0f}x vs RK4, {rec['BDF (sparse Jacobian)'][0] / base:.0f}x vs BDF"))
    for Ne, rece in e3b:
        eu = rece["Euler"]
        b = rece["exact-inj. Strang"][0]
        if eu[1] >= FAITHFUL:
            L.append((f"Benefit: wall-clock vs Euler at N={Ne}, r=100", eu[0] / b > 1.0, f"splitting is {eu[0] / b:.1f}x faster"))
        else:
            L.append((f"Benefit: wall-clock vs Euler at N={Ne}, r=100", None,
                      f"Euler (step tuned at N=100) is unfaithful here: per-spin {eu[1]:.3f}, cut/ref {eu[3]:.3f}"))
    for N0, rec0 in e3b:
        if "LSODA (dense Jacobian)" in rec0:
            ls = rec0["LSODA (dense Jacobian)"][0] / rec0["exact-inj. Strang"][0]
            L.append((f"Benefit: wall-clock vs LSODA at N={N0}, r=100", ls > 1.0,
                      f"splitting is {ls:.1f}x {'faster' if ls > 1 else 'SLOWER (LSODA wins)'}"))
    i = e3c["ksh"].index(0.25)
    L.append(("Noise: splitting beats Euler-Maruyama on stationary-law bias at Ks*h=0.25",
              abs(e3c["st"][i]) < abs(e3c["em"][i]), f"{100 * e3c['st'][i]:+.1f}% vs {100 * e3c['em'][i]:+.1f}%"))
    j = e3c["ksh"].index(1.0)
    L.append(("Noise: splitting keeps correct thermal width at Ks*h=1", abs(e3c["st"][j]) < 0.1, f"{100 * e3c['st'][j]:+.0f}% error"))
    print("\n" + "=" * 100)
    print("VERDICT (computed from the measurements above)")
    print("=" * 100)
    for name, flag, note in L:
        print(f"  [{ok(flag)}] {name:<68s} {note}")
    print("=" * 100)
    return L


def make_figure(path, e1, e2, e3a, e3b, e3c, lines):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 4, figsize=(21, 9.5))
    a = ax.ravel()

    rho = np.array(e1["rho"])
    a[0].loglog(rho, e1["e0"], "o-", label="ignore coupling (exact injection only)")
    a[0].loglog(rho, e1["e1"], "s-", label="+ closed-form 1st-order correction")
    a[0].loglog(rho, 0.07 * rho, "k:", lw=1)
    a[0].loglog(rho, 0.015 * rho ** 2, "k:", lw=1)
    a[0].set(xlabel=r"$\rho=K\lambda_{max}/K_s$", ylabel="median phase error", title="(1) coupling: degradation is perturbative")
    a[0].legend(fontsize=8)

    ms = e2["basin"]["ms"]
    xx = np.linspace(-3, 3, 200)
    s25 = 0.25 / math.sqrt(4 * PI)
    a[1].plot(xx, norm.cdf(xx), "k--", lw=1, label="Gaussian (linearised)")
    a[1].plot(xx, [committor(0.5 + v * s25, 0.25) for v in xx], "-", color="C1", label=r"exact committor, $\sigma=0.25$")
    a[1].step([-3, 0, 3], [0, 1, 1], where="post", color="gray", lw=1, label="noise-free: README (a step)")
    for r, mk in zip(e2["basin"]["rows"], ("o", "s")):
        a[1].plot(ms, r[2], mk, label=rf"simulation, $\sigma={r[0]}$")
    a[1].set(xlabel=r"start offset from saddle, $m=(x_0-1/2)/s$,  $s=\sigma/\sqrt{4\pi}$", ylabel="P(settles in the right-hand integer)",
             title="(2a) noise blurs the basin boundary")
    a[1].legend(fontsize=8)

    f = e2["floor"]
    a[2].semilogy(f["t"], f["sim"], "o-", label="noisy simulation")
    a[2].semilogy(f["t"], f["det"], "k--", label="README formula (noise-free)")
    a[2].axhline(f["floor"], color="C3", ls=":", label="exact stationary floor")
    a[2].set(xlabel="t", ylabel=r"unresolved fraction  $P(|x-\mathrm{round}\,x|>\delta)$", title=r"(2b) $\sigma=0.1,\ \delta=0.1$")
    a[2].legend(fontsize=8)

    e, qs, qp, det_t, samples = e2["restime"]["rows"][0]
    srt = np.sort(samples[~np.isnan(samples)])
    a[3].plot(srt, np.arange(1, srt.size + 1) / srt.size, label="noisy simulation")
    a[3].axvline(det_t, color="k", ls="--", label=f"noise-free t* = {det_t:.2f}")
    a[3].axvline(qp[1], color="C1", ls=":", label="predicted median")
    a[3].set(xlabel="time to come within delta of an integer", ylabel="CDF", xlim=(0, 1.6), title=r"(2c) start at $\varepsilon=10^{-6}$ from the saddle")
    a[3].legend(fontsize=8)

    rs = e3a["r"]
    for key, lab, mk in (("euler", "Euler", "v"), ("rk4", "RK4", "o"), ("strang", "exact-injection Strang", "s")):
        a[4].loglog(rs, [v if v else np.nan for v in e3a[key]], mk + "-", label=lab)
    a[4].loglog(rs, [v[0] for v in e3a["lsoda"]], "k^--", label="LSODA RHS calls (+ Jacobians)")
    a[4].set(xlabel=r"injection stiffness $r=K_{s,max}/(K\lambda_{max})$", ylabel="coupling evaluations for 95% identical spins",
             title="(3a) where the exact flow pays off")
    a[4].legend(fontsize=8)

    Ns = [n for n, _ in e3b]
    for name, mk in (("exact-inj. Strang", "s"), ("Euler", "v"), ("RK4", "o"), ("BDF (sparse Jacobian)", "d"), ("LSODA (dense Jacobian)", "^")):
        pts = [(n, rec[name][0]) for n, rec in e3b if name in rec and rec[name][1] >= 0.999]   # unfaithful runs are not timing-comparable
        if pts:
            a[5].loglog(*zip(*pts), mk + "-", label=name)
    a[5].set(xlabel="number of oscillators N", ylabel="ms per run", title="(3b) wall-clock at strong injection (r = 100)")
    a[5].legend(fontsize=8)

    a[6].semilogx(e3c["ksh"], [100 * v for v in e3c["em"]], "o-", label="Euler-Maruyama")
    a[6].semilogx(e3c["ksh"], [100 * v for v in e3c["st"]], "s-", label="exact-injection splitting")
    a[6].axhline(0, color="k", lw=0.8)
    a[6].set(xlabel=r"$K_s h$", ylabel="error of thermal width (%)", ylim=(-100, 100), title="(3c) with noise the advantage is bias, not step size")
    a[6].legend(fontsize=8)

    a[7].axis("off")
    import textwrap
    txt = "\n".join(textwrap.fill(f"{'~' if fl is None else ('+' if fl else '-')} {nm}: {note}", width=88,
                                  subsequent_indent="    ") for nm, fl, note in lines[1:])
    a[7].text(0, 1, "Measured verdicts  (+ confirmed, - not confirmed, ~ not comparable)\n\n" + txt, va="top", fontsize=6.6, family="monospace")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    print(f"\nfigure saved to {path}")


# =========================================================================== #
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true", help="smaller samples / fewer sizes (about 30 s)")
    ap.add_argument("--no-plot", action="store_true", help="skip the figure")
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    t0 = time.time()
    print(f"OIM exact-injection study   (seed={args.seed}, {'quick' if args.quick else 'full'} mode)")

    e0 = experiment0_exactness(rng)
    e1 = experiment1_coupling(rng, args.quick)
    e2 = experiment2_noise(rng, args.quick)
    e3a = experiment3a_cost_vs_stiffness(rng, args.quick)
    h_rk4 = e3a["h"].get(("rk4", 100)) or 0.005
    h_euler = e3a["h"].get(("euler", 100)) or 0.005
    e3b = experiment3b_wallclock(rng, args.quick, h_rk4, h_euler)
    e3c = experiment3c_noise_bias(rng, args.quick)
    lines = verdict(e0, e1, e2, e3a, e3b, e3c)
    if not args.no_plot:
        try:
            make_figure(os.path.join(args.outdir, "oim_exact_flow_study.png"), e1, e2, e3a, e3b, e3c, lines)
        except ImportError:
            print("matplotlib not installed: figure skipped")
    print(f"total time {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
