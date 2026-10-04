#!/usr/bin/env python3
"""
discreteness.py  --  "Discreteness as the late-time limit of a smooth flow"

An intuitive, runnable visualisation of the gradient flow

        dx/dt = -sin(2 pi x),        V(x) = -cos(2 pi x) / (2 pi),

whose stable fixed points are the integers and whose unstable ones are the
half-integers.  The flow is solved exactly: inside the cell around the integer n,

        R_t(x) = n + (1/pi) * arctan( e^{-2 pi t} * tan(pi (x - n)) ).

At every finite time R_t is a smooth, strictly increasing bijection of the line
(and R_{s+t} = R_s o R_t).  As t -> infinity it becomes "round to the nearest
integer".  Discreteness is the late-time limit of a perfectly smooth process.

What you see (five linked panels, one shared time t)
----------------------------------------------------
  1  Marbles on the potential V(x): they roll downhill into the integer wells.
  2  The same marbles wrapped on a circle z = exp(2 pi i x): the integers are ONE
     point (z = +1), the half-integers are the opposite point (z = -1).
  3  The map x -> R_t(x): starts as the identity, sharpens into a staircase.
  4  Space-time: trajectories funnel into the integers; points started extremely
     close to a half-integer linger, then peel off (logarithmic delay).
  5  Distance to the nearest integer on a log axis: every curve ends up parallel
     to exp(-2 pi t).  Circles mark the predicted resolution time
     t* = ln(1 / (pi^2 * eps * delta)) / (2 pi).

Usage
-----
  python discreteness.py                 interactive explorer (slider + Play)
  python discreteness.py --filmstrip     save a still "filmstrip" PNG
  python discreteness.py --gif           save an animated GIF
  python discreteness.py --snapshot 0.4  save one dashboard frame at t = 0.4
  python discreteness.py --selftest      numerically verify the maths
  python discreteness.py --all           selftest + filmstrip + snapshot + gif

Options: --tmax 3.0  --marbles 150  --outdir .  --dpi 110  --frames 70  --fps 20
Keys in the interactive window: space = play/pause, left/right = step.

Requirements: Python >= 3.9, numpy, matplotlib (>= 3.5).  GIF export uses Pillow,
which is installed together with matplotlib.
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

TWO_PI = 2.0 * math.pi

# --------------------------------------------------------------------------- #
# 1. The mathematics                                                          #
# --------------------------------------------------------------------------- #


def potential(x):
    """V(x) = -cos(2 pi x) / (2 pi)."""
    return -np.cos(TWO_PI * np.asarray(x, dtype=float)) / TWO_PI


def flow(t, x):
    """Exact time-t map R_t(x) of dx/dt = -sin(2 pi x).  Broadcasts over t and x.

    Cell form with the half-angle written as atan2, which is numerically stable
    for large t.  Exact half-integers are the unstable fixed points: they are
    returned unchanged (floating-point sin(pi) is not exactly zero, so they must
    be treated explicitly).
    """
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    k = np.exp(-TWO_PI * t)                      # contraction factor
    n = np.round(x)
    y = x - n                                    # offset inside the cell, |y| <= 1/2
    out = n + np.arctan2(k * np.sin(np.pi * y), np.cos(np.pi * y)) / np.pi
    return np.where(np.abs(y) == 0.5, x, out)


def flow_closed_form(t, x):
    """The global closed form  x - arctan((1-k) sin 2pi x / D) / pi  (for checks)."""
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    k = np.exp(-TWO_PI * t)
    c, s = np.cos(TWO_PI * x), np.sin(TWO_PI * x)
    D = (1.0 + c) + k * (1.0 - c)
    return x - np.arctan((1.0 - k) * s / D) / np.pi


def _rk4(x0, T, n=4000):
    """Independent check: integrate dx/dt = -sin(2 pi x) with classic RK4."""
    x = np.array(x0, dtype=float)
    h = np.asarray(T, dtype=float) / n
    def f(u): return -np.sin(TWO_PI * u)
    for _ in range(n):
        k1 = f(x)
        k2 = f(x + 0.5 * h * k1)
        k3 = f(x + 0.5 * h * k2)
        k4 = f(x + h * k3)
        x = x + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x


def selftest(verbose=True):
    """Check the exact solution against independent computations."""
    rng = np.random.default_rng(0)
    N = 600
    x = rng.uniform(-3, 3, N)
    t = rng.uniform(0.0, 1.5, N)
    s = rng.uniform(0.0, 1.5, N)

    results = {}
    results["cell form  vs  global closed form"] = np.max(
        np.abs(flow(t, x) - flow_closed_form(t, x)))
    results["flow  vs  RK4 integration of the ODE"] = np.max(
        np.abs(flow(t, x) - _rk4(x, t)))
    results["group law  R_{s+t} = R_s o R_t"] = np.max(
        np.abs(flow(s + t, x) - flow(s, flow(t, x))))

    z = np.exp(1j * TWO_PI * x)
    r = np.tanh(np.pi * t)
    results["Moebius form  (z+r)/(1+rz)"] = np.max(
        np.abs(np.exp(1j * TWO_PI * flow(t, x)) - (z + r) / (1 + r * z)))

    # slopes by central differences; kept to moderate t where the map is not yet a near-step
    h = 1e-6
    tm = rng.uniform(0.0, 0.6, N)
    k = np.exp(-TWO_PI * tm)
    E = (1 + np.cos(TWO_PI * x)) + k**2 * (1 - np.cos(TWO_PI * x))
    slope_fd = (flow(tm, x + h) - flow(tm, x - h)) / (2 * h)
    results["slope  dR/dx = 2k/E  (relative)"] = np.max(
        np.abs(slope_fd / (2 * k / E) - 1))

    n = np.arange(-3, 4, dtype=float)
    tt = 0.7
    sl_int = (flow(tt, n + h) - flow(tt, n - h)) / (2 * h)
    sl_half = (flow(tt, n + 0.5 + h) - flow(tt, n + 0.5 - h)) / (2 * h)
    results["slope at integers  = e^{-2 pi t}"] = np.max(
        np.abs(sl_int - math.exp(-TWO_PI * tt)))
    results["slope at half-integers = e^{+2 pi t}"] = np.max(
        np.abs(sl_half / math.exp(TWO_PI * tt) - 1))

    xs = np.sort(rng.uniform(-2, 2, 4000))
    results["monotone (min increment, want > 0)"] = - \
        min(0.0, np.min(np.diff(flow(1.3, xs))))

    ok = all(v < 1e-6 for v in results.values())
    if verbose:
        print("self-test of the exact solution")
        print("-" * 58)
        for name, v in results.items():
            print(f"  {name:<42s} {v:9.2e}")
        print("-" * 58)
        print("  PASS" if ok else "  FAIL")
    return ok


# --------------------------------------------------------------------------- #
# 2. Shared styling and helpers                                               #
# --------------------------------------------------------------------------- #

X_LO, X_HI = -1.5, 1.5                          # displayed window: three cells
INK = "#1d1d1f"
GREY = "#8d99ae"
RED = "#c1121f"
CELL_COLORS = ["#2a9d8f", "#e76f51", "#457b9d"]  # teal, orange, blue (cycled)


def cell_color(n):
    return CELL_COLORS[int(n) % 3]


def make_marbles(n):
    """n starting points, uniformly spread over the window, never exactly on a
    half-integer.  Returns positions and the integer each one is destined for."""
    x0 = X_LO + (np.arange(n) + 0.5) * (X_HI - X_LO) / n
    return x0, np.round(x0).astype(int)


def shade_cells(ax, orient):
    """Light background tint per basin of attraction (cell around each integer)."""
    for n in (-1, 0, 1):
        span = ax.axvspan if orient == "v" else ax.axhspan
        span(n - 0.5, n + 0.5, color=cell_color(n), alpha=0.09, lw=0, zorder=0)


def snapped_fraction(xt, tol=0.02):
    return float(np.mean(np.abs(xt - np.round(xt)) < tol))


def style_axes(ax, title):
    ax.set_title(title, loc="left", fontsize=10.5,
                 fontweight="bold", color=INK, pad=6)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=8, colors="#444")


# --------------------------------------------------------------------------- #
# 3. Panels (each has update(t))                                              #
# --------------------------------------------------------------------------- #


class LandscapePanel:
    """Marbles rolling down the potential."""

    def __init__(self, ax, x0, dest, compact=False):
        self.ax, self.x0 = ax, x0
        xs = np.linspace(X_LO, X_HI, 801)
        shade_cells(ax, "v")
        ax.plot(xs, potential(xs), color=INK, lw=2.2, zorder=2)
        for n in (-1, 0, 1):
            ax.plot([n], [potential(n)], "o", ms=8 if compact else 10,
                    mfc="white", mec=INK, mew=1.6, zorder=3)
        for h in (-0.5, 0.5):
            ax.plot([h], [potential(h)], "^",
                    ms=8 if compact else 10, color=RED, zorder=3)
            ax.axvline(h, color=RED, lw=0.8, ls=(
                0, (3, 3)), alpha=0.5, zorder=1)
        self.sc = ax.scatter(
            x0, potential(x0), s=14 if compact else 30, c=[cell_color(n) for n in dest],
            edgecolors="white", linewidths=0.5, zorder=5,
        )
        ax.set_xlim(X_LO, X_HI)
        ax.set_ylim(-0.215, 0.235)
        ax.set_yticks([])
        ax.set_xticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
        ax.set_xticklabels(["", "-1", "-½", "0", "½", "1", ""])
        if not compact:
            ax.set_xlabel("position x", fontsize=9)
            ax.text(0.0, -0.205, "stable: integers", ha="center",
                    va="bottom", fontsize=8, color=INK)
            ax.text(0.5, 0.185, "unstable: half-integers",
                    ha="center", va="bottom", fontsize=8, color=RED)
            ax.text(-0.5, 0.185, "unstable", ha="center",
                    va="bottom", fontsize=8, color=RED)
        for sp in ("left", "top", "right"):
            ax.spines[sp].set_visible(False)

    def update(self, t):
        xt = flow(t, self.x0)
        self.sc.set_offsets(np.column_stack([xt, potential(xt)]))


class CirclePanel:
    """The same marbles on the circle z = exp(2 pi i x): one ring per cell."""

    def __init__(self, ax, x0, dest):
        self.x0 = x0
        self.radius = 1.0 + 0.085 * (dest + 1)
        th = np.linspace(0, TWO_PI, 400)
        for n in (-1, 0, 1):
            rr = 1.0 + 0.085 * (n + 1)
            ax.plot(rr * np.cos(th), rr * np.sin(th),
                    color=cell_color(n), lw=0.9, alpha=0.35, zorder=1)
        # the two special points, marked just outside the rings so marbles never hide them
        ax.plot([1.31], [0.0], "o", ms=9, mfc="white",
                mec=INK, mew=1.8, zorder=3, clip_on=False)
        ax.plot([-1.31], [0.0], "^", ms=9, color=RED, zorder=3, clip_on=False)
        ax.text(1.31, -0.17, "+1", ha="center",
                va="top", fontsize=8.5, color=INK)
        ax.text(-1.31, -0.17, "-1", ha="center",
                va="top", fontsize=8.5, color=RED)
        ax.text(0.0, 0.10, r"$z=e^{2\pi i x}$",
                ha="center", va="center", fontsize=12, color=INK)
        ax.text(0.0, -0.17, "one ring per cell", ha="center",
                va="center", fontsize=7.5, color="#555")
        ax.text(0.0, -1.58, "○  z = +1 : every integer (attracting)",
                ha="center", va="top", fontsize=8, color=INK)
        ax.text(0.0, -1.80, "▲  z = −1 : every half-integer (repelling)",
                ha="center", va="top", fontsize=8, color=RED)
        self.sc = ax.scatter(
            self.radius * np.cos(TWO_PI * x0), self.radius * np.sin(TWO_PI * x0), s=22,
            c=[cell_color(n) for n in dest], edgecolors="white", linewidths=0.4, zorder=5,
        )
        ax.set_aspect("equal")
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(-2.0, 1.4)
        ax.axis("off")

    def update(self, t):
        ang = TWO_PI * flow(t, self.x0)
        self.sc.set_offsets(np.column_stack(
            [self.radius * np.cos(ang), self.radius * np.sin(ang)]))


class MapPanel:
    """The map x -> R_t(x): identity at t = 0, staircase as t -> infinity."""

    def __init__(self, ax, compact=False, ghost_times=(0.05, 0.1, 0.2, 0.4, 0.8, 1.6)):
        N = 1500
        self.xs = X_LO + (np.arange(N) + 0.5) * (X_HI - X_LO) / N
        shade_cells(ax, "v")
        ax.plot([X_LO, X_HI], [X_LO, X_HI], color=GREY, lw=1, ls=":", zorder=1)
        for n in (-1, 0, 1):
            ax.hlines(n, n - 0.5, n + 0.5, colors=INK,
                      linestyles=(0, (4, 3)), lw=1.3, zorder=2)
        for g in ghost_times:
            ax.plot(self.xs, flow(g, self.xs), color=GREY,
                    lw=0.9, alpha=0.45, zorder=2)
        (self.line,) = ax.plot(self.xs, self.xs, color=RED, lw=2.4, zorder=4)
        for n in (-1, 0, 1):
            ax.plot([n], [n], "o", ms=6, mfc="white",
                    mec=INK, mew=1.3, zorder=5)
        for h in (-0.5, 0.5):
            ax.plot([h], [h], "^", ms=6, color=RED, zorder=5)
        ax.set_xlim(X_LO, X_HI)
        ax.set_ylim(X_LO, X_HI)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([-1, 0, 1])
        ax.set_yticks([-1, 0, 1])
        if not compact:
            ax.set_xlabel("starting point x", fontsize=9)
            ax.set_ylabel("position after time t:  R_t(x)", fontsize=9)
            ax.text(
                0.98, 0.03, "", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=8, color=INK, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.8),
            )
            self.note = ax.texts[-1]
            ax.text(0.03, 0.97, "dotted: t = 0   grey: earlier t\ndashed: t → ∞ (rounding)",
                    transform=ax.transAxes, ha="left", va="top", fontsize=7.5, color="#555")
        else:
            self.note = None

    def update(self, t):
        self.line.set_ydata(flow(t, self.xs))
        if self.note is not None:
            k = math.exp(-TWO_PI * t)
            self.note.set_text(
                f"slope at integers  {k:.3g}\nslope at half-integers  {1 / k:.3g}")


class SpaceTimePanel:
    """Trajectories x(t): basins funnel into the integers."""

    EPS_EXP = (2, 4, 6, 8)

    def __init__(self, ax, tmax, x0, dest):
        self.tmax, self.x0 = tmax, x0
        tg = np.linspace(0.0, tmax, 700)
        shade_cells(ax, "h")
        for h in (-0.5, 0.5):
            ax.axhline(h, color=RED, lw=1, ls=(0, (4, 3)), alpha=0.8, zorder=1)
        for n in (-1, 0, 1):
            ax.axhline(n, color=cell_color(n), lw=0.8, alpha=0.6, zorder=1)
        base, base_dest = make_marbles(40)
        traj = flow(tg[None, :], base[:, None])
        for xi, di in zip(traj, base_dest):
            ax.plot(tg, xi, color=cell_color(di), lw=0.8, alpha=0.55, zorder=2)
        for e in self.EPS_EXP:                         # starts hugging the unstable point
            for sgn in (-1, +1):
                xs0 = 0.5 + sgn * 10.0 ** (-e)
                ax.plot(tg, flow(tg, xs0), color=cell_color(
                    round(xs0)), lw=1.7, alpha=0.95, zorder=3)
        self.vline = ax.axvline(0.0, color=INK, lw=1.4, zorder=6)
        self.sc = ax.scatter(
            np.zeros_like(x0), x0, s=9, c=[cell_color(n) for n in dest],
            edgecolors="white", linewidths=0.3, zorder=7,
        )
        ax.set_xlim(0.0, tmax)
        ax.set_ylim(X_LO, X_HI)
        ax.set_yticks([-1.5, -1, -0.5, 0, 0.5, 1, 1.5])
        ax.set_yticklabels(["", "-1", "-½", "0", "½", "1", ""])
        ax.set_xlabel(
            "time t      (dashed red: unstable half-integers)", fontsize=9)

    def update(self, t):
        self.vline.set_xdata([t, t])
        self.sc.set_offsets(np.column_stack(
            [np.full_like(self.x0, t), flow(t, self.x0)]))


class DecayPanel:
    """Distance to the nearest integer vs time (log axis)."""

    GENERIC = (0.45, 0.30, 0.15)
    EPS_EXP = (2, 4, 6, 8)
    DELTA = 0.1

    def __init__(self, ax, tmax, plt):
        tg = np.linspace(0.0, tmax, 900)
        starts, labels, colors = [], [], []
        for g in self.GENERIC:
            starts.append(g)
            labels.append(f"x₀ = {g:.2f}")
            colors.append("#6c757d")
        cmap = plt.cm.plasma(np.linspace(0.05, 0.78, len(self.EPS_EXP)))
        for e, c in zip(self.EPS_EXP, cmap):
            starts.append(0.5 - 10.0 ** (-e))
            labels.append(f"x₀ = ½ − 1e-{e}")
            colors.append(c)
        self.starts = np.array(starts)
        for x, lab, c in zip(starts, labels, colors):
            ax.plot(tg, flow(tg, x), color=c, lw=1.6, label=lab, zorder=3)
        ax.plot(tg, np.exp(-TWO_PI * tg) / np.pi, color=INK, lw=1.1, ls="--", zorder=2,
                label=r"$\propto e^{-2\pi t}$")
        ax.axhline(self.DELTA, color=GREY, lw=0.9, ls=":", zorder=1)
        ax.text(tmax * 0.985, self.DELTA * 1.25,
                f"δ = {self.DELTA}", ha="right", va="bottom", fontsize=8, color="#555")
        first = True
        for e in self.EPS_EXP:                          # predicted resolution times
            ts = math.log(1.0 / (math.pi**2 * 10.0 **
                          (-e) * self.DELTA)) / TWO_PI
            if ts <= tmax:
                ax.plot([ts], [self.DELTA], "o", ms=7, mfc="none", mec=INK, mew=1.3, zorder=6,
                        label="predicted t*" if first else None)
                first = False
        self.vline = ax.axvline(0.0, color=INK, lw=1.4, zorder=5)
        self.sc = ax.scatter(np.zeros(len(starts)), self.starts, s=26, c=colors, edgecolors="white",
                             linewidths=0.6, zorder=7)
        ax.set_yscale("log")
        floor = 10.0 ** (math.floor(math.log10(math.exp(-TWO_PI * tmax) / math.pi)) - 1)
        ax.set_ylim(floor, 1.0)
        ax.set_xlim(0.0, tmax)
        ax.set_xlabel("time t", fontsize=9)
        ax.set_ylabel("distance to nearest integer", fontsize=9)
        ax.legend(loc="lower left", fontsize=6.8, frameon=True, framealpha=0.9, ncol=1, borderpad=0.4,
                  labelspacing=0.25, handlelength=1.6)

    def update(self, t):
        self.vline.set_xdata([t, t])
        self.sc.set_offsets(np.column_stack(
            [np.full(len(self.starts), t), flow(t, self.starts)]))


# --------------------------------------------------------------------------- #
# 4. The dashboard (interactive / GIF / snapshot)                             #
# --------------------------------------------------------------------------- #


class Dashboard:
    def __init__(self, plt, tmax=3.0, n_marbles=150, interactive=False, figsize=(15.0, 8.4)):
        from matplotlib.gridspec import GridSpec

        self.plt, self.tmax = plt, tmax
        self.x0, self.dest = make_marbles(n_marbles)
        self.fig = plt.figure(figsize=figsize, facecolor="white")
        gs = GridSpec(
            2, 3, figure=self.fig, height_ratios=[1.0, 1.1], wspace=0.26, hspace=0.36,
            left=0.055, right=0.985, top=0.855, bottom=0.17 if interactive else 0.075,
        )
        axA = self.fig.add_subplot(gs[0, :2])
        axB = self.fig.add_subplot(gs[0, 2])
        axC = self.fig.add_subplot(gs[1, 0])
        axD = self.fig.add_subplot(gs[1, 1])
        axE = self.fig.add_subplot(gs[1, 2])

        style_axes(axA, "1 · Marbles roll downhill on  V(x) = −cos(2πx)/2π")
        style_axes(axB, "2 · The same flow wrapped on a circle")
        style_axes(axC, "3 · The map R_t: identity → staircase")
        style_axes(axD, "4 · Space–time: basins funnel into integers")
        style_axes(axE, "5 · Distance to the nearest integer")

        self.panels = [
            LandscapePanel(axA, self.x0, self.dest),
            CirclePanel(axB, self.x0, self.dest),
            MapPanel(axC),
            SpaceTimePanel(axD, tmax, self.x0, self.dest),
            DecayPanel(axE, tmax, plt),
        ]
        self.fig.text(0.055, 0.955, "Discreteness as the late-time limit of a smooth flow",
                      fontsize=16, fontweight="bold", color=INK, ha="left", va="center")
        self.fig.text(
            0.055, 0.915,
            r"$\dot x=-\sin(2\pi x)$:  integers are stable, half-integers unstable.   "
            r"Exact solution  $R_t(x)=n+\frac{1}{\pi}\arctan(e^{-2\pi t}\tan\pi(x-n))$",
            fontsize=10.5, color="#444", ha="left", va="center",
        )
        self.clock = self.fig.text(0.985, 0.955, "", fontsize=13, family="monospace", color=INK,
                                   ha="right", va="center")
        self.stat = self.fig.text(
            0.985, 0.915, "", fontsize=10, color="#444", ha="right", va="center")
        self.update(0.0)

    def update(self, t):
        t = float(t)
        for p in self.panels:
            p.update(t)
        k = math.exp(-TWO_PI * t)
        self.clock.set_text(f"t = {t:5.3f}   k = e^(-2πt) = {k:.3g}")
        frac = snapped_fraction(flow(t, self.x0))
        self.stat.set_text(
            f"marbles within 0.02 of an integer: {100 * frac:3.0f}%")

    # -- interactive ------------------------------------------------------- #
    def run_interactive(self):
        from matplotlib.widgets import Button, Slider

        fig, tmax = self.fig, self.tmax
        slider = Slider(fig.add_axes([0.10, 0.065, 0.62, 0.03]), "time t", 0.0, tmax, valinit=0.0,
                        color="#e76f51")
        button = Button(fig.add_axes([0.76, 0.05, 0.075, 0.06]), "Play")
        state = {"playing": False}
        timer = fig.canvas.new_timer(interval=40)
        step = tmax / 160.0

        def on_slide(val):
            self.update(val)
            fig.canvas.draw_idle()

        def stop():
            timer.stop()
            state["playing"] = False
            button.label.set_text("Play")

        def tick():
            v = slider.val + step
            if v >= tmax:
                v = tmax
                stop()
            slider.set_val(v)

        def toggle(_event=None):
            if state["playing"]:
                stop()
            else:
                if slider.val >= tmax - 1e-9:
                    slider.set_val(0.0)
                state["playing"] = True
                button.label.set_text("Pause")
                timer.start()

        def on_key(event):
            if event.key == " ":
                toggle()
            elif event.key in ("left", "right"):
                stop()
                d = (-1 if event.key == "left" else 1) * tmax / 100.0
                slider.set_val(min(tmax, max(0.0, slider.val + d)))

        slider.on_changed(on_slide)
        button.on_clicked(toggle)
        timer.add_callback(tick)
        fig.canvas.mpl_connect("key_press_event", on_key)
        # prevent garbage collection
        self._keepalive = (slider, button, timer)
        self.plt.show()


# --------------------------------------------------------------------------- #
# 5. Static outputs                                                           #
# --------------------------------------------------------------------------- #


def make_snapshot(plt, path, t, tmax, n_marbles, dpi):
    dash = Dashboard(plt, tmax, n_marbles, interactive=False)
    dash.update(t)
    dash.fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(dash.fig)
    print(f"saved {path}")


def make_filmstrip(plt, path, tmax, n_marbles, dpi):
    times = sorted(
        {t for t in (0.0, 0.05, 0.15, 0.4, 1.0) if t < tmax} | {tmax})
    m = len(times)
    x0, dest = make_marbles(n_marbles)
    fig, axes = plt.subplots(
        2, m, figsize=(3.1 * m, 7.2), facecolor="white",
        gridspec_kw=dict(height_ratios=[1.0, 1.35], hspace=0.30, wspace=0.10,
                         left=0.06, right=0.99, top=0.82, bottom=0.07),
    )
    for j, t in enumerate(times):
        axL, axM = axes[0, j], axes[1, j]
        L = LandscapePanel(axL, x0, dest, compact=True)
        M = MapPanel(axM, compact=True)
        L.update(t)
        M.update(t)
        k = math.exp(-TWO_PI * t)
        frac = snapped_fraction(flow(t, x0))
        axL.set_title(f"t = {t:g}   (k = {k:.2g})\n{100 * frac:.0f}% within 0.02 of an integer",
                      fontsize=9.5, color=INK, pad=6)
        if j > 0:
            axM.tick_params(labelleft=False)
        axM.set_xlabel("start x", fontsize=8)
        for ax in (axL, axM):
            ax.tick_params(labelsize=8)
    axes[1, 0].set_ylabel("position after time t:  R_t(x)", fontsize=9)
    axes[0, 0].set_ylabel("marbles on V(x)", fontsize=9)
    fig.text(0.06, 0.945, "Discreteness as the late-time limit of a smooth flow", fontsize=16,
             fontweight="bold", color=INK, ha="left", va="center")
    fig.text(0.06, 0.905,
             r"$\dot x=-\sin(2\pi x)$.  Top: marbles roll into the integer wells (circles); half-integers (triangles) "
             r"are unstable.  Bottom: the smooth map $R_t$ sharpens into rounding (dashed).",
             fontsize=10, color="#444", ha="left", va="center")
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)
    print(f"saved {path}")


def make_gif(plt, path, tmax, n_marbles, frames, fps, dpi):
    from matplotlib.animation import PillowWriter

    dash = Dashboard(plt, tmax, n_marbles,
                     interactive=False, figsize=(13.5, 7.6))
    u = np.linspace(0.0, 1.0, frames)
    # linger on the fast early phase
    times = tmax * u**1.7
    sequence = [times[0]] * 8 + list(times) + [times[-1]] * 14
    writer = PillowWriter(fps=fps)
    print(f"rendering {len(sequence)} frames ...")
    with writer.saving(dash.fig, path, dpi):
        for i, t in enumerate(sequence):
            dash.update(t)
            writer.grab_frame()
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(sequence)}")
    plt.close(dash.fig)
    print(f"saved {path}")


# --------------------------------------------------------------------------- #
# 6. Command line                                                             #
# --------------------------------------------------------------------------- #


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Visualise discreteness as the late-time limit of the flow dx/dt = -sin(2 pi x).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--filmstrip", action="store_true",
                    help="save a still filmstrip PNG")
    ap.add_argument("--gif", action="store_true", help="save an animated GIF")
    ap.add_argument("--snapshot", type=float, metavar="T",
                    help="save one dashboard frame at time T")
    ap.add_argument("--selftest", action="store_true",
                    help="verify the maths numerically")
    ap.add_argument("--all", action="store_true",
                    help="selftest + filmstrip + snapshot + gif")
    ap.add_argument("--tmax", type=float, default=3.0,
                    help="final time (default 3.0)")
    ap.add_argument("--marbles", type=int, default=150,
                    help="number of marbles (default 150)")
    ap.add_argument("--outdir", default=".",
                    help="output directory (default .)")
    ap.add_argument("--dpi", type=int, default=110,
                    help="PNG dpi (default 110)")
    ap.add_argument("--frames", type=int, default=70,
                    help="GIF animation frames (default 70)")
    ap.add_argument("--fps", type=int, default=20,
                    help="GIF frames per second (default 20)")
    args = ap.parse_args(argv)

    if args.tmax <= 0.5 or args.marbles < 10:
        ap.error("--tmax must be > 0.5 and --marbles >= 10")

    saving = args.filmstrip or args.gif or args.snapshot is not None or args.selftest or args.all
    import matplotlib

    if saving:
        matplotlib.use("Agg")                              # no window needed
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 9, "axes.edgecolor": "#555", "axes.linewidth": 0.8})
    os.makedirs(args.outdir, exist_ok=True)
    def out(name): return os.path.join(args.outdir, name)

    if args.selftest or args.all:
        if not selftest():
            return 1
    if args.filmstrip or args.all:
        make_filmstrip(plt, out("discreteness_filmstrip.png"),
                       args.tmax, args.marbles, args.dpi)
    if args.snapshot is not None or args.all:
        t = args.snapshot if args.snapshot is not None else 0.4
        make_snapshot(plt, out("discreteness_dashboard.png"),
                      t, args.tmax, args.marbles, args.dpi)
    if args.gif or args.all:
        make_gif(plt, out("discreteness.gif"), args.tmax,
                 args.marbles, args.frames, args.fps, 72)
    if saving:
        return 0

    # ---- interactive explorer ------------------------------------------- #
    backend = matplotlib.get_backend().lower()
    headless = backend in {"agg", "pdf", "ps",
                           "svg", "cairo", "template", "pgf"}
    if headless:
        print("No interactive display available (matplotlib backend: %s)." % backend)
        print("Writing a filmstrip and a GIF instead; run with a display for the live explorer.")
        make_filmstrip(plt, out("discreteness_filmstrip.png"),
                       args.tmax, args.marbles, args.dpi)
        make_gif(plt, out("discreteness.gif"), args.tmax,
                 args.marbles, args.frames, args.fps, 72)
        return 0
    Dashboard(plt, args.tmax, args.marbles, interactive=True).run_interactive()
    return 0


if __name__ == "__main__":
    sys.exit(main())
