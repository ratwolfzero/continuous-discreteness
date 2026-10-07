# Exact injection flow for oscillator Ising machines: evaluation and verdict (AI-Assisted)

*Source: `oim_exact_flow_study.py` (final run, seed 2024; seeds 1 and 2 and several earlier runs used to check stability). Companion files: `oim_exact_flow_study.png` (figure) and `oim_exact_flow_results.xlsx` (all timing data).*

---

## 1. The verdict in plain words

The exact flow is correct, and it is useful in one specific situation.

- **How much does the clean picture degrade?** Gently and predictably. Turning on the coupling between oscillators adds an error that grows in direct proportion to the coupling strength, and a closed-form correction cuts that error by about 49 times. Noise does not destroy the picture; it replaces its sharp edges with a blur of known width.
- **Is there a real computational advantage?** Yes, but only when the injection ("snap to 0 or π") term is much stronger than the coupling. In that regime the exact flow lets a simulation take large time steps where standard methods must take tiny ones. Measured gain: roughly **5 to 10 times fewer coupling evaluations** than the best standard explicit method (10 to 25 times at the strongest setting on other random seeds).
- **Where is there no advantage?** When the injection is weak or moderate, which is where most OIM schedules operate. There, plain RK4 is as good or better.

---

## 2. The idea in plain language

An oscillator Ising machine is a network of oscillators. Two things act on each oscillator's phase $\varphi_i$:

1. **Coupling** pulls neighbouring phases together or apart, according to the problem (here: max-cut).
2. **Injection locking** pushes every phase toward $0$ or $\pi$. Those two values are the two spin states.

$$
\dot\varphi_i \;=\; -K\sum_j J_{ij}\,\sin(\varphi_i-\varphi_j)\;-\;K_s(t)\,\sin(2\varphi_i)\;\;[+\;\sigma\,\text{noise}]
$$

Your README solves the second term on its own, exactly. With $n$ the nearest multiple of $\pi$ to the starting phase and $S$ the accumulated injection strength:

$$
\varphi(t)\;=\;n\pi+\arctan\!\Big(e^{-2S}\,\tan(\varphi_0-n\pi)\Big),
\qquad
S=\int_0^t K_s(t')\,dt'
$$

Two properties matter in practice:

- **The schedule is free.** Any ramp $K_s(t)$ only enters through the single number $S$.
- **Steps combine exactly.** Applying the flow for $S_1$ and then for $S_2$ is the same as applying it once for $S_1+S_2$:

$$
\mathcal A_{S_2}\circ\mathcal A_{S_1}\;=\;\mathcal A_{S_1+S_2}
$$

The one number that decides everything below is the **injection stiffness**, how much stronger the injection is than the coupling ($\lambda_{\max}$ is the largest eigenvalue magnitude of the coupling matrix):

$$
r\;=\;\frac{K_{s,\max}}{K\,\lambda_{\max}}
$$

---

## 3. How much does the clean picture degrade?

### 3.1 Non-factorable coupling

The coupling cannot be folded into the exact solution. The reason is mathematical, not a matter of the schedule. For fixed neighbours the coupling varies like $\sin\varphi$ while the injection varies like $\sin 2\varphi$. Combining them produces ever higher harmonics, so no finite closed form exists:

$$
\big[\sin\varphi\,\partial_\varphi,\;\sin 2\varphi\,\partial_\varphi\big]
=\Big(-\tfrac32\sin\varphi+\tfrac12\sin 3\varphi\Big)\partial_\varphi
$$

**What the measurements show.** Let $\rho = K\lambda_{\max}/K_s$ be the coupling-to-injection ratio.

- If you ignore the coupling and just apply the exact injection flow, the error grows linearly with $\rho$.
- If you add a first-order correction built from the flow and its exact derivative, the error grows only with $\rho^2$:

$$
\varphi(t)\;\approx\;\mathcal A_{S(t)}(\varphi_0)\;+\;\int_0^t \mathcal A'_{S(t)-S(s)}\!\big(\mathcal A_{S(s)}(\varphi_0)\big)\;F\!\big(\mathcal A_{S(s)}(\varphi_0)\big)\,ds
$$

$$
\text{error}_{0}\;\propto\;\rho^{1.00},
\qquad
\text{error}_{1}\;\propto\;\rho^{2.00}
$$

At $\rho=0.1$ the correction makes the error about **49 times smaller**, and the fitted exponents were identical on all three seeds. A small fraction of oscillators (about 1% at $\rho=0.4$) end up in the opposite basin from the one the exact flow predicts. I report only their size, because there are too few events to fit a scaling law.

Oscillators very close to the unstable point (within 0.35 rad) were excluded from the error statistics, because the expansion is not uniform there. This is the same non-uniform behaviour your README already describes.

### 3.2 A real-world imperfection that costs nothing: frequency mismatch

Not every imperfection breaks exactness. A constant frequency offset $\Delta$ keeps the problem exactly solvable, because the equation stays in the same family of circle maps (Möbius maps, written with a $2\times2$ matrix):

$$
\dot\varphi=\Delta-K_s\sin 2\varphi
\qquad\Longrightarrow\qquad
\varphi_s=\tfrac12\arcsin\!\frac{\Delta}{K_s}\ \text{(stable)},\quad
\tfrac{\pi}{2}-\varphi_s\ \text{(unstable)}
$$

Locking exists only while $|\Delta|<K_s$. I derived the closed form (it is not in the README) and checked it against high-accuracy integration to about $10^{-12}$ for offsets up to $0.95\,K_s$.

### 3.3 Noise

In the README's units, one oscillator with noise obeys:

$$
dx=-a\sin(2\pi x)\,dt+\sigma\,dW
$$

Noise does not break the picture. It swaps sharp features for blurs of known size.

**The basin boundary becomes a smooth S-curve.** The width of the blur follows from the README's saddle growth rate $e^{2\pi a t}$:

$$
s=\frac{\sigma}{\sqrt{4\pi a}}
$$

The exact probability of settling in the right-hand integer, starting from $x$ between $0$ and $1$, is a 1-D diffusion result:

$$
q(x)=\frac{\displaystyle\int_0^{x}e^{-a\cos(2\pi y)/(\pi\sigma^{2})}\,dy}{\displaystyle\int_0^{1}e^{-a\cos(2\pi y)/(\pi\sigma^{2})}\,dy}
$$

Simulation matched it with $|z|\le 1.9$ over 14 test points. A plain Gaussian approximation is off by up to 0.9%.

**Discreteness has a noise floor.** The fraction of phases not yet within $\delta$ of an integer stops decaying and settles at the exact stationary tail of

$$
p(x)\;\propto\;\exp\!\Big(\frac{a}{\pi\sigma^{2}}\cos 2\pi x\Big)
$$

The simulated late-time value matched this floor ($z=-1.60$). In the crossover region the noise-free README formula underestimates the unresolved fraction by up to about 30%.

**Waiting time is capped.** Without noise, the time to resolve a phase that starts a distance $\varepsilon$ from the unstable point grows without limit as $\varepsilon\to0$:

$$
t^{*}(\varepsilon,\delta)=\frac{1}{2\pi a}\,\ln\frac{1}{\tan(\pi\delta)\,\tan(\pi\varepsilon)}
$$

With noise, replace $\varepsilon$ by $|\varepsilon+\eta|$ with $\eta\sim\mathcal N(0,s^{2})$. This predicted the simulated quantiles to within 6.4%. For a start at $\varepsilon=10^{-6}$ the noise-free time is 2.20, while the noisy median is about 0.63.

---

## 4. Where the speed advantage comes from

Fast injection makes the equations **stiff**: explicit solvers must keep their time step small relative to $1/K_s$ to stay stable, even though nothing interesting happens on that time scale once the phases have snapped to their spins.

The exact flow removes that restriction. The scheme is a *splitting*. Each step applies half the injection exactly, then a standard second-order step (Heun) for the coupling $F$, then the other half of the injection exactly:

$$
\varphi_{n+1}
=\mathcal A_{S(t_n+h/2,\;t_{n+1})}\;\circ\;\mathcal B_h\;\circ\;\mathcal A_{S(t_n,\;t_n+h/2)}\,(\varphi_n)
$$

$$
\mathcal B_h:\;\;\varphi\;\mapsto\;\varphi+\frac h2\Big(F(\varphi)+F\big(\varphi+h\,F(\varphi)\big)\Big)
$$

Because of the group law, the two half-steps of neighbouring steps merge, so each step makes **one** call to the exact flow and **two** coupling evaluations. The time step is then limited only by how fast the coupling changes things, not by $K_s$.

**With noise,** the same idea works, with the noise kick inserted between the exact flows. Its accuracy is then limited by how well the discrete scheme reproduces the true thermal spread. For a single oscillator the exact reference is:

$$
\big\langle 1-\cos 2\varphi\big\rangle \;=\; 1-\frac{I_1(\kappa)}{I_0(\kappa)},
\qquad
\kappa=\frac{K_s}{\sigma^{2}}
$$

---

## 5. Measured results

### 5.1 Cost to get the right answer

Cost is the number of coupling evaluations needed for at least 95% of runs to reproduce the exact reference's final spin pattern ($N=100$ sparse graph, seed 2024). "Best explicit" is the cheaper of Euler and RK4.

| Injection stiffness $r$ | Best explicit | Exact-injection splitting | Splitting advantage |
| --- | --- | --- | --- |
| 1 | 120 (RK4) | 240 | 0.5× (worse) |
| 3 | 240 | 240 | 1.0× (equal) |
| 10 | 300 (Euler) | 600 | 0.5× (worse) |
| 30 | 1200 | 240 | **5×** |
| 100 | 1200 | 120 | **10×** |

Other seeds gave 10 to 25× at $r=100$ and 2.5 to 5× at $r=30$. At $r=10$ the result swung from 2× worse to 10× better depending on the seed, so treat $r\approx3$ to $10$ as a toss-up.

### 5.2 Wall-clock time at $r=100$

Splitting was faster by these ranges across all seven full runs. Timings vary by up to about 2× between runs of identical code, so quote ranges, not single values.

| Compared with | $N=100$ | $N=400$ | $N=1500$ |
| --- | --- | --- | --- |
| RK4 | 15 to 27× (all sizes) | | |
| Euler (where Euler is accurate) | 4.2 to 5.9× (all sizes) | | not comparable |
| BDF (sparse Jacobian) | 5 to 8× | 9 to 16× | **55 to 104×** |
| LSODA (dense Jacobian) | 1.2 to 1.7× | 2.6 to 5.6× | not run |

The "not comparable" entry is deliberate. Euler with the step size tuned at $N=100$ fails at $N=1500$ (about half the spins wrong, cut value 19% below the reference), so timing it there would be misleading.

### 5.3 Quality of the answer

In the runs I examined, per-spin agreement with the reference was at least 99.9% and the cut value was within about 0.05%. Exact agreement of the *whole* spin pattern is a stricter test: it was 100% at $N=100$ and 67 to 100% at $N\ge400$ depending on seed.

### 5.4 With noise: bias, not step size

Error in the thermal width of a noisy oscillator:

| $K_s h$ | Euler–Maruyama | Exact-injection splitting |
| --- | --- | --- |
| 0.1 | about +10% | about 0% |
| 0.25 | +32.8% | −2.5% |
| 0.5 | +88.5% | −11.4% |
| 1.0 | about +450% | −36% |

Splitting has about ten times less bias at moderate steps and never blows up. However, at $K_s h=1$ it underestimates the thermal spread by 36%. With noise, the exact flow improves accuracy and stability but does **not** buy large steps.

---

## 6. Advantages and disadvantages against established methods

### 6.1 Advantages

| Compared with | What the exact-injection method does better |
| --- | --- |
| **Explicit Euler / RK4** (the usual choice for OIM simulation) | No step-size ceiling from the injection strength. 5 to 10× fewer coupling evaluations at $r\ge30$, and 4 to 27× less wall-clock time. The injection step can never push a phase across a basin boundary, so large steps do not cause spurious spin flips. |
| **Euler–Maruyama** (noisy simulation) | About ten times less bias in the thermal spread at moderate steps. Stays bounded where Euler–Maruyama diverges ($K_s h\gtrsim0.5$). |
| **Implicit / stiff solvers** (BDF, LSODA) | Matrix-free: cost per step grows with the number of edges, with no Jacobians and no linear solves. The advantage over sparse BDF grows with graph size (55 to 104× at $N=1500$). |
| **Any numerical method** | Gives closed-form predictions that simulation alone does not: the error law in $\rho$, the noise-blur width $s$, the noise floor, and a schedule-design rule for how long the injection must act. |

### 6.2 Disadvantages

| Issue | Detail |
| --- | --- |
| **Narrow regime** | Real advantage only from about $r\gtrsim30$. Typical OIM schedules ramp $K_s$ to only a few times $K\lambda_{\max}$ ($r\approx1$ to $10$), where RK4 or Euler is equal or better. |
| **Not exact once coupling is on** | The coupling step is still ordinary numerical integration (second-order Heun). Accuracy in the early part of the ramp is limited by it, not by the exact flow. |
| **Noise limits step size** | Needs about $K_s h\lesssim0.25$ to keep the thermal width within a few percent. At $K_s h=1$ the error is −36%. |
| **Large-N exactness** | Whole-pattern agreement with the reference drops with $N$ (67 to 100% at $N\ge400$), though per-spin agreement stays above 99.9%. |
| **Competition at small $N$** | A dense-Jacobian LSODA is within 1.2 to 1.7× of splitting at $N=100$, so there the gain is small. |
| **Not new technique** | Splitting methods with exactly solvable parts, and the exact solution of the injection equation, are classical. What is specific here is the combination and the measurements. |
| **Fragile to model changes** | Exactness relies on a purely sinusoidal injection term. Real oscillators with extra harmonics in their response would leave the exactly solvable family. |

---

## 7. When is the advantage significant?

You understood correctly: in specific cases the advantage is large. Those cases combine four conditions:

1. **Strong injection:** $r\gtrsim30$, for example a long final "pinning" or binarization phase.
2. **Large sparse networks:** where each coupling evaluation is the main cost and Jacobians are expensive.
3. **A fixed-step explicit baseline:** which is the usual way OIMs are simulated.
4. **Low or no noise,** or noise handled with $K_s h\lesssim0.25$.

| Your situation | Recommended approach |
| --- | --- |
| $r\lesssim3$ | RK4 (or an adaptive solver). Exact injection gives no gain. |
| $3\lesssim r\lesssim30$ | Toss-up; test both. |
| $r\gtrsim30$, large sparse graph | **Exact-injection splitting.** |
| $r\gtrsim30$, small graph ($N\lesssim100$) | Splitting or LSODA; they are close. |
| Noise present | Splitting, but keep $K_s h\lesssim0.25$. |
| Predicting noise effects or binarization time | Use the closed-form formulas in section 3.3. |

---

## 8. What the script's 16 checks say

The script tests 16 hypotheses against its own measurements. 13 are confirmed, 2 are not confirmed, and 1 is not comparable. "Confirmed" means the result passed the script's threshold; some thresholds are my own choices.

| Status | Plain-language meaning | Measured |
| --- | --- | --- |
| Confirmed | The exact formulas are right (match numerical integration, steps combine exactly, detuned form works). | error about $10^{-12}$ |
| Confirmed | Ignoring coupling gives an error proportional to $\rho$. | exponent 1.00 |
| Confirmed | The first-order fix gives an error proportional to $\rho^2$. | exponent 2.00; 49× smaller at $\rho=0.1$ |
| Confirmed | Noise blurs the basin boundary exactly as the formula says. | max $\lvert z\rvert=1.88$ |
| Confirmed | Unresolved fraction settles at the exact noise floor. | $z=-1.60$ |
| Confirmed | Noisy resolution time follows the README formula with a random $\varepsilon$. | 6.4% error |
| Confirmed | At least 2× fewer coupling evaluations than the best explicit scheme at strong injection. | 10× at $r=100$; from $r=30$ |
| **Not confirmed** | Any benefit at weak injection. | 0.5× at $r=1$ |
| Confirmed | Faster than RK4 and BDF at $N=1500$. | 26× and 97× (final run) |
| Confirmed | Faster than Euler at $N=100$ and $N=400$. | 4.2× and 5.9× |
| **Not comparable** | Faster than Euler at $N=1500$. Euler is inaccurate there. | per-spin agreement 0.50 |
| Confirmed | Faster than LSODA at $N=100$ and $N=400$. | 1.3× and 3.3× |
| Confirmed | Less bias than Euler–Maruyama at $K_s h=0.25$. | −2.5% vs +32.8% |
| **Not confirmed** | Correct thermal width at $K_s h=1$. | −36% |

---

## 9. What was not tested

- **Optimizer quality.** I measured whether the simulation reproduces the true trajectory, not whether the exact flow helps the machine *find better cuts*. The cut value matched the reference to about 0.05%, which says the method does no harm, not that it helps.
- **Other graph types:** only random sparse graphs, with antiferromagnetic couplings.
- **Matrix-free implicit solvers** (such as CVODE with a Krylov solver). These could narrow the gap at large $N$.
- **Hardware realism:** non-sinusoidal injection response, time-varying frequency offsets.
- **Quantization-aware training** and other uses of the same flow. The mechanism (large steps overshooting into a neighbouring cell) is plausibly the same, but I did not test it.
- **Novelty.** I cannot judge it. The ingredients are classical.

---

## 10. Reproducing the results

```txt
python oim_exact_flow_study.py            # full run, about 1.5 to 2 minutes
python oim_exact_flow_study.py --quick    # about 30 seconds
python oim_exact_flow_study.py --seed 1   # check another seed
```

Requires `numpy`, `scipy` and (for the figure) `matplotlib`. Counts of coupling evaluations are reproducible; wall-clock times depend on the machine and vary by up to about 2× between runs.
