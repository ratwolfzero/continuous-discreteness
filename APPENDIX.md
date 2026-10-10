# Appendix to "Discreteness as the Late-Time Limit of a Smooth Flow"

This appendix accompanies the main note, [README.md](README.md). It contains the proofs of Results 1–4, the properties that the main line of argument does not need, the geometric picture, details on the resolution time, related work, and a precise account of what has and has not been verified.

Notation is that of the note: $k = e^{-2\pi t}$, $R_t$ is the time-$t$ map of $\dot x = -\sin 2\pi x$, $n$ is an integer, and $D$ is the denominator of Section 2.2. Positions and times are real unless stated otherwise.

## Contents

- **A.** The global closed form: positivity, the algebraic core, agreement with the cell form
- **B.** Proofs of Results 1–4
- **C.** Further properties: C.1 fixed-point slopes, C.2 dissipation, C.3 monotone approach, C.4 resolution bounds
- **D.** The limit is not uniform
- **E.** Geometric reading: a hyperbolic dilation of the circle
- **F.** Resolution time: details and numerical check
- **G.** Related work, possible uses, references
- **H.** Scope and verification

---

## A. The global closed form

### A.1 Positivity and bounds

Besides $D$, define the second denominator

$$
E = (1+\cos 2\pi x) + k^{2}\,(1-\cos 2\pi x).
$$

In half-angle form,

$$
D = 2\big(\cos^{2}\pi x + k\sin^{2}\pi x\big), \qquad E = 2\big(\cos^{2}\pi x + k^{2}\sin^{2}\pi x\big).
$$

Both are strictly positive for every real $x$ and every $k>0$, because $\cos^2\pi x$ and $\sin^2\pi x$ are non-negative and never vanish together. For $0<k\le1$, $E=(1+k^2)+(1-k^2)\cos 2\pi x$ is affine in $\cos 2\pi x$ with non-negative slope, so it lies between its values at $\cos 2\pi x=-1$ and $+1$:

$$
2k^{2} \le E \le 2 \qquad (0<k\le1).
$$

With this notation the derivative in Result 2 is $\partial_x R_t = 2k/E$.

### A.2 The algebraic core

$$
D^{2} + \big((1-k)\sin 2\pi x\big)^{2} = 2E.
$$

To verify it, substitute $\sin^{2}2\pi x = (1-\cos 2\pi x)(1+\cos 2\pi x)$ and expand. The left side becomes

$$
(1+\cos 2\pi x)^{2} + (1-\cos 2\pi x)(1+\cos 2\pi x)(1+k^{2}) + k^{2}(1-\cos 2\pi x)^{2},
$$

which factors as

$$
2(1+\cos 2\pi x) + 2k^{2}(1-\cos 2\pi x) = 2E.
$$

Equivalently, with $w=(1-k)\sin 2\pi x/D$, one has $1+w^2 = 2E/D^2$. This is what makes the direct differentiation of the closed form of Section 2.2 collapse to $2k/E$ on the whole line, including the half-integers where the tangent is undefined.

### A.3 The cell form and the closed form agree

Let $x$ lie in the open cell around $n$ and set $T=\tan\pi(x-n)$. Then

$$
\sin 2\pi x = \frac{2T}{1+T^{2}}, \qquad
1+\cos 2\pi x = \frac{2}{1+T^{2}}, \qquad
1-\cos 2\pi x = \frac{2T^{2}}{1+T^{2}},
$$

so that $D = 2(1+kT^2)/(1+T^2)$, $E = 2(1+k^2T^2)/(1+T^2)$ and

$$
q := \frac{(1-k)\sin 2\pi x}{D} = \frac{(1-k)\,T}{1+kT^{2}}.
$$

A direct computation gives $T-q = kT(1+T^2)/(1+kT^2)$ and $1+qT = (1+T^2)/(1+kT^2)$, hence

$$
\tan\!\big(\pi(x-n) - \arctan q\big) = \frac{T-q}{1+qT} = kT.
$$

The angle $\pi(x-n)-\arctan q$ is the closed form's $\pi(R_t(x)-n)$. Treat $k$ as a continuous parameter starting from the value one, where $q=0$ and the angle is $\pi(x-n)\in(-\tfrac\pi2,\tfrac\pi2)$. Its tangent $kT$ is finite for every $k$, so the angle can never reach $\pm\tfrac\pi2$ and stays strictly inside $(-\tfrac\pi2,\tfrac\pi2)$. Hence it equals $\arctan(kT)$, and the closed form and the cell form coincide. The closed form has the advantage of being valid on the cell boundaries too.

*A branch-free alternative.* Both forms, as functions of $t$ for a fixed start $x$, are solutions of $\dot y=-\sin 2\pi y$ with $y(0)=x$: the cell form by the differentiation in B.1, the closed form by direct differentiation, which uses the algebraic core A.2 (and is checked symbolically for all $x$, half-integers included, in `prove_paper.py`). The right-hand side is globally Lipschitz, so by uniqueness (B.3) the two coincide for all real $t$.

---

## B. Proofs of Results 1–4

**B.1 Result 1.** Inside an open cell, with $T=\tan\pi(x-n)$ and $\partial_t k = -2\pi k$, differentiating the cell form gives

$$
\frac{\partial R_t}{\partial t} = \frac{1}{\pi}\cdot\frac{T\,\partial_t k}{1+k^{2}T^{2}} = -\frac{2kT}{1+k^{2}T^{2}}.
$$

On the other hand, the double-angle formula applied to the arctangent gives

$$
\sin\!\big(2\pi R_t\big) = \sin\!\big(2\arctan(kT)\big) = \frac{2kT}{1+k^{2}T^{2}},
$$

which proves the equation of motion inside the cell. At $t=0$, $k=1$ and $R_0(x) = n + \frac1\pi\arctan\tan\pi(x-n) = x$. On the cell boundaries every point is fixed (the numerator of the arctangent argument in the closed form vanishes there), so the equation of motion holds trivially. Using the expressions of A.3 for $E$ and $\sin 2\pi x$, the sine of the transformed angle has the closed form

$$
\sin\!\big(2\pi R_t(x)\big) = \frac{2kT}{1+k^{2}T^{2}} = \frac{2k\,\sin 2\pi x}{E},
$$

valid on every open cell; at a half-integer both sides vanish (the point is fixed and $\sin 2\pi x=0$), so it holds on the whole line. For the linearizing form, the cell form gives $\tan(\pi(R_t(x)-n)) = kT$, and the tangent is periodic with period one, so $n$ drops out: $\tan(\pi R_t(x)) = k\tan(\pi x)$.

**B.2 Result 2.** Inside a cell,

$$
\frac{\partial R_t}{\partial x} = \frac{k\sec^{2}\pi(x-n)}{1+k^{2}\tan^{2}\pi(x-n)} = \frac{k}{\cos^{2}\pi(x-n)+k^{2}\sin^{2}\pi(x-n)} = \frac{2k}{E},
$$

using $\cos^{2}\pi(x-n) = (1+\cos 2\pi x)/2$ and $\sin^{2}\pi(x-n) = (1-\cos 2\pi x)/2$. The formula holds on the cell boundaries by continuity: both sides of $\partial_xR_t=2k/E$ are continuous on the real line (the left side because the closed form is smooth, the right side because $E>0$), and they agree on every open cell. Alternatively, differentiate the closed form directly and use the algebraic core of A.2. Positivity of the derivative follows from $E>0$. Shift-equivariance follows from the periodicity of sine and cosine. A continuous, strictly increasing map that commutes with integer shifts is onto the real line, so $R_t$ is a bijection; it is smooth because the closed form is, and its inverse is $R_{-t}$ by Result 3.

**B.3 Result 3.** Fix the time shift $t$ and the starting point $x$, and consider the two functions of the running time

$$
\sigma \mapsto R_{\sigma+t}(x), \qquad \sigma \mapsto R_{\sigma}\big(R_t(x)\big).
$$

Both solve the same ordinary differential equation $\dot y = -\sin(2\pi y)$, and both take the value $R_t(x)$ at $\sigma=0$. The right-hand side is globally Lipschitz, because its derivative is bounded in absolute value by $2\pi$. By uniqueness for Lipschitz ordinary differential equations (Picard–Lindelöf, or Grönwall's inequality), the two functions coincide for all real $\sigma$. Taking the second time equal to minus the first, and using $R_0=\mathrm{id}$, shows that each map is inverted by the map at the opposite time.

**B.4 Result 4.** For a point in an open cell the tangent $T$ is finite and the contraction factor tends to zero. The arctangent is continuous at zero, so

$$
R_t(x) = n + \frac{1}{\pi}\arctan(kT) \;\longrightarrow\; n .
$$

The statement about the supremum is proved in Appendix D.

---

## C. Further properties

**C.1 Fixed points and their slopes.** Integers and half-integers are fixed at every time, with exact slopes

$$
\partial_x R_t(n) = e^{-2\pi t}, \qquad \partial_x R_t\!\left(n+\tfrac12\right) = e^{+2\pi t}.
$$

Integers contract and half-integers expand, with exactly the rates predicted by the linearization of Section 1.

*Proof.* At an integer the cosine equals one, so $E=2$ and $2k/E = k = e^{-2\pi t}$. At a half-integer the cosine equals minus one, so $E=2k^{2}$ and $2k/E = 1/k = e^{+2\pi t}$. Both kinds of point are fixed because $\sin 2\pi x=0$ there, so the numerator of the arctangent argument vanishes. $\square$

**C.2 Dissipation.** The potential is a Lyapunov function:

$$
\frac{d}{dt}\, V\big(R_t(x)\big) = -\sin^{2}\!\big(2\pi R_t(x)\big) \le 0.
$$

*Proof.* By the chain rule and Result 1, $V'(R_t)\,\partial_t R_t = \sin(2\pi R_t)\cdot\big(-\sin 2\pi R_t\big)$. $\square$

**C.3 Monotone approach to the integers.** The distance to the nearest integer never increases along the flow:

$$
|R_t(x) - n| \le |R_s(x) - n|, \qquad s \le t \ \text{(real, of either sign)}, \quad x \in \left[n-\tfrac12,\; n+\tfrac12\right].
$$

Each closed cell is mapped into itself.

*Proof.* The map is increasing and fixes every half-integer, so at every time it maps each closed cell into itself. On the open cell the cell form gives $|R_t(x)-n| = \frac1\pi\arctan(k\,|T|)$ with $T=\tan\pi(x-n)$. This is nondecreasing in $k$, because $\partial_k\arctan(kT)=T/(1+k^{2}T^{2})$ has the sign of $T$. Since $k=e^{-2\pi t}$ decreases with $t$, the distance is nonincreasing in $t$, for all real times. On the cell boundary the distance is constantly $\tfrac12$. $\square$

**C.4 Resolution bounds.** For non-negative time, no two points are brought together faster than the contraction factor allows:

$$
e^{-2\pi t}\,|x-y| \;\le\; |R_t(x) - R_t(y)| \;\le\; e^{2\pi t}\,|x-y|, \qquad t \ge 0.
$$

*Proof.* For $t\ge0$ the contraction factor lies in $(0,1]$, so A.1 gives $2k^{2} \le E \le 2$. Therefore the slope satisfies

$$
e^{-2\pi t} = k \;\le\; \frac{2k}{E} \;\le\; \frac{1}{k} = e^{2\pi t},
$$

and the mean value theorem gives the stated inequalities. $\square$

---

## D. The limit is not uniform

Result 4 states that the limit holds pointwise on the open cells and uniformly only on compact subsets of them. Because the limit map is discontinuous, uniform convergence on the whole line is impossible: for every finite $t$,

$$
\sup_x \big|R_t(x)-\mathrm{round}(x)\big| = \tfrac{1}{2}.
$$

*Proof.* $R_t$ is continuous with $R_t(n+\tfrac12)=n+\tfrac12$, so $R_t(x)\to n+\tfrac12$ as $x\uparrow n+\tfrac12$, while $\mathrm{round}(x)=n$ on the open cell; the gap therefore tends to $\tfrac12$. $\square$

**The resolved and unresolved fractions.** For $0<\delta<\tfrac12$, the fraction of a uniformly distributed start lying within $\delta$ of an integer at time $t$ is

$$
2R_{-t}(\delta) = \frac{2}{\pi}\arctan\!\frac{\tan\pi\delta}{k}, \qquad k = e^{-2\pi t},
$$

so the unresolved fraction is exactly

$$
1 - 2R_{-t}(\delta) = \frac{2}{\pi}\arctan\frac{k}{\tan\pi\delta}.
$$

*Proof.* Since $R_t$ is increasing and fixes $n$, $|R_t(x)-n|\le\delta$ exactly when $|x-n|\le R_{-t}(\delta)$ (Result 3), so each unit cell contributes length $2R_{-t}(\delta)$, and the cell form gives $R_{-t}(\delta)=\frac1\pi\arctan(\tan(\pi\delta)/k)$. The unresolved fraction follows from $\arctan a+\arctan(1/a)=\pi/2$. $\square$

**Asymptotics.** If $k \ll \delta \ll 1$, equivalently $t \gg \frac{1}{2\pi}\ln\frac1\delta$, the unresolved fraction is

$$
\approx \frac{2\,e^{-2\pi t}}{\pi^{2}\delta},
$$

with relative error close to $\frac13\big[(k/\pi\delta)^{2} + (\pi\delta)^{2}\big]$. This follows from $\arctan z\sim z$ applied to $z=k/\tan\pi\delta$ (which requires $z\to0$, i.e. $k/\delta\to0$) together with $\tan z\sim z$. Smallness of $k$ and of $\delta$ separately is not enough: the approximation requires $k/\tan\pi\delta$ to be small, and for $\delta \lesssim k$ it fails completely (the exact fraction then tends to $1$, not to $0$). In every case the unresolved fraction decays exponentially in $t$ but is positive at every finite time.

---

## E. Geometric reading: a hyperbolic dilation of the circle

Map the line to the unit circle by $z = e^{2\pi i x}$. Then the transformed point is

$$
e^{2\pi i R_t(x)} = \frac{z+r}{1+rz}, \qquad r = \tanh(\pi t).
$$

This is a Möbius transformation of the unit circle that fixes the points plus one and minus one. Its multiplier at plus one is

$$
\frac{1-r}{1+r} = e^{-2\pi t},
$$

so plus one is the attracting point (the image of the integers) and minus one the repelling point (the image of the half-integers). Composition of two such maps corresponds to multiplying the matrices

$$
\begin{pmatrix} 1 & r_s \\ r_s & 1 \end{pmatrix}
\begin{pmatrix} 1 & r_t \\ r_t & 1 \end{pmatrix},
$$

which yields the velocity-addition law

$$
r_{s+t} = \frac{r_s + r_t}{1 + r_s r_t}.
$$

So the group law is **relativistic velocity addition** in rapidity form. The hyperbolic distance in the Poincaré disc (curvature $-1$) from the origin to the point $r$ is

$$
d_{\mathbb{H}}(0,r) = 2\,\mathrm{artanh}(r) = 2\pi t,
$$

so time is, up to a constant, hyperbolic distance travelled. Extended to the unit disc, the map sends the origin to the point $r$. The uniform distribution on the circle is the harmonic measure seen from the origin, so its push-forward is the harmonic measure seen from $r$; its density with respect to $d\theta/2\pi$, at the image angle $\theta = 2\pi R_t(x)$, is the Poisson kernel

$$
P_r(\theta) = \frac{1-r^{2}}{1-2r\cos\theta+r^{2}}.
$$

(Equivalently, $1/\partial_xR_t = E/(2k) = P_r(2\pi R_t(x))$.)

On the circle there is a single attracting point, $+1$, and a single repelling point, $-1$. The integers are the fibre of the covering map $x \mapsto e^{2\pi i x}$ over $+1$, and $R_t$ is the lift to the line of this circle flow. This restates where the lattice comes from (the period of the potential) rather than deriving it. In the limit, every point of the circle except $-1$ is carried to $+1$.

---

## F. Resolution time: details and numerical check

This appendix supplements Section 4 of the note.

**Accuracy of the small-scale form.** Because $\tan z>z$ on $(0,\pi/2)$, the estimate $\frac{1}{2\pi}\ln\frac{1}{\pi^{2}\varepsilon\delta}$ always overshoots the exact time. Expanding $\ln(\tan z/z)=z^2/3+O(z^4)$ gives

$$
t_\star^{\text{approx}} - t_\star \;\approx\; \frac{\pi}{6}\,\big(\varepsilon^{2}+\delta^{2}\big).
$$

**Numerical check.** Exact values from the closed form (confirmed by a direct root-find of the flow) against the small-scale estimate:

| $\varepsilon$ | $\delta$ | exact $t_\star$ | estimate | estimate − exact |
| --- | --- | --- | --- | --- |
| $10^{-3}$ | $10^{-3}$ | 1.834427 | 1.834428 | $1.0\times10^{-6}$ |
| $10^{-6}$ | $10^{-4}$ | 3.300299 | 3.300299 | $5.2\times10^{-9}$ |
| $10^{-9}$ | $10^{-6}$ | 5.132638 | 5.132638 | $5.2\times10^{-13}$ |
| $10^{-2}$ | $10^{-1}$ | 0.729612 | 0.735025 | $5.4\times10^{-3}$ |
| $0.2$ | $0.2$ | 0.101687 | 0.147921 | $4.6\times10^{-2}$ |

The first three rows lie deep in the small-scale regime, where the estimate is accurate to the digits shown. The last two show it degrading as $\varepsilon$ and $\delta$ grow, which is why the exact formula is preferable outside that regime.

**Link to Appendix D.** Fix a time $t$ and invert the formula for $t_\star$: the starts that are still unresolved at that time are exactly those within

$$
\varepsilon_c(t) = \frac1\pi\arctan\!\frac{k}{\tan\pi\delta}
$$

of a half-integer, that is, $t_\star(\varepsilon_c(t),\delta)=t$. Two such ends per cell give the unresolved fraction $1-2R_{-t}(\delta)=2\varepsilon_c(t)$ of Appendix D. Since $\varepsilon_c(t)>0$ at every finite time, no finite time resolves every point.

---

## G. Related work, possible uses, references

**Related work.** The ingredients are known. After rescaling phase and time, the equation of motion is the Adler equation with zero detuning (an overdamped pendulum without applied torque), and its tangent half-angle solution is classical [1]. That the flow of a phase forced at the first harmonic is a Möbius map of the circle is the single-oscillator case of the Möbius-group structure behind the Watanabe–Strogatz reduction [2, 3], and the invariance of the Poisson-kernel family under Möbius maps underlies the Ott–Antonsen ansatz [4]; the relation between that ansatz and the Watanabe–Strogatz reduction is worked out in [6]. Maps of the form arctangent-of-scaled-tangent also appear as "staircase" families in the literature, though no specific source has been traced here. Periodic regularizers have been used to push neural-network weights towards a grid; second-harmonic locking (called subharmonic injection locking in [5]) is used to binarize oscillator Ising machines; and continuous-collapse models describe measurement outcomes as limits of smooth dynamics (these are stochastic, unlike the deterministic flow here). What the note adds is the packaging: the exact group structure, the slope statements, the non-uniformity of the limit made quantitative, and the reading of time as a resolution parameter. This is based on a web search only, not on a systematic literature review, and novelty has not been established.

**Possible uses (speculative).**

- **Quantization-aware training.** Annealing weights onto an integer grid with a schedule that composes exactly, in place of ad hoc sharpness parameters. The exact composition holds for the regularizer flow alone; combined with a loss gradient the two flows do not commute, and the weights nearest a half-integer are exactly those that stay unresolved longest (Section 4 of the note).
- **Oscillator-based Ising machines.** The phase equation of a single oscillator under second-harmonic injection locking is this equation after rescaling the phase and the time, so the closed form gives the exact binarization trajectory in the absence of coupling [5]. If the injection strength $a(t)$ enters as a prefactor, the same closed form holds with $t$ replaced by $\int a(t)\,dt$. With Ising couplings present the dynamics is no longer this flow.
- **Sine-Gordon, Josephson and Frenkel–Kontorova systems.** A solvable toy model for relaxation into quantized states.
- **Measurement and decoherence.** A cartoon in which an outcome is the infinite-time limit of smooth dynamics. The flow is deterministic: the outcome is fixed by the initial cell and there is no Born rule or probability, so this illustrates the limiting mechanism only.
- **Coarse-graining.** A concrete example of a group action whose boundary limit is non-invertible, as in renormalization-style pictures.

**References.**

1. R. Adler, "A study of locking phenomena in oscillators", Proc. IRE 34, 351–357 (1946).
2. S. Watanabe and S. H. Strogatz, "Integrability of a globally coupled oscillator array", Phys. Rev. Lett. 70, 2391 (1993).
3. S. A. Marvel, R. E. Mirollo and S. H. Strogatz, "Identical phase oscillators with global sinusoidal coupling evolve by Möbius group action", Chaos 19, 043104 (2009).
4. E. Ott and T. M. Antonsen, "Low dimensional behavior of large systems of globally coupled oscillators", Chaos 18, 037113 (2008).
5. T. Wang and J. Roychowdhury, "OIM: Oscillator-based Ising machines for solving combinatorial optimisation problems", UCNC 2019, LNCS 11493, pp. 232–256 (2019); arXiv:1903.07163.
6. A. Pikovsky and M. Rosenblum, "Partially integrable dynamics of ensembles of nonidentical oscillators", arXiv:1001.1299.

---

## H. Scope and verification

The labels printed by the three scripts are the section labels of this appendix and the note ("Res1" to "Res4" are Results 1 to 4 of the note), so every check can be traced to the statement it supports.

**Numerical verification (`verify_paper.py`).** All identities stated in the note and the appendix were cross-checked in 50-digit arithmetic on random positions and times of either sign. The checks cover the group law and inverse, the equation of motion, the spatial derivative and the sine formula, the linearizing identity, the slopes at integers and half-integers, the Lyapunov derivative, the Möbius form together with the velocity-addition law, the Poisson-kernel push-forward, the hyperbolic distance, the resolved-fraction formula and the non-uniform limit (Appendix D), the exact resolution time (Section 4) together with its small-scale estimate, the leading-order error of that estimate and the unresolved radius $\varepsilon_c(t)$ (Appendix F), and, on random samples, the resolution bounds (C.4) and the monotone approach (C.3; the latter for times of either sign, over cells $n=-5,\dots,5$). For the exact identities the residuals were below $10^{-30}$ (worst observed in the default run about $4\times10^{-39}$); no violations were found in the inequalities; and the exact resolution time agrees with a direct root-find of the flow to better than $10^{-40}$.

**Cross-check of `discreteness.py`.** The script also compares the double-precision map in `discreteness.py` against the same reference: absolute agreement of order $10^{-16}$ on random samples; in a dedicated sweep of starts between $10^{-16}$ and $10^{-1}$ from a half-integer, on both sides and at times from $t=-3$ to $t=6$, the absolute error stays below $10^{-14}$ (measured: $3\times10^{-16}$) and the relative error of the distance to the integer in the cell $n=0$ below $10^{-13}$ (measured: $2\times10^{-15}$). For $|t|$ up to $1000$ of either sign the map stays finite, keeps every start in its closed cell and leaves integers and half-integers exactly fixed. The bundled `discreteness.py --selftest` is a smaller double-precision check restricted to non-negative times.

**Exact proof of the algebra (`prove_paper.py`).** The numerical checks above support the statements but are not proofs. The algebraic content is proved exactly by `prove_paper.py` (SymPy): the positivity and bounds on $D$ and $E$ (A.1), the algebraic core (A.2), the tangent-coordinate algebra and the branch-free argument of A.3, $\partial_x R_t = 2k/E$ on the whole line including the half-integers, $\partial_t R_t = -\sin(2\pi R_t)$ with $R_0=\mathrm{id}$, the sine formula (B), the slopes (C.1), the linearizing identity (Result 1), the Möbius form with the multiplier, velocity-addition law and Poisson density (E), and the algebra of Section 4 and Appendices D and F. Each identity is shown to hold for every real position and every positive contraction factor, by polynomial division modulo $\sin^2+\cos^2=1$, and a set of deliberately wrong formulas is confirmed to be rejected.

**What is not mechanised.** Three standard theorems that the proofs apply directly: Picard–Lindelöf uniqueness (B.3), the fact that a continuous, strictly increasing map commuting with integer shifts is onto (B.2), and the mean value theorem (C.4). The argument is therefore a hand proof whose algebraic steps are machine-checked in computer algebra; it has not been formalised in a proof assistant.

**Not addressed.** A literature search beyond the informal one mentioned in Appendix G, and the choice of which physical or computational setting, if any, this model fits best.

**Where each statement is checked.**

| Statement | `verify_paper.py` label | `prove_paper.py` label |
| --- | --- | --- |
| A.1 positivity and bounds | A.1 | A.1 |
| A.2 algebraic core | A.2 | A.2 |
| A.3 cell form = closed form | A.3 | A.3 |
| Result 1 (equation of motion, linearizing form) | Res1, B | B, Res1 |
| Result 2 (derivative, shift equivariance) | Res2 | B |
| Result 3 (group law, inverse) | Res3 | standard theorem T1 |
| Result 4 (limit; non-uniformity) | D | Res4, D |
| C.1 fixed points and slopes | C.1 | C.1 |
| C.2 dissipation | C.2 | C.2 |
| C.3 monotone approach | C.3 | C.3 |
| C.4 resolution bounds | C.4 | A.1 (bound), standard theorem T3 |
| D resolved and unresolved fractions | D | D |
| E Möbius form, velocity addition, Poisson kernel | E | E |
| Section 4 exact resolution time | 4 | 4 |
| F small-scale estimate, table, $\varepsilon_c(t)$ | F | F |
