# sympy_check: independent Python reproduction of the paper's analytic results

This directory re-derives, independently of the companion Mathematica notebook and of the
Black Hole Perturbation Toolkit, the analytic calculations and the numbers quoted in

> E. Barausse, *Gravitational radiation from a photon on the light ring: multipole fluxes and
> the logarithmic divergence* (equation numbers as in the submitted
> version of the paper).

Symbolic parts use **sympy**, numerical ones **mpmath** (arbitrary precision) and **scipy/numpy**.

## Contents

| file | what it is |
|---|---|
| `reproduce_paper.py` | the script: one function per section of the paper, each printing lines `quantity -> recomputed -> paper -> deviation -> PASS/FAIL` |
| `report.txt` | the saved output of a run with the current paper text |
| `README.md` | this file |

## How to run

```
cd sympy_check
../venv/bin/python reproduce_paper.py        # ~3-4 minutes on a laptop; writes report.txt
```

The script needs only `sympy`, `mpmath`, `numpy`, `scipy` (all in `../venv`). Exit status is 0 if
every check passes. It reads (but never regenerates) three stored data files of the repository:

* `../photon_big_results.m` -- Schwarzschild Zerilli/Regge-Wheeler run, rows `{l, j, Edot_I, Edot_H, u(0), u'(0)}`, l = 10 ... 12800, j <= 4 (Sec. III D);
* `../kerr_py/kerr_results.json` -- Kerr Teukolsky runs, rows with `a, sign, l, j, FluxI, FluxH, lam`, l <= 800 (Sec. IV, Fig. 2, Table I last column);
* `../timelike_results.m` -- timelike circular orbits, rows `{r0, l, Edot_I/E^2, Edot_H/E^2}` with l = m (Sec. VI A).

Two small integrators are part of the script: `zrw_mode`, an independent Zerilli/Regge-Wheeler
integrator (scipy DOP853 in the tortoise coordinate, iterated-Riccati WKB boundary data), and
`scalar_mode`, the same for a scalar charge. They supply the few quantities that are quoted in the
paper but not stored in the repository (the odd-mode ratios of Sec. V, the exact v(0), v'(0) at
l = 6400 of Sec. III D, the exact scalar flux used in the Chrzanowski-Misner comparison of Sec. IV),
and they are validated against the stored data (and, for the scalar one, against the weak-field
dipole formula) before being used.

## What is checked, section by section

**Abstract and Sec. II (Eqs. 1-4).** Light ring, impact parameter, orbital frequency, Upsilon_t, the
null condition, the Wronskian W = 2 i omega A_inc, the largest l of the stored runs. Eq. (2) and the
jump conditions are tested through the independent integrator: its fluxes at l = 10, 20 agree with
the stored Zerilli/RW run to ~1e-6 and with the stored Teukolsky run (a = 0) to ~1e-6 (the paper's
"Teukolsky = Zerilli/RW" check, which independently tests the flux normalisations c_+-).

**Sec. III A (Eqs. 5-6).** Large-l expansion of the jumps for m = l - j (even and odd sectors, exact
where the paper says exact), the O(l^0) cancellation of Eq. (6) proportional to p.p, its value
-4 pi Y mu^2/(E r0) for a massive particle, the separate A and F contributions, and the leading-order
coefficients with which A and F enter the Zerilli jumps.

**Sec. III B (Eqs. 7-9) and its footnote.** V0, k, V3 of the Zerilli and Regge-Wheeler potentials
expanded about r = 3M in r_*, epsilon and eta_j = j + 1/2 (symbolic series), the reduction to Weber's
equation, omega/(2k)^{1/4} = sqrt(l/2) and omega (2k)^{1/4} = sqrt2 l^{3/2}/27, the cubic coefficient
V3 z^3/(2k)^{5/4} = O(l^{-1/2}) z^3; the footnote: the WKB phase shift -(V3/36)(2/k)^{1/2} x^3, the
criterion |x|^3 << sqrt(k)/V3 proportional to 1/l, the local condition |x| << k/V3 = 9M/2 ~ M, the
parabolic region (sqrt(k)/V3)^{1/3} ~ M l^{-1/3}, i.e. |z| << l^{1/6}; lambda_L = sqrt(2k)/(2 omega) =
Omega, Re omega_QNM = (l + 1/2) Omega, the quasinormal condition eta = i(n + 1/2) of the parabolic
barrier and Im omega_QNM = -(n + 1/2) sqrt(2k)/(2 omega), the distance (j + 1/2) Omega of the mode
from resonance; the fact that D_nu(-e^{-i pi/4} z) solves Weber's equation (mpmath `pcfd`), the
normalisation |C|^2 and the transmission coefficient (0.041, 96% reflected) from a WKB decomposition
of the exact parabolic-cylinder solution at |z| = 40, the same from a direct numerical integration of
Weber's equation, the closed forms of D_nu(0), D_nu'(0), and hence Eqs. (8)-(9); mirror symmetry
|v(0)| = |u(0)|.

**Sec. III C (Eqs. 10-14).** Closed forms of |Y_{l,l-j}(pi/2)|^2 and |d_theta Y|^2 against
`scipy.special.sph_harm_y`, their large-l limits with c_j, d_j (scipy at l = 600, closed form at
l = 1e7, sympy limits), the derivation of Eqs. (11)-(12) from the flux formula (symbolic limit),
the numerical kappa_j, their sum and kappa = 0.063653, the fractions missed by j <= 4 and j <= 3,
Eq. (14), the cross term 1.0390, the cubic-barrier contribution -0.231 (perturbative treatment of
V3 z^3/(2k)^{5/4}: the model equation u'' + (z^2/4 - eta - g z^3) u = 0 is solved numerically on
both sides for several g and the linear coefficient of (|u(0)|^2 - |v(0)|^2)/(|u(0)|^2 + |v(0)|^2)
is extrapolated to g -> 0), sigma_0 = 0.81, sigma_j < 0 for j >= 1, sigma_bar = 0.61, the 12 % at
l = 100 and the horizon fraction 40 % for l_max = 100-140 (with the asymptotic formulae).

**Sec. III D.** From the stored run: the three-term fit (0.0318263), the kappa_j at l = 12800
(relative deviations below 1e-3 for j <= 4, below 1e-4 for kappa_0 and 2e-5 for kappa_1), the
consistency of the stored u(0), u'(0) with the stored flux and with Eqs. (8)-(9), the decomposition
of sigma_0 at l = 6400 into the cross term (1.0389) and the barrier asymmetry (-0.231) with u and v
recomputed by the independent integrator (the stored file contains u but not v), sigma_j and
sigma_bar from the data at l = 12800, 3200.

**Sec. IV (Eqs. 15-21, Table I).** All closed forms as sympy identities in r0 (a = sqrt(r0)(3-r0)/2):
b, the factorisation of R~, Upsilon_t, beta_b = sqrt3 r0, R~''(r0) = 6 r0^2, b^2 - a^2 = 3 r0^2,
k~, Omega_theta = lambda_L (and their derivation from the polar and radial geodesic equations: vertical
epicyclic frequency and Lyapunov exponent), eta_j = j + 1/2, the Schroedinger form of the Teukolsky
equation and the potential Q of Eq. (16), Im Q, its vanishing at the light ring and its O(l^{-1/2}) z
size in Weber's variable, Re Q with the eikonal Lambda, the expansion of the angular equation about the
equator (harmonic oscillator of frequency omega beta_b, hence Lambda and the limits of |S(pi/2)|^2 and
|S'(pi/2)|^2 with c_j, d_j from the Hermite functions), the Lambda expansion against the eigenvalues
stored with the Kerr runs, N = 16/r0^6 and h_inf, e^{2 J_H} and the horizon/infinity ratio = 1 (mpmath
quadrature at seven r0 on both branches, plus the polynomial identity that makes the ratio exactly 1),
Omega > Omega_H, the source expansion (`alphaExpand` transcribed from the notebook, evaluated
symbolically in r0): the two cancellations, Ahat_j = (2j+1) Ahat_0, Bhat independent of j, the closed
forms of Ahat_0 and Bhat, the identity |Bhat/Ahat_0|^2 sqrt(2 k~) beta_b = 2, Eq. (19) and its odd
counterpart reproducing the Schwarzschild kappa_j at a = 0, all columns of Table I, Eq. (21) with its
maximum 6 sqrt3 - 9 at a = 3^{1/4}(3 - sqrt3)/2, 81/112, g(0.98), g(0.995), the Richardson
extrapolation of the stored Kerr runs (better than 1.5e-5), the last column of Table I, the
asymmetries sqrt(l)(I-H)/(I+H) = 0.78, 1.10, 0.51, 0.45 at l = 800 (0.61 at a = 0), and the
Chrzanowski-Misner comparison: their null limit (12/sqrt(pi)) e^{-pi/2} E^2 (r0-M)/[r0^2 (r0+3M) m] with
the spin dependence of Eq. (21); the m-counting of the 1973-74 formulae (BCHM 1973 Eq. (5.4) = 2.00
times the exact single-m scalar flux at l = m = 40, delta' = 0.01, with the scalar integrator validated
in the weak field: 0.955 x M^2/(6 r0^4) at r0 = 100 M, tending to the dipole value as 1 - 5 M/r0);
their constant 0.94 of 2 kappa_0; in the massive regime their scalar formula = BCHM's exact one, their
tensor formula = 4 x (BRTV's exact one per |m|), the (s!)^2 = 4 of their master formula; the WKB
barrier factor eta^{3/2} e^{-pi eta} at eta = 1/2 is 1/(4 x 1.065) of the exact (eta+1/2)^2 F(eta).

**Sec. V (Eqs. 22-23).** The null limits (kappa_0 exactly, kappa_1/4 exactly, via the reflection
formula), kappa = 0.0574 and 0.0580 with the 1973 odd term, and the odd-mode ratios 3.15, 3.48, 3.64,
3.77 at delta = 1e-4 recomputed with the independent integrator (these exact Teukolsky fluxes are not
stored in the repository), with the fits 4.06 / 4.05.

**Sec. VI (Eq. 24) and Conclusions.** gamma^2(delta), the shifted jump and eta for timelike orbits,
Eq. (24) against the stored timelike fluxes over exactly the ranges quoted in the paper (better than 1%
for gamma = 20 and 120 <= l <= 4800, 0.4-4% for gamma = 10 and 30 <= l <= 1200, 2-4% for gamma = 5 and
25 <= l <= 125), the cut-off (l/6 gamma^2)^{3/2} exp(-pi l/6 gamma^2), c_1 = pi/6, 0.127 ln gamma =
kappa ln gamma^2, the orbital period 6 sqrt3 pi M, Delta E = 2 x 2 pi b x kappa = 4.2 (E^2/M) ln l_max,
the M87* estimate (l_max = 2 pi r0/lambda ~ 1e17 for lambda = 1.3 mm and M = 6.5e9 Msun, E/M ~ 1e-79,
Delta E/E ~ 1e-77).

Items that need the Toolkit (MST solutions at l <= 5, the comparison of Lambda and S(pi/2) with the
Toolkit's spheroidal harmonics) are listed as `NOT CHECKED` with the reason.

## Conventions

* G = c = 1 and **M = 1** in the code. Where the paper writes explicit powers of M they are noted in
  the labels (e.g. k = [2 l^2 + O(l)]/(729 M^4) is checked as 729 k/l^2 -> 2).
* E = 1: all fluxes are **per unit E^2**, for a **single m > 0** and a single (l, m) mode, with
  **j = l - m**. The physical flux per multipole is 2 x sum over m > 0 (the -m modes), which is how
  kappa = 2 sum_j kappa_j and the Kerr comparisons (`l (Edot_I + Edot_H)` summed over j <= 3) are built.
* Jumps: [Psi] = Psi(r0+) - Psi(r0-), [d_r Psi] = d_r Psi(r0+) - d_r Psi(r0-), and
  [d_r* Psi] = f(r0) [d_r Psi]. The symbols Y, dY in the jump formulae stand for conj(Y_lm) and
  conj(d_theta Y_lm) at (pi/2, 0), both real; the source coefficients A, C, F follow Eq. (2)
  (Toolkit ReggeWheeler conventions), the jump formulae in terms of A, C, F are transcribed from the
  notebook's `jumpsZRW`.
* delta: the paper's delta, r0 = 3M(1 + delta); Breuer et al.'s delta' is 3 delta (their
  epsilon = 2 eta = 1 + 3 m delta).
* eta_j = j + 1/2; F(eta) = e^{-pi eta/2}/[cosh(pi eta) |Gamma(3/4 + i eta/2)|^2] as in the paper;
  the script also uses the shorthand G(eta) = e^{-pi eta/2}/[cosh(pi eta) |Gamma(1/4 + i eta/2)|^2]
  for the odd sector (not a notation of the paper).
* sigma_j: Edot^{inf,H}_{l,l-j} = (kappa_j/l)[1 +- sigma_j/sqrt(l) + O(1/l)]; sigma_j = cross term
  (from the phase of Eq. 14) + cubic-barrier asymmetry; sigma_bar = sum_j kappa_j sigma_j / sum_j kappa_j.
* Weber model: z = (2k)^{1/4} x, g = V3/[6 (2k)^{5/4}] = g0 l^{-1/2}; "In" = incident from z -> +inf
  (infinity side), "Up" = incident from z -> -inf; |u(0)|^2 and |u'(0)|^2 are quoted for unit incident
  WKB amplitude in units where omega/(2k)^{1/4} = 1 and omega (2k)^{1/4} = 1 respectively.
* Kerr: r0 parametrises the light ring, a = sqrt(r0)(3 - r0)/2 with 1 < r0 < 4 (prograde for
  r0 < 3, retrograde for r0 > 3, negative a); s = -2; the Kerr identities are checked for symbolic r0
  (both branches) and numerically at the seven spins of Table I.
* Scalar field (Chrzanowski-Misner comparison): Box Phi = 4 pi q Int dtau delta^4/sqrt(-g), stress
  tensor (1/4 pi)(dPhi dPhi - g (dPhi)^2/2), radial equation u'' + (omega^2 - V0) u = S0 delta(r* - r0*)
  with V0 = f (l(l+1)/r^2 + 2M/r^3) and S0 = 4 pi conj(Y_lm(pi/2,0))/(u^t r0), flux (1/4 pi) omega^2 |Z|^2,
  per q^2 and per single m.

## Tolerances and the report

Each numeric line shows the relative deviation (absolute when the paper value is 0) and PASS/FAIL
with a tolerance chosen from the number of digits quoted in the paper (e.g. 2e-6 for kappa_0 =
0.0277668) or from the expected size of the neglected corrections (e.g. O(l^-1/2) for comparisons
with finite-l data); statements of the form "better than X" are checked as upper bounds. Symbolic
identities show `0 (symbolic)`. The summary at the end of `report.txt` gives the number of
passed/failed checks, the maximum relative deviation among the passed numeric checks against a quoted
value, and the total runtime.

## Result of the saved run (`report.txt`)

285 checks pass, 0 fail, 16 lines are informational or `NOT CHECKED`; runtime 171 s. Among the
quantitative checks the largest relative deviation from a quoted value is 1.7e-2 (2 sigma_bar/sqrt(l) at l = 100, quoted as 12%: recomputed 12.2%); the
others are at the level of the rounding of the quoted digits (e.g. kappa_j to 1e-6 ... 2e-4, Table I
entries to < 6e-4, the cross term 1.0390 to 5e-5, the cubic term -0.231 to 7e-4, the odd-mode ratios
3.15/3.48/3.64/3.77 to ~1e-3, the Kerr asymmetries at l = 800 to < 1e-2).

Points worth noting:

1. **Timelike comparison.** The script checks exactly the ranges stated in the paper and also prints
   the deviation for every stored row (it reaches 11-15% at the smallest l and, for gamma = 5, beyond
   the cut-off).
2. **Orders of magnitude.** The M87* cut-off l_max ~ 2 pi r0/lambda evaluates to 1.4e17 (paper: ~1e17),
   E/M to 1.3e-79 (paper: ~1e-79) and Delta E/E to 2.2e-77 per orbit (paper: ~1e-77).
3. **sigma_0 at l = 6400.** The stored Schwarzschild run contains u(0), u'(0) but not v(0), v'(0); the
   independent integrator reproduces the stored fluxes and u at l = 6400 to < 1e-5 and gives the cross
   term 1.03893 and the barrier asymmetry -0.2309 (paper: 1.0389 and -0.231); the analytic sigma_j
   (cross + cubic) agree with the stored l = 12800 asymmetries to < 0.25% for every j <= 4.
4. **Sec. V ratios.** The four exact odd-mode Teukolsky ratios are not stored in the repository; the
   independent integrator gives 3.153, 3.479, 3.645, 3.767 (paper: 3.15, 3.48, 3.64, 3.77) and the
   fits 4.062 / 4.047 (paper: 4.06 / 4.05).
5. **Scalar check.** The scalar flux of a charge on a circular orbit at r0 = 100 M, l = m = 1, is
   0.9547 x M^2/(6 r0^4); (1 - ratio) r0/M = 4.26, 4.53, 4.73, 4.86 at r0/M = 50, 100, 200, 400, i.e. the
   ratio tends to 1 as 1 - 5.0 M/r0. BCHM Eq. (5.4) is 1.997 times the exact single-m flux at l = m = 40,
   delta' = 0.01 (2.008 and 1.992 at (100, 0.003) and (400, 0.001)): their power at a given frequency
   collects m and -m.
6. Items that need the Black Hole Perturbation Toolkit (MST solutions at l <= 5, the comparison with
   the Toolkit's spheroidal harmonics) are listed as `NOT CHECKED`; the Lambda expansion is derived
   from the angular equation and checked against the eigenvalues stored with the Kerr runs instead.
