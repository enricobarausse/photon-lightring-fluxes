#!/usr/bin/env python
"""
reproduce_paper.py -- independent Python re-derivation of the analytic results and of the quoted
numbers of

    E. Barausse, "Gravitational radiation from a photon on the light ring: multipole fluxes and the
    logarithmic divergence"   (equation numbers as in the submitted version of the paper)

The symbolic parts use sympy, the numerical ones mpmath (arbitrary precision) and scipy/numpy.
Nothing here depends on Mathematica, on the companion notebook or on the Black Hole Perturbation
Toolkit. The only formulae transcribed from the notebook source (notebook/make_notebook.py) are, as
allowed by the task: the jump conditions as linear combinations of the source coefficients A, C, F
(function jumpsZRW), the circular-orbit Teukolsky source projection alpha_lm of Hughes (2000)
(function alphaExpand), and the Chrzanowski-Misner / Breuer et al. 1973 formulae of Section 7.
Everything else is derived here from the equations written in the paper.

Conventions (identical to the paper unless stated):
  * G = c = 1 and M = 1 in the code; the paper's powers of M are restored by dimensional analysis and
    are noted in the labels.
  * Fluxes are per unit E^2 (E = Killing energy of the photon), for a single m > 0 and a single
    (l, m) mode, with j = l - m. The physical flux per multipole is twice the sum over m > 0.
  * Jumps: [Psi] = Psi(r0+) - Psi(r0-), [d_r Psi] likewise; tortoise-coordinate jump
    [d_r* Psi] = f(r0) [d_r Psi]. The angular factors Y, dY stand for conj(Y_lm) and
    conj(d_theta Y_lm) at (theta, phi) = (pi/2, 0), where both are real.
  * delta: the paper's delta (r0 = 3M(1 + delta)); Breuer et al.'s delta' = 3 delta.
  * eta_j = j + 1/2; F(eta) = e^{-pi eta/2}/(cosh(pi eta) |Gamma(3/4 + i eta/2)|^2);
    G(eta) = e^{-pi eta/2}/(cosh(pi eta) |Gamma(1/4 + i eta/2)|^2)   (G is our shorthand, not the paper's).
  * sigma_j: Edot^{inf,H}_{l,l-j} = (kappa_j/l) [1 +- sigma_j/sqrt(l) + O(1/l)].

Stored numerical data used (never regenerated):
  * ../photon_big_results.m    Schwarzschild run, rows {l, j, Edot_I, Edot_H, u(0), u'(0)}, l = 10..12800
  * ../kerr_py/kerr_results.json   Kerr Teukolsky runs, rows {a, sign, l, j, FluxI, FluxH, lam, ...}, l <= 800
  * ../timelike_results.m      timelike orbits, rows {r0, l, Edot_I/E^2, Edot_H/E^2}, l = m
  * ../toolkit_checks/schwarzschild_toolkit_results.m   Black Hole Perturbation Toolkit run at l <= 5, rows
        {l, m, method, Edot_I, Edot_H} with method "RW-MST" (Zerilli/RW with MST solutions), "Teuk" (the null
        source in the Toolkit's Teukolsky solver), "MP" (the production integrator), "NI" (arbitrary precision)
  * ../schwarzschild/asym_check_results.m   decomposition of sqrt(l)(I-H)/(I+H) at j = 0 into the barrier
        asymmetry A1, the cross term A2 and the O(1/l) remainder A3, from the production integrator, l = 100..6400
  * ../kerr_py/check_odd_BRTV_results.txt   exact Teukolsky fluxes of timelike orbits vs the 1973 formulae
  * ../notebook/toolkit_harmonics.m   Toolkit spheroidal harmonics, rows {a, sign, l, j, Lambda, |S(pi/2)|^2, |S'(pi/2)|^2}
  The last four are optional: a missing file produces NOT CHECKED lines instead of failures. The spectral
  spheroidal-harmonic solver of ../kerr_py/swsh.py (numpy + mpmath; validated against the Toolkit at l <= 6)
  is imported for the Lambda and S(pi/2) asymptotics of Sec. IV.
Numbers recomputed here rather than read from the repository, with a small independent Zerilli/Regge-Wheeler
integrator (scipy DOP853 in the tortoise coordinate, iterated-Riccati WKB boundary data; function zrw_mode),
which is validated against the stored data and against the Toolkit MST fluxes at l <= 5 (5e-8): the four
exact odd-mode ratios of Sec. V (3.15, 3.48, 3.64, 3.77), and the exact v(0), v'(0) at l = 6400 that the
decomposition of sigma_0 at the end of Sec. III D needs (both are also read from the stored files above
and the two are compared). A scalar version of the same integrator (scalar_mode) provides the exact scalar
flux with which Sec. IV establishes that the 1973-74 formulae collect the modes m and -m.

Structure: one function per section of the paper (section_II ... section_VI), each printing lines
    quantity -> recomputed -> paper value -> deviation -> PASS/FAIL
through the Report class; identities print "0 (symbolic)"; informational lines (no quantitative claim in
the paper) are marked "info"; a missing optional data file gives NOT CHECKED with the reason.

Run:  ../venv/bin/python reproduce_paper.py   (writes report.txt next to this file; ~4-5 minutes)
"""
import os, re, sys, json, time, math
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import solve_ivp, quad
from scipy.special import sph_harm_y, gammaln

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
mp.mp.dps = 30
T_START = time.time()

# =============================================================================================
# Report machinery
# =============================================================================================
class Report:
    """Collects checks as (name, recomputed, paper, agreement) lines, section by section."""
    def __init__(self):
        self.lines = []; self.npass = 0; self.nfail = 0; self.ninfo = 0
        self.maxdev = 0.0; self.maxdev_name = ''; self.fails = []

    def out(self, s=''):
        print(s, flush=True); self.lines.append(s)

    def section(self, title):
        self.out(); self.out('=' * 110); self.out(title); self.out('=' * 110)
        self.out(f"{'quantity':<62}{'recomputed':>18}{'paper':>14}{'deviation':>11}  status")

    def note(self, text):
        for line in text.split('\n'):
            self.out('   ' + line)

    @staticmethod
    def fmt(v):
        if isinstance(v, str): return v
        v = float(v)
        if v == 0: return '0'
        if 1e-3 <= abs(v) < 1e5: return f'{v:.8g}'
        return f'{v:.6e}'

    def check(self, name, value, paper, tol, kind='rel', note=''):
        """Numeric check. paper may be a string (as quoted in the paper); deviation is relative
        (kind='rel', default) or absolute (kind='abs'); status PASS if deviation <= tol."""
        v = float(value)
        try: p = float(paper)
        except ValueError: p = float(sp.sympify(paper))
        if kind == 'max':            # PASS if value <= paper (an upper bound quoted in the paper)
            dev = max(v - p, 0.0); ok = v <= p
        else:
            dev = abs(v - p) / abs(p) if (kind == 'rel' and p != 0) else abs(v - p)
            ok = dev <= tol
        if ok: self.npass += 1
        else: self.nfail += 1; self.fails.append((name, v, p, dev))
        # the max-deviation statistic covers the quantitative checks only (tolerance <= 5%), not order-of-magnitude statements
        if kind == 'rel' and p != 0 and ok and tol <= 0.05 and dev > self.maxdev:
            self.maxdev = dev; self.maxdev_name = name
        tag = 'PASS' if ok else 'FAIL'
        ptxt = paper if isinstance(paper, str) else self.fmt(p)
        self.out(f"{name:<62}{self.fmt(v):>18}{ptxt:>14}{dev:>11.1e}  {tag}" + (f'   [{note}]' if note else ''))
        return ok

    def ident(self, name, expr, note=''):
        """Symbolic identity: PASS if expr simplifies to zero (or vanishes numerically at random points)."""
        how = is_zero(expr)
        ok = how is not None
        if ok: self.npass += 1
        else: self.nfail += 1; self.fails.append((name, 'not zero', 0, float('nan')))
        self.out(f"{name:<62}{('0 (' + how + ')') if ok else 'NOT ZERO':>18}{'identity':>14}{'':>11}  {'PASS' if ok else 'FAIL'}"
                 + (f'   [{note}]' if note else ''))
        return ok

    def info(self, name, value, note=''):
        self.ninfo += 1
        self.out(f"{name:<62}{self.fmt(value):>18}{'--':>14}{'':>11}  info" + (f'   [{note}]' if note else ''))

    def unavailable(self, name, why):
        self.ninfo += 1
        self.out(f"{name:<62}{'--':>18}{'--':>14}{'':>11}  NOT CHECKED   [{why}]")


def is_zero(expr):
    """Return a short description if sympy expr is zero (symbolically, or numerically at random points), else None."""
    e = sp.nsimplify(expr) if expr == 0 else expr
    try:
        s = sp.simplify(e)
        if s == 0: return 'symbolic'
        s2 = sp.simplify(sp.radsimp(sp.expand(s)))
        if s2 == 0: return 'symbolic'
    except Exception:
        s = e
    syms = sorted(s.free_symbols, key=lambda x: x.name)
    rng = np.random.default_rng(12345)
    worst = 0
    for _ in range(4):
        vals = {x: sp.Rational(int(rng.integers(11, 39)), 10) for x in syms}   # points in (1.1, 3.9): fine for r0 and l
        num = sp.N(s.subs(vals), 40)
        try: worst = max(worst, abs(complex(num)))
        except Exception: return None
    return f'numeric<{max(worst, 1e-30):.0e}' if worst < 1e-25 else None


rep = Report()

# =============================================================================================
# Shared sympy objects
# =============================================================================================
pi = sp.pi; I = sp.I
l, j, E0, L0, Y, dY = sp.symbols('l j E0 L0 Y dY', positive=True)
m = sp.symbols('m', positive=True)
r, r0 = sp.symbols('r r0', positive=True)
eps = sp.symbols('epsilon', positive=True)   # eps = 1/l for large-l series

def f_schw(rr): return 1 - 2/rr

def series_in_l(expr, n=3):
    """Large-l expansion: substitute l = 1/eps and expand in eps (eps = 1/l) to O(eps^n)."""
    return sp.expand(sp.series(expr.subs(l, 1/eps), eps, 0, n).removeO())

# ---------------------------------------------------------------------------------------------
# Eq. (2) and the jump conditions
# ---------------------------------------------------------------------------------------------
def source_ACF(l, m, r0, E0, L0, Y, dY):
    """Eq. (2): stress-energy coefficients A, C, F in the Toolkit ReggeWheeler conventions, with
    E -> E0 and bE -> L0 (so b^2 E = L0^2/E0). Y, dY = conj(Y_lm), conj(d_theta Y_lm) at (pi/2, 0)."""
    r0, E0, L0 = map(sp.sympify, (r0, E0, L0))   # keep exact arithmetic when integers are passed
    A = -16*pi*(r0 - 2)*E0/r0**3*Y
    C = -16*pi*(r0 - 2)*L0/r0**4/(l*(l + 1))*dY
    F = -16*pi*(r0 - 2)*(L0**2/E0)/r0**5*(l*(l + 1) - 2*m**2)/((l - 1)*l*(l + 1)*(l + 2))*Y
    return A, C, F

def jumps(l, m, r0, E0, L0, Y, dY, even):
    """[Psi], [d_r Psi] at r0 as linear combinations of A, C, F (Sago, Nakano & Sasaki 2003 / Toolkit
    conventions; formulae transcribed from the companion notebook's jumpsZRW, as permitted)."""
    r0, E0, L0 = map(sp.sympify, (r0, E0, L0))
    A, C, F = source_ACF(l, m, r0, E0, L0, Y, dY)
    if even:
        n = (l - 1)*(l + 2); rm2M = r0 - 2; np6M = n*r0 + 6
        term1 = -r0**4*A/2/np6M/rm2M
        dterm1 = -r0**4*A/2/np6M/rm2M*(4/r0 - 1/rm2M - n/np6M)
        coeff = 2/r0/rm2M
        term2 = r0**3/4*(r0**2*n*(n - 2) + r0*(14*n - 36) + 96)*A/(rm2M*np6M)**2 + (n + 2)*r0**2/4*F/rm2M
        return term1, -coeff*term1 - dterm1 + term2
    return (r0**3/(r0 - 2)*C,
            -2*r0**2/(r0 - 2)**2*C + r0**2/(r0 - 2)*C - 3*r0**2/(r0 - 2)*C + r0**3/(r0 - 2)**2*C)

# numeric versions of the jumps (fast evaluation)
_jl, _jm, _jr0, _jE, _jL, _jY, _jdY = sp.symbols('jl jm jr0 jE jL jY jdY')
_jumps_even_num = sp.lambdify((_jl, _jm, _jr0, _jE, _jL, _jY, _jdY), jumps(_jl, _jm, _jr0, _jE, _jL, _jY, _jdY, True), 'numpy')
_jumps_odd_num = sp.lambdify((_jl, _jm, _jr0, _jE, _jL, _jY, _jdY), jumps(_jl, _jm, _jr0, _jE, _jL, _jY, _jdY, False), 'numpy')

def jumps_num(l, m, r0, E0, L0, Yv, dYv):
    even = (l + m) % 2 == 0
    fn = _jumps_even_num if even else _jumps_odd_num
    P, dP = fn(float(l), float(m), r0, E0, L0, Yv, dYv)
    return float(P), float(dP)

def c_pm(l, m):
    """Flux normalization c_+ (Zerilli) or c_- (Regge-Wheeler), text after Eq. (4)."""
    if (l + m) % 2 == 0: return (l - 1)*(l + 2)/(4*math.pi*l*(l + 1))
    return l*(l + 1)/(16*math.pi*(l - 1)*(l + 2))

# ---------------------------------------------------------------------------------------------
# Potentials (Sec. III B)
# ---------------------------------------------------------------------------------------------
def V_RW(l, r): return f_schw(r)*(l*(l + 1)/r**2 - 6/r**3)
def V_Z(l, r):
    lam = (l - 1)*(l + 2)/2
    return f_schw(r)*(2*lam**2*(lam + 1)*r**3 + 6*lam**2*r**2 + 18*lam*r + 18)/(r**3*(lam*r + 3)**2)
def Dst(e): return f_schw(r)*sp.diff(e, r)   # d/dr*

# numeric potentials as functions of x = r - 2 (f = x/(x+2) keeps full precision near the horizon)
def V_RW_num(l, x):
    rr = x + 2; return x/rr*(l*(l + 1)/rr**2 - 6/rr**3)
def V_Z_num(l, x):
    lam = (l - 1)*(l + 2)/2.0; rr = x + 2
    return x/rr*(2*lam**2*(lam + 1)*rr**3 + 6*lam**2*rr**2 + 18*lam*rr + 18)/(rr**3*(lam*rr + 3)**2)

# ---------------------------------------------------------------------------------------------
# Spherical harmonics at the equator
# ---------------------------------------------------------------------------------------------
def Y_eq(l, m):
    """Y_lm(pi/2, 0) (real) from scipy."""
    return float(np.real(sph_harm_y(l, m, np.pi/2, 0.0)))
def dY_eq(l, m):
    """d_theta Y_lm(pi/2, 0) (real) from scipy."""
    val, grad = sph_harm_y(l, m, np.pi/2, 0.0, diff_n=1)
    return float(np.real(grad[0]))
def Ysq_closed(l, j):
    """|Y_{l,l-j}(pi/2,0)|^2 for even j, closed form (log-gamma, stable at large l)."""
    if j % 2: return 0.0
    lg = (math.log((2*l + 1)/(4*math.pi)) + gammaln(j + 1) + gammaln(2*l - j + 1) - 2*l*math.log(2)
          - 2*gammaln(j/2 + 1) - 2*gammaln((2*l - j)/2 + 1))
    return math.exp(lg)
def Y_closed(l, j):
    """Y_{l,l-j}(pi/2,0) for even j including the sign (-1)^{(2l-j)/2}."""
    return (-1)**((2*l - j)//2)*math.sqrt(Ysq_closed(l, j))
def dYsq_closed(l, j):
    """|d_theta Y_{l,l-j}(pi/2,0)|^2 for odd j: j (2l-j+1) |Y_{l,l-j+1}|^2."""
    if j % 2 == 0: return 0.0
    return j*(2*l - j + 1)*Ysq_closed(l, j - 1)
def dY_closed(l, j):
    """d_theta Y_{l,l-j}(pi/2,0), odd j, with the sign of scipy's derivative."""
    # sign fixed empirically against scipy (checked in Sec. III C); magnitude from the closed form
    return math.sqrt(dYsq_closed(l, j))*(-1)**((2*l - j + 1)//2)
def c_j(j): return sp.binomial(j, sp.Rational(j, 2))/sp.Integer(2)**j
def d_j(j): return j*sp.binomial(j - 1, sp.Rational(j - 1, 2))/sp.Integer(2)**(j - 1)

# ---------------------------------------------------------------------------------------------
# Barrier-top functions (mpmath)
# ---------------------------------------------------------------------------------------------
def Fpar(eta): return mp.exp(-mp.pi*eta/2)/(mp.cosh(mp.pi*eta)*abs(mp.gamma(mp.mpf(3)/4 + 1j*eta/2))**2)
def Gpar(eta): return mp.exp(-mp.pi*eta/2)/(mp.cosh(mp.pi*eta)*abs(mp.gamma(mp.mpf(1)/4 + 1j*eta/2))**2)
def kappa_j(j):
    """Eqs. (11)-(12), M = 1."""
    x = mp.mpf(j) + mp.mpf(1)/2
    if j % 2 == 0: return mp.sqrt(mp.pi)/9*(2*j + 1)**2*mp.mpf(c_j(j))*Fpar(x)
    return 8*mp.sqrt(mp.pi)/9*mp.mpf(d_j(j))*Gpar(x)

# =============================================================================================
# Independent Zerilli / Regge-Wheeler integrator (scipy), used for the validations of Sec. II and
# for the timelike odd-mode ratios of Sec. V.  Variable x = r - 2 keeps precision near the horizon.
# =============================================================================================
def rstar_of_r(rr): return rr + 2*math.log(rr/2 - 1)
def x_of_rstar(rs):
    x = 2*math.exp((rs - 2)/2)
    for _ in range(60):
        Fv = x + 2 + 2*math.log(x/2) - rs
        x -= Fv/(1 + 2/x)
    return x

def riccati_Q(V, w, x, sign, niter=3):
    """Iterated Riccati-WKB local wavenumber: Q_{n+1} = sign sqrt(w^2 - V + i f dQ_n/dr), f = x/(x+2)."""
    def Q(k, x):
        if k == 0: return sign*np.sqrt(w*w - V(x) + 0j)
        h = 1e-3*x
        return sign*np.sqrt(w*w - V(x) + 1j*(x/(x + 2))*(Q(k - 1, x + h) - Q(k - 1, x - h))/(2*h))
    return Q(niter, x)

def _integrate(V, w, rs_from, x_from, rs_to, Q0, amp=1.0, rtol=1e-11):
    y0 = np.array([amp, 0.0, (1j*Q0*amp).real, (1j*Q0*amp).imag, x_from])
    def rhs(rs, y):
        x = y[4]; q = w*w - V(x)
        return [y[2], y[3], -q*y[0], -q*y[1], x/(x + 2)]
    sol = solve_ivp(rhs, (rs_from, rs_to), y0, method='DOP853', rtol=rtol, atol=1e-14)
    y = sol.y[:, -1]
    return y[0] + 1j*y[1], y[2] + 1j*y[3]

def zrw_mode(l, m, r0, E0, L0, Omega, rsA=-40.0, rB=40.0):
    """Fluxes (Edot_I, Edot_H) per unit E0^2-normalization of a circular orbit source with energy E0,
    angular momentum L0, radius r0, orbital frequency Omega, mode (l, m): Eq. (4) and the flux formula.
    Returns also u(0), u'(0), v(0), v'(0) (unit incident amplitude) and the jumps P, dP."""
    even = (l + m) % 2 == 0
    V = (lambda x: V_Z_num(l, x)) if even else (lambda x: V_RW_num(l, x))
    w = m*Omega
    xA = x_of_rstar(rsA); rs0 = rstar_of_r(r0)
    QA = riccati_Q(V, w, xA, -1); QB = riccati_Q(V, w, rB - 2, +1)
    # WKB amplitudes: Psi = exp(i int Q dr*), so |Psi(r_B)| = exp(+int_{r_B}^inf Im Q dr*) gives unit amplitude at
    # infinity, and |Psi(r_A)| = exp(-int_{-inf}^{r_A} Im Q dr*) unit amplitude at the horizon (dr* = dx/f, f = x/(x+2))
    ampB = math.exp(quad(lambda x: riccati_Q(V, w, x, +1).imag*(x + 2)/x, rB - 2, np.inf, limit=200)[0])
    ampA = math.exp(-quad(lambda x: riccati_Q(V, w, x, -1).imag*(x + 2)/x, 0.0, xA, limit=200)[0])   # ~1 - O(1e-9)
    PsiIn, dPsiIn = _integrate(V, w, rsA, xA, rs0, QA, ampA)
    PsiUp, dPsiUp = _integrate(V, w, rstar_of_r(rB), rB - 2, rs0, QB, ampB)
    W = PsiIn*dPsiUp - PsiUp*dPsiIn
    Ainc = W/(2j*w)
    P, dPr = jumps_num(l, m, r0, E0, L0, Y_eq(l, m), dY_eq(l, m))
    dP = (1 - 2/r0)*dPr
    ZI = (PsiIn*dP - P*dPsiIn)/W; ZH = (PsiUp*dP - P*dPsiUp)/W
    c = c_pm(l, m)
    return dict(I=c*w*w*abs(ZI)**2, H=c*w*w*abs(ZH)**2, u0=PsiIn/Ainc, Du0=dPsiIn/Ainc,
                v0=PsiUp/Ainc, Dv0=dPsiUp/Ainc, P=P, dP=dP, w=w)

# ---------------------------------------------------------------------------------------------
# Scalar field: the same integrator for a scalar charge q on a circular geodesic (used in Sec. IV to
# establish the m-counting of the 1973-74 scalar formula of Breuer, Chrzanowski, Hughes and Misner).
# Conventions (as in BCHM 1973 and in the companion notebook): field equation
#   Box Phi = 4 pi q Int dtau delta^4(x - z(tau))/sqrt(-g),  stress tensor (1/4 pi)(dPhi dPhi - g (dPhi)^2/2);
# with Phi = sum (u_lm(r)/r) Y_lm e^{-i omega t}, omega = m Omega, the radial function obeys
#   u'' + (omega^2 - V0) u = S0 delta(r* - r0*),   V0 = f (l(l+1)/r^2 + 2M/r^3),   S0 = 4 pi conj(Y_lm(pi/2,0))/(u^t r0),
# the solution is u = u_in(r<) u_up(r>) S0/W, and the flux at infinity is (1/4 pi) omega^2 |Z|^2 with
# Z = u_in(r0) S0/W. Everything is per q^2; a single m > 0.
# ---------------------------------------------------------------------------------------------
def V_scalar_num(l, x):
    rr = x + 2; return x/rr*(l*(l + 1)/rr**2 + 2/rr**3)

def scalar_mode(l, m, r0, rsA=-40.0, rB=None):
    """Flux at infinity, per q^2, of the (l, m) mode of a scalar charge on the circular geodesic of radius r0
    (Schwarzschild, M = 1): u^t = (1 - 3M/r0)^{-1/2}, Omega = r0^{-3/2}. The outgoing data are imposed at
    rB, outside the evanescent region (rB >~ 8/omega), with the Riccati-WKB wavenumber."""
    V = lambda x: V_scalar_num(l, x)
    w = m*r0**-1.5; ut = 1/math.sqrt(1 - 3/r0)
    S0 = 4*math.pi*np.conj(sph_harm_y(l, m, math.pi/2, 0.0))/(ut*r0)
    if rB is None: rB = max(40.0, 3*r0, 8/w)
    xA = x_of_rstar(rsA); rs0 = rstar_of_r(r0)
    QA = riccati_Q(V, w, xA, -1); QB = riccati_Q(V, w, rB - 2, +1)
    # amplitude integrals as in zrw_mode; for the weak-field case (omega ~ 1e-3) the finite-difference Riccati
    # tail is noisy at the 1e-9 level and quad reports a (harmless) convergence warning, which is suppressed
    import warnings
    from scipy.integrate import IntegrationWarning
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', IntegrationWarning)
        ampB = math.exp(quad(lambda x: riccati_Q(V, w, x, +1).imag*(x + 2)/x, rB - 2, np.inf, limit=400)[0])
        ampA = math.exp(-quad(lambda x: riccati_Q(V, w, x, -1).imag*(x + 2)/x, 0.0, xA, limit=400)[0])
    uin, duin = _integrate(V, w, rsA, xA, rs0, QA, ampA)
    uup, duup = _integrate(V, w, rstar_of_r(rB), rB - 2, rs0, QB, ampB)
    W = uin*duup - duin*uup
    return w*w*abs(uin*S0/W)**2/(4*math.pi)

# =============================================================================================
# Stored data
# =============================================================================================
def load_mathematica_list(path):
    txt = open(path).read()
    s = re.sub(r'\s+', '', txt).replace('*^', 'e').replace('I', '1j').replace('{', '[').replace('}', ']')
    return eval(s, {'__builtins__': {}})

def load_mathematica_file(path):
    """Mathematica Put[] output with comments, line continuations, strings, rationals and arbitrary-precision
    numbers (0.0186...`20.): returns (nested lists, list of the precision marks found, as floats)."""
    txt = open(path).read()
    s = re.sub(r'\(\*.*?\*\)', '', txt, flags=re.S)        # comments
    s = re.sub(r'\\\s*\n', '', s)                            # line continuations inside long numbers
    marks = [float(x) for x in re.findall(r'`([\d.]+)', s) if x not in ('', '.')]
    s = re.sub(r'`[\d.]*', '', s)                            # precision marks
    s = re.sub(r'\s+', '', s).replace('*^', 'e').replace('{', '[').replace('}', ']')
    return eval(s, {'__builtins__': {}}), marks

BIG = load_mathematica_list(os.path.join(ROOT, 'photon_big_results.m'))          # {l, j, EdotI, EdotH, u0, Du0}
TIMELIKE = load_mathematica_list(os.path.join(ROOT, 'timelike_results.m'))       # {r0, l, EdotI/E^2, EdotH/E^2}
KERR = json.load(open(os.path.join(ROOT, 'kerr_py', 'kerr_results.json')))       # dicts a, sign, l, j, FluxI, FluxH, lam
# optional stored files (each check that needs one prints NOT CHECKED if the file is absent)
def _optional(path, loader):
    try: return loader(path)
    except (OSError, SyntaxError, ValueError): return None
SCHW_TK = _optional(os.path.join(ROOT, 'toolkit_checks', 'schwarzschild_toolkit_results.m'), load_mathematica_file)   # {l, m, method, EdotI, EdotH}
ASYM = _optional(os.path.join(ROOT, 'schwarzschild', 'asym_check_results.m'), lambda p: load_mathematica_file(p)[0])  # {l, sqrt(l)(I-H)/(I+H), A1, A2, A3}
HARM_TK = _optional(os.path.join(ROOT, 'notebook', 'toolkit_harmonics.m'), lambda p: load_mathematica_file(p)[0])     # {a, sign, l, j, Lambda, S2, dS2}
def _load_brtv(path):
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) == 6 and not line.startswith('#'):
            rows.append((float(p[0]), int(p[1]), p[2], float(p[3]), float(p[4]), float(p[5])))
    return rows or None
BRTV_STORED = _optional(os.path.join(ROOT, 'kerr_py', 'check_odd_BRTV_results.txt'), _load_brtv)   # (r0, m, parity, exact, BRTV, ratio)
try:
    sys.path.insert(0, os.path.join(ROOT, 'kerr_py'))
    from swsh import S_at_equator as SWSH_at_equator      # the spectral spheroidal-harmonic solver of kerr_py (numpy + mpmath)
except Exception:
    SWSH_at_equator = None
def big_row(lv, jv):
    return next(rw for rw in BIG if rw[0] == lv and rw[1] == jv)
def kerr_rows(a, sign, lv, jmax=3):
    return [rw for rw in KERR if rw['a'] == abs(a) and rw['sign'] == sign and rw['l'] == lv and rw['j'] <= jmax]

# =============================================================================================
# Section II
# =============================================================================================
def section_II():
    rep.section('Sec. II  The problem and the formalism  [Eqs. (1)-(4)]')
    # light ring and impact parameter from the null-geodesic potential f/r^2
    rr = sp.symbols('rr', positive=True)
    lr = sp.solve(sp.diff(f_schw(rr)/rr**2, rr), rr)
    rep.check('light-ring radius r0/M (max of f/r^2)', float(lr[0]), '3', 0)
    bsq = (rr**2/f_schw(rr)).subs(rr, 3)
    rep.ident('b^2 = r0^2/f(r0) = 27 M^2  (b = 3 sqrt3 M)', bsq - 27)
    rep.check('Omega M = 1/b', float(1/sp.sqrt(bsq)), float(1/(3*math.sqrt(3))), 1e-15)
    rep.ident('Omega = 1/b is the orbital frequency of the null geodesic: dphi/dt = L f(r0)/(E r0^2) = 1/b at r0 = 3M',
              (3*sp.sqrt(3))*f_schw(sp.Integer(3))/sp.Integer(9) - 1/(3*sp.sqrt(3)))
    rep.ident('p^t = E/f(r0) = 3E, Upsilon_t = r0^2 p^t = 27 M^2 E', 9/f_schw(sp.Integer(3)) - 27)
    rep.info('largest l in the stored runs (abstract: l = 12800 in Schwarzschild, l = 800 in Kerr)',
             f"{max(rw[0] for rw in BIG)}, {max(rw['l'] for rw in KERR)}")
    rep.ident('null condition L^2/r0^3 = E^2/(r0-2M) for L = bE, r0 = 3M', (3*sp.sqrt(3))**2/sp.Integer(27) - 1/(sp.Integer(3) - 2))
    rep.note('Eq. (1)-(2): the source coefficients A, C, F are definitions (transcribed from Eq. (2)); they are tested\n'
             'below through everything built on them: Eq. (5) [Sec. III A], and the fluxes of the independent integrator\n'
             'against the stored Zerilli/RW data (photon_big_results.m) and the stored Teukolsky data (kerr_results.json, a = 0).')
    # Wronskian W = 2 i w A_inc from the asymptotic forms
    w, rs, Ainc, Aref = sp.symbols('omega r_* A_inc A_ref')
    PsiIn = Ainc*sp.exp(-I*w*rs) + Aref*sp.exp(I*w*rs); PsiUp = sp.exp(I*w*rs)
    rep.ident('W = Psi_in d Psi_up - Psi_up d Psi_in = 2 i omega A_inc (at infinity)',
              sp.simplify(PsiIn*sp.diff(PsiUp, rs) - PsiUp*sp.diff(PsiIn, rs) - 2*I*w*Ainc))
    # independent integrator vs stored data (photon source: E0 = 1, L0 = b, Omega = 1/b)
    b = 3*math.sqrt(3)
    worst_big = 0; worst_teuk = 0
    for lv in (10, 20):
        for jv in range(5):
            md = zrw_mode(lv, lv - jv, 3.0, 1.0, b, 1/b)
            rw = big_row(lv, jv)
            worst_big = max(worst_big, abs(md['I']/rw[2] - 1), abs(md['H']/rw[3] - 1))
            if jv <= 3:
                kr = kerr_rows(0.0, 1, lv)[jv]
                worst_teuk = max(worst_teuk, abs(md['I']/kr['FluxI'] - 1), abs(md['H']/kr['FluxH'] - 1))
    rep.check('integrator vs stored Zerilli/RW fluxes, l=10,20, j<=4 (max rel. dev.)', worst_big, '0', 2e-6, 'abs',
              'tests Eq. (2), jumps, Eq. (4), c_pm')
    rep.check('integrator vs stored Teukolsky null-source fluxes (a=0), l=10,20, j<=3', worst_teuk, '0', 1e-5, 'abs',
              'Sec. II last paragraph: Teukolsky = Zerilli/RW')
    # ---- stored Toolkit cross-check at l <= 5 (toolkit_checks/schwarzschild_toolkit_results.m, rows {l, m, method, Edot_I, Edot_H}):
    # "RW-MST" = Zerilli/RW fluxes with the Toolkit's MST solutions and the null source of Sec. II; "Teuk" = the same null
    # stress-energy tensor in the Toolkit's Teukolsky solver (the check quoted at the end of Sec. II and in Sec. V);
    # "MP" = the production integrator of Sec. III D; "NI" = an arbitrary-precision direct integration (not quoted).
    if SCHW_TK is None:
        rep.unavailable('Teukolsky vs Zerilli/RW agreement to 9-20 digits at l<=5 (Toolkit MST solutions)', 'toolkit_checks/schwarzschild_toolkit_results.m not found')
        rep.unavailable('production integrator vs Toolkit: 1e-7 at l = 5, 1e-5 at l = 2 (Sec. III D)', 'toolkit_checks/schwarzschild_toolkit_results.m not found')
    else:
        TK = {(lv, mv, meth): (EI, EH) for lv, mv, meth, EI, EH in SCHW_TK[0]}
        pairs = sorted(set((lv, mv) for lv, mv, meth in TK))
        lvals = sorted(set(lv for lv, mv in pairs))
        rep.info('stored Toolkit run: (l, m) pairs; precision marks (digits) of the stored MST fluxes: min, max',
                 f'{pairs}; {min(SCHW_TK[1]):.1f}, {max(SCHW_TK[1]):.1f}')
        def digits(x, y): return min(-math.log10(abs(x/y - 1)), 20.0) if x != y else 20.0   # stored with 20 digits
        dig = {(lv, mv): min(digits(TK[(lv, mv, 'Teuk')][0], TK[(lv, mv, 'RW-MST')][0]), digits(TK[(lv, mv, 'Teuk')][1], TK[(lv, mv, 'RW-MST')][1]))
               for lv, mv in pairs}
        digmax = {(lv, mv): max(digits(TK[(lv, mv, 'Teuk')][0], TK[(lv, mv, 'RW-MST')][0]), digits(TK[(lv, mv, 'Teuk')][1], TK[(lv, mv, 'RW-MST')][1]))
                  for lv, mv in pairs}
        rep.check(f'Sec. II: Teukolsky vs Zerilli/RW with the Toolkit solutions, l = {lvals}, m = l, l-1: worst agreement in significant digits (paper: "between nine and seventeen")',
                  min(dig.values()), '9', 1e-2, note='per (l, m): ' + ', '.join(f'{k}: {v:.1f}' for k, v in dig.items()))
        rep.info('  ... best agreement in significant digits (the paper\'s "seventeen"; this count caps at the 20 stored digits, the notebook\'s finer count gives 17.3)',
                 f'{max(digmax.values()):.1f}', note='per (l, m), best of the two fluxes: ' + ', '.join(f'{k}: {v:.1f}' for k, v in digmax.items()))
        for lv in lvals:
            w = max(abs(TK[(lv, mv, 'MP')][k]/TK[(lv, mv, 'RW-MST')][k] - 1) for (l2, mv) in pairs if l2 == lv for k in (0, 1))
            if lv == 5: rep.check('Sec. III D: production integrator (MP) vs Toolkit MST at l = 5, both m, both fluxes: max rel. dev. (paper: 10^-7; bound 5e-7)', w, '5e-7', 0, 'max', note=f'{w:.2e}')
            elif lv == 2: rep.check('Sec. III D: production integrator (MP) vs Toolkit MST at l = 2 (paper: 10^-5)', w, '1e-5', 0, 'max', note=f'{w:.2e}')
            else: rep.info(f'  ... production integrator (MP) vs Toolkit MST at l = {lv} (not quoted)', w)
        rep.info('  ... arbitrary-precision direct integration (NI) vs Toolkit MST, max rel. dev. over all stored (l, m) (not quoted)',
                 max(abs(TK[(lv, mv, 'NI')][k]/TK[(lv, mv, 'RW-MST')][k] - 1) for lv, mv in pairs for k in (0, 1)))
        # this script's integrator against the Toolkit MST fluxes (wider domain than the default, so that the WKB data are accurate)
        w = 0
        for lv, mv in pairs:
            md = zrw_mode(lv, mv, 3.0, 1.0, b, 1/b, rsA=-80.0, rB=120.0)
            w = max(w, abs(md['I']/TK[(lv, mv, 'RW-MST')][0] - 1), abs(md['H']/TK[(lv, mv, 'RW-MST')][1] - 1))
        rep.check(f'this script\'s integrator (zrw_mode, r* from -80 M to r = 120 M) vs Toolkit MST fluxes, l = {lvals}: max rel. dev.', w, '0', 1e-7, 'abs')

# =============================================================================================
# Section III A
# =============================================================================================
def section_IIIA():
    rep.section('Sec. III A  The source  [Eqs. (5)-(6)]')
    b = 3*sp.sqrt(3)
    PsiE, dPsiE = jumps(l, m, 3, 1, b, Y, dY, True)
    rep.ident('Eq. (5) even: [Psi] = 8 pi E Y/(l(l+1))  (exact)', PsiE - 8*pi*Y/(l*(l + 1)))
    ser = series_in_l(dPsiE.subs(m, l - j), 3)
    rep.ident('Eq. (5) even: [d_r Psi] = -8 pi (2j+1) E Y/(l M) + O(l^-2)  (coefficient of 1/l)',
              ser.coeff(eps, 1) + 8*pi*(2*j + 1)*Y)
    rep.ident('Eq. (5) even: no O(l^0) term in [d_r Psi]', ser.coeff(eps, 0))
    PsiO, dPsiO = jumps(l, m, 3, 1, b, Y, dY, False)
    rep.ident('Eq. (5) odd: [Psi] = -16 sqrt3 pi E dY/(l(l+1))  (exact)', PsiO + 16*sp.sqrt(3)*pi*dY/(l*(l + 1)))
    rep.ident('Eq. (5) odd: [d_r Psi] = 16 pi E dY/(sqrt3 l(l+1) M)  (exact)', dPsiO - 16*pi*dY/(sp.sqrt(3)*l*(l + 1)))
    rep.check('[d_r* Psi] = f(3M)[d_r Psi] = [d_r Psi]/3', float(f_schw(3)), '1/3', 1e-15)
    # Eq. (6): O(l^0) part of the even derivative jump for generic E0, L0, r0, m = l
    psiGen, dPsiGen = jumps(l, l, r0, E0, L0, Y, 0, True)
    lead = series_in_l(dPsiGen, 1).coeff(eps, 0)
    rep.ident('Eq. (6): O(l^0) [d_r Psi] = 4 pi Y/E (L^2/r0^3 - E^2/(r0-2M))', lead - 4*pi*Y/E0*(L0**2/r0**3 - E0**2/(r0 - 2)))
    rep.ident('Eq. (6) no negative powers of l (no O(l) term)', series_in_l(dPsiGen, 1).coeff(eps, -1))
    massive = lead.subs({L0: r0/sp.sqrt(r0 - 3), E0: (r0 - 2)/sp.sqrt(r0*(r0 - 3))})
    rep.ident('massive particle: O(l^0) term = -4 pi Y mu^2/(E r0)  [per mu: -4 pi Y/(E0 r0)]',
              massive + 4*pi*Y/(((r0 - 2)/sp.sqrt(r0*(r0 - 3)))*r0))
    leadA = lead.subs(L0, 0); leadF = sp.simplify(lead - leadA)
    ph = {r0: 3, E0: 1, L0: b}
    rep.check('photon: A (energy density) contribution / (4 pi Y)', float(leadA.subs(ph)/(4*pi*Y)), '-1', 0)
    rep.check('photon: F (azimuthal pressure) contribution / (4 pi Y)', float(leadF.subs(ph)/(4*pi*Y)), '1', 0)
    rep.ident('photon: A + F contributions cancel (p.p = 0)', (leadA + leadF).subs(ph))
    # Abstract, Introduction, Sec. VI B: "without this cancellation the flux per multipole would grow like l" (a factor l^2
    # with respect to a generic source). With a generic O(l^0) derivative jump J0 Y (e.g. the massive-particle term
    # -4 pi Y mu^2/(E r0)) in place of the photon's -8 pi (2j+1) Y/(l M), the even flux l c_+/4 |u(0)|^2 ([d_r* Psi])^2,
    # with |u(0)|^2 ~ l^{1/2} (Eq. 8) and |Y|^2 ~ l^{1/2}, grows like l.
    J0, Fs, cs = sp.symbols('J_0 F c_j', positive=True)
    u0sq = sp.sqrt(l/2)*pi*Fs/sp.sqrt(2); cplus = (l - 1)*(l + 2)/(4*pi*l*(l + 1)); Ysq = cs*sp.sqrt(l)/(2*pi**sp.Rational(3, 2))
    flux_gen = cplus/4*u0sq*(J0/3)**2*Ysq                       # [d_r* Psi] = J0 Y/3, |Y|^2 -> Ysq
    flux_photon = cplus/4*u0sq*(8*pi*(2*j + 1)/(3*l))**2*Ysq
    rep.ident('without the cancellation (generic O(l^0) jump J0 Y): the flux per multipole grows like l^{+1} (exponent of l)',
              sp.limit(sp.log(flux_gen)/sp.log(l), l, sp.oo) - 1)
    rep.ident('  ... i.e. l^2 times the photon flux: lim flux_generic/(l^2 flux_photon) = [J0/(8 pi (2j+1))]^2',
              sp.limit(flux_gen/(l**2*flux_photon), l, sp.oo) - (J0/(8*pi*(2*j + 1)))**2)
    rep.ident('massive particle: the O(l^0) term -4 pi Y mu^2/(E r0) equals the photon term -8 pi Y/(l M) (j = 0) at l = 2 E^2 r0/(mu^2 M) = 6 gamma^2 (Introduction: the ingredients fail beyond l ~ gamma^2)',
              sp.simplify((2*E0**2*3/(sp.Symbol('mu2'))).subs(E0**2, sp.Symbol('gamma2')*sp.Symbol('mu2')) - 6*sp.Symbol('gamma2')))
    # how A and F enter the Zerilli jumps at leading order (text before Eq. 6)
    A, C, F = source_ACF(l, l, r0, E0, L0, Y, 0)
    coefA = sp.simplify(dPsiGen.subs(L0, 0)/A)
    coefF = sp.simplify((dPsiGen - dPsiGen.subs(L0, 0))/F)
    rep.ident('[Psi] ~ -r0^3 A/(2 l^2 (r0-2M))  (leading)', series_in_l(sp.simplify(psiGen/A), 3).coeff(eps, 2) + r0**3/(2*(r0 - 2)))
    rep.ident('[d_r Psi] ~ r0^3 A/(4 (r0-2M)^2)  (A coefficient, leading)', series_in_l(coefA, 1).coeff(eps, 0) - r0**3/(4*(r0 - 2)**2))
    rep.ident('[d_r Psi] ~ l^2 r0^2 F/(4 (r0-2M))  (F coefficient, leading)', sp.limit(coefF/l**2, l, sp.oo) - r0**2/(4*(r0 - 2)))
    rep.ident('C sources only the Regge-Wheeler sector (even jumps independent of C)', sp.diff(dPsiGen, dY))
    rep.ident('Zerilli jumps independent of dY; RW jumps independent of Y', sp.diff(PsiO, Y) + sp.diff(dPsiO, Y))

# =============================================================================================
# Section III B
# =============================================================================================
def section_IIIB():
    rep.section('Sec. III B  The barrier top  [Eqs. (7)-(9)]')
    for name, V in (('Zerilli', V_Z(l, r)), ('Regge-Wheeler', V_RW(l, r))):
        V0 = V.subs(r, 3); kk = -Dst(Dst(V)).subs(r, 3); V3 = Dst(Dst(Dst(V))).subs(r, 3)
        rep.check(f'{name}: V0 -> [l(l+1) + O(1)]/27: lim 27 V0/l^2', float(sp.limit(27*V0/l**2, l, sp.oo)), '1', 0)
        c0 = sp.limit(27*V0 - l*(l + 1), l, sp.oo)
        rep.check(f'{name}: 27 V0 - l(l+1) -> O(1) constant', float(c0), str(c0), 0, note='finite limit')
        rep.check(f'{name}: k -> [2 l^2 + O(l)]/729: lim 729 k/l^2', float(sp.limit(729*kk/l**2, l, sp.oo)), '2', 0)
        rep.check(f'{name}: V3 -> [4 l^2 + O(l)]/6561: lim 6561 V3/l^2', float(sp.limit(6561*V3/l**2, l, sp.oo)), '4', 0)
        epsilon = (l - j)**2/27 - V0
        rep.check(f'{name}: Eq. (7) 27 eps + (2j+1) l -> O(1) (lim at j=0)', float(sp.limit((27*epsilon + (2*j + 1)*l).subs(j, 0), l, sp.oo)),
                  str(sp.limit((27*epsilon + (2*j + 1)*l).subs(j, 0), l, sp.oo)), 0, note='finite limit')
        eta = -epsilon/sp.sqrt(2*kk)
        ser = series_in_l(eta, 2)
        rep.ident(f'{name}: Eq. (7) eta_j = j + 1/2 + O(1/l)', sp.simplify(ser.coeff(eps, 0) - j - sp.Rational(1, 2)))
        rep.ident(f'{name}: eta_j has no O(l) term', ser.coeff(eps, -1))
        # Secs. II and III B: V = l^2 f/r^2 + O(l) as a function of r (not only its Taylor coefficients at r = 3M)
        rep.ident(f'{name}: V = l^2 f/r^2 + O(l) for every r: lim V/l^2 = f/r^2', sp.simplify(sp.limit(V/l**2, l, sp.oo) - f_schw(r)/r**2))
        rep.ident(f'{name}: (V - l^2 f/r^2)/l -> f/r^2 for every r (the remainder is O(l), not larger)',
                  sp.simplify(sp.limit((V - l**2*f_schw(r)/r**2)/l, l, sp.oo) - f_schw(r)/r**2))
        # the expansion of Sec. III B is about r = 3M, the maximum of the leading term: the linear term dV/dr* at 3M is O(l^0)
        V1 = Dst(V).subs(r, 3)
        rep.ident(f'{name}: linear term dV/dr* at r = 3M is O(l^0): lim = 2/(243 M^3)', sp.limit(V1, l, sp.oo) - sp.Rational(2, 243))
        rep.ident(f'{name}: in Weber\'s variable the linear term V1 z/(2k)^{{3/4}} is O(l^{{-3/2}}) z (negligible): l^{{3/2}} V1/(2k)^{{3/4}} -> 1/sqrt6',
                  sp.simplify(sp.limit(l**sp.Rational(3, 2)*V1/(4*l**2/sp.Integer(729))**sp.Rational(3, 4), l, sp.oo) - 1/sp.sqrt(6)))
        rep.ident(f'{name}: the maximum of V is at r = 3M + M/l^2 + O(l^-3): shift x_max = V1/k in r*, f(3M) V1/k in r: l^2 f V1/k -> 1',
                  sp.limit(l**2*V1/kk*f_schw(3), l, sp.oo) - 1)
        rmax = sp.nsolve(sp.diff(V, r).subs(l, 1000), r, 3)
        rep.check(f'{name}: l^2 (r_max - 3M) at l = 1000 (numerical maximum of V)', float((rmax - 3)*10**6), '1', 2e-3)
    # Weber's equation from the parabolic barrier
    x, z, kk_, ep_ = sp.symbols('x z k epsilon_b', positive=True)
    u = sp.Function('u')
    eta_ = -ep_/sp.sqrt(2*kk_)
    lhs = sp.diff(u(z), z, 2) + (z**2/4 - eta_)*u(z)                      # Weber
    # u_xx + (eps + k x^2/2) u with x = z/(2k)^(1/4):  d/dx = (2k)^(1/4) d/dz
    phys = (2*kk_)**sp.Rational(1, 2)*sp.diff(u(z), z, 2) + (ep_ + kk_*(z/(2*kk_)**sp.Rational(1, 4))**2/2)*u(z)
    rep.ident("Weber: u_xx + (eps + k x^2/2) u = sqrt(2k) [u_zz + (z^2/4 - eta) u], z = (2k)^(1/4) x",
              sp.simplify(phys - sp.sqrt(2*kk_)*lhs))
    # omega/(2k)^(1/4) = sqrt(l/2),  omega (2k)^(1/4) = sqrt2 l^(3/2)/27 at leading order
    om = l/(3*sp.sqrt(3)); k2 = 4*l**2/729
    rep.ident('omega/(2k)^(1/4) = sqrt(l/2)', sp.simplify(om/k2**sp.Rational(1, 4) - sp.sqrt(l/2)))
    rep.ident('omega (2k)^(1/4) = sqrt2 l^(3/2)/(27 M^2)', sp.simplify(om*k2**sp.Rational(1, 4) - sp.sqrt(2)*l**sp.Rational(3, 2)/27))
    # cubic term  V3 z^3/(2k)^(5/4) = O(l^(-1/2)) z^3
    gco = (4*l**2/6561)/(6*k2**sp.Rational(5, 4))
    rep.info('cubic coefficient g = V3/(6 (2k)^(5/4)) = g0 l^(-1/2): g0 (used in Sec. III C)', float(sp.limit(gco*sp.sqrt(l), l, sp.oo)))
    rep.ident('g l^(1/2) is l-independent', sp.simplify(sp.diff(gco*sp.sqrt(l), l)))
    # ---- the parabolic region and the footnote on the cubic phase shift (leading-order k, V3)
    xx, kk3, V3s, w2 = sp.symbols('x k V3 omega2', positive=True)
    V3lead = 4*l**2/sp.Integer(6561)
    p0 = sp.sqrt(kk3/2)*xx                                   # p0 = (k/2)^{1/2} x for |x| >> k^{-1/4}
    shift = sp.integrate(-V3s*xx**3/6/(2*p0), (xx, 0, xx))   # int delta/(2 p0) dx with delta = -V3 x^3/6
    rep.ident('footnote: WKB phase shift int delta/(2 p0) dx = -(V3/36)(2/k)^{1/2} x^3 for p0 = (k/2)^{1/2} x',
              sp.simplify(shift + V3s/36*sp.sqrt(2/kk3)*xx**3))
    klead = k2/2                                             # k2 = 2k (leading order), see above
    rep.ident('footnote: sqrt(k)/V3 = (243 sqrt2/4) M/l  (the criterion |x|^3 << sqrt(k)/V3 is proportional to 1/l)',
              sp.simplify(sp.sqrt(klead)/V3lead - 243*sp.sqrt(2)/(4*l)))
    rep.ident('footnote: local condition |x| << k/V3 = 9M/2, i.e. ~ M, independent of l', sp.simplify(klead/V3lead - sp.Rational(9, 2)))
    rep.ident('parabolic region (sqrt(k)/V3)^{1/3} ~ M l^{-1/3}: l^{1/3} (sqrt(k)/V3)^{1/3} is l-independent',
              sp.simplify(sp.diff(l**sp.Rational(1, 3)*(sp.sqrt(k2)/V3lead)**sp.Rational(1, 3), l)))
    rep.ident('parabolic region in z: (2k)^{1/4} (sqrt(k)/V3)^{1/3} ~ l^{1/6}  (|z| << l^{1/6})',
              sp.simplify(sp.diff(l**sp.Rational(-1, 6)*(2*k2)**sp.Rational(1, 4)*(sp.sqrt(k2)/V3lead)**sp.Rational(1, 3), l)))
    # the hypothesis of the footnote, "for |x| >> k^{-1/4} the quadratic term dominates, p0 ~ (k/2)^{1/2} x": with
    # omega^2 - V0 = eps = -eta sqrt(2k) (Eq. 7) and eta = O(1), the ratio of the constant to the quadratic term is O(1/(sqrt(k) x^2))
    eta_s = sp.symbols('eta', positive=True)
    eps_s = -eta_s*sp.sqrt(2*kk3)
    rep.ident('footnote hypothesis: (omega^2 - V0)/(k x^2/2) = -2 sqrt2 eta/(sqrt(k) x^2), i.e. O(1) only for |x| ~ k^{-1/4}',
              sp.simplify(eps_s/(kk3*xx**2/2) + 2*sp.sqrt(2)*eta_s/(sp.sqrt(kk3)*xx**2)))
    rep.ident('  ... p0^2/[(k/2) x^2] = 1 - 4 eta/z^2 with z = (2k)^{1/4} x: p0 -> (k/2)^{1/2} x for |z| >> 1',
              sp.simplify((eps_s + kk3*xx**2/2)/(kk3*xx**2/2) - (1 - 4*eta_s/((2*kk3)**sp.Rational(1, 4)*xx)**2)))
    # ---- Lyapunov exponent and the eikonal quasinormal frequency of the parabolic barrier
    rep.ident('lambda_L = sqrt(2k)/(2 omega) = Omega = 1/(3 sqrt3 M) at leading order', sp.simplify(sp.sqrt(k2)/(2*om) - 1/(3*sp.sqrt(3))))
    V0Z = V_Z(l, r).subs(r, 3)
    rep.ident('Re omega_QNM = sqrt(V0) = (l + 1/2) Omega + O(1/l)  (Zerilli V0)',
              sp.simplify(sp.limit(sp.sqrt(V0Z) - (l + sp.Rational(1, 2))/(3*sp.sqrt(3)), l, sp.oo)))
    # QNM of the parabolic barrier: outgoing on both sides <=> eta = i (n + 1/2), nu = n; then
    # epsilon = -eta sqrt(2k) and omega^2 = V0 + epsilon give Im omega = -(n + 1/2) sqrt(2k)/(2 omega)
    n_, V0s, ks = sp.symbols('n V0 k', positive=True)
    eta_qnm = I*(n_ + sp.Rational(1, 2))
    rep.ident('QNM condition: nu = -1/2 - i eta = n (integer) for eta = i(n + 1/2)', sp.simplify(-sp.Rational(1, 2) - I*eta_qnm - n_))
    om_qnm = sp.sqrt(V0s - eta_qnm*sp.sqrt(2*ks))
    rep.ident('Im omega_QNM = -(n + 1/2) sqrt(2k)/(2 omega) for the parabolic barrier (first order in sqrt(k)/V0)',
              sp.simplify(sp.series(om_qnm, ks, 0, 1).removeO() - sp.sqrt(V0s) + I*(n_ + sp.Rational(1, 2))*sp.sqrt(2*ks)/(2*sp.sqrt(V0s))))
    rep.ident('n = 0: D_0(-e^{-i pi/4} z) = exp(i z^2/4), outgoing toward both z -> +inf and z -> -inf',
              sp.simplify(sp.exp(-(-sp.exp(-I*pi/4)*z)**2/4) - sp.exp(I*z**2/4)))
    rep.ident('the mode (l, l-j) lies (j + 1/2) Omega below Re omega_QNM = (l + 1/2) Omega', (l - j)/(3*sp.sqrt(3)) - ((l + sp.Rational(1, 2)) - (j + sp.Rational(1, 2)))/(3*sp.sqrt(3)))

    # --- parabolic cylinder functions: D_nu(-e^{-i pi/4} z), nu = -1/2 - i eta, solves Weber's equation
    mp.mp.dps = 30
    def Dsol(eta, zz):
        nu = -mp.mpf(1)/2 - 1j*eta
        return mp.pcfd(nu, -mp.exp(-1j*mp.pi/4)*zz)
    resid = 0
    for eta in (mp.mpf(1)/2, mp.mpf(3)/2):
        for zz in (mp.mpf('-2.3'), mp.mpf('0.7'), mp.mpf('3.1')):
            val = Dsol(eta, zz)
            d2 = mp.diff(lambda t: Dsol(eta, t), zz, 2)
            resid = max(resid, abs(d2 + (zz**2/4 - eta)*val)/abs(val))
    rep.check('D_nu(-e^{-i pi/4} z) solves u_zz + (z^2/4 - eta) u = 0 (max residual)', float(resid), '0', 1e-18, 'abs')
    # matching / connection formula: incident-wave amplitude of the exact solution at large z
    def riccati_model(eta, g, sign, zz, niter=3):
        p2 = lambda t: t*t/4 - eta - g*t**3
        def Q(k, t):
            if k == 0: return sign*mp.sqrt(p2(t))
            h = mp.mpf('1e-4')
            return sign*mp.sqrt(p2(t) + 1j*(Q(k - 1, t + h) - Q(k - 1, t - h))/(2*h))
        return Q(niter, zz)
    Zb = mp.mpf(40)
    worst_C = 0; worst_pure = 0; worst_T = 0
    for eta in (mp.mpf(1)/2, mp.mpf(3)/2, mp.mpf(5)/2):
        uB = Dsol(eta, Zb); duB = mp.diff(lambda t: Dsol(eta, t), Zb)
        Qp = riccati_model(eta, 0, +1, Zb); Qm = riccati_model(eta, 0, -1, Zb)
        Ainc = (duB - 1j*Qp*uB)/(1j*(Qm - Qp)); Aref = (duB - 1j*Qm*uB)/(1j*(Qp - Qm))
        pB = mp.sqrt(Zb**2/4 - eta)
        C2 = 1/(abs(Ainc)**2*pB)          # |C|^2 (2k)^(1/4)/omega  for unit incident WKB amplitude
        worst_C = max(worst_C, abs(C2/(mp.exp(-mp.pi*eta/2)/mp.cosh(mp.pi*eta)) - 1))
        # purity at z -> -infinity: no wave travelling toward +z
        uA = Dsol(eta, -Zb); duA = mp.diff(lambda t: Dsol(eta, t), -Zb)
        Qp2 = riccati_model(eta, 0, +1, -Zb); Qm2 = riccati_model(eta, 0, -1, -Zb)
        Aright = (duA - 1j*Qm2*uA)/(1j*(Qp2 - Qm2)); Aleft = (duA - 1j*Qp2*uA)/(1j*(Qm2 - Qp2))
        worst_pure = max(worst_pure, abs(Aright/Aleft))
        T2 = abs(Aleft)**2*mp.sqrt(Zb**2/4 - eta)/(abs(Ainc)**2*pB)
        worst_T = max(worst_T, abs(T2*(1 + mp.exp(2*mp.pi*eta)) - 1))
    rep.check('|C|^2 = omega e^{-pi eta/2}/[(2k)^{1/4} cosh(pi eta)] from WKB matching of D_nu (max rel. dev., eta=1/2,3/2,5/2)',
              float(worst_C), '0', 2e-6, 'abs', 'z = 40: residual WKB truncation error ~ 6e-7')
    rep.check('D_nu(-e^{-i pi/4} z) purely transmitted at z -> -inf (reflected/transmitted)', float(worst_pure), '0', 2e-6, 'abs')
    rep.check('|T|^2 = 1/(1 + e^{2 pi eta}) from the same matching (max rel. dev.)', float(worst_T), '0', 2e-6, 'abs')
    rep.check('|T|^2 at eta = 1/2', float(1/(1 + mp.exp(mp.pi))), '0.041', 2e-2)
    rep.check('reflected fraction 1 - |T|^2 at eta = 1/2 (96%)', float(1 - 1/(1 + mp.exp(mp.pi))), '0.96', 2e-3)
    # Eqs. (8)-(9): values at z = 0
    worst8 = 0; worst9 = 0
    for eta in (mp.mpf(1)/2, mp.mpf(3)/2, mp.mpf(5)/2, mp.mpf(7)/2):
        nu = -mp.mpf(1)/2 - 1j*eta
        D0 = mp.pcfd(nu, 0); D0f = 2**(nu/2)*mp.sqrt(mp.pi)/mp.gamma((1 - nu)/2)
        dD0 = mp.diff(lambda t: mp.pcfd(nu, t), 0); dD0f = -2**((nu + 1)/2)*mp.sqrt(mp.pi)/mp.gamma(-nu/2)
        worst8 = max(worst8, abs(D0/D0f - 1), abs(dD0/dD0f - 1))
        C2 = mp.exp(-mp.pi*eta/2)/mp.cosh(mp.pi*eta)       # |C|^2 (2k)^(1/4)/omega
        u0sq = C2*abs(D0)**2                              # |u(0)|^2 / (omega/(2k)^(1/4))
        Du0sq = C2*abs(dD0)**2                            # |u'(0)|^2 / (omega (2k)^(1/4))  (|d/dz D_nu(-e^{-i pi/4} z)| = |D_nu'|)
        worst9 = max(worst9, abs(u0sq/(mp.pi*Fpar(eta)/mp.sqrt(2)) - 1),
                     abs(Du0sq/(mp.sqrt(2)*mp.pi*Gpar(eta)) - 1))
    rep.check('D_nu(0), D_nu\'(0) closed forms (max rel. dev. vs mpmath pcfd)', float(worst8), '0', 1e-25, 'abs')
    rep.check('Eqs. (8)-(9): |u(0)|^2 = sqrt(l/2) pi F(eta)/sqrt2, |u\'(0)|^2 = (sqrt2 l^{3/2}/27) sqrt2 pi G(eta)',
              float(worst9), '0', 1e-25, 'abs', 'eta-dependent factors, eta = 1/2..7/2')
    # independent numerical check of the matching by integrating Weber's equation (no pcfd)
    for eta in (0.5, 1.5):
        md = weber_model(eta, 0.0, 'In')
        rep.check(f'|u(0)|^2 numerical ODE / Eq. (8), eta = {eta}', md['u0sq']/float(mp.pi*Fpar(eta)/mp.sqrt(2)), '1', 2e-6)
        rep.check(f'|u\'(0)|^2 numerical ODE / Eq. (9), eta = {eta}', md['Du0sq']/float(mp.sqrt(2)*mp.pi*Gpar(eta)), '1', 2e-6)
        rep.check(f'|T|^2 numerical ODE / (1+e^{{2 pi eta}})^-1, eta = {eta}', md['T2']*(1 + math.exp(2*math.pi*eta)), '1', 2e-6)
    # mirror symmetry of the parabolic barrier: v(z) = u(-z)
    mdU = weber_model(0.5, 0.0, 'Up')
    mdI = weber_model(0.5, 0.0, 'In')
    rep.check('|v(0)|^2 = |u(0)|^2 (parabolic barrier, g = 0)', mdU['u0sq']/mdI['u0sq'], '1', 1e-8)
    rep.check("|v'(0)|^2 = |u'(0)|^2 (parabolic barrier, g = 0)", mdU['Du0sq']/mdI['Du0sq'], '1', 1e-8)
    # d_r* v(0) = -d_r* u(0) (sign, not only modulus; it fixes the opposite signs of the cross term in the two fluxes): the up
    # solution, defined by its own boundary conditions (incident from z -> -inf), has v'(0)/v(0) = -u'(0)/u(0); the ratio does
    # not depend on the normalisation of either solution. (With v(z) = u(-z) the relation is exact.)
    for eta in (0.5, 1.5):
        mU = weber_model(eta, 0.0, 'Up'); mI = weber_model(eta, 0.0, 'In')
        rep.check(f"d_r* v(0) = -d_r* u(0): v'(0)/v(0) = -u'(0)/u(0) on the parabolic barrier (model ODE, complex ratio), eta = {eta}",
                  abs(mU['ratio']/mI['ratio'] + 1), '0', 1e-7, 'abs', note=f"v'(0)/v(0) = {mU['ratio']:.6f}, u'(0)/u(0) = {mI['ratio']:.6f}")


def weber_model(eta, g, side, Z=40.0):
    """Numerical solution of u'' + (z^2/4 - eta - g z^3) u = 0 with unit-amplitude WKB data at the
    transmitted end (side 'In': incident from z -> +inf, transmitted to -inf; 'Up': the reverse).
    Returns |u(0)|^2, |u'(0)|^2 and |T|^2 for unit incident WKB amplitude (units omega/(2k)^{1/4} = 1)."""
    p2 = lambda t: t*t/4 - eta - g*t**3
    def Q(sign, t, k=3):
        if k == 0: return sign*np.sqrt(p2(t) + 0j)
        h = 1e-3
        return sign*np.sqrt(p2(t) + 1j*(Q(sign, t + h, k - 1) - Q(sign, t - h, k - 1))/(2*h))
    if side == 'In': zA, zB, d = -Z, Z, -1
    else: zA, zB, d = Z, -Z, +1
    QA = Q(d, zA)
    y0 = np.array([1.0, 0.0, (1j*QA).real, (1j*QA).imag])
    def rhs(t, y):
        q = p2(t); return [y[2], y[3], -q*y[0], -q*y[1]]
    s1 = solve_ivp(rhs, (zA, 0.0), y0, method='DOP853', rtol=1e-12, atol=1e-14)
    y0v = s1.y[:, -1]
    s2 = solve_ivp(rhs, (0.0, zB), y0v, method='DOP853', rtol=1e-12, atol=1e-14)
    yB = s2.y[:, -1]
    u0 = y0v[0] + 1j*y0v[1]; du0 = y0v[2] + 1j*y0v[3]
    uB = yB[0] + 1j*yB[1]; duB = yB[2] + 1j*yB[3]
    Qp, Qm = Q(+1, zB), Q(-1, zB); pB = math.sqrt(p2(zB))
    Ainc = (duB - 1j*Qp*uB)/(1j*(Qm - Qp)) if side == 'In' else (duB - 1j*Qm*uB)/(1j*(Qp - Qm))
    nrm = abs(Ainc)**2*pB
    return dict(u0sq=abs(u0)**2/nrm, Du0sq=abs(du0)**2/nrm, T2=math.sqrt(p2(zA))/nrm, ratio=du0/u0)

# =============================================================================================
# Section III C
# =============================================================================================
SIGMA = {}   # filled here, used in Sec. III D

def section_IIIC():
    rep.section('Sec. III C  The fluxes  [Eqs. (10)-(14)] and the O(l^-1/2) asymmetry')
    # ---- harmonics: closed forms vs scipy, and the large-l limits with c_j, d_j
    worst = 0
    for lv in range(2, 41):
        for jv in range(0, lv - 1):
            if jv % 2 == 0:
                worst = max(worst, abs(Ysq_closed(lv, jv)/Y_eq(lv, lv - jv)**2 - 1))
            else:
                worst = max(worst, abs(dYsq_closed(lv, jv)/dY_eq(lv, lv - jv)**2 - 1))
    rep.check('|Y_{l,l-j}(pi/2)|^2 and |d_theta Y|^2 closed forms vs scipy sph_harm_y, l<=40', worst, '0', 1e-10, 'abs')
    rep.check('Y_{l,l-j}(pi/2) = 0 for odd j, d_theta Y = 0 for even j (scipy, l<=40)',
              max(abs(Y_eq(lv, lv - jv)) for lv in range(2, 41) for jv in range(1, lv, 2)) +
              max(abs(dY_eq(lv, lv - jv)) for lv in range(2, 41) for jv in range(0, lv, 2)), '0', 1e-12, 'abs')
    for jv in range(0, 13):
        cj = float(c_j(jv)) if jv % 2 == 0 else float(d_j(jv))
        L1 = 600   # largest l at which scipy's sph_harm_y is finite for m = l-4
        if jv % 2 == 0:
            lim_cf = Ysq_closed(10**7, jv)/(cj*math.sqrt(1e7)/(2*math.pi**1.5))
            rep.check(f'|Y_{{l,l-{jv}}}(pi/2)|^2 / [c_{jv} sqrt(l)/(2 pi^3/2)], c_{jv} = {cj:.6g}: closed form at l = 1e7', lim_cf, '1', 2e-6*(1 + jv))
            if jv <= 4:
                lim_sc = Y_eq(L1, L1 - jv)**2/(cj*math.sqrt(L1)/(2*math.pi**1.5))
                rep.check(f'  ... scipy sph_harm_y at l = {L1} (O(1/l) corrections)', lim_sc, '1', 1e-2)
        else:
            lim_cf = dYsq_closed(10**7, jv)/(cj*1e7**1.5/math.pi**1.5)
            rep.check(f'|d_theta Y_{{l,l-{jv}}}|^2 / [d_{jv} l^3/2/pi^3/2], d_{jv} = {cj:.6g}: closed form at l = 1e7', lim_cf, '1', 2e-6*(1 + jv))
            if jv <= 4:
                lim_sc = dY_eq(L1, L1 - jv)**2/(cj*L1**1.5/math.pi**1.5)
                rep.check(f'  ... scipy sph_harm_y at l = {L1} (O(1/l) corrections)', lim_sc, '1', 1e-2)
    # symbolic limits of the closed forms (Stirling): sympy limit of the ratio
    lsym = sp.symbols('l', positive=True)
    for jv in (0, 2, 4):
        Ysq_sym = (2*lsym + 1)/(4*pi)*sp.factorial(jv)*sp.gamma(2*lsym - jv + 1)/(4**lsym*sp.factorial(sp.Rational(jv, 2))**2*sp.gamma((2*lsym - jv)/2 + 1)**2)
        lim = sp.limit(Ysq_sym/(c_j(jv)*sp.sqrt(lsym)/(2*pi**sp.Rational(3, 2))), lsym, sp.oo)
        rep.check(f'sympy limit |Y_{{l,l-{jv}}}|^2/(c_{jv} sqrt(l)/(2 pi^3/2))', float(lim), '1', 0)
    # ---- Eqs. (11)-(12) from the flux formula (symbolic, leading order)
    Fs, Gs, cs, ds = sp.symbols('F G c_j d_j', positive=True)
    u0sq = sp.sqrt(l/2)*pi*Fs/sp.sqrt(2)                          # Eq. (8)
    Du0sq = sp.sqrt(2)*l**sp.Rational(3, 2)/27*sp.sqrt(2)*pi*Gs    # Eq. (9)
    cplus = (l - 1)*(l + 2)/(4*pi*l*(l + 1)); cminus = l*(l + 1)/(16*pi*(l - 1)*(l + 2))
    dPstar_even = -8*pi*(2*j + 1)*Y/(3*l)                           # Eq. (5)/3
    P_odd = -16*sp.sqrt(3)*pi*dY/(l*(l + 1))
    Ysq = cs*sp.sqrt(l)/(2*pi**sp.Rational(3, 2)); dYsq = ds*l**sp.Rational(3, 2)/pi**sp.Rational(3, 2)
    flux_even = (cplus/4*u0sq*dPstar_even**2).subs(Y**2, Ysq)
    flux_odd = (cminus/4*Du0sq*P_odd**2).subs(dY**2, dYsq)
    kap_even = sp.limit(l*flux_even, l, sp.oo); kap_odd = sp.limit(l*flux_odd, l, sp.oo)
    rep.ident('Eq. (11): kappa_j = sqrt(pi)/9 (2j+1)^2 c_j F(eta_j)  (even j, from l c_+/4 |u(0)|^2 [d_r* Psi]^2)',
              kap_even - sp.sqrt(pi)/9*(2*j + 1)**2*cs*Fs)
    rep.ident('Eq. (12): kappa_j = 8 sqrt(pi)/9 d_j G(eta_j)  (odd j, from l c_-/4 |u\'(0)|^2 [Psi]^2)',
              kap_odd - 8*sp.sqrt(pi)/9*ds*Gs)
    rep.ident('even sector: |u(0)|^2 [d_r* Psi]^2 ~ l^-1', sp.limit(sp.log(u0sq*dPstar_even**2/Y**2*Ysq)/sp.log(l), l, sp.oo) + 1)
    rep.ident('odd sector: |u\'(0)|^2 [Psi]^2 ~ l^-1', sp.limit(sp.log(Du0sq*P_odd**2/dY**2*dYsq)/sp.log(l), l, sp.oo) + 1)
    # ---- numerical values
    paper_k = {0: '0.0277668', 1: '0.0037670', 2: '2.7521e-4', 3: '1.6604e-5', 4: '9.23e-7'}
    tol_k = {0: 2e-6, 1: 2e-5, 2: 2e-5, 3: 3e-5, 4: 6e-4}
    for jv in range(5):
        rep.check(f'kappa_{jv}', kappa_j(jv), paper_k[jv], tol_k[jv])
    sumk = mp.fsum(kappa_j(jv) for jv in range(40))
    rep.check('sum_j kappa_j', sumk, '0.0318265', 2e-6)
    rep.check('kappa = 2 sum_j kappa_j  [Eq. (13)]', 2*sumk, '0.063653', 1e-5)
    rep.check('kappa agrees with the 2021 fit 0.064 +- 0.001', 2*sumk, '0.064', 1e-3/0.064)
    rep.info('kappa_{j+1}/kappa_j for j = 0..3 (fall-off ~ e^{-pi j}; e^{-pi} = 0.0432)',
             ', '.join(f'{float(kappa_j(jv + 1)/kappa_j(jv)):.4f}' for jv in range(4)))
    # "The kappa_j fall off roughly like e^{-pi j} (up to powers of j)": with Stirling's |Gamma(x+iy)|^2 ~ 2 pi |y|^{2x-1} e^{-pi|y|},
    # F(eta) ~ e^{-pi eta}/(pi sqrt(eta/2)) and G(eta) ~ e^{-pi eta} sqrt(eta/2)/pi, and with c_j ~ sqrt(2/(pi j)), d_j ~ sqrt(2j/pi),
    # Eqs. (11) and (12) both give kappa_j -> (8/(9 pi)) e^{-pi/2} j e^{-pi j} [1 + O(1/j)], the same for even and odd j.
    mp.mp.dps = 60
    kap_lim = 8/(9*mp.pi)*mp.exp(-mp.pi/2)
    r200 = kappa_j(200)*mp.exp(200*mp.pi)/200; r400 = kappa_j(400)*mp.exp(400*mp.pi)/400
    rep.check('kappa_j e^{pi j}/j -> (8/(9 pi)) e^{-pi/2} = 0.058818 (Richardson in 1/j from j = 200, 400; both parities)', 2*r400 - r200, kap_lim, 1e-4,
              note=f'kappa_j e^(pi j)/j = {float(r200):.6f} (j=200), {float(r400):.6f} (j=400)')
    rep.check('  ... odd j: kappa_401 e^{401 pi}/401 vs the same limit (O(1/j) corrections)', kappa_j(401)*mp.exp(401*mp.pi)/401, kap_lim, 3e-3)
    rep.check('e^{pi} kappa_{j+1}/kappa_j = 1 + 1/j + O(j^-2): [e^pi kappa_401/kappa_400]/(1 + 1/400) - 1', mp.exp(mp.pi)*kappa_j(401)/kappa_j(400)/(1 + mp.mpf(1)/400) - 1, '0', 1e-5, 'abs')
    rep.info('  ... e^{pi} kappa_{j+1}/kappa_j at j = 4, 10, 50 (-> 1: the fall-off is e^{-pi j} up to powers of j)',
             ', '.join(f'{float(mp.exp(mp.pi)*kappa_j(jv + 1)/kappa_j(jv)):.4f}' for jv in (4, 10, 50)))
    mp.mp.dps = 30
    rep.check('fraction of sum_j kappa_j missed by j <= 4 (Sec. III D: 2e-6)', float(mp.fsum(kappa_j(jv) for jv in range(5, 40))/sumk), '2e-6', 0.3)
    rep.check('fraction of sum_j kappa_j missed by j <= 3 (Sec. IV: 3e-5)', float(mp.fsum(kappa_j(jv) for jv in range(4, 40))/sumk), '3e-5', 0.1)
    # ---- Eq. (14): phase ratio
    worst = 0
    for eta in (mp.mpf(1)/2, mp.mpf(3)/2, mp.mpf(5)/2):
        nu = -mp.mpf(1)/2 - 1j*eta
        ratio = mp.pcfd(nu, 0)/(-mp.exp(-1j*mp.pi/4)*mp.diff(lambda t: mp.pcfd(nu, t), 0))   # u(0)/u_z(0)
        formula = mp.exp(1j*mp.pi/4)*mp.gamma(mp.mpf(1)/4 + 1j*eta/2)/(mp.sqrt(2)*mp.gamma(mp.mpf(3)/4 + 1j*eta/2))
        worst = max(worst, abs(ratio/formula - 1))
    rep.check('Eq. (14): (2k)^{1/4} u(0)/u\'(0) = e^{i pi/4} Gamma(1/4+i eta/2)/(sqrt2 Gamma(3/4+i eta/2))', float(worst), '0', 1e-25, 'abs')
    phase0 = mp.arg(mp.exp(1j*mp.pi/4)*mp.gamma(mp.mpf(1)/4 + 1j/4)/mp.gamma(mp.mpf(3)/4 + 1j/4))
    rep.info('phase of u(0)/u\'(0) at eta = 1/2 (radians)', float(phase0))
    # ---- sigma_j: cross term + cubic barrier asymmetry
    g0 = float((4/sp.Integer(6561))/(6*(4/sp.Integer(729))**sp.Rational(5, 4)))      # g = g0 l^(-1/2)
    def cross_term(jv):
        """sqrt(l) x cross-term contribution to (I-H)/(I+H), from the phase of Eq. (14)."""
        eta = mp.mpf(jv) + mp.mpf(1)/2
        phi = mp.arg(mp.exp(1j*mp.pi/4)*mp.gamma(mp.mpf(1)/4 + 1j*eta/2)/mp.gamma(mp.mpf(3)/4 + 1j*eta/2))
        F = Fpar(eta); G = Gpar(eta)
        # Eqs. (8)-(9): |u'(0)|^2/|u(0)|^2 = (2k)^{1/2} (sqrt2 pi G)/(pi F/sqrt2) = (2k)^{1/2} 2G/F, (2k)^{1/4} = sqrt(2l/27)
        if jv % 2 == 0:   # (I-H)/(I+H) = -2 (P/dP) |u'|/|u| cos(phi),  P/dP -> -3/((2j+1) l)
            return float(6/(2*jv + 1)*mp.sqrt(mp.mpf(2)/27)*mp.sqrt(2*G/F)*mp.cos(phi))
        # odd: (I-H)/(I+H) = -2 (dP/P) |u|/|u'| cos(phi), dP/P = -1/9
        return float(mp.mpf(2)/9*mp.sqrt(mp.mpf(27)/2)*mp.sqrt(F/(2*G))*mp.cos(phi))
    def cubic_term(jv):
        """sqrt(l) x barrier-asymmetry contribution: (|u0|^2-|v0|^2)/(|u0|^2+|v0|^2) (even j) or the same
        with u'(0) (odd j), from the cubic model equation, linear coefficient in g times g0."""
        eta = jv + 0.5; key = 'u0sq' if jv % 2 == 0 else 'Du0sq'
        gs = np.array([0.001, 0.002, 0.003, 0.005])*(min(1.0, 4.5/eta))**1.5; asym = []   # linear regime: g (2 sqrt(eta))^3 small
        for g in gs:
            a = weber_model(eta, g, 'In')[key]; b = weber_model(eta, g, 'Up')[key]
            asym.append((a - b)/(a + b)/g)
        co = np.polyfit(gs, asym, 2)          # quadratic in g, extrapolate to g -> 0
        return co[-1]*g0, asym
    cross0 = cross_term(0)
    rep.check('cross term of sigma_0 (analytic, from Eq. 14)', cross0, '1.0390', 1e-4)
    rep.check('cross term of sigma_0 (paper text: +1.04)', cross0, '1.04', 2e-3)
    cub0, asym0 = cubic_term(0)
    rep.check('cubic-barrier contribution to sigma_0 (model equation, g -> 0 extrapolation)', cub0, '-0.231', 5e-3,
              note='asym/g at g=0.001..0.005: ' + ', '.join(f'{a:.4f}' for a in asym0))
    sig0 = cross0 + cub0
    rep.check('sigma_0 = cross + cubic', sig0, '0.81', 1e-2)
    SIGMA[0] = sig0; SIGMA[('cross', 0)] = cross0; SIGMA[('cubic', 0)] = cub0
    for jv in range(1, 9):
        cr = cross_term(jv); cu, _ = cubic_term(jv); SIGMA[jv] = cr + cu; SIGMA[('cross', jv)] = cr; SIGMA[('cubic', jv)] = cu
        rep.check(f'sigma_{jv} = cross ({cr:+.4f}) + cubic ({cu:+.4f}) is negative (paper: sigma_j < 0 for j >= 1)',
                  float(SIGMA[jv] < 0), '1', 0, note=f'sigma_{jv} = {SIGMA[jv]:+.4f}' + ('' if jv <= 4 else ' (j > 4: beyond the stored data)'))
    sigbar = sum(float(kappa_j(jv))*SIGMA[jv] for jv in range(5))/sum(float(kappa_j(jv)) for jv in range(5))
    rep.check('sigma_bar = sum_j kappa_j sigma_j / sum_j kappa_j (analytic, j<=4)', sigbar, '0.61', 2e-2)
    SIGMA['bar'] = sigbar
    rep.check('2 sigma_bar/sqrt(l) at l = 100 (12%)', 2*0.61/10, '0.12', 2e-2)
    def horizon_fraction(lmax, sb=0.61):
        ls = np.arange(2, lmax + 1)
        return np.sum((1 - sb/np.sqrt(ls))/ls)/np.sum(2.0/ls)
    rep.check('horizon fraction, asymptotic formulae, l_max = 100 (paper: 40% for l_max = 100-140)', horizon_fraction(100), '0.40', 1.5e-2)
    rep.check('horizon fraction, asymptotic formulae, l_max = 140 (paper: 40% for l_max = 100-140)', horizon_fraction(140), '0.40', 1.5e-2)
    rep.info('horizon fraction with the asymptotic formulae at l_max = 1e3, 1e4, 1e6 (approaches 1/2 only logarithmically)',
             ', '.join(f'{horizon_fraction(n):.3f}' for n in (1000, 10000, 1000000)))
    # Secs. III C and VII: "the fraction approaches 1/2 only logarithmically". With the asymptotic formulae,
    # f_H = 1/2 - (sigma_bar/2) sum_l l^{-3/2} / sum_l l^{-1}: the numerator converges to zeta(3/2) - 1 while the denominator
    # is the harmonic sum ~ ln l_max, so (1/2 - f_H) sum_{l<=l_max} 1/l -> sigma_bar (zeta(3/2) - 1)/2 = 0.4918.
    rep.ident('sum_{l>=2} l^{-3/2} = zeta(3/2) - 1', sp.summation(l**sp.Rational(-3, 2), (l, 2, sp.oo)) - (sp.zeta(sp.Rational(3, 2)) - 1))
    Hsum = lambda n: float(np.sum(1/np.arange(2, n + 1.0)))
    hf_const = 0.61*(float(mp.zeta(1.5)) - 1)/2
    rep.check('(1/2 - horizon fraction) x sum_{l<=l_max} 1/l -> sigma_bar (zeta(3/2)-1)/2 = 0.4918 at l_max = 1e6 (1/ln law; O(l_max^-1/2) corrections)',
              (0.5 - horizon_fraction(10**6))*Hsum(10**6), hf_const, 3e-3,
              note='(1/2 - f_H) ln l_max at l_max = 1e3, 1e4, 1e6: ' + ', '.join(f'{(0.5 - horizon_fraction(n))*math.log(n):.3f}' for n in (1000, 10000, 1000000)))

# =============================================================================================
# Section III D  (stored Schwarzschild data)
# =============================================================================================
def section_IIID():
    rep.section('Sec. III D  Numerical check: stored run photon_big_results.m (l = 10..12800, j <= 4)')
    rep.note('The raw fluxes are NOT recomputed (production run); the quoted numbers are re-derived from the stored rows.')
    ls = sorted(set(rw[0] for rw in BIG))
    mean = {lv: sum(rw[2] + rw[3] for rw in BIG if rw[0] == lv)/2 for lv in ls}   # single-m mean flux, j<=4
    # three-term fit of l x mean over 400 <= l <= 12800 in powers of l^-1/2 (the l^-1/2 term cancels in the mean)
    Lfit = np.array([lv for lv in ls if lv >= 400], float)
    yfit = np.array([lv*mean[lv] for lv in Lfit])
    Amat = np.vstack([np.ones_like(Lfit), 1/Lfit, Lfit**-1.5]).T
    co = np.linalg.lstsq(Amat, yfit, rcond=None)[0]
    rep.check('three-term fit {1, 1/l, l^-3/2} of l x mean flux, 400<=l<=12800: sum_j kappa_j', co[0], '0.0318263', 2e-6)
    rep.check('  ... against the analytic sum_j kappa_j = 0.0318265', co[0], float(mp.fsum(kappa_j(jv) for jv in range(40))), 1e-5)
    # "At l = 12800 the individual kappa_j, j <= 4, are reproduced to better than 1e-3 in relative terms, and
    #  kappa_0 and kappa_1 to 1e-4 and 2e-5" (mean of the two fluxes, whose O(l^-1/2) terms cancel)
    for jv in range(5):
        rw = big_row(12800, jv)
        ratio = 12800*(rw[2] + rw[3])/2/float(kappa_j(jv))
        rep.check(f'l=12800: |l x mean flux / kappa_{jv} - 1| (paper: better than 1e-3, relative, for j <= 4)', abs(ratio - 1), '1e-3', 0, 'max',
                  note=f'ratio = {ratio:.7f}')
    rep.check('l=12800: |l x mean flux / kappa_0 - 1| <= 1e-4 (paper: kappa_0 to 1e-4)', abs(12800*(big_row(12800, 0)[2] + big_row(12800, 0)[3])/2/float(kappa_j(0)) - 1), '1e-4', 0, 'max')
    rw = big_row(12800, 1)
    rep.check('l=12800: |l x mean flux / kappa_1 - 1| <= 2e-5 (paper: kappa_1 to 2e-5; Sec. V: Eq. (12) to 2e-5)', abs(12800*(rw[2] + rw[3])/2/float(kappa_j(1)) - 1), '2e-5', 0, 'max')
    # consistency: stored u(0), u'(0) + our jumps reproduce the stored Edot_I  (checks the stored normalisation u = Psi_in/A_inc)
    worst = 0
    for lv in (400, 6400, 12800):
        for jv in range(5):
            rw = big_row(lv, jv); mm = lv - jv
            Yv = Y_closed(lv, jv) if jv % 2 == 0 else 0.0
            dYv = dY_closed(lv, jv) if jv % 2 else 0.0
            P, dPr = jumps_num(lv, mm, 3.0, 1.0, 3*math.sqrt(3), Yv, dYv); dP = dPr/3
            EI = c_pm(lv, mm)/4*abs(rw[4]*dP - P*rw[5])**2
            worst = max(worst, abs(EI/rw[2] - 1))
    rep.check('Edot_I = (c_pm/4)|u(0)[d_r* Psi] - [Psi] u\'(0)|^2 from stored u(0), u\'(0) (l=400,6400,12800)', worst, '0', 1e-9, 'abs')
    # |u(0)|^2 and |u'(0)|^2 of the numerical solution vs Eqs. (8)-(9)
    for lv in (400, 6400, 12800):
        rw = big_row(lv, 0)
        rep.check(f'l={lv}, j=0: |u(0)|^2 / Eq. (8)', abs(rw[4])**2/(math.sqrt(lv/2)*float(mp.pi*Fpar(0.5)/mp.sqrt(2))), '1', 2/math.sqrt(lv),
                  note='tolerance 2/sqrt(l): O(l^-1/2) cubic corrections')
        rw1 = big_row(lv, 1)
        rep.check(f'l={lv}, j=1: |u\'(0)|^2 / Eq. (9)', abs(rw1[5])**2/(math.sqrt(2)*lv**1.5/27*float(mp.sqrt(2)*mp.pi*Gpar(1.5))), '1', 2/math.sqrt(lv),
                  note='tolerance 2/sqrt(l): O(l^-1/2) cubic corrections')
    # asymmetry decomposition at l = 6400, j = 0 ("the largest l for which we computed v as well"): the stored
    # file has u(0), u'(0) but not v(0), v'(0), so u and v are recomputed with the independent integrator
    # (~20 s) and checked against the stored u. With I0 = |u dP - P u'|^2, H0 = |v dP - P v'|^2, den = (|u|^2+|v|^2) dP^2:
    #   sqrt(l)(I0-H0)/(I0+H0) = A1 (barrier asymmetry, |u|^2-|v|^2) + A2 (cross term) + A3 (O(1/l), |u'|^2-|v'|^2)
    lv = 6400; rw = big_row(lv, 0)
    md = zrw_mode(lv, lv, 3.0, 1.0, 3*math.sqrt(3), 1/(3*math.sqrt(3)))
    rep.check('l=6400, j=0: independent integrator vs stored Edot_I, Edot_H, u(0), u\'(0) (max rel. dev.)',
              max(abs(md['I']/rw[2] - 1), abs(md['H']/rw[3] - 1), abs(md['u0']/rw[4] - 1), abs(md['Du0']/rw[5] - 1)), '0', 2e-5, 'abs')
    u0, Du0, v0, Dv0, P, dP = md['u0'], md['Du0'], md['v0'], md['Dv0'], md['P'], md['dP']
    I0 = abs(u0*dP - P*Du0)**2; H0 = abs(v0*dP - P*Dv0)**2; den = (abs(u0)**2 + abs(v0)**2)*dP**2; sl = math.sqrt(lv)
    A1 = sl*(abs(u0)**2 - abs(v0)**2)*dP**2/den
    A2 = sl*(-2*((u0*np.conj(Du0)).real - (v0*np.conj(Dv0)).real)*dP*P)/den
    A3 = sl*(abs(Du0)**2 - abs(Dv0)**2)*P**2/den
    tot = math.sqrt(lv)*(rw[2] - rw[3])/(rw[2] + rw[3])
    rep.check('l=6400, j=0: cross term (exact u and v from the independent integrator; paper: 1.0389)', A2, '1.0389', 2e-4,
              note=f'O(1/l) remainder A3 = {A3:.1e}')
    rep.check('l=6400, j=0: barrier asymmetry |u(0)|^2 vs |v(0)|^2 (paper: -0.231)', A1, '-0.231', 2e-3)
    rep.check('l=6400, j=0: sqrt(l)(I-H)/(I+H) from the stored fluxes = cross + asymmetry (paper: 0.81 = 1.039 - 0.231)', tot, A1 + A2 + A3, 1e-3,
              note='the two integrators differ by ~1e-5 in the fluxes, i.e. by ~3e-4 in the asymmetry')
    rep.check('l=6400, j=0: cross term vs the analytic 1.0390 of Eq. (14)', A2, '1.0390', 2e-4)
    rat_uv = (v0/Dv0)/(u0/Du0)
    rep.check("l=6400, j=0: (v(0)/v'(0))/(u(0)/u'(0)) = -1 + O(l^-1/2) from the exact solutions (d_r* v(0) = -d_r* u(0), Sec. III B)",
              abs(rat_uv + 1), '0', 1e-2, 'abs', note=f'ratio = {rat_uv:.5f}')
    # the same decomposition from the production integrator, stored in schwarzschild/asym_check_results.m (rows {l, sqrt(l)(I-H)/(I+H), A1, A2, A3})
    if ASYM is None:
        rep.unavailable('stored decomposition of sigma_0 (asym_check_results.m): cross term 1.0389 and barrier asymmetry -0.231 at l = 6400', 'schwarzschild/asym_check_results.m not found')
    else:
        AS = {int(rw[0]): rw for rw in ASYM}
        rep.check('stored (production integrator) l=6400: cross term A2 (paper: 1.0389)', AS[6400][3], '1.0389', 1e-4)
        rep.check('stored l=6400: barrier asymmetry A1 (paper: -0.231)', AS[6400][2], '-0.231', 2e-3)
        rep.check('stored l=6400: A1 + A2 + A3 = sqrt(l)(I-H)/(I+H) (= sigma_0 + O(1/l); paper: 0.81)', AS[6400][1], '0.81', 3e-3,
                  note=f'A1 + A2 + A3 - total = {AS[6400][2] + AS[6400][3] + AS[6400][4] - AS[6400][1]:.1e}')
        rep.check('stored vs this script\'s integrator at l = 6400: cross term A2', AS[6400][3], A2, 1e-4)
        rep.check('stored vs this script\'s integrator at l = 6400: barrier asymmetry A1', AS[6400][2], A1, 3e-3)
        rep.check('stored A2(l) -> 1.0390 with O(1/l) corrections: Richardson in 1/l from l = 1600, 6400 vs the analytic 1.0389513',
                  (4*AS[6400][3] - AS[1600][3])/3, SIGMA[('cross', 0)], 2e-5, note='A2 at l = 100, 400, 1600, 6400: ' + ', '.join(f'{AS[n][3]:.5f}' for n in (100, 400, 1600, 6400)))
        rep.check('stored A1(l) Richardson-extrapolated in 1/l from l = 1600, 6400 vs the cubic-barrier model (-0.231)',
                  (4*AS[6400][2] - AS[1600][2])/3, SIGMA[('cubic', 0)], 4e-3, note='A1 at l = 100, 400, 1600, 6400: ' + ', '.join(f'{AS[n][2]:.5f}' for n in (100, 400, 1600, 6400)))
    # sigma_j and sigma_bar from the stored data at l = 12800 vs the analytic values of Sec. III C
    for jv in range(5):
        rw = big_row(12800, jv)
        rep.check(f'l=12800: sqrt(l)(I-H)/(I+H) for j={jv} vs analytic sigma_{jv} = {SIGMA[jv]:+.4f}', math.sqrt(12800)*(rw[2] - rw[3])/(rw[2] + rw[3]), SIGMA[jv], 3e-2)
    for lv in (3200, 12800):
        SI = sum(rw[2] for rw in BIG if rw[0] == lv); SH = sum(rw[3] for rw in BIG if rw[0] == lv)
        rep.check(f'l={lv}: m-summed sqrt(l)(I-H)/(I+H) (sigma_bar = 0.61)', math.sqrt(lv)*(SI - SH)/(SI + SH), '0.61', 2e-2)
    SI = sum(rw[2] for rw in BIG if rw[0] == 100); SH = sum(rw[3] for rw in BIG if rw[0] == 100)
    rep.check('l=100: (I-H)/mean, stored data (paper: ~12%)', (SI - SH)/((SI + SH)/2), '0.12', 0.25)
    rep.info('l x mean flux / kappa at l = 100, 400, 12800 (Fig. 1: mean on the dashed line from l ~ 100)',
             ', '.join(f'{lv*mean[lv]/float(mp.fsum(kappa_j(jv) for jv in range(40))):.5f}' for lv in (100, 400, 12800)))
    # Fig. 1 text: "their mean sits on the dashed line already at l ~ 100"; Sec. III C: the sum of the two fluxes approaches its
    # limit with relative corrections O(1/l), the l^-1/2 term cancelling (the data are summed over j <= 4, hence sum_{j<=4} kappa_j)
    sumk4 = float(mp.fsum(kappa_j(jv) for jv in range(5)))
    dev_mean = {lv: lv*mean[lv]/sumk4 - 1 for lv in ls}
    rep.check('|l x mean flux / sum_{j<=4} kappa_j - 1| <= 1% for every stored l >= 100 (the mean on the dashed line from l ~ 100)',
              max(abs(dev_mean[lv]) for lv in ls if lv >= 100), '0.01', 0, 'max', note=f'l = 100: {dev_mean[100]:+.5f}')
    Lf = np.array([lv for lv in ls if lv >= 100], float); yf = np.array([dev_mean[lv] for lv in Lf])
    cof = np.linalg.lstsq(np.vstack([Lf**-0.5, 1/Lf, Lf**-1.5]).T, yf, rcond=None)[0]
    rep.check('mean flux: fit of l x mean/kappa - 1 with {l^-1/2, 1/l, l^-3/2} over l >= 100: |coefficient of l^-1/2| negligible vs sigma_bar = 0.61 (bound 2e-3)',
              abs(cof[0]), '2e-3', 0, 'max', note=f'coefficients: {cof[0]:+.1e} l^-1/2, {cof[1]:+.3f}/l, {cof[2]:+.3f} l^-3/2')
    lrel = {lv: lv*dev_mean[lv] for lv in ls if lv >= 50}
    rep.check('mean flux: l (l x mean/kappa - 1) is bounded, i.e. the relative corrections are O(1/l): max/min over 50 <= l <= 12800',
              max(lrel.values())/min(lrel.values()), '1', 0.1, note=', '.join(f'{lv}: {v:.3f}' for lv, v in lrel.items()))
    rep.note('The production integrator is validated against the Toolkit to 1e-7 at l = 5 (1e-5 at l = 2) in the Sec. II block above (stored Toolkit run).')

# =============================================================================================
# Section IV  Kerr
# =============================================================================================
KERR_CLOSED = {}

def section_IV():
    rep.section('Sec. IV  Kerr  [Eqs. (15)-(21), Table I]   (M = 1; r0 parametrizes the light ring, a = sqrt(r0)(3-r0)/2)')
    s = -2
    aS = sp.sqrt(r0)*(3 - r0)/2
    Delta = lambda rr, a: rr**2 - 2*rr + a**2
    psi, xW, zW = sp.symbols('psi x z', real=True)      # polar offset theta - pi/2, radial offset r - r0, Weber's z
    rep.ident('light-ring equation r0 (r0-3)^2 = 4 a^2 with a = sqrt(r0)(3-r0)/2', r0*(r0 - 3)**2 - 4*aS**2)
    worst = 0
    for rv in np.linspace(1.01, 3.99, 25):
        av = math.sqrt(rv)*(3 - rv)/2
        worst = max(worst, abs(2*(1 + math.cos(2/3*math.acos(-av)))/rv - 1))
    rep.check('r0(a) = 2{1 + cos[(2/3) arccos(-a)]} inverts a(r0) on 1 < r0 < 4 (both branches)', worst, '0', 1e-12, 'abs')
    bb = sp.symbols('b', positive=True)
    Rt = (r**2 + aS**2 - aS*bb)**2 - Delta(r, aS)*(bb - aS)**2
    bS = sp.sqrt(r0)*(r0 + 3)/2
    rep.ident('Eq. (15): R~(r0) = 0 for b = sqrt(r0)(r0+3)/2', Rt.subs(r, r0).subs(bb, bS))
    rep.ident("Eq. (15): R~'(r0) = 0 for b = sqrt(r0)(r0+3)/2", sp.diff(Rt, r).subs(r, r0).subs(bb, bS))
    roots = sp.solve(sp.diff(Rt, r).subs(r, r0), bb)
    rep.info('roots in b of R~\'(r0) = 0: b = sqrt(r0)(r0+3)/2 and a negative root -sqrt(r0)(r0^2+3)/(2(r0-1)) (not a light ring)', str([sp.factor(x) for x in roots]))
    RtS = Rt.subs(bb, bS)
    rep.ident('Eq. (15): R~ = (r-r0)^2 r (r+2 r0)', sp.expand(RtS - (r - r0)**2*r*(r + 2*r0)))
    D0 = sp.simplify(Delta(r0, aS)); P0 = r0**2 + aS**2 - aS*bS
    rep.ident('Delta_0 = r0 (r0-1)^2/4', D0 - r0*(r0 - 1)**2/4)
    rep.ident('P_0 = 2 r0 Delta_0/(r0-1)', P0 - 2*r0*D0/(r0 - 1))
    UpsS = (r0**2 + aS**2)*P0/D0 + aS*(bS - aS)
    rep.ident('Upsilon_t^ = (r0^2+a^2) P0/Delta_0 + a(b-a) = r0^2 (r0+3)/(r0-1)', UpsS - r0**2*(r0 + 3)/(r0 - 1))
    rep.ident('Upsilon_t^(a=0) = 27', UpsS.subs(r0, 3) - 27)
    # "The orbital frequency is Omega = 1/b": from the equatorial geodesic equations with E = 1, L = b,
    # Sigma dphi/dlambda = (b - a) + a P/Delta and Sigma dt/dlambda = a (b - a) + (r^2 + a^2) P/Delta = Upsilon_t^
    dphidt = ((bS - aS) + aS*P0/D0)/UpsS
    rep.ident('Omega = 1/b: dphi/dt = [(b-a) + a P0/Delta_0]/Upsilon_t^ = 1/b on the Kerr light ring (geodesic equations, every r0)', sp.simplify(dphidt - 1/bS))
    # "it is the Kerr version of the null condition L^2/r0^3 = E^2/(r0-2M)": at a = 0, R~ = r^4 - Delta b^2 with b = L/E
    Lsym, Esym_ = sp.symbols('L E', positive=True)
    Rt_a0 = r0**4 - (r0**2 - 2*r0)*(Lsym/Esym_)**2
    rep.ident('R~(r0) at a = 0 equals r0^4 (r0-2M)/E^2 [E^2/(r0-2M) - L^2/r0^3]: R~(r0) = 0 is the null condition of Eq. (6)',
              sp.simplify(Rt_a0 - r0**4*(r0 - 2)/Esym_**2*(Esym_**2/(r0 - 2) - Lsym**2/r0**3)))
    rep.ident('b^2 - a^2 = 3 r0^2  (beta_b = sqrt3 r0)', bS**2 - aS**2 - 3*r0**2)
    Rpp = sp.diff(RtS, r, 2).subs(r, r0)
    rep.ident("R~''(r0) = 6 r0^2", Rpp - 6*r0**2)
    ktS = D0**2*Rpp/(r0**2 + aS**2)**4
    rep.ident('k~ = Delta_0^2 R~\'\'(r0)/(r0^2+a^2)^4 = 3 r0^4 (r0-1)^4/[8 (r0^2+a^2)^4]', ktS - 3*r0**4*(r0 - 1)**4/(8*(r0**2 + aS**2)**4))
    DstK = lambda e: Delta(r, aS)/(r**2 + aS**2)*sp.diff(e, r)
    rep.ident('k~ = d^2/dr*^2 [R~/(r^2+a^2)^2] at r0 (curvature of Re Q/omega^2)', DstK(DstK(RtS/(r**2 + aS**2)**2)).subs(r, r0) - ktS)
    betab = sp.sqrt(3)*r0
    rep.ident('Omega_theta = lambda_L: beta_b/Ups = sqrt(R~\'\'/(2 Ups^2))', betab/UpsS - sp.sqrt(Rpp/(2*UpsS**2)))
    # Omega_theta and lambda_L from the geodesic equations (E = 1): the polar equation Sigma^2 (dtheta/dlambda)^2 =
    # Theta = Q + cos^2(theta) (a^2 - b^2/sin^2(theta)) with Q = 0, expanded about the equator (theta = pi/2 + psi),
    # is a harmonic oscillator of frequency beta_b/Sigma in lambda, i.e. beta_b/Upsilon_t in t (dt/dlambda = Ups/Sigma);
    # the radial equation Sigma^2 (dr/dlambda)^2 = R~ ~ (1/2) R~''(r0) x^2 gives the growth rate sqrt(R~''/2)/Sigma in
    # lambda, i.e. lambda_L = sqrt(R~''/(2 Ups^2)) in t.
    Theta = sp.cos(sp.pi/2 + psi)**2*(aS**2 - bS**2/sp.sin(sp.pi/2 + psi)**2)
    coef2 = sp.series(Theta, psi, 0, 4).removeO().coeff(psi, 2)
    rep.ident('polar equation about the equator: Theta = -(b^2 - a^2) psi^2 + O(psi^4) = -beta_b^2 psi^2',
              sp.simplify(sp.series(Theta, psi, 0, 4).removeO() + betab**2*psi**2))
    rep.ident('Omega_theta = sqrt(-coef)/Upsilon_t^ = beta_b/Upsilon_t^ is the vertical epicyclic frequency (dt/dlambda = Ups/Sigma; psi\'\' = coef psi/Sigma^2)',
              sp.simplify(sp.sqrt(-coef2)/UpsS - betab/UpsS))
    rep.ident('radial equation about r0: R~ = (R~\'\'(r0)/2) x^2 + O(x^3) with R~\'\'(r0)/2 = 3 r0^2',
              sp.simplify(sp.series(RtS.subs(r, r0 + xW), xW, 0, 3).removeO() - 3*r0**2*xW**2))
    rep.ident('Schwarzschild: lambda_L = sqrt(R~\'\'/(2 Ups^2)) = 1/(3 sqrt3 M) = Omega', sp.sqrt(Rpp/(2*UpsS**2)).subs(r0, 3) - 1/(3*sp.sqrt(3)))
    # eta_j from Re Q = omega^2 R~/(r^2+a^2)^2 - (2j+1) omega beta_b Delta/(r^2+a^2)^2:  eps = Re Q(r0), k = omega^2 k~
    w = sp.symbols('omega', positive=True)
    epsQ = -(2*j + 1)*w*betab*D0/(r0**2 + aS**2)**2
    rep.ident('Eq. (17): eta_j = -eps/sqrt(2k) = j + 1/2 for every r0', sp.simplify(-epsQ/sp.sqrt(2*w**2*ktS) - j - sp.Rational(1, 2)))
    # ---- Eq. (16): Teukolsky equation in Schroedinger form, and Im Q
    a, mm, lam = sp.symbols('a m lambda', real=True)
    y0, y1, y2 = sp.symbols('y0 y1 y2')          # Y(r0), Y'(r0), Y''(r0) as independent symbols
    Yf = sp.Function('Y')
    Kf = (r**2 + a**2)*w - a*mm
    Dl = r**2 - 2*r + a**2
    teuk = lambda R: Dl**(-s)*sp.diff(Dl**(s + 1)*sp.diff(R, r), r) + ((Kf**2 - 2*I*s*(r - 1)*Kf)/Dl + 4*I*s*w*r - lam)*R
    Gf = s*(r - 1)/(r**2 + a**2) + r*Dl/(r**2 + a**2)**2
    DstG = lambda e: Dl/(r**2 + a**2)*sp.diff(e, r)
    Qexpr = (Kf**2 - 2*I*s*(r - 1)*Kf + Dl*(4*I*s*w*r - lam))/(r**2 + a**2)**2 - Gf**2 - DstG(Gf)
    tosym = lambda e: e.subs(sp.Derivative(Yf(r), (r, 2)), y2).subs(sp.Derivative(Yf(r), r), y1).subs(Yf(r), y0)
    schro = tosym(DstG(DstG(Yf(r))) + Qexpr*Yf(r))
    ypp = sp.solve(schro, y2)[0]
    resid = tosym(teuk(Dl**sp.Rational(-s, 2)*(r**2 + a**2)**sp.Rational(-1, 2)*Yf(r))).subs(y2, ypp)
    # the residual is a rational function of r, a, omega, m, lambda (times y0, y1): test it exactly at random rational points
    rng = np.random.default_rng(7); worst = 0
    for _ in range(6):
        pt = {r: sp.Rational(int(rng.integers(30, 90)), 10), a: sp.Rational(int(rng.integers(-9, 10)), 10), w: sp.Rational(int(rng.integers(1, 50)), 7),
              mm: int(rng.integers(1, 40)), lam: sp.Rational(int(rng.integers(1, 900)), 3), y0: sp.Rational(int(rng.integers(1, 9)), 5), y1: sp.Rational(int(rng.integers(1, 9)), 3)}
        worst = max(worst, abs(complex(sp.N(resid.subs(pt), 50))))
    rep.check('Eq. (16): R = Delta^{-s/2}(r^2+a^2)^{-1/2} Y turns the Teukolsky equation into Y_r*r* + Q Y = 0 (residual at 6 random points, 50 digits)',
              worst, '0', 1e-40, 'abs')
    P_ = r**2 + a**2 - a*bb
    Qres = Qexpr.subs(mm, w*bb)
    ImQ_target = s*w*(-2*(r - 1)*P_ + 4*r*Dl)/(r**2 + a**2)**2
    # Im Q: (Q - conj(Q))/(2i) with all parameters real; evaluated exactly at random points
    worst = 0
    for _ in range(6):
        pt = {r: sp.Rational(int(rng.integers(30, 90)), 10), a: sp.Rational(int(rng.integers(-9, 10)), 10), w: sp.Rational(int(rng.integers(1, 50)), 7),
              bb: sp.Rational(int(rng.integers(20, 70)), 10), lam: sp.Rational(int(rng.integers(1, 900)), 3)}
        qv = complex(sp.N(Qres.subs(pt), 50)); tv = float(sp.N(ImQ_target.subs(pt), 50))
        worst = max(worst, abs(qv.imag - tv))
    rep.check('Im Q = s omega [-2(r-M)P + 4 r Delta]/(r^2+a^2)^2 for m = omega b (6 random points, 50 digits)', worst, '0', 1e-40, 'abs')
    rep.ident('Im Q = 0 at the light ring: -2(r0-M)P0 + 4 r0 Delta_0 = 0 for every spin', -2*(r0 - 1)*P0 + 4*r0*D0)
    dImQ = sp.diff((-2*(r - 1)*(r**2 + aS**2 - aS*bS) + 4*r*Delta(r, aS)), r).subs(r, r0)
    rep.ident('Im Q odd in x at leading order: d/dr[-2(r-M)P + 4 r Delta] at r0 = 3 r0 (r0-M) (nonzero for r0 > M)', sp.simplify(dImQ - 3*r0*(r0 - 1)))
    # Im Q = O(omega x/M^2) against the curvature term O(omega^2 x^2/M^4): in Weber's variable z = (2k)^{1/4} x,
    # k = omega^2 k~, the equation divided by sqrt(2k) has Im Q/sqrt(2k) = O(omega^{-1/2}) z
    f0K = D0/(r0**2 + aS**2)                                   # dr/dr* at r0: x_r = f0 x
    ImQ_lin = s*w*3*r0*(r0 - 1)*f0K/(r0**2 + aS**2)**2         # Im Q = ImQ_lin x + O(x^2)
    weber_term = (ImQ_lin*xW/sp.sqrt(2*w**2*ktS)).subs(xW, zW/(2*w**2*ktS)**sp.Rational(1, 4))
    rep.ident('Im Q in Weber\'s equation is an O(l^{-1/2}) z perturbation: omega^{1/2} x [Im Q/sqrt(2k)] is omega-independent',
              sp.simplify(sp.diff(sp.sqrt(w)*weber_term, w)))
    rep.ident('  ... and it is odd in z (coefficient of z only; the curvature term is z^2/4)',
              sp.simplify(weber_term + weber_term.subs(zW, -zW)))
    # Re Q with the eikonal eigenvalue
    j_, bbs, lam0 = sp.symbols('j beta_b lambda_0', real=True)
    ReQlead = sp.simplify((((r**2 + a**2)*w - a*mm)**2 - Dl*lam)/(r**2 + a**2)**2).subs({mm: w*bb, lam: (w*bb - a*w)**2 + (2*j_ + 1)*w*bbs + lam0})
    rep.ident('Re Q = omega^2 R~/(r^2+a^2)^2 - (2j+1) omega beta_b Delta/(r^2+a^2)^2 + O(1) with Lambda = (m-a omega)^2 + (2j+1) omega beta_b + O(1)',
              sp.simplify(ReQlead - (w**2*((r**2 + a**2 - a*bb)**2 - Dl*(bb - a)**2) - (2*j_ + 1)*w*bbs*Dl)/(r**2 + a**2)**2 + lam0*Dl/(r**2 + a**2)**2))
    # Lambda expansion against the eigenvalues stored with the Kerr runs
    worst_abs = 0; worst_rel = 0; worst_conv = 0
    for av, sg in ((0.5, 1), (0.5, -1), (0.9, 1), (0.9, -1)):
        aa = sg*av; rv = 2*(1 + math.cos(2/3*math.acos(-aa))); bv = math.sqrt(rv)*(rv + 3)/2; bbv = math.sqrt(3)*rv
        rem = {}
        for rw in KERR:
            if rw['a'] == av and rw['sign'] == sg and rw['l'] >= 100:
                mv = rw['l'] - rw['j']; wv = mv/bv
                diff = rw['lam'] - ((mv - aa*wv)**2 + (2*rw['j'] + 1)*wv*bbv)
                rem[(rw['j'], rw['l'])] = diff
                worst_abs = max(worst_abs, abs(diff)); worst_rel = max(worst_rel, abs(diff)/rw['lam'])
        for jv in range(4): worst_conv = max(worst_conv, abs(rem[(jv, 800)] - rem[(jv, 400)]))
    rep.check('Lambda - (m-a omega)^2 - (2j+1) omega beta_b = O(1): bounded (max |.| < 20) over stored a=+-0.5,+-0.9, l=100..800, j<=3 eigenvalues',
              worst_abs, '0', 20, 'abs', note=f'max |remainder| = {worst_abs:.2f} (j = 3), relative to Lambda: {worst_rel:.1e}')
    rep.check('  ... and l-independent: max |remainder(800) - remainder(400)|', worst_conv, '0', 0.1, 'abs')
    # The expansion of the angular equation about the equator (the origin of Lambda and of |S(pi/2)|^2). In the
    # Toolkit convention the angular equation is (1/sin) d/dtheta (sin dS/dtheta) + [g^2 cos^2 - m^2/sin^2 - 2 g s cos
    # - 2 m s cos/sin^2 - s^2 cot^2 + s + A] S = 0, g = a omega, A = Lambda - g^2 + 2 m g. With theta = pi/2 + psi the
    # bracket is (A + s - m^2) + 2 s (g + m) psi - (m^2 - g^2 + s^2) psi^2 + O(psi^3): a harmonic oscillator of
    # frequency sqrt(m^2 - g^2) = omega beta_b at large omega, whose j-th eigenvalue (2j+1) omega beta_b gives
    # A = m^2 + (2j+1) omega beta_b + O(1), i.e. Lambda = (m - a omega)^2 + (2j+1) omega beta_b + O(1).
    gS, As_ = sp.symbols('g A', real=True)
    th = sp.pi/2 + psi
    bracket = gS**2*sp.cos(th)**2 - mm**2/sp.sin(th)**2 - 2*gS*s*sp.cos(th) - 2*mm*s*sp.cos(th)/sp.sin(th)**2 - s**2*sp.cot(th)**2 + s + As_
    ser_b = sp.series(bracket, psi, 0, 3).removeO()
    rep.ident('angular equation about the equator: bracket = (A + s - m^2) + 2 s (g + m) psi - (m^2 - g^2 + s^2) psi^2 + O(psi^3)',
              sp.simplify(ser_b - ((As_ + s - mm**2) + 2*s*(gS + mm)*psi - (mm**2 - gS**2 + s**2)*psi**2)))
    rep.ident('oscillator frequency^2 = m^2 - a^2 omega^2 = omega^2 beta_b^2 for m = omega b', sp.simplify((w*bb)**2 - (a*w)**2 - w**2*(bb**2 - a**2)))
    rep.ident('j-th eigenvalue (2j+1) omega beta_b: A = m^2 - s + (2j+1) omega beta_b + O(1) gives Lambda = (m - a omega)^2 + (2j+1) omega beta_b + O(1)',
              sp.simplify(((w*bb)**2 - s + (2*j_ + 1)*w*bbs) + (a*w)**2 - 2*(w*bb)*(a*w) - ((w*bb - a*w)**2 + (2*j_ + 1)*w*bbs) + s))
    # the harmonic coefficients c_j, d_j of Sec. III C from the normalised Hermite functions h_j(x) of the oscillator:
    # |S(pi/2)|^2 = h_j(0)^2 sqrt(omega beta_b)/(2 pi) = c_j sqrt(omega beta_b)/(2 pi^{3/2}), with h_j(0)^2 sqrt(pi) = c_j,
    # and |dS/dtheta(pi/2)|^2 = h_j'(0)^2 (omega beta_b)^{3/2}/(2 pi) = d_j (omega beta_b)^{3/2}/pi^{3/2}, h_j'(0)^2 sqrt(pi)/2 = d_j
    xh = sp.symbols('x_h', real=True)
    worst_h = 0
    for jv in range(0, 13):
        hj = sp.hermite(jv, xh)*sp.exp(-xh**2/2)/sp.sqrt(2**jv*sp.factorial(jv)*sp.sqrt(sp.pi))
        if jv % 2 == 0: worst_h = max(worst_h, abs(sp.simplify(hj.subs(xh, 0)**2*sp.sqrt(sp.pi) - c_j(jv))))
        else: worst_h = max(worst_h, abs(sp.simplify(sp.diff(hj, xh).subs(xh, 0)**2*sp.sqrt(sp.pi)/2 - d_j(jv))))
    rep.check('|S(pi/2)|^2 -> c_j sqrt(omega beta_b)/(2 pi^{3/2}) and |S\'|^2 -> d_j (omega beta_b)^{3/2}/pi^{3/2}: c_j = sqrt(pi) h_j(0)^2, d_j = sqrt(pi) h_j\'(0)^2/2 (j <= 12)',
              float(worst_h), '0', 0, 'abs', note='normalised Hermite functions h_j; at a = 0, omega beta_b = l (the spherical limits of Sec. III C)')
    # ---- amplitude integrals, Eq. (18): numerics with mpmath at several r0, and the closed forms
    mp.mp.dps = 30
    def lr_num(rv):
        rv = mp.mpf(rv); av = mp.sqrt(rv)*(3 - rv)/2; bv = mp.sqrt(rv)*(rv + 3)/2
        return dict(r0=rv, a=av, b=bv, D0=rv*(rv - 1)**2/4, P0=rv**2 + av**2 - av*bv, Ups=rv**2*(rv + 3)/(rv - 1), betab=mp.sqrt(3)*rv,
                    kt=3*rv**4*(rv - 1)**4/(8*(rv**2 + av**2)**4), rp=1 + mp.sqrt(1 - av**2))
    def Jinf(d):
        av, bv, rv = d['a'], d['b'], d['r0']
        Dl_ = lambda rr: rr**2 - 2*rr + av**2
        P_ = lambda rr: rr**2 + av**2 - av*bv
        q_drs = lambda rr: s*(-2*(rr - 1)*P_(rr) + 4*rr*Dl_(rr))/(2*Dl_(rr)*(rr - rv)*mp.sqrt(rr*(rr + 2*rv)))
        integrand = lambda rr: q_drs(rr) - s/rr*(rr**2 + av**2)/Dl_(rr)
        return mp.quad(integrand, [rv, rv + 1, rv + 10, mp.inf])
    def hinf(d):
        av, rv = d['a'], d['r0']
        return mp.quad(lambda rr: (rr**2 + av**2)/(rr*(rr**2 - 2*rr + av**2)) - 1/rr, [rv, rv + 10, mp.inf])
    def JH(d):
        av, bv, rv = d['a'], d['b'], d['r0']
        Dl_ = lambda rr: rr**2 - 2*rr + av**2
        P_ = lambda rr: rr**2 + av**2 - av*bv
        integrand = lambda rr: s/(2*Dl_(rr))*((-2*(rr - 1)*P_(rr) + 4*rr*Dl_(rr))/((rv - rr)*mp.sqrt(rr*(rr + 2*rv))) + 2*(rr - 1))
        return mp.quad(integrand, [d['rp'], (d['rp'] + rv)/2, rv])
    worstN = 0; worsth = 0; worstJH = 0; worstHI = 0; worstOm = 0
    # the seven spins of Table I (r0 from the exact inversion), plus one more prograde point
    r0_table = [2*(1 + mp.cos(mp.mpf(2)/3*mp.acos(-mp.mpf(av)))) for av in ('0', '0.5', '0.9', '0.99', '-0.5', '-0.9', '-0.99')]
    for rv in r0_table + [mp.mpf('1.3')]:
        d = lr_num(rv)
        J0 = Jinf(d); h0 = hinf(d)
        N0 = d['r0']**(2*s)*mp.exp(-2*s*h0)*mp.exp(-2*J0)
        worstN = max(worstN, abs(N0/(16/d['r0']**6) - 1))
        rp, rm = d['rp'], 1 - mp.sqrt(1 - d['a']**2)
        hcl = 2/(rp - rm)*mp.log((d['r0'] - rm)/(d['r0'] - rp))   # h_inf = int 2/Delta dr (r- = 0 at a = 0)
        worsth = max(worsth, abs(h0/hcl - 1))
        JHv = JH(d)
        e2JH = (d['r0']**2 + d['r0'] + 4*mp.sqrt(4 - d['r0']) - 8)**4/d['r0']**8
        worstJH = max(worstJH, abs(mp.exp(2*JHv)/e2JH - 1))
        # horizon/infinity ratio built from the Wronskian at the horizon, J_H and the Teukolsky-Press factor (large-omega limit)
        OmH = d['a']/(2*rp); kap = 1 - d['b']*OmH
        NH = (2*rp)*d['D0']**s*mp.exp(-2*JHv)/((2*rp)**2*kap**2)
        alphaH = 256*(2*rp)**5*kap**5/(d['b'] - d['a'])**8
        worstHI = max(worstHI, abs(alphaH*kap*NH/N0 - 1))
        worstOm = min(worstOm, 1/d['b'] - OmH) if worstOm else 1/d['b'] - OmH
    rep.check('Eq. (18): N = r0^{2s} e^{-2s h_inf} e^{-2 J_inf} = 16/r0^6 (mpmath quadrature, the 7 spins of Table I and r0 = 1.3)', float(worstN), '0', 1e-20, 'abs')
    rep.check('h_inf = int 2/Delta dr (elementary logarithm) vs quadrature', float(worsth), '0', 1e-20, 'abs')
    rep.check('e^{2 J_H} = [r0^2 + r0 + 4 sqrt(4-r0) - 8]^4/r0^8 vs quadrature (same 8 values of r0)', float(worstJH), '0', 1e-18, 'abs')
    rep.check('horizon/infinity flux ratio = 1 (alpha_H kappa N_H/N, numerically, same 8 values of r0)', float(worstHI), '0', 1e-18, 'abs')
    rep.check('footnote: the ratio is 1 "to better than 1e-17 at all the spins of Table I" (30-digit quadrature at a = 0, +-0.5, +-0.9, +-0.99)', float(worstHI), '1e-17', 0, 'max')
    rep.ident('horizon/infinity ratio = 1 symbolically: 4(2 r+ - a b) = (r0-1)(r0^2 + r0 + 4 sqrt(4-r0) - 8) (reduces the ratio to 1 given e^{2J_H})',
              sp.expand((4*(2*(1 + (r0 - 1)*sp.sqrt(4 - r0)/2) - aS*bS) - (r0 - 1)*(r0**2 + r0 + 4*sp.sqrt(4 - r0) - 8))))
    rep.ident('sqrt(1 - a^2) = (r0-1) sqrt(4-r0)/2 (used above; r+ = 1 + sqrt(1-a^2))', sp.simplify(1 - aS**2 - (r0 - 1)**2*(4 - r0)/4))
    # no superradiance: Omega > Omega_H
    worstOm = min(float(1/sp.sqrt(rv)/(rv + 3)*2 - sp.sqrt(rv)*(3 - rv)/2/(2*(1 + (rv - 1)*sp.sqrt(4 - rv)/2))) for rv in [sp.Rational(k, 100) for k in range(101, 400)])
    rep.check('Omega - Omega_H > 0 on 1 < r0 < 4 (min over a grid)', float(worstOm > 0), '1', 0, note=f'min = {worstOm:.3e}; -> 0 only as r0 -> 1 (a -> 1)')
    rep.ident('Omega - Omega_H = (r0-1)[...]: vanishes at r0 = 1', (1/bS - aS/(2*(1 + (r0 - 1)*sp.sqrt(4 - r0)/2))).subs(r0, 1))
    # symbolic proof: with s = sqrt(4 - r0) (0 < s < sqrt3 on 1 < r0 < 4), Omega - Omega_H = (2 r+ - a b)/(2 b r+) and 2 r+ - a b is a
    # polynomial in s whose real roots are all outside the interval except s = sqrt3 (r0 = M, a = M), where it vanishes
    sS = sp.symbols('s', positive=True)
    r0s = 4 - sS**2; aSs = sp.sqrt(r0s)*(3 - r0s)/2; bSs = sp.sqrt(r0s)*(r0s + 3)/2; rps = 1 + (r0s - 1)*sS/2
    polyOm = sp.expand(2*rps - aSs*bSs)
    rep.ident('Omega > Omega_H symbolically: 2 r+ - a b = (2-s)^2 (s+1)(s+3)(3-s^2)/4 with s = sqrt(4 - r0) (r+ = 1 + (r0-1) s/2)',
              sp.expand(polyOm - (2 - sS)**2*(sS + 1)*(sS + 3)*(3 - sS**2)/4))
    rootsOm = sp.real_roots(sp.Poly(polyOm, sS))
    rep.check('  ... no real root of 2 r+ - a b in 0 < s < sqrt3 (1 < r0 < 4) and 2 r+ - a b > 0 at s = 1: Omega > Omega_H for every light ring, -> 0 only as a -> M',
              float(bool(all(not (0 < rt < sp.sqrt(3)) for rt in rootsOm)) and bool(polyOm.subs(sS, 1) > 0)), '1', 0, note=f'real roots in s: {[sp.N(x, 5) for x in rootsOm]}')
    # ---- source expansion (alphaExpand transcribed from the notebook), closed forms for Ahat_0, Bhat
    ae0 = alpha_expand(aS, bS, betab, D0, P0, 0)
    ae1 = alpha_expand(aS, bS, betab, D0, P0, 1)
    rep.ident('source: O(omega^2) coefficient vanishes for every r0 (even j = 0)', ae0['w2'])
    rep.ident('source: O(omega^{3/2}) coefficient vanishes for every r0 (even j = 0)', ae0['w32'])
    rep.ident('source: O(omega^2) and O(omega^{3/2}) coefficients vanish (odd j = 1)', ae1['w2'] + ae1['w32'])
    AhS = sp.simplify(ae0['Ahat']); BhS = sp.simplify(ae1['Bhat'])
    rep.ident('Ahat_0 = sqrt3 r0^4/(4 M sqrt(r0^2+a^2))', AhS - sp.sqrt(3)*r0**4/(4*sp.sqrt(r0**2 + aS**2)))
    rep.ident('Bhat = -i r0^{5/2} sqrt(r0^2+a^2)/(2 sqrt(M) (r0-M))', BhS + I*r0**sp.Rational(5, 2)*sp.sqrt(r0**2 + aS**2)/(2*(r0 - 1)))
    rep.ident('Schwarzschild: Ahat_0 = (9/4) sqrt27 M^2 = (3/4) r0 b', AhS.subs(r0, 3) - sp.Rational(9, 4)*sp.sqrt(27))
    rep.ident('Schwarzschild: Bhat = -i Ahat_0', (BhS + I*AhS).subs(r0, 3))
    ae2 = alpha_expand(aS, bS, betab, D0, P0, 2); ae3 = alpha_expand(aS, bS, betab, D0, P0, 3)
    rep.ident('Ahat_2/Ahat_0 = 5 = 2j+1', sp.simplify(ae2['Ahat']/ae0['Ahat']) - 5)
    rep.ident('Bhat independent of j (j = 3 vs j = 1)', sp.simplify(ae3['Bhat']/ae1['Bhat']) - 1)
    rep.ident('source cancellations also at j = 2, 3', ae2['w2'] + ae2['w32'] + ae3['w2'] + ae3['w32'])
    worstA = 0; worstB = 0; worstC = 0
    for jv in range(4, 13):
        ae = alpha_expand(aS, bS, betab, D0, P0, jv)
        worstC = max(worstC, abs(complex(sp.N((ae['w2'] + ae['w32']).subs(r0, sp.Rational(23, 10)), 30))))
        if jv % 2 == 0: worstA = max(worstA, abs(complex(sp.N((ae['Ahat']/ae0['Ahat'] - (2*jv + 1)).subs(r0, sp.Rational(23, 10)), 30))),
                                     abs(complex(sp.N((ae['Ahat']/ae0['Ahat'] - (2*jv + 1)).subs(r0, sp.Rational(37, 10)), 30))))
        else: worstB = max(worstB, abs(complex(sp.N((ae['Bhat']/ae1['Bhat'] - 1).subs(r0, sp.Rational(23, 10)), 30))),
                           abs(complex(sp.N((ae['Bhat']/ae1['Bhat'] - 1).subs(r0, sp.Rational(37, 10)), 30))))
    rep.check('Ahat_j = (2j+1) Ahat_0 for j = 4, 6, ..., 12 (exact expressions evaluated at r0 = 2.3 and 3.7 to 30 digits)', worstA, '0', 1e-25, 'abs')
    rep.check('Bhat independent of j for j = 5, 7, ..., 11 (same evaluation)', worstB, '0', 1e-25, 'abs')
    rep.check('source cancellations (omega^2, omega^{3/2}) for j = 4 ... 12', worstC, '0', 1e-25, 'abs')
    # "The null geodesic enters alpha_lm only through E, L = bE and Upsilon_t, so that alpha_lm propto E^2": with p^mu -> E p^mu
    # (L = bE) the source projection is homogeneous of degree 2 in E (Upsilon_t -> E Upsilon_t divides Z, so Z propto E, fluxes propto E^2)
    aG, bG, bwG, EG_ = sp.symbols('a b b_w E', positive=True)
    g0 = alpha_expand_general(aG, bG, r0, 0); gE = alpha_expand_general(aG, bG, r0, 0, Esym=EG_)
    rep.ident('alpha_lm propto E^2: alpha(E p^mu)/alpha(p^mu) = E^2 (symbolic in a, b, r0; even j)', sp.simplify(gE['alpha']/g0['alpha'] - EG_**2))
    g1 = alpha_expand_general(aG, bG, r0, 1); g1E = alpha_expand_general(aG, bG, r0, 1, Esym=EG_)
    rep.ident('alpha_lm propto E^2 (odd j)', sp.simplify(g1E['alpha']/g1['alpha'] - EG_**2))
    # "both cancellations are consequences of the null condition R~(r0) = 0": the same expansion with the orbit data free
    # (a, b = L/E and r0 independent, so that R~(r0) = P0^2 - Delta_0 (b-a)^2 is not zero), keeping omega = m/b
    RtG = (r0**2 + aG**2 - aG*bG)**2 - (r0**2 - 2*r0 + aG**2)*(bG - aG)**2
    DG = r0**2 - 2*r0 + aG**2
    rep.ident('orbit data free, omega = m/b: the O(omega^2) coefficient of alpha is -S Y(0) (b-a)^2 R~(r0)/[8 sqrt(r0^2+a^2) Delta_0], i.e. propto the null condition (even j)',
              sp.simplify(g0['w2'] + g0['S']*g0['Y0']*(bG - aG)**2*RtG/(8*sp.sqrt(r0**2 + aG**2)*DG)))
    rep.ident('orbit data free, omega = m/b: the O(omega^{3/2}) coefficient vanishes identically (even j)', g0['w32'])
    rep.ident('orbit data free: both coefficients vanish identically for odd j (S(pi/2) = O(1/omega) S\'(pi/2))', g1['w2'] + g1['w32'])
    # ... and the frequency condition omega = m/b, i.e. Omega = 1/b, is itself equivalent to the null condition: with the wave
    # frequency decoupled from the orbit (m = omega b_w) the O(omega^{3/2}) term is propto (b - b_w), and for a circular orbit
    # with E = 1, L = b at r0 the geodesic equations give dphi/dt - 1/b = -R~(r0)/[r0 b (r0^3 + a^2 r0 + 2a^2 - 2ab)]
    gw = alpha_expand_general(aG, bG, r0, 0, bwave=bwG)
    rep.ident('wave frequency decoupled (m = omega b_w): O(omega^{3/2}) coefficient = i S y1 r0^2 (a-b) sqrt(r0^2+a^2) (b - b_w)/(4 Delta_0): vanishes iff Omega = 1/b',
              sp.simplify(gw['w32'] - I*gw['S']*gw['y1']*r0**2*(aG - bG)*sp.sqrt(r0**2 + aG**2)*(bG - bwG)/(4*DG)))
    PG = r0**2 + aG**2 - aG*bG
    dphidtG = ((bG - aG) + aG*PG/DG)/(aG*(bG - aG) + (r0**2 + aG**2)*PG/DG)
    rep.ident('circular orbit with E = 1, L = b at r0: dphi/dt - 1/b = -R~(r0)/[r0 b (r0^3 + a^2 r0 + 2 a^2 - 2 a b)]: Omega = 1/b iff R~(r0) = 0',
              sp.simplify(dphidtG - 1/bG + RtG/(r0*bG*(r0**3 + aG**2*r0 + 2*aG**2 - 2*aG*bG))))
    rep.ident('identity |Bhat/Ahat_0|^2 sqrt(2 k~) beta_b = 2 for every spin',
              sp.simplify(sp.expand_complex(sp.Abs(BhS)**2/sp.Abs(AhS)**2)*sp.sqrt(2*ktS)*betab - 2))
    # ---- Eq. (19) and the odd counterpart; Table I
    NS = 16/r0**6
    def kappa_kerr(d, jv):
        """Eq. (19) (even j) and the odd analogue from |u'(0)|^2 and Bhat, with the closed forms (mpmath)."""
        x = mp.mpf(jv) + mp.mpf(1)/2
        Ah = mp.sqrt(3)*d['r0']**4/(4*mp.sqrt(d['r0']**2 + d['a']**2))*(2*jv + 1)
        Bh = d['r0']**mp.mpf(2.5)*mp.sqrt(d['r0']**2 + d['a']**2)/(2*(d['r0'] - 1))
        N0 = 16/d['r0']**6
        if jv % 2 == 0:
            return mp.sqrt(2*mp.pi)*d['b']*mp.mpf(c_j(jv))*mp.sqrt(d['betab'])*N0*Fpar(x)*Ah**2/((2*d['kt'])**mp.mpf(0.25)*d['Ups']**2)
        Du = (2*d['kt'])**mp.mpf(0.25)*N0*mp.sqrt(2)*mp.pi*Gpar(x)
        return 4*mp.pi*Bh**2*mp.mpf(d_j(jv))*d['betab']**mp.mpf(1.5)/mp.pi**mp.mpf(1.5)*Du*d['b']/d['Ups']**2
    d0 = lr_num(3)
    for jv in range(5):
        rep.check(f'Eq. (19)/odd analogue at a = 0 reproduces kappa_{jv} of Eqs. (11)-(12)', kappa_kerr(d0, jv)/kappa_j(jv), '1', 1e-20)
    k00 = kappa_kerr(d0, 0)
    table = [(0, 3, 5.196, 27, 7.407e-2, 11.69, 2.195e-2, 0.02777, 1), (0.5, 2.347, 4.096, 21.87, 3.408e-2, 5.477, 9.566e-2, 0.03428, 1.235),
             (0.9, 1.558, 2.844, 19.83, 1.948e-3, 1.418, 1.119, 0.03781, 1.362), (0.99, 1.168, 2.252, 33.89, 1.825e-5, 0.5258, 6.313, 0.02212, 0.797),
             (-0.5, 3.532, 6.138, 32.18, 9.149e-2, 18.89, 8.240e-3, 0.02329, 0.839), (-0.9, 3.910, 6.832, 36.31, 9.360e-2, 25.23, 4.476e-3, 0.02065, 0.744),
             (-0.99, 3.991, 6.983, 37.23, 9.317e-2, 26.72, 3.959e-3, 0.02014, 0.725)]
    cols = ['r0', 'b', 'Ups_t', 'k~', '|Ahat_0|', 'N', 'kappa_0(a)', 'g(a)']
    worst_tab = 0
    for row in table:
        av = row[0]; rv = 2*(1 + mp.cos(mp.mpf(2)/3*mp.acos(-mp.mpf(av)))); d = lr_num(rv)
        vals = [d['r0'], d['b'], d['Ups'], d['kt'], mp.sqrt(3)*d['r0']**4/(4*mp.sqrt(d['r0']**2 + d['a']**2)), 16/d['r0']**6,
                kappa_kerr(d, 0), kappa_kerr(d, 0)/k00]
        devs = [abs(float(v)/p - 1) for v, p in zip(vals, row[1:])]
        worst_tab = max(worst_tab, max(devs))
        rep.check(f'Table I, a = {av}: max rel. dev. over {", ".join(cols)}', max(devs), '0', 6e-4, 'abs',
                  note=' '.join(f'{c}={float(v):.5g}' for c, v in zip(cols, vals)))
        KERR_CLOSED[av] = d
    gS_closed = 27*(r0 - 1)/(r0**2*(r0 + 3))
    kap0S = bS*sp.sqrt(betab)*NS*AhS**2/((2*ktS)**sp.Rational(1, 4)*UpsS**2)     # kappa_0(a) up to spin-independent factors
    rep.ident('Eq. (21): g(a) = kappa_0(a)/kappa_0(0) = 27 M^2 (r0-M)/[r0^2 (r0+3M)]', sp.simplify(kap0S/kap0S.subs(r0, 3) - gS_closed))
    rep.ident('g(0) = 1', gS_closed.subs(r0, 3) - 1)
    rep.ident("g'(r0) = 0 at r0 = sqrt3 M", sp.diff(gS_closed, r0).subs(r0, sp.sqrt(3)))
    rep.ident('max g = 6 sqrt3 - 9', sp.simplify(gS_closed.subs(r0, sp.sqrt(3)) - (6*sp.sqrt(3) - 9)))
    rep.check('max g = 6 sqrt3 - 9 = 1.392', float(6*math.sqrt(3) - 9), '1.392', 5e-4)
    rep.ident('a/M at the maximum = 3^{1/4}(3 - sqrt3)/2', sp.simplify(aS.subs(r0, sp.sqrt(3)) - 3**sp.Rational(1, 4)*(3 - sp.sqrt(3))/2))
    rep.check('a/M at the maximum = 0.834', float(aS.subs(r0, sp.sqrt(3))), '0.834', 1e-3)
    rep.check('g -> 0 as r0 -> M (extremal prograde)', float(gS_closed.subs(r0, 1)), '0', 0, 'abs')
    rep.ident('g(r0 = 4M) = 81/112 (extremal retrograde)', gS_closed.subs(r0, 4) - sp.Rational(81, 112))
    rep.check('81/112 = 0.723', 81/112, '0.723', 1e-3)
    gnum = lambda av: float(gS_closed.subs(r0, 2*(1 + sp.cos(sp.Rational(2, 3)*sp.acos(-sp.nsimplify(av))))))
    rep.check('g(a = 0.98) = 0.99', gnum(0.98), '0.99', 6e-3)
    rep.check('g(a = 0.995) = 0.62', gnum(0.995), '0.62', 1e-2)
    rep.check('g(a = 0.83) close to the maximum (paper: maximum near a/M = 0.83)', gnum(0.83)/(6*math.sqrt(3) - 9), '1', 1e-3)
    # k~ -> 0 as r0 -> 1 (barrier flattens)
    rep.check('k~ -> 0 as r0 -> M (barrier flattens): k~(r0=1.001)', float(ktS.subs(r0, sp.Rational(1001, 1000))), '0', 1e-11, 'abs')
    # Table I caption: at a = 0.99 (prograde) "the barrier is so flat that the asymptotic regime sets in only at l >> 800". The two
    # odd O(l^-1/2) perturbations of Weber's equation are the cubic term g z^3, g = V3/[6 (2k)^{5/4}] with V3 = -omega^2 d^3/dr*^3
    # [R~/(r^2+a^2)^2] at r0 (Re Q = omega^2 R~/(r^2+a^2)^2 at leading order), and the Im Q term i s c1 omega^{-1/2} z/(2k~)^{3/4}
    # with Im Q = s omega c1 x (c1 = 3 r0 (r0-M) f0/(r0^2+a^2)^2, see above); both scale as omega^{-1/2} = (b/l)^{1/2} times a
    # spin-dependent coefficient. The l at which the a = 0.99 coefficients equal their a = 0 (or a = 0.9) values at l = 800 is
    # 800 x (coefficient ratio)^2 (the corrections to the mean flux are the squares, O(1/l)).
    R3 = DstK(DstK(DstK(RtS/(r**2 + aS**2)**2))).subs(r, r0)
    gcub = -R3/(6*(2*ktS)**sp.Rational(5, 4))*sp.sqrt(bS)            # cubic coefficient x l^{1/2} (omega = l/b at j = 0)
    cimq = s*3*r0*(r0 - 1)*f0K/(r0**2 + aS**2)**2/(2*ktS)**sp.Rational(3, 4)*sp.sqrt(bS)   # Im Q coefficient x l^{1/2}
    rep.ident('Kerr cubic coefficient at a = 0 reduces to the Schwarzschild g0 l^{-1/2} of Sec. III B (|g0| = 0.0680)',
              sp.simplify(sp.Abs(gcub.subs(r0, 3)) - (4/sp.Integer(6561))/(6*(4/sp.Integer(729))**sp.Rational(5, 4))))
    coefs = {}
    for av in (0, 0.5, 0.9, 0.99, -0.5, -0.9, -0.99):
        rv = sp.Float(2*(1 + mp.cos(mp.mpf(2)/3*mp.acos(-mp.mpf(av)))), 30)
        coefs[av] = (abs(float(gcub.subs(r0, rv))), abs(float(cimq.subs(r0, rv))))
    rep.info('O(l^-1/2) coefficients (x l^{1/2}) of the cubic and Im Q terms at a = 0, 0.5, 0.9, 0.99, -0.5, -0.9, -0.99',
             '; '.join(f'{av}: {c[0]:.3f}, {c[1]:.3f}' for av, c in coefs.items()))
    lstar = {nm: 800*(coefs[0.99][k]/coefs[0][k])**2 for k, nm in ((0, 'cubic'), (1, 'ImQ'))}
    lstar9 = {nm: 800*(coefs[0.99][k]/coefs[0.9][k])**2 for k, nm in ((0, 'cubic'), (1, 'ImQ'))}
    rep.check('a = 0.99: l at which the cubic-term coefficient equals its a = 0 value at l = 800 is >> 800 (bound: > 8000)',
              float(lstar['cubic'] > 8000), '1', 0, note=f"l* = {lstar['cubic']:.2e} (vs a = 0.9 at l = 800: {lstar9['cubic']:.2e})")
    rep.check('a = 0.99: the same for the Im Q term (bound: > 1600)', float(lstar['ImQ'] > 1600), '1', 0,
              note=f"l* = {lstar['ImQ']:.2e} (vs a = 0.9 at l = 800: {lstar9['ImQ']:.2e})")
    # ---- stored Kerr Teukolsky runs: Richardson extrapolation and Table I last column
    rep.note('Stored Kerr runs (kerr_py/kerr_results.json, l <= 800) are NOT recomputed; the quoted comparisons are re-derived from them.')
    tabnum = {0: 1.0012, 0.5: 1.0014, 0.9: 1.0017, -0.5: 1.0011, -0.9: 1.0010}
    worst_rich = 0
    for av in (0, 0.5, -0.5, 0.9, -0.9):
        sg = 1 if av >= 0 else -1
        d = KERR_CLOSED[av]
        kA = 2*mp.fsum(kappa_kerr(d, jv) for jv in range(4))
        fl = {lv: lv*sum(rw['FluxI'] + rw['FluxH'] for rw in kerr_rows(av, sg, lv)) for lv in (400, 800)}
        rich = 2*fl[800] - fl[400]
        worst_rich = max(worst_rich, abs(rich/float(kA) - 1))
        rep.check(f'Table I last column, a = {av}: num./an. at l = 800', fl[800]/float(kA), str(tabnum[av]), 1e-4,
                  note=f'Richardson (2 f800 - f400)/analytic - 1 = {rich/float(kA) - 1:+.2e}')
    rep.check('Richardson extrapolation vs 2 g(a) sum_{j<=3} kappa_j: max |ratio - 1| (paper: better than 1.5e-5 in relative terms)', worst_rich, '1.5e-5', 0, 'max')
    # the O(l^-1/2) asymmetry in Kerr, from the stored runs at l = 800 (m-summed, j <= 3)
    asym_paper = {0.5: '0.78', 0.9: '1.10', -0.5: '0.51', -0.9: '0.45', 0: '0.61'}
    for av in (0.5, 0.9, -0.5, -0.9, 0):
        sg = 1 if av >= 0 else -1
        SI = sum(rw['FluxI'] for rw in kerr_rows(av, sg, 800)); SH = sum(rw['FluxH'] for rw in kerr_rows(av, sg, 800))
        rep.check(f'sqrt(l)(I-H)/(I+H) at l = 800, a = {av} (paper: {asym_paper[av]}{"; the asymptotic sigma_bar of Sec. III C" if av == 0 else ""})',
                  math.sqrt(800)*(SI - SH)/(SI + SH), asym_paper[av], 1e-2)
    # ---- Chrzanowski-Misner null limit (their formulae transcribed as in the notebook, Sec. 7)
    rr_, m_, g2 = sp.symbols('r m gamma2', positive=True)
    mcrit = (2*sp.sqrt(3)/pi)*(rr_ + 3)/sp.sqrt(rr_)*g2                                   # CM Eq. (2.35b)
    PCM = 2*sp.sqrt(pi)*(1/g2)*(rr_ - 1)*sp.sqrt(3*rr_)/(rr_**2*(rr_ + 3)**2)*(mcrit/m_)*sp.exp(-pi/2)   # CM Eq. (4.28), eps = 1, k = 0, per E^2
    nullCM = sp.simplify(m_*PCM)
    rep.ident('Chrzanowski-Misner null limit: m P_m/E^2 = (12/sqrt(pi)) e^{-pi/2} (r0-M)/[r0^2 (r0+3M)]',
              sp.simplify(nullCM - 12/sp.sqrt(pi)*sp.exp(-pi/2)*(rr_ - 1)/(rr_**2*(rr_ + 3))))
    rep.ident('its spin dependence coincides with Eq. (21)', sp.simplify(nullCM/nullCM.subs(rr_, 3) - gS_closed.subs(r0, rr_)))
    ratioCM = float(nullCM.subs(rr_, 3))/float(kappa_j(0))
    rep.info('CM constant / kappa_0 (single-m kappa_0)', ratioCM)
    # "Their power at a given frequency collects the modes m and -m (we checked this on the scalar formula of
    #  Ref. [BCHM], which their master formula reproduces, against an exact scalar flux)": BCHM 1973 Eq. (5.4),
    #  P_m = (1/(27 pi^{5/2})) (m/m_crit) e^{-pi eps/4} |Gamma(1/4 + i eps/4)|^2 with m_crit = 4/(pi delta'),
    #  eps = 1 + (4/pi) m/m_crit = 1 + m delta', for a scalar charge on the circular orbit r0 = (3 + delta') M
    #  (their delta' = 3 x the delta of this paper), per q^2, against the exact single-m scalar flux of scalar_mode.
    def bchm54(mv, dB):
        mcrit = 4/(mp.pi*dB); e = 1 + (4/mp.pi)*mv/mcrit
        return (1/(27*mp.pi**mp.mpf(2.5)))*(mv/mcrit)*mp.exp(-mp.pi*e/4)*abs(mp.gamma(mp.mpf(1)/4 + 1j*e/4))**2
    # validation of the scalar integrator in the weak field: a scalar charge on a circular orbit of radius r0 >> M
    # radiates the dipole power q^2 acc^2/3 = q^2 M^2/(3 r0^4) (acc = M/r0^2), shared equally by m = +1 and m = -1,
    # so the single-m flux at l = m = 1 tends to M^2/(6 r0^4), with O(M/r0) corrections
    wf = {rv: scalar_mode(1, 1, rv)*6*rv**4 for rv in (50.0, 100.0, 200.0, 400.0)}
    rep.check('scalar integrator, weak field: single-m flux at l = m = 1, r0 = 100 M, over M^2/(6 r0^4) = 0.955', wf[100.0], '0.955', 1e-3)
    rep.info('  ... tends to 1 as 1 - c M/r0 with c -> 5.0: (1 - ratio) r0/M at r0/M = 50, 100, 200, 400',
             ', '.join(f'{(1 - v)*rv:.2f}' for rv, v in wf.items()), note='Richardson in 1/r0 from 200 and 400: ' + f'{2*(1 - wf[400.0])*400 - (1 - wf[200.0])*200:.2f}')
    num40 = scalar_mode(40, 40, 3.01)
    rep.check('BCHM Eq. (5.4) / exact single-m scalar flux, l = m = 40, delta\' = 0.01 (eps = 1.4): 2.00, i.e. their power collects m and -m',
              float(bchm54(40, mp.mpf('0.01')))/num40, '2.00', 5e-3, note=f'exact single-m flux {num40:.6e} q^2, BCHM {float(bchm54(40, mp.mpf("0.01"))):.6e}')
    others = {(100, '0.003'): None, (400, '0.001'): None}
    for (mv, dB) in others: others[(mv, dB)] = float(bchm54(mv, mp.mpf(dB)))/scalar_mode(mv, mv, 3 + float(dB))
    rep.info('  ... the same ratio at (l, delta\') = (100, 0.003), (400, 0.001)  [O(m^-1/2) deviations from 2]',
             ', '.join(f'{v:.3f}' for v in others.values()))
    rep.check('CM constant / (2 kappa_0) = 0.94 (paper: "their constant is 0.94 times ours", per |m|; Introduction: agrees to 6%)', ratioCM/2, '0.94', 2e-3,
              note=f'recomputed {ratioCM/2:.4f}')
    # the origin of the factor: CM master formula (5.1)-(5.2) vs BCHM scalar (5.3) and BRTV tensor Eq. (22)
    mp.mp.dps = 30
    def epsCM(mv, d): return 1 + 3*mv*d
    def mcrit0(d): return 4/(3*mp.pi*d)
    def PBCHM(mv, d):   # BCHM 1973 Eq. (5.3), exact scalar, per mu^2
        e = epsCM(mv, d); return (1/(27*mp.pi**mp.mpf(2.5)))*(mv/mcrit0(d))*mp.exp(-mp.pi*e/4)*abs(mp.gamma(mp.mpf(1)/4 + 1j*e/4))**2
    ACM3 = 2*(3 - 1)*mp.sqrt(9)/(mp.sqrt(mp.pi)*9*36)   # CM Eq. (5.2) at r = 3
    def PCMs(sv, mv, d):  # CM Eq. (5.1), k = 0 term
        e = epsCM(mv, d); return ACM3*mp.factorial(sv)**2*(4*mv/(e*mp.pi*mcrit0(d)))**(1 - sv)*mp.sqrt(e)*mp.exp(-mp.pi*e/2)
    def Peven(mv, d):     # Eq. (22), per mu^2
        eta = mp.mpf(1)/2 + 3*mv*d/2
        return mp.exp(-mp.pi*eta/2)/(162*mp.pi**mp.mpf(1.5)*mv*d)*abs((eta + mp.mpf(1)/2)*mp.gamma(mp.mpf(1)/4 + 1j*eta/2) + mp.sqrt(2)*(1 - 1j)/mp.sqrt(3*mv)*mp.gamma(mp.mpf(3)/4 + 1j*eta/2))**2
    dd = mp.mpf('1e-3'); mv = mp.mpf(10)**6
    rep.check('massive regime (eta >> 1): CM scalar / BCHM exact scalar -> 1 (m delta = 1000)', PCMs(0, mv, dd)/PBCHM(mv, dd), '1', 1e-3)
    rep.check('massive regime: CM tensor / BRTV exact tensor (single m) -> 8, i.e. 4 x the exact per-|m| flux (paper: factor 4)', PCMs(2, mv, dd)/Peven(mv, dd), '8', 1e-3)
    rep.check('massive regime: CM tensor/scalar = (s!)^2 = 4 (-> 4 as eps -> inf: ratio = 4 (eps/(eps-1))^2)', PCMs(2, mv, dd)/PCMs(0, mv, dd), '4', 1e-3)
    rep.check('massive regime: exact BRTV (single m) / BCHM (per |m|) -> 1/2, i.e. equal per |m|', Peven(mv, dd)/PBCHM(mv, dd), '0.5', 3e-3)
    dd = mp.mpf('1e-20'); mv = mp.mpf(10)**12
    rep.check('null limit: CM tensor / BRTV even (-> 1.877)', PCMs(2, mv, dd)/Peven(mv, dd), '1.877', 5e-4)
    wkb_exact = (mp.mpf(1)/4)*(2/mp.pi)*mp.exp(-mp.pi/2)/Fpar(mp.mpf(1)/2)
    rep.check('WKB barrier factor eta^{3/2} e^{-pi eta} / exact (eta+1/2)^2 F(eta) at eta = 1/2: 0.2347 (=1/4.26)', wkb_exact, '0.2347', 5e-4)
    rep.check('1/0.2347 = 4 x 1.065 (paper)', 1/wkb_exact, '4.26', 2e-3)
    rep.check('4 x 0.2347 = 0.94 (paper: the Stirling factor that remains per |m|)', 4*wkb_exact, '0.9386', 5e-4)
    rep.check('8 x 0.2347 = 1.877 (per single m)', 8*wkb_exact, '1.877', 5e-4)
    # ---- Lambda = (m - a omega)^2 + (2j+1) omega beta_b + O(1), |S(pi/2)|^2 -> c_j sqrt(omega beta_b)/(2 pi^{3/2}) [and |S'(pi/2)|^2 ->
    # d_j (omega beta_b)^{3/2}/pi^{3/2}, S'/S -> (2j+1) s (a+b)/beta_b, S/S' -> -s (a+b)/(omega beta_b^2), the inputs of alpha_expand]
    # "checked against the Toolkit's harmonics": here against the spectral spheroidal-harmonic solver of kerr_py (swsh.py, itself
    # validated against the Toolkit at l <= 6 and the source of the eigenvalues stored with the Kerr runs), and against the stored
    # Toolkit harmonics (notebook/toolkit_harmonics.m) when that file exists. Retrograde orbits: a -> -|a|, m > 0.
    def harm_targets(aa, lv, jv):
        rv = 2*(1 + math.cos(2/3*math.acos(-aa))); bv = math.sqrt(rv)*(rv + 3)/2; bbv = math.sqrt(3)*rv
        mv = lv - jv; wv = mv/bv
        lam_t = (mv - aa*wv)**2 + (2*jv + 1)*wv*bbv
        S2_t = float(c_j(jv))*math.sqrt(wv*bbv)/(2*math.pi**1.5) if jv % 2 == 0 else None
        dS2_t = float(d_j(jv))*(wv*bbv)**1.5/math.pi**1.5 if jv % 2 else None
        ratio_t = (2*jv + 1)*s*(aa + bv)/bbv if jv % 2 == 0 else -s*(aa + bv)/(wv*bbv**2)   # S'/S (even j) or S/S' (odd j)
        return dict(m=mv, w=wv, lam=lam_t, S2=S2_t, dS2=dS2_t, ratio=ratio_t)
    spins = ((0.5, 1), (0.5, -1), (0.9, 1), (0.9, -1)); lvals = (100, 200, 400, 800)
    if SWSH_at_equator is None:
        rep.unavailable('Lambda and |S(pi/2)|^2 asymptotics vs the spectral spheroidal harmonics of kerr_py', 'kerr_py/swsh.py could not be imported')
    else:
        res = {}
        for av, sg in spins:
            for lv in lvals:
                for jv in range(4):
                    aa = sg*av; tg = harm_targets(aa, lv, jv)
                    lam_s, S0, dS0, d2S0 = SWSH_at_equator(s, lv, tg['m'], aa*tg['w'])
                    rat = S0**2/tg['S2'] if jv % 2 == 0 else dS0**2/tg['dS2']
                    ratio_in = (dS0/S0)/tg['ratio'] if jv % 2 == 0 else (S0/dS0)/tg['ratio']
                    res[(aa, lv, jv)] = dict(rem=lam_s - tg['lam'], rat=rat, ratio_in=ratio_in, lam=lam_s)
        stored_lam = {(rw['a']*rw['sign'], rw['l'], rw['j']): rw['lam'] for rw in KERR if rw['l'] >= 100 and rw['j'] <= 3}
        rep.check('spectral solver reproduces the eigenvalues stored with the Kerr runs (a = +-0.5, +-0.9, l = 100..800, j <= 3): max |dLambda|/Lambda',
                  max(abs(v['lam']/stored_lam[k] - 1) for k, v in res.items() if k in stored_lam), '0', 1e-9, 'abs')
        rep.check('Lambda - (m - a omega)^2 - (2j+1) omega beta_b = O(1) (spectral solver, a = +-0.5, +-0.9, l = 100..800, j <= 3): max |remainder| (bound 20)',
                  max(abs(v['rem']) for v in res.values()), '20', 0, 'max', note='j = 0 at l = 800: ' + ', '.join(f'{sg*av}: {res[(sg*av, 800, 0)]["rem"]:+.2f}' for av, sg in spins))
        rep.check('  ... and l-independent: max |remainder(800) - remainder(400)|',
                  max(abs(res[(sg*av, 800, jv)]['rem'] - res[(sg*av, 400, jv)]['rem']) for av, sg in spins for jv in range(4)), '0', 0.1, 'abs')
        w800 = max(abs(res[(sg*av, 800, jv)]['rat'] - 1) for av, sg in spins for jv in range(4))
        wrich = max(abs(2*res[(sg*av, 800, jv)]['rat'] - res[(sg*av, 400, jv)]['rat'] - 1) for av, sg in spins for jv in range(4))
        rep.check('|S(pi/2)|^2 / [c_j sqrt(omega beta_b)/(2 pi^{3/2})] (even j) and |S\'(pi/2)|^2 / [d_j (omega beta_b)^{3/2}/pi^{3/2}] (odd j) -> 1: max |ratio - 1| at l = 800 (O(1/l) corrections; bound 0.07)',
                  w800, '0.07', 0, 'max', note='largest for j = 3 at a = 0.9; j = 0: ' + ', '.join(f'{sg*av}: {res[(sg*av, 800, 0)]["rat"] - 1:+.1e}' for av, sg in spins))
        rep.check('  ... Richardson extrapolation in 1/l from l = 400, 800: max |ratio - 1| (all spins, j <= 3)', wrich, '0', 5e-3, 'abs')
        rep.check('  ... l (ratio - 1) is l-independent: max |l (ratio-1)|_800 / |l (ratio-1)|_400 - 1|',
                  max(abs((800*(res[(sg*av, 800, jv)]['rat'] - 1))/(400*(res[(sg*av, 400, jv)]['rat'] - 1)) - 1) for av, sg in spins for jv in range(4)), '0', 0.1, 'abs')
        wi800 = max(abs(res[(sg*av, 800, jv)]['ratio_in'] - 1) for av, sg in spins for jv in range(4))
        wirich = max(abs(2*res[(sg*av, 800, jv)]['ratio_in'] - res[(sg*av, 400, jv)]['ratio_in'] - 1) for av, sg in spins for jv in range(4))
        rep.check('inputs of the source expansion: S\'/S -> (2j+1) s (a+b)/beta_b (even j), S/S\' -> -s (a+b)/(omega beta_b^2) (odd j): max |ratio - 1| at l = 800 (bound 0.02)',
                  wi800, '0.02', 0, 'max')
        rep.check('  ... Richardson extrapolation in 1/l from l = 400, 800', wirich, '0', 3e-3, 'abs')
    if HARM_TK is None:
        rep.unavailable('Lambda and |S(pi/2)|^2 asymptotics against the stored Toolkit harmonics (notebook/toolkit_harmonics.m)', 'file not found (not yet produced)')
    else:
        # rows {a, sign, l, j, Lambda, S2, dS2} with S2 = |S(pi/2)|^2, dS2 = |S'(pi/2)|^2 (S normalised to unity on the sphere); the stored
        # Lambda and |S|^2 do not depend on whether the retrograde case is represented as (a, -m, -omega) or (-a, m, omega)
        TKH = {}
        for rw in HARM_TK:
            av, sg, lv, jv = float(rw[0]), int(rw[1]), int(rw[2]), int(rw[3])
            TKH[(sg*av, lv, jv)] = (float(rw[4]), float(rw[5]), float(rw[6]))
        keys = sorted(TKH)
        rem = {k: TKH[k][0] - harm_targets(k[0], k[1], k[2])['lam'] for k in keys}
        rep.info('stored Toolkit harmonics: (a, l, j) rows', f'{len(keys)} rows, a = {sorted(set(k[0] for k in keys))}, l = {sorted(set(k[1] for k in keys))}, j <= {max(k[2] for k in keys)}')
        rep.check('Toolkit: Lambda - (m - a omega)^2 - (2j+1) omega beta_b = O(1): max |remainder| (bound 20)', max(abs(v) for v in rem.values()), '20', 0, 'max',
                  note='j = 0: ' + ', '.join(f'{k[0]},l={k[1]}: {rem[k]:+.2f}' for k in keys if k[2] == 0))
        lmax_t = max(k[1] for k in keys); lhalf = max(k[1] for k in keys if k[1] <= lmax_t/2) if any(k[1] <= lmax_t/2 for k in keys) else None
        if lhalf:
            rep.check(f'Toolkit: the remainder is l-independent: max |rem({lmax_t}) - rem({lhalf})|',
                      max(abs(rem[(a_, lmax_t, j_)] - rem[(a_, lhalf, j_)]) for (a_, l_, j_) in keys if l_ == lmax_t and (a_, lhalf, j_) in rem), '0', 0.1, 'abs')
        ratT = {}
        for k in keys:
            tg = harm_targets(*k)
            ratT[k] = TKH[k][1]/tg['S2'] if k[2] % 2 == 0 else TKH[k][2]/tg['dS2']
        rep.check(f'Toolkit: |S(pi/2)|^2 -> c_j sqrt(omega beta_b)/(2 pi^3/2) (even j), |S\'|^2 -> d_j (omega beta_b)^3/2/pi^3/2 (odd j): max |ratio - 1| at l = {lmax_t} (bound 0.07)',
                  max(abs(ratT[k] - 1) for k in keys if k[1] == lmax_t), '0.07', 0, 'max')
        if lhalf:
            rep.check(f'  ... Richardson extrapolation in 1/l from l = {lhalf}, {lmax_t}: max |ratio - 1|',
                      max(abs((lmax_t*ratT[(a_, lmax_t, j_)] - lhalf*ratT[(a_, lhalf, j_)])/(lmax_t - lhalf) - 1) for (a_, l_, j_) in keys if l_ == lmax_t and (a_, lhalf, j_) in ratT), '0', 5e-3, 'abs')
        if SWSH_at_equator is not None:
            wl = 0; wS = 0
            for k in keys:
                tg = harm_targets(*k)
                lam_s, S0, dS0, d2S0 = SWSH_at_equator(s, k[1], tg['m'], k[0]*tg['w'])
                wl = max(wl, abs(lam_s/TKH[k][0] - 1)); wS = max(wS, abs(S0**2/TKH[k][1] - 1) if k[2] % 2 == 0 else abs(dS0**2/TKH[k][2] - 1))
            rep.check('Toolkit vs spectral solver: eigenvalues (max rel. dev.)', wl, '0', 1e-8, 'abs')
            rep.check('Toolkit vs spectral solver: |S(pi/2)|^2 (even j), |S\'(pi/2)|^2 (odd j) (max rel. dev.)', wS, '0', 1e-6, 'abs')


def alpha_expand(a, b, bb, D0, P0, jv, s=-2):
    """Large-omega expansion of the circular-orbit Teukolsky source projection alpha_lm (Hughes 2000, Toolkit
    conventions) for the null source, transcribed from the companion notebook (alphaExpand). Unknowns: Y(0) and
    Y'(0) = omega^{1/2} y1 at the light ring; S0 = S(pi/2) with S'(pi/2)/S0 = (2j+1) s (a+b)/beta_b (even j),
    S0/S'(pi/2) = -s (a+b)/(omega beta_b^2) (odd j); S'' from the angular equation. Returns the coefficients of
    omega^2, omega^{3/2} (must vanish), Ahat_j (coefficient of omega S0 Y(0)) and Bhat (of omega^{1/2} S0' y1)."""
    w, u, SS, YY0, yy1, lam0, rr = sp.symbols('omega u S Y0 y1 lambda_0 r', positive=True)
    m = w*b; c = a*w; lam = (m - c)**2 + (2*jv + 1)*w*bb + lam0; A = lam - c**2 + 2*m*c; Kt = w*P0; L = -m + c
    if jv % 2 == 0: S0 = SS; dS0 = (2*jv + 1)*s*(a + b)/bb*SS
    else: dS0 = SS; S0 = -s*(a + b)/bb**2*SS/w
    d2S0 = (m**2 - s - A)*S0
    L2S = dS0 + L*S0; L1L2S = d2S0 + 2*L*dS0 - 2*S0 + L**2*S0
    rho = -1/r0; rhob = -1/r0; Sig = r0**2
    Ann0 = -rho**-2*rhob**-1*(sp.sqrt(2)*D0)**-2*(rho**-1*L1L2S + 2*I*a*L2S)
    Anmb0 = rho**-3*(sp.sqrt(2)*D0)**-1*((2*rho - I*Kt/D0)*L2S)
    Anmb1 = -rho**-3*(sp.sqrt(2)*D0)**-1*L2S
    Ambmb0 = (Kt**2*S0*rhob)/(4*D0**2*rho**3) + (I*Kt*S0*(1 - r0 + D0*rho)*rhob)/(2*D0**2*rho**3) + (I*r0*S0*rhob*w)/(2*D0*rho**3)
    Ambmb1 = -rho**-3*rhob*S0/2*(I*Kt/D0 - rho)
    Ambmb2 = -rho**-3*rhob*S0/4
    rc = P0/(2*Sig); tc = I*(b - a)/(sp.sqrt(2)*r0)
    Cnn, Cnmb, Cmbmb = rc**2, rc*tc, tc**2
    pref = (rr**2 - 2*rr + a**2)*(rr**2 + a**2)**sp.Rational(-1, 2)
    pref0 = pref.subs(rr, r0); dpref0 = sp.diff(pref, rr).subs(rr, r0)
    R = pref0*YY0; dR = dpref0*YY0 + pref0*((r0**2 + a**2)/D0)*u*yy1
    d2R = (-(-lam + 2*I*r0*s*2*w + (-2*I*(-1 + r0)*s*(-a*m + (a**2 + r0**2)*w) + (-a*m + (a**2 + r0**2)*w)**2)/D0)*R - (-2 + 2*r0)*(1 + s)*dR)/D0
    alpha = (Ann0*Cnn + Anmb0*Cnmb + Ambmb0*Cmbmb)*R - (Anmb1*Cnmb + Ambmb1*Cmbmb)*dR + Ambmb2*Cmbmb*d2R
    poly = sp.expand(alpha.subs(w, u**2)*u**8)
    # robust power-by-power collection in u (products like u**11/u**2 are reduced term by term)
    by_power = {}
    for term in sp.Add.make_args(poly):
        cft, ex = term.as_coeff_exponent(u)
        by_power[ex] = by_power.get(ex, 0) + cft
    co = lambda p: sp.expand(by_power.get(sp.Integer(2*p + 8), sp.Integer(0)))
    return dict(w2=co(2), w32=co(sp.Rational(3, 2)), Ahat=co(1).coeff(SS*YY0), Bhat=co(sp.Rational(1, 2)).coeff(SS*yy1))


def alpha_expand_general(a, b, r0, jv, s=-2, Esym=1, bwave=None):
    """The same expansion as alpha_expand, but with the orbit data free: spin a, impact parameter b = L/E and radius r0 are
    independent symbols (the orbit need not be a null geodesic, so R~(r0) = P0^2 - Delta_0 (b-a)^2 is not zero), the particle
    momentum is multiplied by Esym (p^mu -> E p^mu, L = bE), and the wave frequency is omega = m/b_w with b_w = b unless given
    (b_w = 1/Omega; Omega = 1/b for the light ring). The angular relations S'/S follow from the equatorial oscillator with the
    wave's m and omega, as in alpha_expand. Returns the coefficients of omega^2, omega^{3/2}, omega, omega^{1/2}, the full
    alpha and the symbols S, Y0, y1 used."""
    w, u, SS, YY0, yy1, lam0, rr = sp.symbols('omega u S Y0 y1 lambda_0 r', positive=True)
    if bwave is None: bwave = b
    D0 = r0**2 - 2*r0 + a**2; P0 = r0**2 + a**2 - a*b; Pw = r0**2 + a**2 - a*bwave; bb = sp.sqrt(bwave**2 - a**2)
    m = w*bwave; c = a*w; lam = (m - c)**2 + (2*jv + 1)*w*bb + lam0; A = lam - c**2 + 2*m*c; Kt = w*Pw; L = -m + c
    if jv % 2 == 0: S0 = SS; dS0 = (2*jv + 1)*s*(a + bwave)/bb*SS
    else: dS0 = SS; S0 = -s*(a + bwave)/bb**2*SS/w
    d2S0 = (m**2 - s - A)*S0
    L2S = dS0 + L*S0; L1L2S = d2S0 + 2*L*dS0 - 2*S0 + L**2*S0
    rho = -1/r0; rhob = -1/r0; Sig = r0**2
    Ann0 = -rho**-2*rhob**-1*(sp.sqrt(2)*D0)**-2*(rho**-1*L1L2S + 2*I*a*L2S)
    Anmb0 = rho**-3*(sp.sqrt(2)*D0)**-1*((2*rho - I*Kt/D0)*L2S)
    Anmb1 = -rho**-3*(sp.sqrt(2)*D0)**-1*L2S
    Ambmb0 = (Kt**2*S0*rhob)/(4*D0**2*rho**3) + (I*Kt*S0*(1 - r0 + D0*rho)*rhob)/(2*D0**2*rho**3) + (I*r0*S0*rhob*w)/(2*D0*rho**3)
    Ambmb1 = -rho**-3*rhob*S0/2*(I*Kt/D0 - rho)
    Ambmb2 = -rho**-3*rhob*S0/4
    rc = Esym*P0/(2*Sig); tc = I*Esym*(b - a)/(sp.sqrt(2)*r0)       # particle momentum on the tetrad: propto E (r^2+a^2) - a L and L - a E
    Cnn, Cnmb, Cmbmb = rc**2, rc*tc, tc**2
    pref = (rr**2 - 2*rr + a**2)*(rr**2 + a**2)**sp.Rational(-1, 2)
    pref0 = pref.subs(rr, r0); dpref0 = sp.diff(pref, rr).subs(rr, r0)
    R = pref0*YY0; dR = dpref0*YY0 + pref0*((r0**2 + a**2)/D0)*u*yy1
    d2R = (-(-lam + 2*I*r0*s*2*w + (-2*I*(-1 + r0)*s*(-a*m + (a**2 + r0**2)*w) + (-a*m + (a**2 + r0**2)*w)**2)/D0)*R - (-2 + 2*r0)*(1 + s)*dR)/D0
    alpha = (Ann0*Cnn + Anmb0*Cnmb + Ambmb0*Cmbmb)*R - (Anmb1*Cnmb + Ambmb1*Cmbmb)*dR + Ambmb2*Cmbmb*d2R
    poly = sp.expand(alpha.subs(w, u**2)*u**8)
    by_power = {}
    for term in sp.Add.make_args(poly):
        cft, ex = term.as_coeff_exponent(u)
        by_power[ex] = by_power.get(ex, 0) + cft
    co = lambda p: sp.expand(by_power.get(sp.Integer(2*p + 8), sp.Integer(0)))
    return dict(w2=co(2), w32=co(sp.Rational(3, 2)), w1=co(1), w12=co(sp.Rational(1, 2)), alpha=alpha, S=SS, Y0=YY0, y1=yy1)

# =============================================================================================
# Section V  Comparison with the 1973 formulae
# =============================================================================================
def section_V():
    rep.section('Sec. V  Comparison with the geodesic synchrotron radiation literature  [Eqs. (22)-(23)]')
    mp.mp.dps = 30
    def Peven(mv, d):
        eta = mp.mpf(1)/2 + 3*mv*d/2
        return mp.exp(-mp.pi*eta/2)/(162*mp.pi**mp.mpf(1.5)*mv*d)*abs((eta + mp.mpf(1)/2)*mp.gamma(mp.mpf(1)/4 + 1j*eta/2) + mp.sqrt(2)*(1 - 1j)/mp.sqrt(3*mv)*mp.gamma(mp.mpf(3)/4 + 1j*eta/2))**2
    def Podd(mv, d):
        eta = mp.mpf(3)/2 + 3*mv*d/2
        return mp.exp(-mp.pi*eta/2)/(162*mp.pi**mp.mpf(1.5)*mv*d)*abs(mp.sqrt(2)*mp.gamma(mp.mpf(3)/4 + 1j*eta/2) + (1 + 1j)/(2*mp.sqrt(3*mv))*mp.gamma(mp.mpf(1)/4 + 1j*eta/2))**2
    # null limit: multiply by mu^2/E^2 = 9 delta and by m, delta -> 0, m -> infinity
    limEven = 9*mp.exp(-mp.pi/4)/(162*mp.pi**mp.mpf(1.5))*abs(mp.gamma(mp.mpf(1)/4 + 1j/4))**2
    limOdd = 9*mp.exp(-3*mp.pi/4)/(162*mp.pi**mp.mpf(1.5))*2*abs(mp.gamma(mp.mpf(3)/4 + 3j/4))**2
    rep.check('Eq. (22) null limit 9 delta m P^even -> kappa_0 exactly (closed-form limit / kappa_0)', limEven/kappa_j(0), '1', 1e-25)
    rep.check('Eq. (22) null limit, direct evaluation at delta = 1e-20, m = 1e12', 9*mp.mpf('1e-20')*mp.mpf(10)**12*Peven(mp.mpf(10)**12, mp.mpf('1e-20'))/kappa_j(0), '1', 1e-5)
    rep.check('Eq. (23) null limit 9 delta m P^odd -> kappa_1/4 exactly (closed-form limit / kappa_1)', limOdd/kappa_j(1), '0.25', 1e-25)
    rep.check('Eq. (23) null limit, direct evaluation at delta = 1e-20, m = 1e12', 9*mp.mpf('1e-20')*mp.mpf(10)**12*Podd(mp.mpf(10)**12, mp.mpf('1e-20'))/kappa_j(1), '0.25', 1e-5)
    # the exactness rests on the reflection formula |Gamma(1/4+ix)|^2 |Gamma(3/4+ix)|^2 = 2 pi^2/cosh(2 pi x)
    xs = sp.symbols('x', positive=True)
    rep.ident('reflection identity |Gamma(1/4+ix) Gamma(3/4+ix)|^2 = 2 pi^2/cosh(2 pi x) (from Gamma(z)Gamma(1-z) = pi/sin(pi z))',
              sp.simplify(sp.expand_complex(sp.Abs(pi/sp.sin(pi*(sp.Rational(1, 4) + I*xs)))**2) - 2*pi**2/sp.cosh(2*pi*xs)))
    # "the second term inside the modulus [of Eq. (22)] is the m^{-1/2} cross term of the flux at infinity": expanding the modulus,
    # |A + B|^2 = |A|^2 (1 + 2 Re(B/A) + ...) with A = (eta+1/2) Gamma(1/4+i eta/2), B = sqrt2 (1-i) Gamma(3/4+i eta/2)/sqrt(3m), the relative
    # O(m^-1/2) term of Eq. (22) at eta = 1/2 is m^{-1/2} x 2 Re[sqrt2 (1-i) Gamma(3/4+i/4)/(sqrt3 Gamma(1/4+i/4))], to be compared with the
    # cross term of sigma_0 (from the phase of Eq. 14, Sec. III C), which enters the flux at infinity with the + sign
    cross22 = 2*mp.re(mp.sqrt(2)*(1 - 1j)*mp.gamma(mp.mpf(3)/4 + 1j/4)/(mp.sqrt(3)*mp.gamma(mp.mpf(1)/4 + 1j/4)))
    rep.check('Eq. (22): relative m^{-1/2} term 2 Re[sqrt2 (1-i) Gamma(3/4+i/4)/(sqrt3 Gamma(1/4+i/4))] vs the cross term of sigma_0 from Eq. (14) (sign included)',
              cross22, SIGMA[('cross', 0)], 1e-8, note=f'{float(cross22):.7f}')
    rep.check('  ... = 1.0390 (Sec. III C)', cross22, '1.0390', 1e-4)
    cross23 = 2*mp.re((1 + 1j)*mp.gamma(mp.mpf(1)/4 + 3j/4)/(2*mp.sqrt(6)*mp.gamma(mp.mpf(3)/4 + 3j/4)))
    rep.check('odd analogue (not stated in the paper): relative m^{-1/2} term of Eq. (23), 2 Re[(1+i) Gamma(1/4+3i/4)/(2 sqrt6 Gamma(3/4+3i/4))], vs the cross term of sigma_1',
              cross23, SIGMA[('cross', 1)], 1e-8, note=f'{float(cross23):.7f}: the 1973 odd formula has the right cross term, only its normalisation is off by 4')
    rep.info('kappa with the 1973 odd term, j <= 1 only: 2(kappa_0 + kappa_1/4) (= 0.0574; this number is no longer quoted in the paper)', 2*(kappa_j(0) + kappa_j(1)/4))
    rep.check('kappa with all odd kappa_j divided by 4: 0.0580', 2*mp.fsum(kappa_j(jv)/(1 if jv % 2 == 0 else 4) for jv in range(40)), '0.0580', 1e-3)
    rep.check('our kappa, 0.0637', 2*mp.fsum(kappa_j(jv) for jv in range(40)), '0.0637', 1e-3)
    rep.check('0.0580 is outside the 2021 fit 0.064 +- 0.001, 0.0637 is inside', float(abs(0.0580 - 0.064) > 0.001 and abs(0.0637 - 0.064) <= 0.001), '1', 0)
    # Sec. VI A: "Eq. (24) is the leading term of Eq. (22)" for every eta: with mu^2/E^2 -> 9 delta and the cross term dropped, l x Eq. (22)/E^2
    # -> e^{-pi eta/2} (eta+1/2)^2 |Gamma(1/4+i eta/2)|^2/(18 pi^{3/2}), and l x Eq. (24)/E^2 = kappa_0 (eta+1/2)^2 F(eta)/F(1/2); the two agree
    # identically because of the reflection identity (F(eta) = e^{-pi eta/2} |Gamma(1/4+i eta/2)|^2/(2 pi^2))
    def lP22_lead(eta): return mp.exp(-mp.pi*eta/2)*(eta + mp.mpf(1)/2)**2*abs(mp.gamma(mp.mpf(1)/4 + 1j*eta/2))**2/(18*mp.pi**mp.mpf(1.5))
    def lEq24(eta): return kappa_j(0)*(eta + mp.mpf(1)/2)**2*Fpar(eta)/Fpar(mp.mpf(1)/2)
    rep.check('Eq. (24) = leading term of Eq. (22) (cross term dropped, mu^2/E^2 -> 9 delta) at eta = 1/2, 1, 2, 5: max |ratio - 1|',
              max(abs(lP22_lead(mp.mpf(e))/lEq24(mp.mpf(e)) - 1) for e in ('0.5', '1', '2', '5')), '0', 1e-25, 'abs')
    etaS = sp.symbols('eta', positive=True)
    Frefl = sp.exp(-pi*etaS/2)*sp.Abs(sp.gamma(sp.Rational(1, 4) + I*etaS/2))**2/(2*pi**2)        # F(eta) via the reflection identity
    kappa0S = sp.sqrt(pi)/9*Frefl.subs(etaS, sp.Rational(1, 2))
    rep.ident('  ... symbolically: e^{-pi eta/2} |Gamma(1/4+i eta/2)|^2/(18 pi^{3/2}) = kappa_0 F(eta)/F(1/2) with F from the reflection identity',
              sp.simplify(sp.exp(-pi*etaS/2)*sp.Abs(sp.gamma(sp.Rational(1, 4) + I*etaS/2))**2/(18*pi**sp.Rational(3, 2)) - kappa0S*Frefl/Frefl.subs(etaS, sp.Rational(1, 2))))
    # the odd ratio in the 1973 regime: independent integrator, delta = 1e-4 (r0 = 3.0003), m = 40..400
    rep.note('Exact timelike fluxes of Sec. V (ratios 3.15, 3.48, 3.64, 3.77) are not stored in the repository: recomputed here with the\n'
             'independent Zerilli/RW integrator (per mu^2, flux at infinity), E0 = (r0-2)/sqrt(r0(r0-3)), L0 = r0/sqrt(r0-3), Omega = r0^-3/2.')
    delta = 1e-4; r0v = 3*(1 + delta)
    E0v = (r0v - 2)/math.sqrt(r0v*(r0v - 3)); L0v = r0v/math.sqrt(r0v - 3); Om = r0v**-1.5
    paper_ratios = {40: 3.15, 100: 3.48, 200: 3.64, 400: 3.77}
    ratios = {}; even_ratios = {}
    for mv in (40, 100, 200, 400):
        md_odd = zrw_mode(mv + 1, mv, r0v, E0v, L0v, Om)
        md_even = zrw_mode(mv, mv, r0v, E0v, L0v, Om)
        ratios[mv] = md_odd['I']/float(Podd(mv, delta)); even_ratios[mv] = md_even['I']/float(Peven(mv, delta))
        rep.check(f'exact odd flux / Eq. (23), delta = 1e-4, m = {mv}', ratios[mv], str(paper_ratios[mv]), 5e-3,
                  note=f'even: exact/Eq.(22) = {even_ratios[mv]:.4f}')
    mvs = np.array([40., 100., 200., 400.]); yv = np.array([ratios[int(x)] for x in mvs])
    A3 = np.vstack([np.ones(4), mvs**-0.5, 1/mvs]).T; A2 = np.vstack([np.ones(4), mvs**-0.5]).T
    rho3 = np.linalg.lstsq(A3, yv, rcond=None)[0][0]; rho2 = np.linalg.lstsq(A2, yv, rcond=None)[0][0]
    rep.check('three-term fit rho_inf + rho_1 m^-1/2 + rho_2/m of the recomputed ratios: rho_inf = 4.06', rho3, '4.06', 1.5e-2)
    rep.check('two-term fit: rho_inf = 4.05', rho2, '4.05', 1.5e-2)
    yq = np.array([3.15, 3.48, 3.64, 3.77])
    rep.check('same fits on the four quoted values (3.15, 3.48, 3.64, 3.77): three-term', np.linalg.lstsq(A3, yq, rcond=None)[0][0], '4.06', 5e-3)
    rep.check('same fits on the four quoted values: two-term', np.linalg.lstsq(A2, yq, rcond=None)[0][0], '4.05', 5e-3)
    rep.check('even formula accurate to O(m^-1/2): |exact/Eq.(22) - 1| x sqrt(m) bounded (m = 400)', abs(even_ratios[400] - 1)*math.sqrt(400), '0', 2.0, 'abs',
              note='values: ' + ', '.join(f'm={k}: {v:.4f}' for k, v in even_ratios.items()))
    # the same ratios from the stored Teukolsky run of kerr_py/check_odd_BRTV.py (columns r0, m, parity, exact, BRTV, ratio; r0 = 3.0003 is delta = 1e-4)
    if BRTV_STORED is None:
        rep.unavailable('stored Teukolsky odd-mode ratios 3.15, 3.48, 3.64, 3.77 (check_odd_BRTV_results.txt)', 'kerr_py/check_odd_BRTV_results.txt not found')
    else:
        st = {(rw[1], rw[2]): rw[5] for rw in BRTV_STORED if abs(rw[0] - 3.0003) < 1e-9}
        for mv in (40, 100, 200, 400):
            rep.check(f'stored Teukolsky: exact odd flux / Eq. (23), delta = 1e-4, m = {mv} (paper: {paper_ratios[mv]})', st[(mv, 'odd')], str(paper_ratios[mv]), 5e-3,
                      note=f'this script (Zerilli/RW integrator): {ratios[mv]:.4f}; even: {st[(mv, "even")]:.4f}')
        rep.check('stored Teukolsky vs this script\'s Zerilli/RW ratios: max rel. dev. (m = 40..400, odd)', max(abs(st[(mv, 'odd')]/ratios[mv] - 1) for mv in (40, 100, 200, 400)), '0', 2e-3, 'abs')
        ys = np.array([st[(int(x), 'odd')] for x in mvs])
        rep.check('stored: three-term fit rho_inf = 4.06', np.linalg.lstsq(A3, ys, rcond=None)[0][0], '4.06', 5e-3)
        rep.check('stored: two-term fit rho_inf = 4.05', np.linalg.lstsq(A2, ys, rcond=None)[0][0], '4.05', 5e-3)
        rep.check('stored: even formula accurate to O(m^-1/2): |exact/Eq.(22) - 1| sqrt(m) at m = 400 bounded', abs(st[(400, 'even')] - 1)*20, '0', 2.0, 'abs')
        oth = {(rw[1], rw[2]): rw[5] for rw in BRTV_STORED if abs(rw[0] - 3.003) < 1e-9}
        if oth: rep.info('  ... stored odd ratios at delta = 1e-3 (not quoted), m = 40, 100, 200, 400', ', '.join(f'{oth[(mv, "odd")]:.3f}' for mv in (40, 100, 200, 400) if (mv, 'odd') in oth))

# =============================================================================================
# Section VI  Discussion
# =============================================================================================
def section_VI():
    rep.section('Sec. VI  Discussion: timelike orbits [Eq. (24)], cut-off, M87* estimate')
    dd, sig = sp.symbols('delta sigma', positive=True)
    r0t = 3*(1 + dd); L0t = r0t/sp.sqrt(r0t - 3); E0t = (r0t - 2)/sp.sqrt(r0t*(r0t - 3))
    rep.ident('gamma^2 = E^2/mu^2 = (1+3 delta)^2/[9 delta (1+delta)]', sp.simplify(E0t**2 - (1 + 3*dd)**2/(9*dd*(1 + dd))))
    rep.ident('gamma = (9 delta)^-1/2 at leading order (Introduction)', sp.limit(E0t**2*9*dd, dd, 0) - 1)
    dPsit = jumps(l, l, r0t, E0t, L0t, Y, 0, True)[1]/E0t
    ser = series_in_l(dPsit.subs(dd, sig*eps), 2)        # delta = sigma/l
    rep.ident('[d_r Psi]/E = -8 pi Y (1 + 3 l delta/2)/(l M) for m = l, r0 = 3M(1+delta)', sp.simplify(ser.coeff(eps, 1) + 8*pi*Y*(1 + 3*sig/2)))
    rep.ident('  ... no O(l^0) term', ser.coeff(eps, 0))
    V0t = V_Z(l, r).subs(r, 3); kt = -Dst(Dst(V_Z(l, r))).subs(r, 3)
    etat = series_in_l(((V0t - l**2/r0t**3)/sp.sqrt(2*kt)).subs(dd, sig*eps), 1)
    rep.ident('eta = 1/2 + 3 l delta/2 for m = l', sp.simplify(etat.coeff(eps, 0) - sp.Rational(1, 2) - 3*sig/2))
    # Eq. (24) assembled: flux = (c+/4)|u(0)|^2 [d_r* Psi]^2 with eta shifted -> kappa_0/l (eta+1/2)^2 F(eta)/F(1/2)
    Fs, etas = sp.symbols('F eta', positive=True)
    rep.ident('Eq. (24): (1 + 3 l delta/2)^2 = (eta + 1/2)^2', sp.simplify((1 + 3*sig/2)**2 - (sp.Rational(1, 2) + 3*sig/2 + sp.Rational(1, 2))**2))
    # Eq. (24) vs stored timelike fluxes (per E^2, l = m), gamma = 5, 10, 20
    def Edot24(lv, d):
        eta = mp.mpf(1)/2 + 3*lv*d/2
        return kappa_j(0)/lv*(eta + mp.mpf(1)/2)**2*Fpar(eta)/Fpar(mp.mpf(1)/2)
    worst = 0; worstI = 0; txt = []
    for rw in TIMELIKE:
        rv, lv, EI, EH = rw; d = (rv - 3)/3; gam = math.sqrt((1 + 3*d)**2/(9*d*(1 + d)))
        if lv > 400: continue
        pred = float(Edot24(lv, mp.mpf(d)))
        dev = abs((EI + EH)/2/pred - 1); worst = max(worst, dev); worstI = max(worstI, abs(EI/pred - 1))
        txt.append(f'gamma={gam:.1f} l={lv}: mean/Eq.(24)-1={((EI + EH)/2/pred - 1):+.3f}, I/Eq.(24)-1={(EI/pred - 1):+.3f}')
    rep.info('Eq. (24) vs stored timelike fluxes (mean of I and H, gamma = 5, 10, 20, all stored l <= 400): max |ratio - 1|', worst, 'the paper quotes the agreement per gamma and l range, see the lines below')
    sub = [abs((rw[2] + rw[3])/2/float(Edot24(rw[1], mp.mpf((rw[0] - 3)/3))) - 1) for rw in TIMELIKE if 50 <= rw[1] <= 125 or (rw[1] <= 400 and (rw[0] - 3)/3 < 2e-3 and rw[1] >= 100)]
    rep.check('  ... restricted to 50 <= l <= 125 (gamma = 5), 100 <= l <= 300 (gamma = 10), 40 <= l <= 400 (gamma = 20): max |ratio - 1|', max(sub), '0.03', 0, 'max')
    for g0, lo, hi, tol in ((20.0, 120, 4800, 0.01), (10.0, 30, 1200, 0.04), (5.0, 25, 125, 0.04)):
        devs = [abs((rw[2] + rw[3])/2/float(Edot24(rw[1], mp.mpf((rw[0] - 3)/3))) - 1) for rw in TIMELIKE if lo <= rw[1] <= hi and abs(math.sqrt((1 + 3*(rw[0] - 3)/3)**2/(9*(rw[0] - 3)/3*(1 + (rw[0] - 3)/3))) - g0) < 0.5]
        rep.check(f'  ... paper statement: gamma = {g0:.0f}, {lo} <= l <= {hi}: max |mean/Eq.(24) - 1| <= {tol}', max(devs), str(tol), 0, 'max')
    rep.info('  ... with the flux to infinity alone (not the mean) the deviation is', worstI, 'the O(l^-1/2) asymmetry is not in Eq. (24)')
    rep.note('\n'.join(txt))
    # cut-off
    etas_ = sp.symbols('eta', positive=True)
    rep.check('F(eta) ~ e^{-pi eta}/(pi sqrt(eta/2)): ratio at eta = 1000', float(Fpar(1000)/(mp.exp(-1000*mp.pi)/(mp.pi*mp.sqrt(500)))), '1', 1e-3)
    rep.ident('3 l delta/2 = l/(6 gamma^2) at leading order (delta = 1/(9 gamma^2))', sp.simplify(3*l*(1/(9*sp.Symbol('gamma2')))/2 - l/(6*sp.Symbol('gamma2'))))
    gm2 = sp.symbols('gamma2', positive=True)
    cut = (etas_ + sp.Rational(1, 2))**2*sp.exp(-pi*etas_)/(pi*sp.sqrt(etas_/2))
    rep.ident('(eta+1/2)^2 F(eta) ~ eta^{3/2} e^{-pi eta} up to a constant: exponent of eta = 3/2',
              sp.limit(sp.log(cut*sp.exp(pi*etas_))/sp.log(etas_), etas_, sp.oo) - sp.Rational(3, 2))
    rep.check('c_1 = pi/6 = 0.5236 (coefficient of l/gamma^2 in the exponent)', math.pi/6, '0.5236', 1e-4)
    rep.check('c_1 = pi/6 vs the 2021 fitted value 0.42 +- 0.02 (paper: fitted value is smaller)', float(math.pi/6 > 0.44), '1', 0)
    sumk = mp.fsum(kappa_j(jv) for jv in range(40))
    rep.check('4 sum_j kappa_j = 0.127 (coefficient of ln gamma in each flux)', 4*sumk, '0.127', 3e-3)
    rep.check('  ... consistent with the 2021 slope k_1 = 0.12 +- 0.01', float(abs(float(4*sumk) - 0.12) <= 0.01), '1', 0)
    # "Summing Eq. (24) and its j >= 1 analogues over l and m up to the cut-off at l ~ gamma^2, each of the two fluxes is kappa ln gamma^2 + const = 0.127 ln gamma".
    # The m-sum needs the timelike formula for j > 0, obtained as Eq. (24): for m = l - j the barrier parameter is
    # eta_j = j + 1/2 + 3 l delta/2, the even derivative jump is -8 pi Y (2j+1 + 3 l delta/2)/(l M) (its O(l^0) massive term adds to
    # the O(1/l) term) and the odd jump [Psi] has no O(l delta) correction, so
    #   Edot_{l,l-j} = (kappa_j/l) [(2j+1 + 3 l delta/2)/(2j+1)]^2 F(eta_j)/F(j+1/2)   (even j),   (kappa_j/l) G(eta_j)/G(j+1/2)   (odd j)
    dPj = jumps(l, l - j, r0t, E0t, L0t, Y, 0, True)[1]/E0t
    serj = series_in_l(dPj.subs(dd, sig*eps), 2)
    rep.ident('timelike, m = l - j (even j): [d_r Psi]/E = -8 pi Y (2j+1 + 3 l delta/2)/(l M), no O(l^0) term',
              sp.simplify(serj.coeff(eps, 1) + 8*pi*Y*(2*j + 1 + 3*sig/2)) + serj.coeff(eps, 0))
    etaj = series_in_l(((V0t - (l - j)**2/r0t**3)/sp.sqrt(2*kt)).subs(dd, sig*eps), 1)
    rep.ident('timelike, m = l - j: eta_j = j + 1/2 + 3 l delta/2', sp.simplify(etaj.coeff(eps, 0) - j - sp.Rational(1, 2) - 3*sig/2))
    Pj_odd = jumps(l, l - j, r0t, E0t, L0t, 0, dY, False)[0]/E0t
    serodd = series_in_l((Pj_odd/Pj_odd.subs(dd, 0)).subs(dd, sig*eps), 1)
    rep.ident('timelike, odd j: [Psi](delta)/[Psi](0) = 1 + O(delta), no O(l delta) term (the mass enters the odd modes only through eta_j)',
              sp.simplify(serodd.coeff(eps, 0) - 1))
    from scipy.special import loggamma
    def logcosh(x): x = np.abs(x); return x + np.log1p(np.exp(-2*x)) - np.log(2)
    def Fnp(eta): return np.exp(-np.pi*eta/2 - logcosh(np.pi*eta) - 2*np.real(loggamma(0.75 + 0.5j*eta)))      # overflow-safe F, G
    def Gnp(eta): return np.exp(-np.pi*eta/2 - logcosh(np.pi*eta) - 2*np.real(loggamma(0.25 + 0.5j*eta)))
    kjf = [float(kappa_j(jv)) for jv in range(5)]; F0 = [float(Fpar(jv + 0.5)) for jv in range(5)]; G0 = [float(Gpar(jv + 0.5)) for jv in range(5)]
    def summand(lv, jv, d):
        eta = jv + 0.5 + 1.5*lv*d
        if jv % 2 == 0: return kjf[jv]/lv*((2*jv + 1 + 1.5*lv*d)/(2*jv + 1))**2*Fnp(eta)/F0[jv]
        return kjf[jv]/lv*Gnp(eta)/G0[jv]
    def S_total(gam, jmax=4, L1=20000):
        """2 x sum_{j<=jmax} sum_{l>=2} Edot_{l,l-j} (both signs of m) for r0 = 3M(1+delta), delta = 1/(9 gamma^2): exact sum to l = L1,
        Euler-Maclaurin integral beyond (the summand is smooth there), cut where eta_j = 60 (e^{-60 pi} ~ 1e-82)."""
        d = 1/(9*gam**2); tot = 0; lcut = 59.5/(1.5*d)
        for jv in range(jmax + 1):
            tot += np.sum(summand(np.arange(2, min(L1, int(lcut)) + 1.0), jv, d))
            if lcut > L1: tot += quad(lambda x: summand(x, jv, d), L1 + 0.5, lcut, limit=500)[0]
        return 2*tot
    Stot = {g: S_total(g) for g in (10, 30, 100, 300, 1000)}
    rep.check('each flux summed over l (and m, j <= 4) for r0 = 3M(1+delta), delta = 1/(9 gamma^2): slope in ln gamma between gamma = 100 and 1000 = 2 kappa = 0.1273',
              (Stot[1000] - Stot[100])/math.log(10), float(4*sumk), 1e-3,
              note='slopes 10->30, 30->100, 100->300, 300->1000: ' + ', '.join(f'{(Stot[g2] - Stot[g1])/math.log(g2/g1):.5f}' for g1, g2 in ((10, 30), (30, 100), (100, 300), (300, 1000))))
    rep.check('  ... = 0.127 ln gamma (paper)', (Stot[1000] - Stot[100])/math.log(10), '0.127', 3e-3,
              note='the constant, S - kappa ln gamma^2, at gamma = 100, 1000: ' + ', '.join(f'{Stot[g] - float(2*sumk)*math.log(g**2):+.5f}' for g in (100, 1000)) + ' (happens to be ~1e-4)')
    S0t = {g: S_total(g, jmax=0) for g in (100, 1000)}
    rep.check('  ... with Eq. (24) alone (j = 0, both signs of m) the slope would be 4 kappa_0 = 0.111: the 0.127 needs all j', (S0t[1000] - S0t[100])/math.log(10), float(4*kappa_j(0)), 1e-3)
    # "A pure exponential in l/gamma^2, fitted to l Edot over a finite range of l, necessarily returns a smaller effective coefficient, which explains the value
    # 0.42 +- 0.02": the local coefficient -d ln(l Edot_ll)/d(l/gamma^2) of Eq. (24) is below pi/6 for every finite l/gamma^2 (it tends to
    # pi/6 from below), and fits of ln(l Edot) = a - c_1 l/gamma^2 over the ranges of the 2021 fit (l <~ 400, gamma <~ 18) give c_1 ~ 0.35-0.43
    def log_lEdot24(lv, gam):      # ln of l Edot_ll/(kappa_0 E^2/M^2) from Eq. (24), without underflow
        eta = 0.5 + lv/(6*gam**2)
        return 2*np.log(eta + 0.5) - np.pi*eta/2 - logcosh(np.pi*eta) - 2*np.real(loggamma(0.75 + 0.5j*eta)) - math.log(F0[0])
    xs = np.array([0.25, 0.5, 1, 2, 4, 8, 16, 32, 100, 1000, 1e4]); h = 1e-5
    ceff = np.array([-(log_lEdot24((x + h)*100, 10.0) - log_lEdot24((x - h)*100, 10.0))/(2*h) for x in xs])
    rep.check('local coefficient -d ln(l Edot_ll)/d(l/gamma^2) of Eq. (24) is below pi/6 for every x = l/gamma^2 in [0.25, 1e4] and increases monotonically to pi/6',
              float(np.all(ceff < math.pi/6) and np.all(np.diff(ceff) > 0) and abs(ceff[-1] - math.pi/6) < 1e-3), '1', 0,
              note='c_eff at x = 0.25, 1, 4, 16, 100, 1e4: ' + ', '.join(f'{c:.3f}' for c in ceff[[0, 2, 4, 6, 8, 10]]))
    fits = {}
    for gam, lo, hi in ((5, 25, 400), (10, 30, 400), (10, 100, 400), (18, 100, 400)):
        ls_ = np.arange(lo, hi + 1.0); fits[(gam, lo, hi)] = -np.polyfit(ls_/gam**2, log_lEdot24(ls_, gam), 1)[0]
    rep.check('pure-exponential fits of l Eq. (24), ln(l Edot) = a - c_1 l/gamma^2, over l <= 400 at gamma = 5, 10, 18 (the 2021 ranges): all c_1 < pi/6',
              float(all(v < math.pi/6 for v in fits.values())), '1', 0, note='; '.join(f'gamma={k[0]}, l={k[1]}..{k[2]}: {v:.3f}' for k, v in fits.items()))
    rep.check('  ... and of the order of the 2021 value 0.42: mean of the four fits', float(np.mean(list(fits.values()))), '0.42', 0.1)
    # the same fit on the stored exact fluxes at infinity (timelike_results.m, l <= 400, as in the 2021 fit of the Teukolsky fluxes)
    sfits = {}
    for rv in sorted(set(rw[0] for rw in TIMELIKE)):
        d = (rv - 3)/3; gam = math.sqrt((1 + 3*d)**2/(9*d*(1 + d)))
        rows = [rw for rw in TIMELIKE if rw[0] == rv and rw[1] <= 400]
        ls_ = np.array([rw[1] for rw in rows], float)
        sfits[gam] = (-np.polyfit(ls_/gam**2, np.log([rw[1]*rw[2] for rw in rows]), 1)[0], int(ls_.min()), int(ls_.max()))
    rep.check('pure-exponential fit of the stored exact fluxes to infinity (gamma = 5, 10, 20; stored l <= 400): c_1 within 0.42 +- 0.03 for each gamma',
              float(all(abs(v[0] - 0.42) <= 0.03 for v in sfits.values())), '1', 0, note='; '.join(f'gamma={g:.1f} (l = {v[1]}..{v[2]}): {v[0]:.3f}' for g, v in sfits.items()))
    # Conclusions
    kap = 2*sumk; b = 3*mp.sqrt(3)
    rep.check('Delta E per orbit = 2 kappa (2 pi b) E^2/M ln l_max: coefficient 4.2', 2*kap*2*mp.pi*b, '4.2', 1.5e-2, note=f'exact: {float(2*kap*2*mp.pi*b):.4f}')
    # M87*: 1.3 mm photons, 6.5e9 solar masses (cgs)
    Gc, cc, hh, Msun = 6.6743e-8, 2.99792458e10, 6.62607015e-27, 1.98841e33
    Eph = hh*cc/0.13; EG = Gc*Eph/cc**4; MG = Gc*6.5e9*Msun/cc**2
    rep.ident('orbital period 2 pi b = 6 sqrt3 pi M', 2*pi*3*sp.sqrt(3) - 6*sp.sqrt(3)*pi)
    rep.check('M87*: E/M for a 1.3 mm photon and M = 6.5e9 Msun (paper: E/M ~ 1e-79)', EG/MG, '1e-79', 0.5,
              note=f'recomputed {EG/MG:.2e}')
    lmax = 2*math.pi*3*MG/0.13
    hb, lam_, c_ = sp.symbols('hbar lambda c', positive=True)
    rep.ident('l_max ~ E r0/hbar = 2 pi r0/lambda for a photon of wavelength lambda (E = 2 pi hbar c/lambda, c = 1)',
              sp.simplify((2*pi*hb*c_/lam_)*r0/(hb*c_) - 2*pi*r0/lam_))
    rep.check('M87*: l_max ~ E r0/hbar = 2 pi r0/lambda (paper: ~1e17)', lmax, '1e17', 0.9,
              note='order-of-magnitude statement: recomputed 1.4e17, i.e. log10 = 17.1')
    rep.check('M87*: Delta E/E = 4.2 (E/M) ln l_max ~ 1e-77 per orbit', 4.2*EG/MG*math.log(lmax), '1e-77', 2.0,
              note=f'recomputed {4.2*EG/MG*math.log(lmax):.1e}; with l_max = 1e18: {4.2*EG/MG*math.log(1e18):.1e}')
    # Sec. VI B: "the flux is proportional to ln(E M/m_Pl^2)": with G, c, hbar restored, l_max ~ E r0/(hbar c) with r0 = 3 G M/c^2 and
    # m_Pl^2 = hbar c/G gives l_max = 3 (E/c^2) M/m_Pl^2, so ln l_max = ln(E M/m_Pl^2) + ln 3 in G = c = 1 units
    Gs, cs_, hbs, Es, Ms = sp.symbols('G c hbar E M', positive=True)
    mPl2 = hbs*cs_/Gs
    rep.ident('l_max = E r0/(hbar c) with r0 = 3 G M/c^2 equals 3 (E/c^2) M/m_Pl^2, m_Pl^2 = hbar c/G', sp.simplify(Es*(3*Gs*Ms/cs_**2)/(hbs*cs_) - 3*(Es/cs_**2)*Ms/mPl2))
    rep.ident('ln l_max = ln(E M/m_Pl^2) + ln 3 (G = c = 1): the flux is proportional to ln(E M/m_Pl^2) up to an additive constant',
              sp.expand_log(sp.log(3*Es*Ms/mPl2.subs({hbs: hbs, cs_: 1, Gs: 1})), force=True) - sp.log(Es*Ms/(hbs)) - sp.log(3))

# =============================================================================================
def main():
    rep.out('Reproduction of the analytic results and quoted numbers of photon_lightring.tex (Python: sympy/mpmath/scipy)')
    rep.out(f'sympy {sp.__version__}, mpmath {mp.__version__}, numpy {np.__version__}; units G = c = M = 1, E = 1; fluxes per E^2, single m > 0, j = l - m')
    for sec in (section_II, section_IIIA, section_IIIB, section_IIIC, section_IIID, section_IV, section_V, section_VI):
        t0 = time.time(); sec(); rep.note(f'[section time {time.time() - t0:.1f} s]')
    rep.out(); rep.out('=' * 110)
    rep.out(f'SUMMARY: {rep.npass} PASS, {rep.nfail} FAIL, {rep.ninfo} info/not-checked lines; total time {time.time() - T_START:.1f} s')
    rep.out(f'maximum relative deviation among passed quantitative checks (tolerance <= 5%) against a quoted paper value: {rep.maxdev:.2e}  ({rep.maxdev_name})')
    if rep.fails:
        rep.out('FAILED CHECKS:')
        for nm, v, p, dv in rep.fails: rep.out(f'   {nm}: recomputed {v}, paper {p}, deviation {dv}')
    with open(os.path.join(HERE, 'report.txt'), 'w') as fh:
        fh.write('\n'.join(rep.lines) + '\n')
    print('report written to', os.path.join(HERE, 'report.txt'))
    return rep.nfail

if __name__ == '__main__':
    sys.exit(0 if main() == 0 else 1)
