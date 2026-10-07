"""Sasaki-Nakamura "up" solution of the s = -2 Teukolsky equation (Kerr, M = 1).

Used by kerrflux.py to obtain R_up and dR_up/dr at the light ring: the Teukolsky equation
is not integrated directly from large r because its outgoing solution grows like r^3 and
the ingoing one is contaminated; the Sasaki-Nakamura (SN) transformation gives a
short-range equation whose two solutions have unit modulus at infinity.

Conventions: those of the Toolkit's Teukolsky/Kernel/SasakiNakamura.m. The SN function X
satisfies X_** - F X_* - U X = 0 (derivatives wrt r*). Z = X exp(-int F/2 dr*) satisfies
Z'' + Qsn Z = 0 with Qsn = -U - F^2/4 - F_*/2, which is what the Riccati-WKB boundary
data are built for. R_up, normalised to unit transmission (R -> r^3 e^{i w r*} at infinity),
follows from X_up -> -(c0/(4 w^2)) e^{i w r*}, with c0 the leading coefficient of the SN
"eta" function. All expressions are built with sympy and compiled with lambdify.

Rup_at(a, w, m, lam, r0, rB) integrates X inwards from r = rB (outgoing WKB data, with the
phase and the F/2 amplitude integrals taken from infinity to rB) with scipy's DOP853
(rtol 1e-12, atol 1e-14) and converts X, X', X'' at r0 to R_up and R_up'.
Running this file directly checks the asymptotic normalisation at a = 0.
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

r, a, w, m, lam = sp.symbols('r a w m lam')
M = 1
Dl = r**2 - 2*M*r + a**2
f = Dl/(r**2 + a**2)
# Sasaki-Nakamura eta function and potential (s = -2)
c0 = -12*sp.I*w*M + lam*(lam + 2) - 12*a*w*(a*w - m)
c1 = 8*sp.I*a*(3*a*w - lam*(a*w - m))
c2 = -24*sp.I*a*M*(a*w - m) + 12*a**2*(1 - 2*(a*w - m)**2)
c3 = 24*sp.I*a**3*(a*w - m) - 24*M*a**2
c4 = 12*a**4
eta = c0 + c1/r + c2/r**2 + c3/r**3 + c4/r**4
K = (r**2 + a**2)*w - m*a
V = -((K**2 + 4*sp.I*(r - M)*K)/Dl) + 8*sp.I*w*r + lam
beta = 2*Dl*(-sp.I*K + r - M - 2*Dl/r)
alpha = -sp.I*K*beta/Dl**2 + 3*sp.I*sp.diff(K, r) + lam + 6*Dl/r**2
U1 = V + Dl**2/beta*(sp.diff(2*alpha + sp.diff(beta, r)/Dl, r) - sp.diff(eta, r)/eta*(alpha + sp.diff(beta, r)/Dl))
G = -((2*(r - M))/(r**2 + a**2)) + (r*Dl)/(r**2 + a**2)**2
F = sp.diff(eta, r)/eta*Dl/(r**2 + a**2)
U = ((Dl*U1)/(r**2 + a**2)**2 + G**2 + (Dl*sp.diff(G, r))/(r**2 + a**2) - F*G)
Qsn = -U - F**2/4 - f*sp.diff(F, r)/2
args = (r, a, w, m, lam)
Ff = sp.lambdify(args, F, 'numpy', cse=True)
Uf = sp.lambdify(args, U, 'numpy', cse=True)
Qf = sp.lambdify(args, Qsn, 'numpy', cse=True)
# Riccati-iterated WKB wavenumber for Z (outgoing): q_{n+1} = sqrt(Qsn + i f q_n')
qw = {}
_q = sp.sqrt(Qsn)
qw[0] = sp.lambdify(args, _q, 'numpy', cse=True)
for n in (1, 2):
    _q = sp.sqrt(Qsn + sp.I*f*sp.diff(_q, r))
    qw[n] = sp.lambdify(args, _q, 'numpy', cse=True)
# inverse transformation X -> R (s = -2): R = [(alpha + beta'/Delta) chi - (beta/Delta) chi']/eta,
# chi = X Delta/sqrt(r^2 + a^2); X, dX, d2X are X and its r derivatives at the evaluation point
X, dX, d2X = sp.symbols('X dX d2X')
Xr = sp.Function('Xr')(r)
chi_r = Xr*Dl/sp.sqrt(r**2 + a**2)
Rup_expr = 1/eta*((alpha + sp.diff(beta, r)/Dl)*chi_r - beta/Dl*sp.diff(chi_r, r))
dRup_expr = sp.diff(Rup_expr, r)
subs = {sp.Derivative(Xr, (r, 2)): d2X, sp.Derivative(Xr, r): dX, Xr: X}
Rup_f = sp.lambdify(args + (X, dX, d2X), Rup_expr.subs(subs), 'numpy', cse=True)
dRup_f = sp.lambdify(args + (X, dX, d2X), dRup_expr.subs(subs), 'numpy', cse=True)
c0f = sp.lambdify((a, w, m, lam), c0, 'numpy')


def rstar_fn(a_):
    """Tortoise coordinate r*(r) for spin a_."""
    rp = 1 + np.sqrt(1 - a_*a_); rm = 1 - np.sqrt(1 - a_*a_)
    if a_ == 0:
        return lambda rr: rr + 2*np.log((rr - 2)/2)
    return lambda rr: rr + 2*rp/(rp - rm)*np.log((rr - rp)/2) - 2*rm/(rp - rm)*np.log((rr - rm)/2)


def Rup_at(a_, w_, m_, lam_, r0, rB=40.0, niter=2):
    """R_up and dR_up/dr at r0, normalised so that R_up -> r^3 e^{i w r*} at infinity."""
    rs = rstar_fn(a_)
    fr = lambda rr: (rr*rr - 2*rr + a_*a_)/(rr*rr + a_*a_)
    q = lambda rr: qw[niter](rr, a_, w_, m_, lam_)
    # Z normalisation: Z -> e^{i w r*}; Z(rB) = e^{i w r*_B} exp(i int_inf^rB (q - w) dr*)
    g = lambda rr: 1j*(q(rr) - w_)/fr(rr)
    Iq = quad(lambda rr: g(rr).real, rB, np.inf, limit=400)[0] + 1j*quad(lambda rr: g(rr).imag, rB, np.inf, limit=400)[0]
    # X = exp(int_inf^r F/2 dr*) Z
    hF = lambda rr: 0.5*Ff(rr, a_, w_, m_, lam_)/fr(rr)
    IF = quad(lambda rr: hF(rr).real, rB, np.inf, limit=400)[0] + 1j*quad(lambda rr: hF(rr).imag, rB, np.inf, limit=400)[0]
    ZB = np.exp(1j*w_*rs(rB))*np.exp(-Iq)
    dZB = 1j*q(rB)*ZB
    XB = np.exp(-IF)*ZB
    dXB = np.exp(-IF)*(dZB + 0.5*Ff(rB, a_, w_, m_, lam_)*ZB)   # X_* = (F/2) X + e^{..} Z_*

    def rhs(rst, y):
        Xv, dXv, rr = y
        rr = rr.real
        return [dXv, Ff(rr, a_, w_, m_, lam_)*dXv + Uf(rr, a_, w_, m_, lam_)*Xv, fr(rr)]

    sol = solve_ivp(rhs, [rs(rB), rs(r0)], [XB, dXB, rB + 0j], method='DOP853', rtol=1e-12, atol=1e-14)
    X0, dXs0 = sol.y[0, -1], sol.y[1, -1]
    dX0 = dXs0/fr(r0)                      # d/dr
    # d2X/dr2 from the SN equation in r: f^2 X'' + f(f' - F) X' - U X = 0
    fp = (fr(r0 + 1e-6) - fr(r0 - 1e-6))/2e-6
    d2X0 = (Uf(r0, a_, w_, m_, lam_)*X0 - fr(r0)*(fp - Ff(r0, a_, w_, m_, lam_))*dX0)/fr(r0)**2
    scale = -c0f(a_, w_, m_, lam_)/(4*w_**2)
    Rv = Rup_f(r0, a_, w_, m_, lam_, X0, dX0, d2X0)*scale
    dRv = dRup_f(r0, a_, w_, m_, lam_, X0, dX0, d2X0)*scale
    return Rv, dRv


if __name__ == "__main__":
    # check of the asymptotic normalisation: R_up/(r^3 e^{iwr*}) -> 1 at large r (a = 0, l = m = 2)
    a_, w_, m_, lam_ = 0.0, 2/np.sqrt(27), 2, 4.0
    rs = rstar_fn(a_)
    for rr in [200.0, 400.0, 800.0]:
        Rv, dRv = Rup_at(a_, w_, m_, lam_, rr, rB=1600.0)
        print(rr, Rv/(rr**3*np.exp(1j*w_*rs(rr))))
