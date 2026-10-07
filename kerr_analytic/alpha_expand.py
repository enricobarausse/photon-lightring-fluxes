"""Analytic large-l coefficients kappa_j(a) for a photon on the Kerr light ring.

Supports Sec. IV of "Gravitational radiation from a photon on the light ring" (E. Barausse):
Eq. (kappakerr), its odd counterpart, the spin factor g(a) = kappa_0(a)/kappa_0(0) and the
ingredients listed in Table I (r0, b, Upsilon_t, ktilde, |Ahat_0|, N, kappa_0, g). The driver
run_kerr.py writes kerr_kappa_analytic.json (and the g(a) grids of the other JSON files).

What is computed (M = 1, E = 1, L = b, single m = l - j > 0; a < 0 = retrograde orbit)
---------------------------------------------------------------------------------------
* lr_data(a): the circular null geodesic (r0, b, Omega, Upsilon_t, beta_b = sqrt(b^2 - a^2)).
* alpha_expand(a, j): the large-omega expansion of the Toolkit's s = -2 circular-orbit source
  projection alpha (ConvolveSourcePointParticleCircular, theta0 = pi/2) for the null geodesic,
  with the wave at the light ring, Y(0) = YY0, and its r* derivative, omega^(1/2) yy1, kept as
  unknowns and the spheroidal harmonic at the equator reduced to one symbol SS through the
  large-l relations S' = sig1 S (even j) or S = sig0 S'/omega (odd j) and
  S'' = (m^2 - s - A) S. The numeric data (a, r0, b, ...) are mpmath 40-digit floats, omega is
  a sympy symbol, and the result is a dictionary {power of omega: coefficient}. The omega^2
  and omega^(3/2) terms vanish for a null source (the cancellations of Sec. IV), leaving
  alpha = E^2 Ahat_j omega S0 Y(0) + ... for even j and alpha = E^2 Bhat omega^(1/2) S0' y1
  + ... for odd j.
* J_inf, h_inf, J_hor: the WKB amplitude integrals of Eq. (N) and of the horizon side;
  ktilde(d): the barrier curvature of the real part of the potential.
* kappa_kerr(a, j): assembles kappa_j(a) for the flux to infinity, Eq. (kappakerr) (even j)
  and its odd analogue, the ratio of the horizon to the infinity normalisation (ratioHI,
  which equals 1), the coefficient a_j = (2j+1) beta_b Delta_0/((r0^2+a^2)^2 sqrt(2 ktilde))
  and N_0 (the amplitude factor N of the paper). F(x) is the barrier-top function of
  Sec. III C, cY and cdY the harmonic coefficients c_j and d_j.

Running this file directly reproduces the Schwarzschild kappa_j of Sec. III C for j <= 3
(a = 0: 0.0277668, 0.0037670, 2.7521e-4, 1.6604e-5). Cost: about a second per (a, j).
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40
s = -2


def lr_data(a):
    """Circular null geodesic of Kerr for spin a (mpmath): r0, b, Omega, Upsilon_t, Delta_0, P_0, beta_b."""
    a = mp.mpf(a)
    r0 = 2*(1 + mp.cos(mp.mpf(2)/3*mp.acos(-a)))
    D0 = r0**2 - 2*r0 + a**2
    P0 = 2*r0*D0/(r0 - 1)
    b = mp.sqrt(27) if a == 0 else (r0**2 + a**2 - P0)/a
    Om = 1/b
    Ups = (r0**2 + a**2)*P0/D0 + a*(b - a)
    betab = mp.sqrt(b**2 - a**2)
    # the radial potential and its derivative must vanish at r0
    P = lambda r: (r**2 + a**2) - a*b
    Rt = lambda r: P(r)**2 - (r**2 - 2*r + a**2)*(b - a)**2
    assert abs(Rt(r0)) < mp.mpf(10)**(-30) and abs(mp.diff(Rt, r0)) < mp.mpf(10)**(-25)
    return dict(a=a, r0=r0, b=b, Omega=Om, Ups=Ups, D0=D0, P0=P0, betab=betab)


def alpha_expand(a, j):
    """Expansion of the source projection alpha in powers of omega^(1/2); see the module docstring."""
    d = lr_data(a)
    a_, r0, b, bb, D0, P0 = [sp.Float(str(d[k]), 40) for k in ('a', 'r0', 'b', 'betab', 'D0', 'P0')]
    w, SS, YY0, yy1, lam0 = sp.symbols('w SS YY0 yy1 lam0')
    I = sp.I
    m = w*b; c = a_*w
    lam = (m - c)**2 + (2*j + 1)*w*bb + lam0      # separation constant at large l
    A = lam - c**2 + 2*m*c
    Kt = w*P0
    L = -m + c
    if j % 2 == 0:
        sig1 = (2*j + 1)*s*(a_ + b)/bb
        S0 = SS; dS0 = sig1*SS
    else:
        sig0 = -s*(a_ + b)/bb**2
        dS0 = SS; S0 = sig0*SS/w
    d2S0 = (m**2 - s - A)*S0
    L2S = dS0 + L*S0
    L1L2S = d2S0 + 2*L*dS0 - 2*S0 + L**2*S0
    rho = -1/r0; rhob = -1/r0; Sig = r0**2
    Ann0 = -rho**(-2)*rhob**(-1)*(sp.sqrt(2)*D0)**(-2)*(rho**(-1)*L1L2S + 2*I*a_*L2S)
    Anmb0 = rho**(-3)*(sp.sqrt(2)*D0)**(-1)*((2*rho - I*Kt/D0)*L2S)
    Anmb1 = -rho**(-3)*(sp.sqrt(2)*D0)**(-1)*L2S
    Ambmb0 = (Kt**2*S0*rhob)/(4*D0**2*rho**3) + (I*Kt*S0*(1 - r0 + D0*rho)*rhob)/(2*D0**2*rho**3) + (I*r0*S0*rhob*w)/(2*D0*rho**3)
    Ambmb1 = -rho**(-3)*rhob*S0/2*(I*Kt/D0 - rho)
    Ambmb2 = -rho**(-3)*rhob*S0/4
    rc = P0/(2*Sig); tc = I*(b - a_)/(sp.sqrt(2)*r0)
    Cnn, Cnmb, Cmbmb = rc**2, rc*tc, tc**2
    r = sp.symbols('r')
    pref = (r**2 - 2*r + a_**2)*(r**2 + a_**2)**sp.Rational(-1, 2)     # R = pref Y
    dpref = sp.diff(pref, r)
    pref0 = pref.subs(r, r0); dpref0 = dpref.subs(r, r0)
    R = pref0*YY0
    dR = dpref0*YY0 + pref0*((r0**2 + a_**2)/D0)*sp.sqrt(w)*yy1
    d2R = (-(-lam + 2*I*r0*s*2*w + (-2*I*(-1 + r0)*s*(-a_*m + (a_**2 + r0**2)*w) + (-a_*m + (a_**2 + r0**2)*w)**2)/(a_**2 - 2*r0 + r0**2))*R
           - (-2 + 2*r0)*(1 + s)*dR)/(a_**2 - 2*r0 + r0**2)
    alpha = (Ann0*Cnn + Anmb0*Cnmb + Ambmb0*Cmbmb)*R - (Anmb1*Cnmb + Ambmb1*Cmbmb)*dR + Ambmb2*Cmbmb*d2R
    # expand in powers of sqrt(omega): substitute omega = u^2
    u = sp.symbols('u', positive=True)
    expr = sp.expand(alpha.subs(w, u**2)*u**8)
    poly = sp.Poly(expr, u)
    coeffs = {}
    for (k,), cf in poly.terms():
        coeffs[sp.Rational(k - 8, 2)] = sp.expand(cf)
    return d, coeffs


def chop(expr, tol=1e-25):
    """Drop numerical noise below tol from a sympy expression."""
    return expr.xreplace({n: 0 for n in expr.atoms(sp.Float) if abs(n) < tol})


# --- WKB amplitude integrals (Eq. (N) of the paper) and barrier curvature ---
def qdrs(a, b, rr):
    """q dr*/dr: imaginary part of sqrt(Q) along the resonant mode, Eq. (N)."""
    Dl = rr**2 - 2*rr + a**2; P = (rr**2 + a**2) - a*b; Rt = P**2 - Dl*(b - a)**2
    return s*(-2*(rr - 1)*P + 4*rr*Dl)/(2*Dl*mp.sqrt(Rt))


def drsdr(a, rr):
    return (rr**2 + a**2)/(rr**2 - 2*rr + a**2)


def J_inf(d):
    a, b, r0 = d['a'], d['b'], d['r0']
    f = lambda rr: qdrs(a, b, rr) - s/rr*drsdr(a, rr)
    return mp.quad(f, [r0, r0 + 1, r0 + 10, mp.inf])


def h_inf(d):
    a, r0 = d['a'], d['r0']
    f = lambda rr: (rr**2 + a**2)/(rr*(rr**2 - 2*rr + a**2)) - 1/rr
    return mp.quad(f, [r0, r0 + 10, mp.inf])


def J_hor(d):
    a, b, r0 = d['a'], d['b'], d['r0']
    rp = 1 + mp.sqrt(1 - a**2)
    f = lambda rr: qdrs(a, b, rr) + s*(rr - 1)/(rr**2 + a**2)*drsdr(a, rr)
    return mp.quad(f, [rp, rp + (r0 - rp)/10, r0])


def ktilde(d):
    """Barrier curvature ktilde = d^2/dr*^2 [Rtilde/(r^2+a^2)^2] at r0 (k = omega^2 ktilde)."""
    a, b, r0 = d['a'], d['b'], d['r0']
    Dl = lambda rr: rr**2 - 2*rr + a**2
    Rt = lambda rr: ((rr**2 + a**2) - a*b)**2 - Dl(rr)*(b - a)**2
    g = lambda rr: Rt(rr)/(rr**2 + a**2)**2
    Dst = lambda fn: (lambda rr: Dl(rr)/(rr**2 + a**2)*mp.diff(fn, rr))
    return Dst(Dst(g))(r0)


def F(x):
    """Barrier-top function F(eta) = e^{-pi eta/2}/(cosh(pi eta) |Gamma(3/4 + i eta/2)|^2)."""
    return mp.exp(-mp.pi*x/2)/(mp.cosh(mp.pi*x)*abs(mp.gamma(mp.mpf(3)/4 + 1j*x/2))**2)


def cY(j):
    """c_j = binomial(j, j/2) 2^-j (even j)."""
    return mp.binomial(j, j//2)/mp.mpf(2)**j


def cdY(j):
    """d_j = j binomial(j-1, (j-1)/2) 2^(1-j) (odd j)."""
    return j*mp.binomial(j - 1, (j - 1)//2)/mp.mpf(2)**(j - 1)


def kappa_kerr(a, j, verbose=True):
    """kappa_j(a) and its ingredients; see the module docstring. Returns a dictionary with
    Ahat, Bhat (source coefficients), J0, JH, hinf (amplitude integrals), N0, NH, ktilde,
    kappaI (the coefficient of E^2/(M^2 l) for the flux to infinity), ratioHI, a_j, data."""
    d, co = alpha_expand(a, j)
    SS, YY0, yy1 = sp.symbols('SS YY0 yy1')

    def coef(p, sym1, sym2):
        e = co.get(sp.Rational(p), sp.Integer(0))
        return complex(chop(sp.expand(e)).coeff(sym1).coeff(sym2))

    lead2 = chop(co.get(sp.Rational(2), sp.Integer(0))); lead32 = chop(co.get(sp.Rational(3, 2), sp.Integer(0)))
    Ahat = coef(1, SS, YY0); Bhat = coef(sp.Rational(1, 2), SS, yy1)
    if verbose:
        print(f"a={a} j={j}: omega^2 coeff = {lead2},  omega^(3/2) coeff = {lead32}")
        print(f"   omega^1 coeff = {chop(co.get(sp.Rational(1), 0))}")
        print(f"   omega^(1/2) coeff = {chop(co.get(sp.Rational(1,2), 0))}")
    a_, r0, b, bb, Ups, D0 = d['a'], d['r0'], d['b'], d['betab'], d['Ups'], d['D0']
    kt = ktilde(d); J0 = J_inf(d); JH = J_hor(d); hinf = h_inf(d)
    N0 = r0**(2*s)*mp.exp(-2*s*hinf)*mp.exp(-2*J0)                 # amplitude factor N, Eq. (N)
    rp = 1 + mp.sqrt(1 - a_**2); OmH = a_/(2*rp); kap = 1 - b*OmH
    NH = (2*rp)*D0**s*mp.exp(-2*JH)/((2*rp)**2*kap**2)
    alphaH = 256*(2*rp)**5*kap*kap**2*kap**2/(b - a_)**8           # Teukolsky-Press factor at large l
    ratioHI = alphaH*kap*NH/N0
    x = j + mp.mpf(1)/2                                              # eta_j
    if j % 2 == 0:
        kI = mp.sqrt(2*mp.pi)*b*cY(j)*mp.sqrt(bb)*N0*F(x)*abs(Ahat)**2/((2*kt)**mp.mpf(0.25)*Ups**2)
    else:
        Du = (2*kt)**mp.mpf(0.25)*N0*mp.sqrt(2)*mp.pi*mp.exp(-mp.pi*x/2)/(mp.cosh(mp.pi*x)*abs(mp.gamma(mp.mpf(1)/4 + 1j*x/2))**2)
        kI = 4*mp.pi*abs(Bhat)**2*cdY(j)*bb**mp.mpf(1.5)/mp.pi**mp.mpf(1.5)*Du*b/Ups**2
    aj = (2*j + 1)*bb*D0/((r0**2 + a_**2)**2*mp.sqrt(2*kt))
    return dict(a=a, j=j, Ahat=Ahat, Bhat=Bhat, J0=J0, JH=JH, hinf=hinf, N0=N0, NH=NH, ktilde=kt, kappaI=kI, ratioHI=ratioHI, a_j=aj, data=d)


if __name__ == '__main__':
    # Schwarzschild check: the Zerilli / Regge-Wheeler kappa_j of Sec. III C
    ref = {0: 0.0277667579, 1: 0.0037669678, 2: 0.00027521046, 3: 1.660379e-5}
    for j in [0, 1, 2, 3]:
        res = kappa_kerr(0, j)
        print(f"a=0 j={j}: Ahat={res['Ahat']:.8g} Bhat={res['Bhat']:.8g} J0={mp.nstr(res['J0'],10)} JH={mp.nstr(res['JH'],10)} hinf={mp.nstr(res['hinf'],10)} a_j={mp.nstr(res['a_j'],10)}")
        print(f"   kappa_I = {mp.nstr(res['kappaI'],10)}   (Zerilli {ref[j]})   ratio H/I = {mp.nstr(res['ratioHI'],10)}")
