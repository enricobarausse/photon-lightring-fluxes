"""Teukolsky (s = -2) fluxes from a photon on the equatorial light ring of a Kerr black hole.

Production code for Sec. IV and Fig. 2 of "Gravitational radiation from a photon on the
light ring" (E. Barausse): run_production.py drives it and writes kerr_results.json.

Problem and conventions
-----------------------
* Units M = 1. The photon has energy E = 1 and angular momentum L = b, with b the impact
  parameter of the circular null geodesic of radius r0(a) = 2 (1 + cos(2/3 arccos(-a)));
  fluxes are therefore per E^2, i.e. Edot M^2/E^2.
* A single m > 0 per mode, j = l - m; the frequency is omega = m Omega with Omega = 1/b.
  Prograde orbits are a > 0; retrograde orbits are represented by a -> -a with b > 0 and
  m > 0 (equivalent to a > 0, m < 0: the spheroidal harmonic is evaluated at pi/2 in that
  representation throughout). The flux of the -m mode is equal, so the physical flux per l
  is twice the sum over j.
* Source, homogeneous solutions, Wronskian and fluxes follow the conventions of the
  Teukolsky package of the Black Hole Perturbation Toolkit: R_in -> Delta^{-s} e^{-i kappa r*}
  at the horizon (unit transmission), R_up -> r^3 e^{i omega r*} at infinity (unit
  transmission), W = Delta^{s+1} (R_in R_up' - R_up R_in') = 2 i omega B_inc,
  Z = -8 pi alpha/(W Upsilon_t) with alpha the circular-orbit source projection of
  Hughes (2000) in the Toolkit's form (alpha_lm of Sec. IV of the paper), and
  Edot_I = |Z_I|^2/(4 pi omega^2), Edot_H = alpha_H |Z_H|^2/(4 pi omega^2) with alpha_H the
  Teukolsky-Press horizon factor.

Method
------
The radial Teukolsky equation is brought to the Schroedinger form Y'' + Q Y = 0 in the
tortoise coordinate r* with R = Delta (r^2 + a^2)^(-1/2) Y (Teukolsky and Press 1974; the
complex potential Q is built symbolically with sympy and compiled with lambdify). R_in is
obtained by integrating Y from r = r_+ + 1e-6 up to r0 with scipy's eighth-order Dormand-Prince
method (DOP853, rtol 1e-12, atol 1e-14), starting from the ingoing Riccati-iterated WKB
solution normalised to the Toolkit's unit transmission at the horizon. R_up is obtained in
sn.py from the Sasaki-Nakamura equation, integrated inwards from r = rB (40 or 60) with the
outgoing WKB data. The spin-weighted spheroidal harmonic, its separation constant and its
theta derivatives at the equator come from swsh.py. The source alpha and the horizon
factor alpha_H are evaluated in Mode.alpha and Mode.alphaH.

Usage
-----
    from kerrflux import Mode
    md = Mode(a, l, j, rB=60.0)      # a < 0 for retrograde orbits
    fluxI, fluxH = md.solve()        # per E^2, single m
Running this file directly checks a = 0, l <= 5 against the Toolkit (MST) values stored in
the dictionary below (agreement to ~1e-8). Cost: a second per mode at l ~ 10, of the order of
ten seconds per mode at l = 800.
"""
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp, quad

from swsh import S_at_equator

s = -2   # spin weight


def lr_data(a):
    """Circular null geodesic of Kerr (M = 1) for spin a (a < 0: retrograde).

    Returns r0, the impact parameter b = L/E, Omega = 1/b, Upsilon_t = Sigma dt/dlambda
    (per unit E), Delta and P = r^2 + a^2 - a b at r0, and a check that the radial potential
    Rtilde = P^2 - Delta (b - a)^2 and its derivative vanish at r0 (both should be ~0).
    """
    r0 = 2.0*(1.0 + np.cos(2.0/3.0*np.arccos(-a)))
    D0 = r0*r0 - 2*r0 + a*a
    P0 = 2*r0*D0/(r0 - 1)
    b = np.sqrt(27.0) if a == 0 else (r0*r0 + a*a - P0)/a
    Om = 1.0/b
    Ups = (r0*r0 + a*a)*P0/D0 + a*(b - a)
    P = lambda r: r*r + a*a - a*b
    Rt = lambda r: P(r)**2 - (r*r - 2*r + a*a)*(b - a)**2
    chk = (Rt(r0), (Rt(r0 + 1e-6) - Rt(r0 - 1e-6))/2e-6)
    return dict(a=a, r0=r0, b=b, Omega=Om, Ups=Ups, D0=D0, P0=P0, check=chk)


# --- Schroedinger-form potential Q(r) and the Riccati-WKB wavenumbers, built symbolically ---
_r, _a, _w, _m, _lam = sp.symbols('r a w m lam')
_Dl = _r**2 - 2*_r + _a**2
_K = (_r**2 + _a**2)*_w - _a*_m
_G = s*(_r - 1)/(_r**2 + _a**2) + _r*_Dl/(_r**2 + _a**2)**2
_Q = (_K**2 - 2*sp.I*s*(_r - 1)*_K + _Dl*(4*sp.I*s*_w*_r - _lam))/(_r**2 + _a**2)**2 - _G**2 - _Dl/(_r**2 + _a**2)*sp.diff(_G, _r)
_Qf = sp.lambdify((_r, _a, _w, _m, _lam), _Q, 'numpy')
_dQf = sp.lambdify((_r, _a, _w, _m, _lam), sp.diff(_Q, _r), 'numpy')
_f = _Dl/(_r**2 + _a**2)          # dr/dr*
_qwkb = {}                        # (sign, iterations) -> q(r) with Y ~ exp(i int q dr*)
for _sign in (+1, -1):
    _q = _sign*sp.sqrt(_Q)
    _qwkb[(_sign, 0)] = sp.lambdify((_r, _a, _w, _m, _lam), _q, 'numpy', cse=True)
    for _n in (1, 2):
        _q = _sign*sp.sqrt(_Q + sp.I*_f*sp.diff(_q, _r))
        _qwkb[(_sign, _n)] = sp.lambdify((_r, _a, _w, _m, _lam), _q, 'numpy', cse=True)


class Mode:
    """One (l, m = l - j) mode of the photon source at the light ring of spin a."""

    def __init__(self, a, l, j, rsA=-40.0, rB=40.0, N=None):
        self.a = a; self.l = l; self.j = j
        d = lr_data(a); self.d = d
        self.m = l - j
        self.w = self.m*d['Omega']
        self.gam = a*self.w                       # spheroidicity a omega
        self.lam, self.S0, self.dS0, self.d2S0 = S_at_equator(s, l, self.m, self.gam, N)
        self.rp = 1 + np.sqrt(1 - a*a); self.rm = 1 - np.sqrt(1 - a*a)
        self.rsA = rsA; self.rB = rB

    def rstar(self, r):
        a = self.a; rp, rm = self.rp, self.rm
        if a == 0:
            return r + 2*np.log((r - 2)/2)
        return r + 2*rp/(rp - rm)*np.log((r - rp)/2) - 2*rm/(rp - rm)*np.log((r - rm)/2)

    def f(self, r):
        """dr/dr* = Delta/(r^2 + a^2)."""
        return (r*r - 2*r + self.a**2)/(r*r + self.a**2)

    def Q(self, r):
        return _Qf(r, self.a, self.w, self.m, self.lam)

    def dQ(self, r):
        return _dQf(r, self.a, self.w, self.m, self.lam)

    def q_wkb(self, r, sign, niter=2):
        return _qwkb[(sign, niter)](r, self.a, self.w, self.m, self.lam)

    def solve(self):
        """Integrate R_in and R_up to r0 and return (Edot_I, Edot_H) per E^2 for this mode."""
        from sn import Rup_at
        a, w, m = self.a, self.w, self.m
        r0 = self.d['r0']
        rA = self.rp + 1e-6
        self.rsA = self.rstar(rA)
        rs0 = self.rstar(r0)

        def rhs(rs, y):
            Y, dY, r = y
            r = r.real
            return [dY, -self.Q(r)*Y, self.f(r)]

        qinA = self.q_wkb(rA, -1)
        kap = w - m*(a/(2*self.rp))             # omega - m Omega_H
        DA = (rA - self.rp)*(rA - self.rm)
        # horizon normalisation: R_in -> Delta^{-s} e^{-i kappa r*}, i.e.
        # Y_in -> (rp^2+a^2)^{1/2} Delta^{-s/2} e^{-i kappa r*} exp(int_{rp}^{r} g dr), with
        # g = [i (q_in + kappa)(r^2+a^2) + s (r-1)]/Delta (finite at the horizon)
        g = lambda r: (1j*(self.q_wkb(r, -1) + kap)*(r*r + a*a) + s*(r - 1))/((r - self.rp)*(r - self.rm))
        Ig = quad(lambda r: g(r).real, self.rp + 1e-13, rA, limit=200)[0] + 1j*quad(lambda r: g(r).imag, self.rp + 1e-13, rA, limit=200)[0]
        YA = np.sqrt(self.rp**2 + a*a)*DA**(-s/2)*np.exp(-1j*kap*self.rsA)*np.exp(Ig)
        self.Ig = Ig
        solIn = solve_ivp(rhs, [self.rsA, rs0], [YA, 1j*qinA*YA, rA + 0j], method='DOP853', rtol=1e-12, atol=1e-14)
        Yin0, dYin0 = solIn.y[0, -1], solIn.y[1, -1]
        # back to the Teukolsky function R = Delta (r^2+a^2)^{-1/2} Y and its r derivative
        pref = lambda r: (r*r - 2*r + a*a)*(r*r + a*a)**(-0.5)
        dpref = (pref(r0 + 1e-6) - pref(r0 - 1e-6))/2e-6
        Rin = pref(r0)*Yin0; dRin = dpref*Yin0 + np.sqrt(r0*r0 + a*a)*dYin0
        Rup, dRup = Rup_at(a, w, m, self.lam, r0, rB=self.rB)
        D0 = self.d['D0']
        W = D0**(s + 1)*(Rin*dRup - Rup*dRin)
        self.W = W; self.Binc = W/(2j*w)
        self.Rin, self.dRin, self.Rup, self.dRup = Rin, dRin, Rup, dRup
        ZI = -8*np.pi*self.alpha(Rin, dRin)/(W*self.d['Ups'])
        ZH = -8*np.pi*self.alpha(Rup, dRup)/(W*self.d['Ups'])
        self.ZI, self.ZH = ZI, ZH
        self.fluxI = abs(ZI)**2/(4*np.pi*w*w)
        self.fluxH = self.alphaH()*abs(ZH)**2/(4*np.pi*w*w)
        return self.fluxI, self.fluxH

    def alpha(self, R, dR):
        """Source projection alpha_lm (Toolkit ConvolveSourcePointParticleCircular, s = -2)
        for the null geodesic (E = 1, L = b), evaluated on the radial function R, R' at r0;
        R'' is eliminated with the Teukolsky equation."""
        a, m, w, lam = self.a, self.m, self.w, self.lam
        r0 = self.d['r0']; b = self.d['b']
        S0, dS0, d2S0 = self.S0, self.dS0, self.d2S0
        D = r0*r0 - 2*r0 + a*a
        Kt = (r0*r0 + a*a)*w - m*a
        d2R = (-(-lam + 2j*r0*s*2*w + (-2j*(-1 + r0)*s*(-a*m + (a*a + r0*r0)*w) + (-a*m + (a*a + r0*r0)*w)**2)/D)*R - (-2 + 2*r0)*(1 + s)*dR)/D
        th0 = np.pi/2
        L1 = -m/np.sin(th0) + a*w*np.sin(th0) + np.cos(th0)/np.sin(th0)
        L2 = -m/np.sin(th0) + a*w*np.sin(th0) + 2*np.cos(th0)/np.sin(th0)
        L2S = dS0 + L2*S0
        L2p = m*np.cos(th0)/np.sin(th0)**2 + a*w*np.cos(th0) - 2/np.sin(th0)**2
        L1Sp = d2S0 + L1*dS0
        L1L2S = L1Sp + L2p*S0 + L2*dS0 + L1*L2*S0
        rho = -1/(r0 - 1j*a*np.cos(th0)); rhob = -1/(r0 + 1j*a*np.cos(th0)); Sig = 1/(rho*rhob)
        Ann0 = -rho**(-2)*rhob**(-1)*(np.sqrt(2)*D)**(-2)*(rho**(-1)*L1L2S + 3j*a*np.sin(th0)*L1*S0 + 3j*a*np.cos(th0)*S0 + 2j*a*np.sin(th0)*dS0 - 1j*a*np.sin(th0)*L2*S0)
        Anmb0 = rho**(-3)*(np.sqrt(2)*D)**(-1)*((rho + rhob - 1j*Kt/D)*L2S + (rho - rhob)*a*np.sin(th0)*Kt/D*S0)
        Anmb1 = -rho**(-3)*(np.sqrt(2)*D)**(-1)*(L2S + 1j*(rho - rhob)*a*np.sin(th0)*S0)
        Ambmb0 = (Kt**2*S0*rhob)/(4*D*D*rho**3) + (1j*Kt*S0*(1 - r0 + D*rho)*rhob)/(2*D*D*rho**3) + (1j*r0*S0*rhob*w)/(2*D*rho**3)
        Ambmb1 = -rho**(-3)*rhob*S0/2*(1j*Kt/D - rho)
        Ambmb2 = -rho**(-3)*rhob*S0/4
        E = 1.0; L = b
        rcomp = (E*(r0*r0 + a*a) - a*L)/(2*Sig)
        tcomp = rho*(1j*np.sin(th0)*(a*E - L/np.sin(th0)**2))/np.sqrt(2)
        Cnn, Cnmb, Cmbmb = rcomp**2, rcomp*tcomp, tcomp**2
        return ((Ann0*Cnn + Anmb0*Cnmb + Ambmb0*Cmbmb)*R - (Anmb1*Cnmb + Ambmb1*Cmbmb)*dR + Ambmb2*Cmbmb*d2R)

    def alphaH(self):
        """Teukolsky-Press factor converting |Z_H|^2 into the horizon energy flux."""
        a, m, w, lam = self.a, self.m, self.w, self.lam
        rh = self.rp; OmH = a/(2*rh); kap = w - m*OmH; eps = np.sqrt(1 - a*a)/(4*rh)
        C2 = ((lam + 2)**2 + 4*a*m*w - 4*a*a*w*w)*(lam**2 + 36*m*a*w - 36*a*a*w*w) + (2*lam + 3)*(96*a*a*w*w - 48*m*a*w) + 144*w*w*(1 - a*a)
        return 256*(2*rh)**5*kap*(kap**2 + 4*eps**2)*(kap**2 + 16*eps**2)*w**3/C2


# Toolkit (MST, 40 digits) reference values of l*Edot_I, l*Edot_H at a = 0 (from toolkit_reference.txt)
_REFERENCE_A0 = {(2, 0): (3.72835158e-2, 1.79280159e-2), (2, 1): (9.98603250e-4, 5.84556670e-3),
                 (3, 0): (4.48227094e-2, 2.21949963e-2), (3, 1): (1.70061841e-3, 5.37761036e-3),
                 (5, 0): (4.28889133e-2, 2.31575681e-2), (5, 1): (2.33432273e-3, 4.93422191e-3)}

if __name__ == "__main__":
    for (l, j), (rI, rH) in _REFERENCE_A0.items():
        t = time.time(); md = Mode(0.0, l, j, rB=60.0); fI, fH = md.solve()
        print(f"a=0 l={l} j={j}: l*I={l*fI:.8e} (ref {rI:.8e}, {l*fI/rI-1:+.1e})  l*H={l*fH:.8e} (ref {rH:.8e}, {l*fH/rH-1:+.1e})  t={time.time()-t:.1f}s")
