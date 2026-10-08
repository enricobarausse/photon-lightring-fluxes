"""Exact Teukolsky fluxes for timelike circular orbits near the light ring versus the
1973 formulas of Breuer, Ruffini, Tiomno and Vishveshwara (BRTV), Eqs. (14)-(15).

Supports Sec. V ("Comparison with the geodesic synchrotron radiation literature") of
"Gravitational radiation from a photon on the light ring" (E. Barausse): the ratio of the
exact odd flux (l = m + 1) to the BRTV odd formula, 3.15, 3.48, 3.64, 3.77 for
m = 40, 100, 200, 400 at delta = 1e-4, tending to 4, while the even formula (l = m) is
accurate to O(m^-1/2).

Setup: Schwarzschild (a = 0), massive particle on the circular orbit r0 = 3 (1 + delta)
(delta_here = r0 - 3 = 3 M delta_paper), m = 40 ... 400, flux to infinity. kerrflux.Mode is
reused with the timelike constants of motion (E and L per unit mass mu, Omega = r0^-3/2,
Upsilon_t = r0^2 u^t/E); since the source there is per E^2 of the null geodesic, the result
is multiplied by E^2 to get the flux per mu^2 of the BRTV formulas. brtv(m, dl, parity)
is Eqs. (14)-(15) of BRTV in the paper's notation (eta = 1/2 + m dl/2, eta' = eta + 1,
in units of mu^2/M^2).

    python check_odd_BRTV.py          (about a minute; m = 400 dominates)

Writes check_odd_BRTV_results.txt next to this script (columns r0, m, parity, exact flux per mu^2,
BRTV flux, ratio) and prints the two fits of the odd ratios at delta = 1e-4 quoted in the paper:
rho_inf + rho_1 m^-1/2 + rho_2 m^-1 (rho_inf = 4.06) and rho_inf + rho_1 m^-1/2 (4.05). The stored
copy is read by the companion notebook and by sympy_check/reproduce_paper.py.
"""
import os
import numpy as np
import mpmath as mp

from kerrflux import Mode
from swsh import S_at_equator


class TimelikeMode(Mode):
    """Mode of a massive particle on the circular orbit of radius r0 in Schwarzschild."""

    def __init__(self, r0, l, j, **kw):
        super().__init__(0.0, l, j, **kw)
        f0 = 1 - 2/r0
        E = f0/np.sqrt(1 - 3/r0); L = np.sqrt(r0)/np.sqrt(1 - 3/r0)     # per unit mu
        ut = 1/np.sqrt(1 - 3/r0)
        self.d = dict(self.d); self.d.update(r0=r0, b=L/E, Omega=r0**-1.5, Ups=r0**2*ut/E, D0=r0*r0 - 2*r0)
        self.w = self.m*self.d['Omega']; self.E = E; self.gamma2 = E**2
        self.lam, self.S0, self.dS0, self.d2S0 = S_at_equator(-2, l, self.m, 0.0, None)


def brtv(m, dl, parity):
    """BRTV 1973 power in the dominant even (l = m) or odd (l = m + 1) mode, per mu^2."""
    if parity == 'even':
        eps = 1 + m*dl
        A = (1 + m*dl/2)*mp.gamma(0.25 + 0.25j*eps) + mp.sqrt(2)*(1 - 1j)/mp.sqrt(3*m)*mp.gamma(0.75 + 0.25j*eps)
    else:
        eps = 3 + m*dl
        A = mp.sqrt(2)*mp.gamma(0.75 + 0.25j*eps) + (1 + 1j)/(2*mp.sqrt(3*m))*mp.gamma(0.25 + 0.25j*eps)
    return float(mp.exp(-mp.pi*eps/4)/(54*mp.pi**1.5*m*dl)*abs(A)**2)


if __name__ == "__main__":
    rows = []
    for r0 in [3.003, 3.0003]:
        dl = r0 - 3
        for m in [40, 100, 200, 400]:
            out = []
            for j, par in [(0, 'even'), (1, 'odd')]:
                md = TimelikeMode(r0, m + j, j); fI, fH = md.solve()
                fI_mu2 = fI*md.E**2          # convert per-E^2 to per-mu^2
                pB = brtv(m, dl, par)
                rows.append((r0, m, par, fI_mu2, pB, fI_mu2/pB))
                out.append(f"{par}: exact={fI_mu2:.5e} BRTV={pB:.5e} ratio={fI_mu2/pB:.4f}")
            print(f"r0={r0} m={m} m*delta={m*dl:.3f} | " + " | ".join(out), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, 'check_odd_BRTV_results.txt'), 'w') as f:
        f.write("# exact Teukolsky flux to infinity (per mu^2) of a massive particle on r0 = 3(1 + delta_paper),\n"
                "# l = m (even) or l = m + 1 (odd), versus BRTV 1973 Eqs. (14)-(15); check_odd_BRTV.py\n"
                "# r0  m  parity  exact  BRTV  ratio\n")
        for r0, m, par, ex, pB, ra in rows:
            f.write(f"{r0} {m} {par} {ex:.10e} {pB:.10e} {ra:.6f}\n")
    # fits of the odd ratios at delta_paper = 1e-4 (r0 = 3.0003), as quoted in Sec. V of the paper
    ms = np.array([m for r0, m, par, *_ in rows if r0 == 3.0003 and par == 'odd'], float)
    rs = np.array([ra for r0, m, par, ex, pB, ra in rows if r0 == 3.0003 and par == 'odd'])
    c3 = np.linalg.lstsq(np.c_[np.ones_like(ms), ms**-0.5, 1/ms], rs, rcond=None)[0]
    c2 = np.linalg.lstsq(np.c_[np.ones_like(ms), ms**-0.5], rs, rcond=None)[0]
    print(f"odd ratios at delta = 1e-4: {', '.join(f'{r:.4f}' for r in rs)}")
    print(f"fit rho_inf + rho_1 m^-1/2 + rho_2 m^-1: rho_inf = {c3[0]:.3f}   (paper: 4.06)")
    print(f"fit rho_inf + rho_1 m^-1/2:              rho_inf = {c2[0]:.3f}   (paper: 4.05)")
    print("results written to check_odd_BRTV_results.txt")
