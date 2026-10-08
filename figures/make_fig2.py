"""Fig. 2 of "Gravitational radiation from a photon on the light ring" (E. Barausse):
Kerr fluxes per multipole (top) and the spin factor g(a) (bottom), from the stored data.

Inputs:    ../kerr_py/kerr_results.json (the run of kerr_py/run_production.py: Teukolsky fluxes per
           E^2, single m = l - j, for a = 0, +-0.5, +-0.9, l = 10 ... 800, j <= 3; M = 1) and
           ../kerr_analytic/kerr_kappa_analytic.json (the analytic kappa_j(a) of
           kerr_analytic/run_kerr.py).
Plotted:   top panel: l (Edot_I + Edot_H)/2 summed over the four dominant m (j <= 3) and over +-m
           (a factor 2), i.e. l (Edot_I^l + Edot_H^l) M^2/(2 E^2) in the notation of the paper,
           with the analytic value with the same truncation,
           2 g(a) sum_{j<=3} kappa_j = 2 sum_{j<=3} kappa_j(a), dashed. Bottom: g(a) of
           Eq. (21), with the five spins of the top panel marked.
Outputs:   fig2.pdf
Run:       python make_fig2.py      (run from the figures/ directory, or from anywhere)
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'cm'})

HERE = os.path.dirname(os.path.abspath(__file__))
NUMERICS = os.path.join(HERE, '..', 'kerr_py', 'kerr_results.json')
ANALYTIC = os.path.join(HERE, '..', 'kerr_analytic', 'kerr_kappa_analytic.json')
# (|a|, sign, colour, label); negative sign = retrograde
SPINS = [(0.9, 1, 'tab:red', r'$a/M=0.9$ prograde'), (0.5, 1, 'tab:orange', r'$a/M=0.5$ prograde'), (0.0, 1, 'k', r'$a=0$'),
         (0.5, -1, 'tab:blue', r'$a/M=0.5$ retrograde'), (0.9, -1, 'tab:green', r'$a/M=0.9$ retrograde')]
ELLS = [10, 20, 50, 100, 200, 400, 800]
JMAX = 3


def r0(a):
    """Light-ring radius, M = 1 (a < 0: retrograde)."""
    return 2*(1 + np.cos(2/3*np.arccos(-a)))


def g_closed(a):
    """Spin factor g(a) of Eq. (21)."""
    r = r0(a)
    return 27*(r - 1)/(r**2*(r + 3))


def make_figure(outfile=os.path.join(HERE, 'fig2.pdf')):
    res = json.load(open(NUMERICS))
    ana = json.load(open(ANALYTIC))
    kappa = {(float(row['a']), row['j']): row['kappaI'] for row in ana}   # (a, j) -> kappa_j(a), single m

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.4, 4.9), gridspec_kw={'height_ratios': [1.15, 1], 'hspace': 0.42})
    for a, sign, col, lab in SPINS:
        y = []
        for l in ELLS:
            rows = [r for r in res if r['a'] == a and r['sign'] == sign and r['l'] == l and r['j'] <= JMAX]
            y.append(l*sum(r['FluxI'] + r['FluxH'] for r in rows))       # = 2 l (I + H)/2 summed over j
        ax1.plot(ELLS, y, '-o', color=col, ms=3.5, lw=1, label=lab)
        kA = 2*sum(kappa[(sign*a, j)] for j in range(JMAX + 1))          # 2 g(a) sum_{j<=3} kappa_j
        ax1.plot([8, 1200], [kA, kA], '--', color=col, lw=0.9)
    ax1.set_xscale('log'); ax1.set_xlim(8, 1200); ax1.set_ylim(0.04, 0.10)
    ax1.set_xlabel(r'$\ell$'); ax1.set_ylabel(r'$\ell\,(\dot E^\infty_\ell+\dot E^H_\ell)\,M^2/(2E^2)$')
    ax1.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=7, frameon=False, handlelength=1.6, columnspacing=1.0)

    agrid = np.linspace(-0.999, 0.9995, 2000)
    ax2.plot(agrid, g_closed(agrid), 'k-', lw=1)
    ax2.axhline(1, color='gray', lw=0.6, ls=':')
    for a, sign, col, lab in SPINS:
        ax2.plot([sign*a], [g_closed(sign*a)], 'o', color=col, ms=4)
    ax2.set_xlim(-1, 1); ax2.set_ylim(0.55, 1.45)
    ax2.set_xlabel(r'$a/M$'); ax2.set_ylabel(r'$g(a)$')
    fig.savefig(outfile, bbox_inches='tight')
    print(outfile, 'written')


if __name__ == "__main__":
    make_figure()
