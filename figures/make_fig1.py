"""Fig. 1 of "Gravitational radiation from a photon on the light ring" (E. Barausse):
flux per multipole from a photon on the Schwarzschild light ring, from the stored data.

Input: ../photon_big_results.m (the run of schwarzschild/run_big.wls), rows
{l, j, Edot_I, Edot_H, u(0), u'(0)} for a single m = l - j, per E^2 (M = 1).
Plotted: l times the flux summed over the five dominant m (j <= 4) and over +-m
(a factor 2), to infinity (circles), into the horizon (squares) and their mean (diamonds);
the analytic coefficient kappa = 2 sum_j kappa_j = 0.063653 of Eq. (kappa) (dashed),
computed here from Eqs. (keven)-(kodd); and kappa (1 +- sigmabar/sqrt(l)) with
sigmabar = 0.610, the m-summed O(l^-1/2) asymmetry of Sec. III C (dotted).

    python make_fig1.py      (run from the figures/ directory, or from anywhere) -> fig1.pdf
"""
import os
import re

import numpy as np
from scipy.special import gamma, binom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 9, 'font.family': 'serif', 'mathtext.fontset': 'cm'})

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'photon_big_results.m')
SIGMA_BAR = 0.610       # m-summed asymmetry coefficient, Sec. III C
JMAX = 4                # the five dominant m


def read_results(path=DATA):
    """Parse the Mathematica list {{l, j, FluxI, FluxH, u0, Du0}, ...}; returns {l: {j: (FluxI, FluxH)}}."""
    text = open(path).read().replace('\n', '').replace(' ', '').replace('*^', 'e')
    data = {}
    for l, j, fI, fH in re.findall(r'\{(\d+),(\d+),([-+0-9.e]+),([-+0-9.e]+),', text):
        data.setdefault(int(l), {})[int(j)] = (float(fI), float(fH))
    return data


def kappa_j(j):
    """Analytic kappa_j of Eqs. (keven)-(kodd), eta_j = j + 1/2."""
    eta = j + 0.5
    if j % 2 == 0:
        cj = binom(j, j/2)*2.0**(-j)
        return np.sqrt(np.pi)/9*(2*j + 1)**2*cj*np.exp(-np.pi*eta/2)/(np.cosh(np.pi*eta)*abs(gamma(0.75 + 0.5j*eta))**2)
    dj = j*binom(j - 1, (j - 1)/2)*2.0**(1 - j)
    return 8*np.sqrt(np.pi)/9*dj*np.exp(-np.pi*eta/2)/(np.cosh(np.pi*eta)*abs(gamma(0.25 + 0.5j*eta))**2)


def make_figure(outfile=os.path.join(HERE, 'fig1.pdf')):
    data = read_results()
    ells = np.array(sorted(data))
    yI = np.array([2*l*sum(data[l][j][0] for j in range(JMAX + 1)) for l in ells])
    yH = np.array([2*l*sum(data[l][j][1] for j in range(JMAX + 1)) for l in ells])
    kappa = 2*sum(kappa_j(j) for j in range(0, 13))
    print(f"kappa = 2 sum_j kappa_j = {kappa:.6f}")

    fig, ax = plt.subplots(figsize=(3.4, 2.55))
    lgrid = np.logspace(np.log10(8), np.log10(2e4), 300)
    ax.plot(lgrid, kappa*(1 + SIGMA_BAR/np.sqrt(lgrid)), ':', color='red', lw=1)
    ax.plot(lgrid, kappa*(1 - SIGMA_BAR/np.sqrt(lgrid)), ':', color='blue', lw=1)
    ax.axhline(kappa, color='k', ls='--', lw=0.8)
    ax.plot(ells, yI, 'o', color='red', ms=4, label='infinity')
    ax.plot(ells, yH, 's', color='blue', ms=4, label='horizon')
    ax.plot(ells, (yI + yH)/2, 'D', color='k', ms=3.5, label='mean')
    ax.set_xscale('log'); ax.set_xlim(8, 2e4); ax.set_ylim(0.045, 0.085)
    ax.set_xlabel(r'$\ell$'); ax.set_ylabel(r'$\ell\,\dot E_\ell\,M^2/E^2$')
    ax.legend(loc='upper right', frameon=False, handlelength=1.0)
    fig.savefig(outfile, bbox_inches='tight')
    print(outfile, 'written')


if __name__ == "__main__":
    make_figure()
