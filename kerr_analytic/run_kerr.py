"""Driver for alpha_expand.py: the analytic Kerr coefficients of Sec. IV and Table I.

    python run_kerr.py                 # kerr_kappa_analytic.json: a = 0, +-0.5, +-0.9, +-0.99, j = 0..3
    python run_kerr.py --spins 0.5     # one spin, printed only (nothing written)
    python run_kerr.py --grids         # also kerr_g_finegrid.json and kerr_g_curve.json

Outputs (all in this directory; a < 0 denotes retrograde orbits):
* kerr_kappa_analytic.json: a list of records {a, j, Ahat, Bhat, kappaI, ratioHI, a_j, J0, JH,
  hinf, ktilde, imJ} with kappaI = kappa_j(a) (flux to infinity, single m, per E^2/(M^2 l)),
  Ahat and Bhat the source coefficients as [Re, Im], ratioHI the horizon/infinity ratio
  (equal to 1), and the amplitude integrals. Read by figures/make_fig2.py (the dashed lines
  2 g(a) sum_{j<=3} kappa_j of Fig. 2) and used for Table I.
* kerr_g_finegrid.json: for a = 0.7, 0.8, 0.95, 0.97, 0.98, 0.995, the j = 0 ingredients
  {kappa0, g, N0, ktilde, Ahat, r0} (the values quoted in the text after Eq. (gclosed)).
* kerr_g_curve.json: [a, g(a)] on a grid of spins, from kappa_0(a)/kappa_0(0).

Cost: about a second per (a, j).
"""
import argparse
import json
import os

import mpmath as mp

from alpha_expand import kappa_kerr, lr_data

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SPINS = ['0', '0.5', '0.9', '0.99', '-0.5', '-0.9', '-0.99']
FINEGRID_SPINS = ['0.7', '0.8', '0.95', '0.97', '0.98', '0.995']
GCURVE_SPINS = ([-0.99, -0.96, -0.93] + [x/100 for x in range(-90, 90, 5)] + [0.9]
                + [x/1000 for x in range(905, 1000, 5)])


def kappa_table(spins, verbose=True):
    """kappa_j(a) for j = 0..3 at each spin; returns the JSON records."""
    out = []
    for a in spins:
        tot = mp.mpf(0); rows = []
        for j in range(0, 4):
            r = kappa_kerr(mp.mpf(a), j, verbose=(verbose and j == 0))
            tot += r['kappaI']
            rows.append(dict(a=a, j=j, Ahat=[r['Ahat'].real, r['Ahat'].imag], Bhat=[r['Bhat'].real, r['Bhat'].imag],
                             kappaI=float(mp.re(r['kappaI'])), ratioHI=float(mp.re(r['ratioHI'])), a_j=float(mp.re(r['a_j'])),
                             J0=float(mp.re(r['J0'])), JH=float(mp.re(r['JH'])), hinf=float(mp.re(r['hinf'])),
                             ktilde=float(mp.re(r['ktilde'])), imJ=[float(mp.im(r['J0'])), float(mp.im(r['JH']))]))
            print(f"a={a} j={j}: Ahat={r['Ahat']:.8g} Bhat={r['Bhat']:.8g}  a_j={mp.nstr(r['a_j'],8)}  kappa_I={mp.nstr(r['kappaI'],10)}  H/I={mp.nstr(r['ratioHI'],10)}  ImJ0,ImJH={mp.nstr(mp.im(r['J0']),3)},{mp.nstr(mp.im(r['JH']),3)}")
        d = lr_data(mp.mpf(a))
        print(f"==> a={a}: r0={mp.nstr(d['r0'],8)} b={mp.nstr(d['b'],8)} Omega={mp.nstr(d['Omega'],8)} Ups={mp.nstr(d['Ups'],8)}  sum_j kappa_j (single m) = {mp.nstr(tot,8)}  physical (x2) = {mp.nstr(2*tot,8)}")
        out += rows
    return out


def g_grids():
    """The j = 0 ingredients on the fine grid of prograde spins and g(a) on the wide grid."""
    re = lambda x: float(mp.re(x))      # the mpmath results carry a zero imaginary part for a != 0
    k00 = re(kappa_kerr(mp.mpf('0'), 0, verbose=False)['kappaI'])
    fine = {}
    for a in FINEGRID_SPINS:
        r = kappa_kerr(mp.mpf(a), 0, verbose=False)
        fine[a] = dict(kappa0=re(r['kappaI']), g=re(r['kappaI'])/k00, N0=re(r['N0']), ktilde=re(r['ktilde']),
                       Ahat=abs(r['Ahat']), r0=re(r['data']['r0']))
        print(f"a={a}: " + "  ".join(f"{k}={v:.6g}" for k, v in fine[a].items()))
    curve = []
    for a in GCURVE_SPINS:
        r = kappa_kerr(mp.mpf(str(a)), 0, verbose=False)     # decimal string: exact spin at 40 digits
        curve.append([a, re(r['kappaI'])/k00])
    return fine, curve


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--spins', nargs='+', default=None, help='spins as strings (negative = retrograde); default: the seven of Table I, written to kerr_kappa_analytic.json')
    p.add_argument('--grids', action='store_true', help='also regenerate kerr_g_finegrid.json and kerr_g_curve.json')
    args = p.parse_args()
    if args.spins is None:
        out = kappa_table(DEFAULT_SPINS)
        json.dump(out, open(os.path.join(HERE, 'kerr_kappa_analytic.json'), 'w'), indent=1)
    else:
        kappa_table(args.spins)
    if args.grids:
        fine, curve = g_grids()
        json.dump(fine, open(os.path.join(HERE, 'kerr_g_finegrid.json'), 'w'), indent=1)
        json.dump(curve, open(os.path.join(HERE, 'kerr_g_curve.json'), 'w'))
