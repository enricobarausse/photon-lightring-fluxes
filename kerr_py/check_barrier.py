"""Shape of the real part of the Kerr potential along the resonant modes: one minimum, two turning points.

Paper:     Sec. IV of "Gravitational radiation from a photon on the light ring" (E. Barausse), the
           sentence after "In Kerr this is no longer true": Re Q, the analogue of omega^2 - V, has a
           single minimum, at the light ring, and is positive elsewhere, so that outside the parabolic
           region the WKB waves propagate to the horizon and to infinity without mixing (the two
           steps of Sec. III B); the turning points z^2 = 4 eta_j lie inside the parabolic region.
Computes:  for the spins of Table I (a = 0, +-1/2, +-9/10, +-0.99), l = 100 and 400, j = 0 and 3,
           the real part of Q (the potential of Y'' + Q Y = 0 in the tortoise coordinate, kerrflux.py)
           on a uniform grid in r* from -60 M to r = 500 M: the number of local minima, the number of
           zeros (turning points), the position of the minimum relative to the light ring, the turning
           points in the parabolic variable z = (2k)^(1/4) x with k = d^2 Re Q/dr*^2 at r0 (Weber form
           Q = Q0 + k x^2/2), eta_num = -Re Q(r0)/sqrt(2k) against eta_j = j + 1/2 (reached with O(1/l)
           corrections, slowest for j = 3 at a = 0.99 prograde: 5.39 at l = 100, 3.79 at l = 400, against
           3.5), and Re Q at the two ends of the grid.
Inputs:    none (kerrflux.Mode builds Q and the spheroidal eigenvalue).
Outputs:   check_barrier_results.txt next to this script (one row per case), read by the companion
           notebook and by sympy_check/reproduce_paper.py.
Run:       python check_barrier.py          (about a minute; the l = 400 harmonics dominate)
"""
import os
import warnings

import numpy as np

from kerrflux import Mode

HERE = os.path.dirname(os.path.abspath(__file__))
SPINS = [0.0, 0.5, -0.5, 0.9, -0.9, 0.99, -0.99]
ELLS = [100, 400]
JS = [0, 3]


def r_of_rstar(md, rs):
    """Invert r*(r) by bisection (vectorised)."""
    lo = np.full_like(rs, md.rp + 1e-12)
    hi = np.full_like(rs, 1000.0)
    for _ in range(90):
        mid = 0.5*(lo + hi)
        below = md.rstar(mid) < rs
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return 0.5*(lo + hi)


def barrier_shape(a, l, j, rs_min=-60.0, r_max=500.0, n=40000):
    md = Mode(a, l, j)
    r0 = md.d['r0']
    rs0 = md.rstar(r0)
    rs = np.linspace(rs_min, md.rstar(r_max), n)
    r = r_of_rstar(md, rs)
    ReQ = np.real(md.Q(r))
    dsign = np.diff(np.sign(np.diff(ReQ)))
    n_min = int(np.sum(dsign > 0))
    zeros = np.where(np.diff(np.sign(ReQ)) != 0)[0]
    n_zero = len(zeros)
    r_min = r[np.argmin(ReQ)]
    # curvature in r* at the light ring: k = d^2 Re Q / dr*^2 (Weber form Q = Q0 + k x^2/2)
    h = 1e-2
    Qx = lambda x: np.real(md.Q(r_of_rstar(md, np.array([rs0 + x]))))[0]
    k = (Qx(h) - 2*Qx(0.0) + Qx(-h))/h**2
    z_tp = [(rs[i] - rs0)*(2*k)**0.25 for i in zeros]
    eta_num = -Qx(0.0)/np.sqrt(2*k)
    return dict(a=a, l=l, j=j, n_min=n_min, n_zero=n_zero, rmin_over_r0=r_min/r0, z_tp=z_tp,
                eta_num=eta_num, eta_j=j + 0.5, ReQ_left=ReQ[0], ReQ_right=ReQ[-1])


if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    rows = []
    print("   a     l  j  minima zeros  r_min/r0   z_turning            eta_num  j+1/2   ReQ(r*=-60)  ReQ(r=500)")
    for a in SPINS:
        for l in ELLS:
            for j in JS:
                d = barrier_shape(a, l, j)
                rows.append(d)
                print(f"{a:+5.2f} {l:5d} {j:d}   {d['n_min']:d}     {d['n_zero']:d}   {d['rmin_over_r0']:.6f}  "
                      f"{d['z_tp'][0]:+7.3f} {d['z_tp'][-1]:+7.3f}   {d['eta_num']:.4f}  {d['eta_j']:.1f}   "
                      f"{d['ReQ_left']:.3e}  {d['ReQ_right']:.3e}", flush=True)
    with open(os.path.join(HERE, 'check_barrier_results.txt'), 'w') as f:
        f.write("# Re Q of the Kerr radial equation Y'' + Q Y = 0 along the resonant modes (check_barrier.py):\n"
                "# a sign l j n_minima n_zeros rmin/r0 z_turn_1 z_turn_2 eta_num eta_j ReQ(r*=-60M) ReQ(r=500M)\n")
        for d in rows:
            f.write(f"{abs(d['a'])} {1 if d['a'] >= 0 else -1} {d['l']} {d['j']} {d['n_min']} {d['n_zero']} "
                    f"{d['rmin_over_r0']:.6f} {d['z_tp'][0]:.4f} {d['z_tp'][-1]:.4f} {d['eta_num']:.5f} {d['eta_j']:.1f} "
                    f"{d['ReQ_left']:.6e} {d['ReQ_right']:.6e}\n")
    shape_ok = all(d['n_min'] == 1 and d['n_zero'] == 2 for d in rows)
    # eta_num -> j + 1/2 with O(1/l) corrections (slowest at a = 0.99 prograde, where the O(1) remainder of
    # the eigenvalue expansion is largest): the deviation must shrink by more than a factor 2.5 from l = 100
    # to l = 400 whenever it exceeds 0.02 at l = 100
    conv_ok = True
    for a in SPINS:
        for j in JS:
            dev = {d['l']: abs(d['eta_num'] - d['eta_j']) for d in rows if d['a'] == a and d['j'] == j}
            if dev[100] > 0.02 and dev[400] > dev[100]/2.5:
                conv_ok = False
            if dev[400] > 0.3:
                conv_ok = False
    print("single minimum and two turning points in every case:", shape_ok)
    print("eta_num -> j + 1/2 with O(1/l) corrections (deviation at l = 400 below 0.3 and shrinking by > 2.5x from l = 100):", conv_ok)
    print("results written to check_barrier_results.txt")
