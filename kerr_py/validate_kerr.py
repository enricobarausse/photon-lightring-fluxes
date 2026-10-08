"""Validation of kerrflux.py against the Black Hole Perturbation Toolkit at l <= 5.

Computes:  compares the fluxes of Mode.solve() with the Teukolsky-package (MST, 40 digits) values
           stored in toolkit_reference.txt for every row with l <= 5. The Python fluxes agree with
           the Toolkit to 2e-6 ... 2e-4 at a = 0 and to 1e-4 ... 3e-3 for a != 0 at these small l (flux at
           infinity, limited by the WKB boundary data; the horizon flux to better than 1e-6); the agreement at l = 10 ... 50 with the Toolkit
           integrator of the companion notebook is better, and the production values at l >= 100
           are validated through the analytic limit (Sec. IV). The Toolkit cross-check of the null
           source quoted in Secs. II and V of the paper (Zerilli/Regge-Wheeler vs Teukolsky to 9-17
           digits at l <= 5) was done with the Toolkit's own solvers
           (toolkit_checks/schwarzschild_toolkit.wls), not with this code.
Inputs:    toolkit_reference.txt (columns a, sign, l, j, l*Edot_I, l*Edot_H; a = 0, 1/2, 9/10,
           sign = +1 prograde / -1 retrograde, single m = l - j). The reference values were
           produced by toolkit_checks/kerr_toolkit.wls.
Run:       python validate_kerr.py
Cost:      about ten seconds.
"""
import os

from kerrflux import Mode

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(HERE, 'toolkit_reference.txt')


def read_reference(path=REFERENCE, lmax=5):
    rows = []
    for line in open(path):
        p = line.split()
        if len(p) == 6 and p[0] in ('0', '1/2', '9/10'):
            a = eval(p[0]); rows.append((a, int(p[1]), int(p[2]), int(p[3]), float(p[4]), float(p[5])))
    return [r for r in rows if r[2] <= lmax]


if __name__ == "__main__":
    print(" a   sign l j |   l*I (py)      ref        rel  |   l*H (py)      ref        rel")
    for (a, sgn, l, j, rI, rH) in read_reference():
        aa = a*sgn   # retrograde represented by a -> -a, m > 0
        md = Mode(aa, l, j, rB=60.0); fI, fH = md.solve()
        print(f"{a:4.2f} {sgn:+d} {l} {j} | {l*fI:.6e} {rI:.6e} {l*fI/rI-1:+.1e} | {l*fH:.6e} {rH:.6e} {l*fH/rH-1:+.1e}")
