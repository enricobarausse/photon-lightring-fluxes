"""Production run of the Kerr Teukolsky fluxes: writes kerr_results.json.

Supports Sec. IV, Table I (last column) and Fig. 2 of "Gravitational radiation from a
photon on the light ring" (E. Barausse): the numerical fluxes per multipole for
a = 0, +-0.5, +-0.9 (negative sign = retrograde), l = 10, 20, 50, 100, 200, 400, 800 and
j = l - m = 0 ... 3, and the Richardson extrapolation from l = 400 and 800 quoted in the text.

Each record of the JSON list is {a, sign, l, j, FluxI, FluxH, lam, omega, kdom, ctail}:
a = |a| and sign = +1/-1 (prograde/retrograde), FluxI and FluxH the fluxes to infinity and
into the horizon per E^2 for the single mode m = l - j > 0 (M = 1; the physical flux per l
is twice the sum over j), lam the separation constant, omega = m Omega, and two diagnostics
of the spheroidal expansion (index of the dominant spherical coefficient, which must equal
j, and the largest of the last five coefficients, which must be negligible). The outer
boundary of the up integration is rB = 60 for l < 50 and rB = 40 otherwise.

    python run_production.py                    # the full run (about ten minutes)
    python run_production.py --cases 0.5 --ls 10 --jmax 0 --out test.json   # one mode

The file is rewritten after every mode, so a partial run is usable. A sanity check against
the Toolkit at l <= 5 is validate_kerr.py; summarize.py tabulates the results.
"""
import argparse
import json
import os
import time

import numpy as np

from kerrflux import Mode
from swsh import spheroidal

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CASES = [0.0, 0.5, -0.5, 0.9, -0.9]
DEFAULT_LS = [10, 20, 50, 100, 200, 400, 800]


def run(cases, ls, jmax, outfile):
    out = []
    t0 = time.time()
    for aa in cases:
        a, sgn = abs(aa), (-1 if aa < 0 else 1)
        for l in ls:
            for j in range(jmax + 1):
                t = time.time()
                try:
                    md = Mode(aa, l, j, rB=60.0 if l < 50 else 40.0)
                    # spectral convergence/ordering check of the spheroidal harmonic
                    lam, A, lv, c = spheroidal(-2, l, l - j, md.gam)
                    kdom = int(np.argmax(np.abs(c))); tail = float(np.abs(c[-5:]).max())
                    fI, fH = md.solve()
                    out.append(dict(a=a, sign=sgn, l=l, j=j, FluxI=fI, FluxH=fH, lam=md.lam, omega=md.w, kdom=kdom, ctail=tail))
                    print(f"a={a} sgn={sgn:+d} l={l} j={j}: l*I={l*fI:.8e} l*H={l*fH:.8e}  kdom={kdom} (expect {l-max(l-j,2)})  ctail={tail:.1e}  t={time.time()-t:.1f}s", flush=True)
                except Exception as e:
                    print(f"a={a} sgn={sgn:+d} l={l} j={j}: FAILED {e!r}", flush=True)
                json.dump(out, open(outfile, 'w'), indent=1)
    print("DONE", time.time() - t0, flush=True)
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('--cases', type=float, nargs='+', default=DEFAULT_CASES, help='spins (negative = retrograde)')
    p.add_argument('--ls', type=int, nargs='+', default=DEFAULT_LS, help='multipoles l')
    p.add_argument('--jmax', type=int, default=3, help='largest j = l - m (default 3)')
    p.add_argument('--out', default=os.path.join(HERE, 'kerr_results.json'), help='output JSON file')
    args = p.parse_args()
    run(args.cases, args.ls, args.jmax, args.out)
