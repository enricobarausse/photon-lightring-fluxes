"""Tabulate kerr_results.json: fluxes per multipole for each spin, summed over j.

    python summarize.py [kerr_results.json]

For each case prints, per l, l*Edot_I and l*Edot_H of the j = 0 mode, the sums over the
stored j (single m, per E^2; multiply by 2 for +-m), their mean, and the asymmetry
sqrt(l) (I - H)/(I + H).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def summarize(path):
    res = json.load(open(path))
    cases = sorted(set((r['a'], r['sign']) for r in res))
    for (a, sgn) in cases:
        print(f"\n=== a={a} {'prograde' if sgn>0 else 'retrograde'}  (single m, units E^2; multiply by 2 for +-m) ===")
        print("  l   l*Edot_I(j=0)  l*Edot_H(j=0) |  l*sum_j Edot_I   l*sum_j Edot_H   mean    sqrt(l)(I-H)/(I+H)")
        for l in sorted(set(r['l'] for r in res if r['a'] == a and r['sign'] == sgn)):
            rows = [r for r in res if r['a'] == a and r['sign'] == sgn and r['l'] == l]
            j0 = [r for r in rows if r['j'] == 0][0]
            sI = sum(r['FluxI'] for r in rows); sH = sum(r['FluxH'] for r in rows)
            print(f"{l:5d}  {l*j0['FluxI']:.6e}  {l*j0['FluxH']:.6e} |  {l*sI:.6e}  {l*sH:.6e}  {l*(sI+sH)/2:.6e}  {np.sqrt(l)*(sI-sH)/(sI+sH):+.4f}")


if __name__ == "__main__":
    summarize(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'kerr_results.json'))
