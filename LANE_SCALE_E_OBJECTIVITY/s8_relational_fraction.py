"""S8 -- THE RELATIONAL FRACTION.

    Phi(F) = 1 - ( sum_i chi(R_i : F) ) / chi(whole register : F)

the fraction of everything the environment holds about the k-record register that is NOT
accounted for by the records one at a time.

WHY THIS IS THE RIGHT SCALAR, and what it is not:
  * it is a single number, not a count, so C-35's objection to counts does not apply;
  * it needs no threshold;
  * at k = 1 it is identically 0 -- there is nothing but the record;
  * in the SEPARATE control (k independent carriers on k disjoint baths) the joint Holevo is
    EXACTLY additive, so Phi = 0 at every k.  That is the D-15 control and it is printed in the
    same table;
  * Phi > 0 therefore says exactly: the environment's information about the register is not
    decomposable into per-record parts.

  * CAVEAT CARRIED IN THE TABLE, not buried: chi_ALL saturates at the bath's own capacity, and
    sum_i chi_i falls with k partly because k records share one site's Holevo capacity (C-36,
    conventional).  Phi -> 1 is therefore driven from BOTH ends and the raw chi_ALL and sum_i are
    printed beside it so the reader can see which end is moving.
"""
import sys, io, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")


def phi(k, nq, kind, lam, times=TIMES, seed=7):
    W = weights(kind, k, nq, seed=seed) if k > 1 else np.ones((1, nq))
    B = Broadcast(k, nq, W, lam, times=times)
    Stau = thermal_entropies(nq)
    specs = [(f"R{i}", [(0.5, {i: +1}), (0.5, {i: -1})]) for i in range(k)]
    subs = [(j,) for j in range(nq)] + [tuple(range(nq))]
    acc = {s: np.zeros(k + 1) for s in subs}
    for ti in range(len(times)):
        for s in subs:
            r = chi_batch(B, list(s), ti, specs)
            acc[s] += np.array([r[f"R{i}"] for i in range(k)] +
                               [max(r['_Savg'] - float(Stau[list(s)].sum()), 0.0)])
    for s in subs: acc[s] /= len(times)
    wh = acc[tuple(range(nq))]
    site = np.array([acc[(j,)] for j in range(nq)])
    out = dict(sum_wh=float(wh[:k].sum()), all_wh=float(wh[k]),
               sum_site=float(site[:, :k].sum(axis=1).mean()), all_site=float(site[:, k].mean()))
    out['phi_wh'] = 1 - out['sum_wh'] / out['all_wh'] if out['all_wh'] > 1e-12 else float('nan')
    out['phi_site'] = 1 - out['sum_site'] / out['all_site'] if out['all_site'] > 1e-12 else float('nan')
    return out


P("=" * 132)
P("S8  THE RELATIONAL FRACTION  Phi = 1 - (sum_i chi_i) / chi_ALL")
P("=" * 132)
P("")
P("k = 1 uses a single coupled record (the other record of [[4,2,2]] left uncoupled); k > 1 uses")
P("[[k+2,k,2]] with all k records coupled.  25 times in [1,13].  'sep' is the D-15 CONTROL:")
P("k independent carriers on k disjoint baths, where chi_ALL is additive by construction.")
P("")
for nq in (6, 4):
    for lam in (0.4, 0.8, 1.2):
        P("=" * 132)
        P(f"nq = {nq}, lam = {lam}")
        P(f"{'k':>3} | {'Phi whole crd':>13} {'Phi whole sym':>13} {'Phi whole SEP':>13} | "
          f"{'Phi site crd':>12} {'Phi site sym':>12} {'Phi site SEP':>12} | "
          f"{'chiALL wh crd':>13} {'sum_i wh crd':>12} | {'chiALL wh sep':>13} {'sum_i wh sep':>12}")
        P("-" * 132)
        for k in range(1, 11):
            r = {kind: phi(k, nq, kind, lam) for kind in ('crowded', 'sym', 'separate')}
            P(f"{k:>3} | {r['crowded']['phi_wh']:>13.6f} {r['sym']['phi_wh']:>13.6f} "
              f"{r['separate']['phi_wh']:>13.2e} | {r['crowded']['phi_site']:>12.6f} "
              f"{r['sym']['phi_site']:>12.6f} {r['separate']['phi_site']:>12.2e} | "
              f"{r['crowded']['all_wh']:>13.5f} {r['crowded']['sum_wh']:>12.5f} | "
              f"{r['separate']['all_wh']:>13.5f} {r['separate']['sum_wh']:>12.5f}")
        P("-" * 132)
        P("")

P("READ: fill from the numbers above.  Phi = 0 at k = 1 by definition and Phi = 0 at every k in")
P("      the SEPARATE control; whatever Phi does in the crowded and sym columns is the effect of")
P("      putting several records on the same environment.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s8_relational_fraction.txt",
     "w").write(OUT.getvalue())
