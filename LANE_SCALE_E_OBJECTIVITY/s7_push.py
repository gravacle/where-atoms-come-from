"""S7 -- PUSH THE k-LADDER UNTIL SOMETHING STOPS IT, AND SAY WHAT STOPPED IT.

The Walsh route (s4_helpers.group_scan_chi_fast, checked against the explicit-mask route in s6)
costs O(2^k * k * d^2) per time instead of O(4^k * d^2), so the wall moves from k ~ 8 to memory:
the array of 2^k conditional fragment states, 2^k * d^2 complex numbers.  This script climbs the
ladder with a memory budget, prints every rung it reaches, and names the rung that stopped it.

Carrier at each rung is [[k+2, k, 2]], whose records are verified in s1 -- by the matrix route to
n = 10 and by the F_2 route (validated against the matrix route where both run) to n = 16.
"""
import sys, io, gc, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY")
from common import *
from s4_helpers import fwht_axis0, group_scan_chi_fast

OUT = io.StringIO()
def P(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True); OUT.write(s + "\n")

LAM = 0.8
BUDGET_GB = 0.5          # cap on the 2^k * d^2 conditional-state array (machine is shared)


def ent2(m):
    """Exact entropy of a stack of 2x2 density matrices, no eigensolver."""
    tr = np.real(m[..., 0, 0] + m[..., 1, 1])
    det = np.real(m[..., 0, 0] * m[..., 1, 1] - m[..., 0, 1] * m[..., 1, 0])
    disc = np.maximum(tr * tr - 4 * det, 0.0)
    s = np.sqrt(disc)
    a = np.clip((tr + s) / 2, 1e-300, None); b = np.clip((tr - s) / 2, 1e-300, None)
    return -(a * np.log2(a) + b * np.log2(b))


def scan(k, nq, W, lam, frags, times=TIMES, chunk=512):
    """Memory-lean Walsh group scan: no 2^k-sized concatenations."""
    B = Broadcast(k, nq, W, lam, times=times)
    nP = 2 ** k
    out = np.zeros((len(frags), nP))
    for fi, F in enumerate(frags):
        acc = np.zeros(nP)
        for ti in range(len(times)):
            rhoF = kron_sites(B.rho, F, ti)
            d = rhoF.shape[-1]
            Wc = fwht_axis0(rhoF.reshape(nP, d * d), k).reshape(nP, d, d) / nP
            del rhoF
            base = Wc[0].copy()
            S0 = float(ent2(base)) if d == 2 else float(entropies(base[None])[0])
            for c0 in range(1, nP, chunk):
                c1 = min(c0 + chunk, nP)
                blk = Wc[c0:c1]
                if d == 2:
                    sp = ent2(base[None] + blk); sm = ent2(base[None] - blk)
                else:
                    sp = entropies(base[None] + blk); sm = entropies(base[None] - blk)
                acc[c0:c1] += np.maximum(S0 - 0.5 * (sp + sm), 0.0)
            del Wc
            gc.collect()
        out[fi] = acc / len(times)
    return out


P("=" * 140)
P("S7  THE k-LADDER, PUSHED")
P("=" * 140)
P("")
P(f"memory budget for the 2^k x d x d conditional-state array: {BUDGET_GB} GB; lam = {LAM};")
P("25 times in [1,13]; bath energies as in common.py; beta = 2.0.")
P("")
P("SELF-CHECK: the memory-lean scanner must reproduce group_scan_chi_fast exactly.")
_d = 0.0
for _k, _nq in ((3, 4), (5, 6)):
    _W = weights('crowded', _k, _nq); _f = [[0], list(range(_nq))]
    _d = max(_d, float(np.abs(scan(_k, _nq, _W, LAM, _f, times=TIMES[:6])
                              - group_scan_chi_fast(_k, _nq, _W, LAM, _f, times=TIMES[:6])).max()))
P(f"            max |diff| = {_d:.3e} -> " + ("PASSED" if _d < 1e-10 else "FAILED, DRAW NO CONCLUSION"))
P("")

for nq in (6, 4):
    P("=" * 140)
    P(f"LADDER at nq = {nq}.  'crd' = generic crowded, 'sym' = equal magnitudes, 'sep' = CONTROL")
    P("=" * 140)
    P(f"{'k':>3} {'carrier':>14} {'geom':>5} | {'chi_i site':>10} {'chi_i whole':>11} {'E':>7} | "
      f"{'chi_d1 site':>11} {'chi_* site':>10} {'exc site':>9} {'d* site':>7} {'Nbet site':>9} | "
      f"{'chi_d1 wh':>9} {'chi_* wh':>9} {'exc wh':>8} {'d* wh':>6} {'ratio wh':>8}")
    P("-" * 152)
    stopped = None
    for k in (2, 4, 6, 8, 10, 12, 14, 16, 18):
        d_whole = 2 ** nq
        need = (2 ** k) * d_whole * d_whole * 16 / 1e9
        if need > BUDGET_GB:
            stopped = (k, need)
            P(f"{k:>3} {'[[%d,%d,2]]'%(k+2,k):>14} {'--':>5} | STOPPED: the 2^k x {d_whole} x {d_whole} "
              f"conditional-state array would be {need:.2f} GB > {BUDGET_GB} GB budget")
            break
        frags = [[j] for j in range(nq)] + [list(range(nq))]
        depth = np.array([bin(m).count('1') for m in range(2 ** k)])
        base = [1 << i for i in range(k)]
        for kind in ('crowded', 'sym', 'separate'):
            W = weights(kind, k, nq)
            g = scan(k, nq, W, LAM, frags)
            cs = g[:nq]; cw = g[nq]
            chi_site = float(np.mean([cs[:, m].mean() for m in base]))
            chi_wh = float(np.mean([cw[m] for m in base]))
            E = float(np.mean([nq * cs[:, m].mean() / cw[m] for m in base if cw[m] > 1e-12]))
            d1 = g[:, depth == 1].max(axis=1); st = g[:, 1:].max(axis=1)
            ds = depth[1 + g[:, 1:].argmax(axis=1)]
            nb = np.array([int((g[q, 1:] > d1[q] + 1e-12).sum()) for q in range(g.shape[0])])
            rat = st[nq] / d1[nq] if d1[nq] > 1e-12 else float('nan')
            P(f"{k:>3} {'[[%d,%d,2]]'%(k+2,k):>14} {kind[:3]:>5} | {chi_site:>10.6f} {chi_wh:>11.6f} "
              f"{E:>7.4f} | {d1[:nq].mean():>11.6f} {st[:nq].mean():>10.6f} "
              f"{(st-d1)[:nq].mean():>9.6f} {ds[:nq].mean():>7.3f} {nb[:nq].mean():>9.2f} | "
              f"{d1[nq]:>9.6f} {st[nq]:>9.6f} {st[nq]-d1[nq]:>8.6f} {ds[nq]:>6.1f} {rat:>8.4f}")
            del g
            gc.collect()
        P("-" * 152)
    if stopped:
        P(f"LADDER STOPPED at k = {stopped[0]} (nq = {nq}): {stopped[1]:.2f} GB required.")
    P("")

# ---------------------------------------------------------------- single-site only, pushed further
P("=" * 140)
P("SINGLE-SITE-ONLY LADDER.  A one-qubit fragment has d = 2, so the array is 2^k x 2 x 2 and the")
P("wall moves much further out.  This is the fragment that matters for objectivity anyway: it is")
P("the individual observer.")
P("=" * 140)
P(f"{'k':>3} {'carrier':>14} {'geom':>5} | {'chi_d1':>9} {'chi_star':>9} {'EXCESS':>9} {'depth*':>7} "
  f"{'Nbetter':>9} {'ratio':>8}")
P("-" * 92)
stopped = None
for k in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
    need = (2 ** k) * 4 * 16 / 1e9
    if need > BUDGET_GB:
        stopped = (k, need); break
    frags = [[j] for j in range(6)]
    depth = np.array([bin(m).count('1') for m in range(2 ** k)])
    for kind in ('crowded', 'sym', 'separate'):
        W = weights(kind, k, 6)
        g = scan(k, 6, W, LAM, frags, chunk=4096)
        d1 = g[:, depth == 1].max(axis=1); st = g[:, 1:].max(axis=1)
        ds = depth[1 + g[:, 1:].argmax(axis=1)]
        nb = np.array([int((g[q, 1:] > d1[q] + 1e-12).sum()) for q in range(g.shape[0])])
        rat = float(np.mean([st[q] / d1[q] if d1[q] > 1e-12 else np.nan for q in range(len(d1))]))
        P(f"{k:>3} {'[[%d,%d,2]]'%(k+2,k):>14} {kind[:3]:>5} | {d1.mean():>9.6f} {st.mean():>9.6f} "
          f"{(st-d1).mean():>9.6f} {ds.mean():>7.3f} {nb.mean():>9.2f} {rat:>8.4f}")
        del g; gc.collect()
    P("-" * 92)
if stopped:
    P(f"SINGLE-SITE LADDER STOPPED at k = {stopped[0]}: {stopped[1]:.2f} GB required for the")
    P(f"2^k x 2 x 2 array (and 2^{stopped[0]} entropy evaluations per time per site).")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_E_OBJECTIVITY/s7_push.txt", "w").write(OUT.getvalue())
